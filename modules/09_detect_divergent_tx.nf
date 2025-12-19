// ============================================================================
// detect_divergent_tx.nf — Statistical Divergent Transcription Detection
// ============================================================================
//
// Purpose:
//   Detects divergent transcription events from strand-specific 3' PRO-seq data
//   using Gaussian Mixture Model with FDR-controlled statistical filtering
//
// Features:
//   • Statistical approach with auto-calibrated thresholds
//   • Gaussian Mixture Model for scoring regions
//   • FDR-controlled filtering (Benjamini-Hochberg-like)
//   • Bayesian balance scoring with Beta-Binomial model
//   • Local background estimation and signal-to-noise ratios
//   • Feature extraction for machine learning-like classification
//
// Algorithm Overview:
//   1. Load strand-specific bedGraphs (pos/neg)
//   2. Auto-calibrate thresholds from background distribution
//      - threshold = 95th percentile of background
//      - sum_thr = 10x threshold
//   3. Call peak blocks on each strand (threshold + merge)
//   4. Pair peaks using relaxed initial criteria
//      - Distance window (default 1000bp)
//      - Optional balance filter (default: disabled for max sensitivity)
//   5. Extract features for each paired region:
//      - Total signal and strand-specific sums
//      - Bayesian balance score (Beta-Binomial model)
//      - Local background and signal-to-background ratio
//      - Region width and signal density
//   6. Fit 2-component Gaussian Mixture Model on features
//      - Identify "positive" (true divergent TX) component
//      - Compute posterior probabilities
//   7. Apply FDR control to filter regions
//   8. Output high-confidence divergent transcription sites
//
// Key Differences from earlier threshold-based approach:
//   • Uses statistics instead of hard thresholds
//   • Auto-calibration eliminates manual threshold tuning
//   • FDR control provides objective quality measure
//   • Confidence scores enable downstream prioritization
//   • More robust to varying coverage depths
//
// Track Type Recommendation:
//   • Primary alignments (unique mappers) for clean detection
//   • Avoids multimapper inflation in repetitive regions
//
// Inputs:
//   tuple(sample_id, pos_bg, neg_bg, condition, timepoint, replicate)
//
// Outputs:
//   ${params.output_dir}/06_divergent_tx/${sample_id}/
//     ├── divergent_transcription.bed       — Detected regions (BED5)
//     ├── divergent_transcription_qc.txt    — QC report with statistics
//     ├── divergent_summary.tsv             — Detection statistics
//     ├── README_divergent.txt              — Documentation
//     └── divergent.log                     — Processing log
//
// Output Format (BED5):
//   Column 1: chromosome       — Chromosome name
//   Column 2: start            — Region start (0-based)
//   Column 3: end              — Region end (exclusive)
//   Column 4: total_signal     — Combined signal from both strands
//   Column 5: confidence_score — Posterior probability from GMM (0-1)
//
// Parameters (params.advanced.*):
//   divergent_threshold    : Per-bin signal threshold (default: auto-calibrate)
//   divergent_sum_thr      : Minimum peak total signal (default: auto-calibrate)
//   divergent_fdr          : False discovery rate threshold (default: 0.05)
//   divergent_nt_window    : Max edge-to-edge gap for pairing (default: 1000bp)
//   divergent_balance      : Min balance ratio for initial pairing (default: 0.0)
//   divergent_bin_gap      : Max gap within peaks (default: 100bp)
//
// Dependencies:
//   • Python 3.7+
//   • numpy, pandas (required)
//   • scikit-learn, scipy (required for GMM and statistics)
//
// ============================================================================

nextflow.enable.dsl = 2

process detect_divergent_tx {

  tag        { sample_id }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/06_divergent_tx/${sample_id}",
             mode: 'copy',
             overwrite: true

  conda (params.conda_divergent ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(pos_bg), path(neg_bg),
          val(condition), val(timepoint), val(replicate)
    
    // Explicit parameter inputs for better cache control
    val threshold
    val sum_thr
    val fdr
    val nt_window
    val balance
    val bin_gap

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path("divergent_transcription.bed"),
          val(condition), val(timepoint), val(replicate),
          emit: bed

    path "divergent_transcription_qc.txt", optional: true, emit: qc_txt
    path "divergent_summary.tsv",                          emit: summary
    path "README_divergent.txt",                           emit: readme
    path "divergent.log",                                  emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -eu
  export LC_ALL=C

  # Redirect all output to log file
  exec > >(tee -a divergent.log) 2>&1
  
  # Trap SIGPIPE to avoid exit code 141
  trap '' PIPE

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "DIVERGENT | START | sample=!{sample_id} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="!{sample_id}"
  CONDITION="!{condition}"
  TIMEPOINT="!{timepoint}"
  REPLICATE="!{replicate}"
  THREADS=!{task.cpus}

  POS_BG="!{pos_bg}"
  NEG_BG="!{neg_bg}"
  
  DETECTOR_SCRIPT="!{projectDir}/bin/detect_divergent_transcription.py"

  # Detection parameters (passed as process inputs for cache control)
  THRESHOLD="!{threshold}"
  SUM_THR="!{sum_thr}"
  FDR=!{fdr}
  NT_WINDOW=!{nt_window}
  BALANCE=!{balance}
  BIN_GAP=!{bin_gap}

  # Feature flags
  DO_QC=$([[ "!{params.advanced?.divergent_qc}" == "false" ]] && echo 0 || echo 1)

  echo "DIVERGENT | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "DIVERGENT | CONFIG | Condition: ${CONDITION}"
  echo "DIVERGENT | CONFIG | Timepoint: ${TIMEPOINT}"
  echo "DIVERGENT | CONFIG | Replicate: ${REPLICATE}"
  echo "DIVERGENT | CONFIG | Threads: ${THREADS}"
  echo "DIVERGENT | CONFIG | Positive bedGraph: ${POS_BG}"
  echo "DIVERGENT | CONFIG | Negative bedGraph: ${NEG_BG}"
  echo ""
  echo "DIVERGENT | CONFIG | Detection Parameters (statistical):"
  echo "DIVERGENT | CONFIG |   Algorithm: Gaussian Mixture Model with FDR control"
  echo "DIVERGENT | CONFIG |   Threshold: ${THRESHOLD} (auto = 95th percentile)"
  echo "DIVERGENT | CONFIG |   Sum threshold: ${SUM_THR} (auto = 10x threshold)"
  echo "DIVERGENT | CONFIG |   FDR threshold: ${FDR}"
  echo "DIVERGENT | CONFIG |   Pairing window: ${NT_WINDOW} bp"
  echo "DIVERGENT | CONFIG |   Balance filter: ${BALANCE} (0.0 = disabled for max sensitivity)"
  echo "DIVERGENT | CONFIG |   Bin gap: ${BIN_GAP} bp"
  echo ""
  echo "DIVERGENT | CONFIG | Feature Extraction:"
  echo "DIVERGENT | CONFIG |   • Bayesian balance scoring"
  echo "DIVERGENT | CONFIG |   • Local background estimation"
  echo "DIVERGENT | CONFIG |   • Signal-to-background ratios"
  echo "DIVERGENT | CONFIG |   • Region width and density"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "DIVERGENT | VALIDATE | Checking inputs..."

  # Check Python script exists
  if [[ ! -f "${DETECTOR_SCRIPT}" ]]; then
    echo "DIVERGENT | ERROR | Detector script not found: ${DETECTOR_SCRIPT}"
    exit 1
  fi
  echo "DIVERGENT | VALIDATE | Detector script: ${DETECTOR_SCRIPT}"

  # Check input bedGraphs
  if [[ ! -s "${POS_BG}" ]]; then
    echo "DIVERGENT | ERROR | Positive bedGraph missing or empty: ${POS_BG}"
    exit 1
  fi

  if [[ ! -s "${NEG_BG}" ]]; then
    echo "DIVERGENT | ERROR | Negative bedGraph missing or empty: ${NEG_BG}"
    exit 1
  fi

  POS_SIZE=$(stat -c%s "${POS_BG}" 2>/dev/null || stat -f%z "${POS_BG}" 2>/dev/null || echo "unknown")
  NEG_SIZE=$(stat -c%s "${NEG_BG}" 2>/dev/null || stat -f%z "${NEG_BG}" 2>/dev/null || echo "unknown")
  POS_LINES=$(wc -l < "${POS_BG}" | tr -d ' ')
  NEG_LINES=$(wc -l < "${NEG_BG}" | tr -d ' ')

  echo "DIVERGENT | VALIDATE | Positive bedGraph: ${POS_SIZE} bytes (${POS_LINES} lines)"
  echo "DIVERGENT | VALIDATE | Negative bedGraph: ${NEG_SIZE} bytes (${NEG_LINES} lines)"

  # Validate tools and dependencies
  TOOLS_OK=1
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_VERSION=$(python3 --version 2>&1 || echo "unknown")
    echo "DIVERGENT | VALIDATE | Python: ${PYTHON_VERSION}"
    
    # Check Python dependencies (including sklearn and scipy for GMM)
    echo "DIVERGENT | VALIDATE | Checking Python dependencies..."
    python3 -c "import numpy, pandas, sklearn.mixture, scipy.stats" 2>/dev/null || {
      echo "DIVERGENT | ERROR | Missing Python dependencies (numpy, pandas, scikit-learn, scipy)"
      echo "DIVERGENT | ERROR | Install with: pip install numpy pandas scikit-learn scipy"
      TOOLS_OK=0
    }
  else
    echo "DIVERGENT | ERROR | python3 not found in PATH"
    TOOLS_OK=0
  fi

  if [[ ${TOOLS_OK} -eq 0 ]]; then
    echo "DIVERGENT | ERROR | Missing required tools or dependencies"
    exit 1
  fi

  ###########################################################################
  # 3) RUN STATISTICAL DIVERGENT TRANSCRIPTION DETECTOR
  ###########################################################################

  echo "DIVERGENT | DETECT | Running statistical detector..."
  echo "DIVERGENT | DETECT | Algorithm: Gaussian Mixture Model with FDR control"
  echo "DIVERGENT | DETECT | This may take several minutes for large genomes..."

  # Build command-line arguments
  THRESHOLD_ARG=""
  if [[ "${THRESHOLD}" != "auto" ]]; then
    THRESHOLD_ARG="--threshold ${THRESHOLD}"
  fi

  SUM_THR_ARG=""
  if [[ "${SUM_THR}" != "auto" ]]; then
    SUM_THR_ARG="--sum-thr ${SUM_THR}"
  fi

  QC_ARG=""
  [[ ${DO_QC} -eq 0 ]] && QC_ARG="--no-report"

  # Run detector (capture exit code)
  DETECT_START=$(date +%s)
  set +e
  python3 "${DETECTOR_SCRIPT}" \
    --sample       "${SAMPLE_ID}" \
    --pos          "${POS_BG}" \
    --neg          "${NEG_BG}" \
    --out          "divergent_transcription.bed" \
    --fdr          "${FDR}" \
    --nt-window    "${NT_WINDOW}" \
    --balance      "${BALANCE}" \
    --bin-gap      "${BIN_GAP}" \
    --ncores       "${THREADS}" \
    --write-summary "divergent_summary.tsv" \
    ${THRESHOLD_ARG} \
    ${SUM_THR_ARG} \
    ${QC_ARG}
  
  DETECT_RC=$?
  set -e
  
  DETECT_END=$(date +%s)
  DETECT_TIME=$((DETECT_END - DETECT_START))

  echo "DIVERGENT | DETECT | Detection completed in ${DETECT_TIME}s (exit code: ${DETECT_RC})"

  ###########################################################################
  # 4) HANDLE DETECTOR EXIT CODES
  ###########################################################################

  if [[ ${DETECT_RC} -ne 0 ]]; then
    echo "DIVERGENT | ERROR | Detector failed with exit code ${DETECT_RC}"
    echo "DIVERGENT | ERROR | Check logs above for details"
    exit ${DETECT_RC}
  fi

  echo "DIVERGENT | STATUS | Detection successful"

  # Ensure output files exist
  if [[ ! -e divergent_transcription.bed ]]; then
    echo "DIVERGENT | WARNING | Output BED missing, creating empty file"
    : > divergent_transcription.bed
  fi
  
  if [[ ! -e divergent_summary.tsv ]]; then
    echo "DIVERGENT | WARNING | Summary missing, creating default"
    cat > divergent_summary.tsv <<SUMMARYEOF
sample	n_pos_pk	n_neg_pk	n_pairs_raw	n_dt	wall_s
${SAMPLE_ID}	0	0	0	0	${DETECT_TIME}
SUMMARYEOF
  fi

  ###########################################################################
  # 5) PARSE DETECTION RESULTS
  ###########################################################################

  echo "DIVERGENT | RESULTS | Parsing detection results..."

  # Count detected regions
  DT_COUNT=$(grep -v '^#' divergent_transcription.bed 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  BED_SIZE=$(stat -c%s divergent_transcription.bed 2>/dev/null || stat -f%z divergent_transcription.bed 2>/dev/null || echo "unknown")

  echo "DIVERGENT | RESULTS | Detected regions: ${DT_COUNT}"
  echo "DIVERGENT | RESULTS | Output BED size: ${BED_SIZE} bytes"

  # Parse summary statistics if available
  if [[ -s divergent_summary.tsv && ${DT_COUNT} -gt 0 ]]; then
    # Extract statistics from summary (skip header)
    SUMMARY_LINE=$(tail -1 divergent_summary.tsv)
    
    N_POS_PK=$(echo "${SUMMARY_LINE}" | cut -f2 || echo "unknown")
    N_NEG_PK=$(echo "${SUMMARY_LINE}" | cut -f3 || echo "unknown")
    N_PAIRS_RAW=$(echo "${SUMMARY_LINE}" | cut -f4 || echo "unknown")
    WALL_TIME=$(echo "${SUMMARY_LINE}" | cut -f6 || echo "unknown")
    
    echo "DIVERGENT | RESULTS | Positive peaks: ${N_POS_PK}"
    echo "DIVERGENT | RESULTS | Negative peaks: ${N_NEG_PK}"
    echo "DIVERGENT | RESULTS | Raw pairs: ${N_PAIRS_RAW}"
    echo "DIVERGENT | RESULTS | Final regions (FDR-filtered): ${DT_COUNT}"
    echo "DIVERGENT | RESULTS | Retention rate: $(awk "BEGIN {printf \"%.1f%%\", (${DT_COUNT}/${N_PAIRS_RAW})*100}" 2>/dev/null || echo "N/A")"
    echo "DIVERGENT | RESULTS | Processing time: ${WALL_TIME}s"
  fi

  # Parse confidence scores if BED has 5 columns
  if [[ ${DT_COUNT} -gt 0 ]]; then
    FIRST_LINE=$(grep -v '^#' divergent_transcription.bed | head -1)
    COL_COUNT=$(echo "${FIRST_LINE}" | awk '{print NF}')
    
    if [[ ${COL_COUNT} -eq 5 ]]; then
      echo "DIVERGENT | RESULTS | Output format: BED5 (with confidence scores)"
      
      # Calculate score statistics
      SCORE_MIN=$(awk '{print $5}' divergent_transcription.bed | sort -n | head -1)
      SCORE_MAX=$(awk '{print $5}' divergent_transcription.bed | sort -n | tail -1)
      SCORE_MEAN=$(awk '{sum+=$5; n++} END {if(n>0) printf "%.4f", sum/n; else print "N/A"}' divergent_transcription.bed)
      
      echo "DIVERGENT | RESULTS | Confidence scores:"
      echo "DIVERGENT | RESULTS |   Min:  ${SCORE_MIN}"
      echo "DIVERGENT | RESULTS |   Mean: ${SCORE_MEAN}"
      echo "DIVERGENT | RESULTS |   Max:  ${SCORE_MAX}"
    fi
  fi

  ###########################################################################
  # 6) QC VALIDATION
  ###########################################################################

  if [[ ${DO_QC} -eq 1 ]]; then
    if [[ -s divergent_transcription_qc.txt ]]; then
      QC_SIZE=$(stat -c%s divergent_transcription_qc.txt 2>/dev/null || stat -f%z divergent_transcription_qc.txt 2>/dev/null || echo "unknown")
      echo "DIVERGENT | QC | QC report generated: ${QC_SIZE} bytes"
      echo "DIVERGENT | QC | Report includes: thresholds, peak counts, score distributions, feature summaries"
    else
      echo "DIVERGENT | QC | WARNING: QC report not generated (check for errors above)"
    fi
  fi

  ###########################################################################
  # 7) CREATE README
  ###########################################################################

  echo "DIVERGENT | README | Creating documentation..."

  cat > README_divergent.txt <<'DOCEOF'
================================================================================
DIVERGENT TRANSCRIPTION — !{sample_id} (Statistical)
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Statistical detection of divergent transcription regions from strand-specific
  3' PRO-seq coverage using Gaussian Mixture Models and FDR control.
  
  Divergent transcription occurs when RNA Polymerase II initiates in both
  directions from a genomic locus, creating paired transcription on opposite
  strands within close proximity. This is a hallmark of active promoters and
  enhancers.

STATISTICAL APPROACH
────────────────────────────────────────────────────────────────────────────
  Major algorithmic upgrade from threshold-based to statistical detection:
  
  1. Auto-Calibration:
     • Automatically determines optimal thresholds from background distribution
     • threshold = 95th percentile of background signal
     • sum_thr = 10x threshold
     • Eliminates need for manual parameter tuning
     • Adapts to varying sequencing depths
  
  2. Feature Extraction:
     • Total signal and strand-specific sums
     • Bayesian balance score using Beta-Binomial model
     • Local background estimation (±5kb windows)
     • Signal-to-background ratios
     • Region width and signal density
  
  3. Gaussian Mixture Model:
     • Fits 2-component GMM on extracted features
     • Identifies "positive" (true divergent TX) and "negative" (noise) components
     • Computes posterior probabilities (0-1 confidence scores)
     • Machine learning-like approach without manual labels
  
  4. FDR Control:
     • Applies Benjamini-Hochberg-like procedure
     • Controls false discovery rate at specified level (default: 5%)
     • Provides objective quality threshold
     • Ranks regions by confidence
  
  Benefits:
    ✓ More robust to varying coverage depths
    ✓ Fewer false positives in low-signal regions
    ✓ Confidence scores enable downstream prioritization
    ✓ No manual threshold tuning required
    ✓ Better balance sensitivity and specificity
    ✓ Statistically principled quality control

ALGORITHM
────────────────────────────────────────────────────────────────────────────
  1. Load strand-specific 3' bedGraphs (positive and negative strands)
  
  2. Auto-calibrate detection thresholds:
     - Sample 100K random bins to estimate background
     - Calculate mean, std, and percentiles
     - Set threshold = 95th percentile
     - Set sum_thr = 10x threshold
  
  3. Call peak blocks on each strand:
     - Filter bins by threshold
     - Merge consecutive bins within bin_gap
     - Keep blocks with total signal ≥ sum_thr
  
  4. Pair peaks with relaxed initial criteria:
     - For each positive peak, find negative peaks within nt_window
     - Calculate overlap and gap between peaks
     - Apply optional balance filter (default: disabled for max sensitivity)
     - Result: candidate paired regions
  
  5. Extract features for each paired region:
     - Query strand-specific signals in region
     - Calculate Bayesian balance score (Beta-Binomial model)
     - Estimate local background (±5kb excluding peak)
     - Compute signal-to-background ratio
     - Calculate region width and signal density
  
  6. Fit Gaussian Mixture Model:
     - Use features: log_total, balance_bayesian, log_snr
     - Fit 2-component GMM with full covariance
     - Identify positive component (higher mean log_total)
     - Compute posterior probabilities for all regions
  
  7. Apply FDR control:
     - Sort regions by posterior probability (descending)
     - Calculate cumulative false positives and FDR
     - Keep regions passing FDR threshold
     - Output with confidence scores

INPUT DATA
────────────────────────────────────────────────────────────────────────────
  Uses: Strand-specific 3' bedGraphs (unnormalized/raw recommended)
  
  Track Type:
    • Primary/Unique mappers (filtered BAM) for clean detection
    • Avoids multimapper inflation in repetitive regions
  
  Why RAW (unnormalized)?
    • Auto-calibration adapts to absolute signal levels
    • Consistent thresholds across varying depths
    • Normalization can distort statistical properties
  
  Strand Convention:
    • Positive strand: Values ≥ 0
    • Negative strand: Values ≤ 0 (negative values)

SAMPLE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Sample:     !{sample_id}
  Condition:  !{condition}
  Timepoint:  !{timepoint}
  Replicate:  !{replicate}

DETECTION PARAMETERS
────────────────────────────────────────────────────────────────────────────
  Threshold:            ${THRESHOLD}
    Per-bin signal minimum (auto = 95th percentile of background)
  
  Sum Threshold:        ${SUM_THR}
    Minimum total signal for peak blocks (auto = 10x threshold)
  
  FDR Threshold:        ${FDR}
    False discovery rate for filtering (e.g., 0.05 = 5% expected false positives)
  
  Pairing Window:       ${NT_WINDOW} bp
    Maximum distance for initial pairing (edge-to-edge or overlapping)
  
  Balance Filter:       ${BALANCE}
    Minimum balance for initial pairing (0.0 = disabled for max sensitivity)
    Note: Final filtering done statistically via GMM
  
  Bin Gap:              ${BIN_GAP} bp
    Maximum gap to merge bins into peak blocks

PROCESSING OPTIONS
────────────────────────────────────────────────────────────────────────────
  Threads:              ${THREADS}
  QC Report:            $([ ${DO_QC} -eq 1 ] && echo "Generated" || echo "Disabled")

FILES
────────────────────────────────────────────────────────────────────────────
  divergent_transcription.bed     — Detected regions (BED5 format)
  divergent_transcription_qc.txt  — QC report with detailed statistics
  divergent_summary.tsv           — Detection summary (TSV)
  README_divergent.txt            — This documentation
  divergent.log                   — Complete processing log

OUTPUT FORMAT (BED5)
────────────────────────────────────────────────────────────────────────────
  Column 1: chromosome       — Chromosome name
  Column 2: start            — Region start (0-based)
  Column 3: end              — Region end (exclusive)
  Column 4: total_signal     — Combined signal from both strands
  Column 5: confidence_score — Posterior probability from GMM (0-1)

Confidence Score Interpretation:
  • Score close to 1.0: High confidence true divergent transcription
  • Score close to 0.5: Uncertain, borderline case
  • Score close to 0.0: Likely noise/artifact (should not appear in output)
  
  All output regions collectively pass FDR control, so even lower-scoring
  regions in the output are expected to be true positives at the specified
  FDR level.

DETECTION RESULTS
────────────────────────────────────────────────────────────────────────────
  Detected regions: ${DT_COUNT}
  $([ -s divergent_summary.tsv ] && [ ${DT_COUNT} -gt 0 ] && cat <<STATS
  Positive peaks:   ${N_POS_PK:-unknown}
  Negative peaks:   ${N_NEG_PK:-unknown}
  Raw pairs:        ${N_PAIRS_RAW:-unknown}
  Final (filtered): ${DT_COUNT}
  Retention rate:   $(awk "BEGIN {printf \"%.1f%%\", (${DT_COUNT}/${N_PAIRS_RAW:-1})*100}" 2>/dev/null || echo "N/A")
  Processing time:  ${WALL_TIME:-unknown}s
STATS
)

QUALITY CONTROL
────────────────────────────────────────────────────────────────────────────

Expected Results:
  • Typical range: 100-10,000 divergent regions per sample
  • Regions cluster near active promoters and enhancers
  • Width typically 200-2000 bp
  • Confidence scores should be well-distributed (avoid all near 1.0)

Interpretation:
  • High count: Active transcription, many promoters/enhancers
  • Low count: Inactive sample, low coverage, or stringent FDR
  • Very high count (>50,000): Check for artifacts or very deep sequencing
  
  Confidence Score Distribution:
    • Mean score ~0.7-0.9: Good separation between signal and noise
    • Mean score >0.95: Possible overfitting or very clean data
    • Mean score <0.6: Possible underfitting or noisy data

Troubleshooting:
  • Zero regions: Lower FDR (e.g., 0.1) or threshold
  • Too many regions: Increase FDR stringency (e.g., 0.01)
  • Low confidence scores: Check input quality, may need more coverage
  • Model convergence issues: Try manual thresholds instead of auto

Parameter Tuning:
  
  More Sensitive (detect more):
    • Increase FDR (e.g., 0.1)
    • Decrease threshold (if not using auto)
    • Decrease sum_thr (if not using auto)
  
  More Specific (higher confidence):
    • Decrease FDR (e.g., 0.01)
    • Use manual thresholds (increase both)
  
  For most applications, auto-calibration with FDR=0.05 works well.

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  These divergent transcription regions can be used for:
  
  1. Promoter Annotation:
     - Identify active promoters via divergent transcription signatures
     - Prioritize by confidence score
  
  2. Enhancer Discovery:
     - Enhancers show divergent transcription
     - Filter for distal regions (not overlapping promoters)
  
  3. Differential Analysis:
     - Compare region counts and scores across conditions
     - Use confidence scores to weight importance
  
  4. Prioritization:
     - Rank regions by confidence_score for follow-up experiments
     - High-score regions more likely to be true functional elements

TECHNICAL NOTES
────────────────────────────────────────────────────────────────────────────
  • Requires: numpy, pandas, scikit-learn, scipy
  • Uses vectorized operations for speed
  • Chromosome-by-chromosome processing
  • Memory-efficient streaming I/O
  • Supports gzipped input bedGraphs

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Date: $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample: !{sample_id}
  Module: detect_divergent_tx (Statistical)
  Detector: detect_divergent_transcription.py v1.0

================================================================================
DOCEOF

  echo "DIVERGENT | README | Documentation created"

  ###########################################################################
  # 8) VALIDATION
  ###########################################################################

  echo "DIVERGENT | VALIDATE | Validating outputs..."

  VALIDATION_OK=1

  # Check BED file exists
  if [[ ! -e divergent_transcription.bed ]]; then
    echo "DIVERGENT | ERROR | Output BED file missing"
    VALIDATION_OK=0
  fi

  # Check summary exists
  if [[ ! -s divergent_summary.tsv ]]; then
    echo "DIVERGENT | ERROR | Summary file missing or empty"
    VALIDATION_OK=0
  fi

  # Validate BED format if non-empty
  if [[ -s divergent_transcription.bed ]]; then
    # Check that we have 5 columns (BED5 with scores)
    FIRST_LINE=$(grep -v '^#' divergent_transcription.bed | head -1 || true)
    if [[ -n "${FIRST_LINE}" ]]; then
      COL_COUNT=$(echo "${FIRST_LINE}" | awk '{print NF}')
      if [[ ${COL_COUNT} -ne 5 ]]; then
        echo "DIVERGENT | WARNING | BED file should have 5 columns, found ${COL_COUNT}"
      fi
      
      # Check coordinates
      START=$(echo "${FIRST_LINE}" | cut -f2)
      END=$(echo "${FIRST_LINE}" | cut -f3)
      if [[ ${START} -ge ${END} ]]; then
        echo "DIVERGENT | ERROR | Invalid BED coordinates: start >= end"
        VALIDATION_OK=0
      fi
      
      # Check score is numeric and in [0,1]
      SCORE=$(echo "${FIRST_LINE}" | cut -f5)
      if ! awk -v s="${SCORE}" 'BEGIN {exit !(s >= 0 && s <= 1)}'; then
        echo "DIVERGENT | WARNING | Confidence score should be in [0,1], found ${SCORE}"
      fi
    fi
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "DIVERGENT | ERROR | Validation failed"
    exit 1
  fi

  echo "DIVERGENT | VALIDATE | All outputs validated"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "DIVERGENT | SUMMARY | Sample: ${SAMPLE_ID}"
  echo "DIVERGENT | SUMMARY | Algorithm: Gaussian Mixture Model + FDR"
  echo "DIVERGENT | SUMMARY | Detected regions: ${DT_COUNT}"
  echo "DIVERGENT | SUMMARY | Output size: ${BED_SIZE} bytes"
  if [[ -s divergent_summary.tsv && ${DT_COUNT} -gt 0 ]]; then
    echo "DIVERGENT | SUMMARY | Positive peaks: ${N_POS_PK}"
    echo "DIVERGENT | SUMMARY | Negative peaks: ${N_NEG_PK}"
    echo "DIVERGENT | SUMMARY | Raw pairs: ${N_PAIRS_RAW}"
    echo "DIVERGENT | SUMMARY | Retention rate: $(awk "BEGIN {printf \"%.1f%%\", (${DT_COUNT}/${N_PAIRS_RAW})*100}" 2>/dev/null || echo "N/A")"
    
    if [[ ${COL_COUNT:-0} -eq 5 ]]; then
      echo "DIVERGENT | SUMMARY | Score range: [${SCORE_MIN}, ${SCORE_MAX}]"
      echo "DIVERGENT | SUMMARY | Mean confidence: ${SCORE_MEAN}"
    fi
  fi
  echo "DIVERGENT | SUMMARY | Detection time: ${DETECT_TIME}s"
  if [[ ${DO_QC} -eq 1 && -s divergent_transcription_qc.txt ]]; then
    echo "DIVERGENT | SUMMARY | QC report: generated"
  fi
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "DIVERGENT | COMPLETE | sample=${SAMPLE_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}
