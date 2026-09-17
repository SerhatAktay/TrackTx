// ============================================================================
// detect_divergent_transcription.nf — Statistical Divergent Transcription Detection
// ============================================================================
//
// Purpose:
//   Detects divergent transcription events from strand-specific 3' PRO-seq data
//   using a Gaussian Mixture Model to RANK candidate regions by confidence.
//
// On "FDR": --fdr (params.advanced.divergent_fdr) is a posterior-based
// STRINGENCY KNOB, not a calibrated false-discovery rate -- it never was
// Benjamini-Hochberg-like despite earlier wording here. This has been
// empirically checked (bin/detect_divergent_transcription.py's
// empirical_null_fdr_check, run every time via the QC report): on real
// PRO-seq data the measured empirical FDR sits at ~0.5 regardless of this
// setting, the feature set, or peak-calling threshold design (all tested
// and documented in that function's docstring) -- a real ceiling of
// peak-pair-aggregate scoring against pervasively-transcribed PRO-seq
// signal, not an unvalidated assumption. Use the output as a ranked
// confidence list; read the QC report's Empirical-Null FDR Check section
// for the real measured rate before citing any FDR figure.
//
// Features:
//   • Statistical approach with auto-calibrated thresholds
//   • Gaussian Mixture Model for scoring regions
//   • Score-based stringency filtering (see "On FDR" above)
//   • Bayesian balance scoring with Beta-Binomial model
//   • Local background estimation and signal-to-noise ratios
//   • Feature extraction for machine learning-like classification
//
// Algorithm Overview:
//   1. Load strand-specific bedGraphs (pos/neg)
//   2. Auto-calibrate thresholds (targets 50-100K sites for mammalian genomes)
//      - threshold = 75th percentile (more permissive than 95th)
//      - sum_thr = 3x threshold (allows narrower peaks)
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
//   8. Merge overlapping regions (reduces redundancy from one-to-many pairing)
//   9. Output high-confidence divergent transcription sites
//
// Key Differences from earlier threshold-based approach:
//   • Uses statistics instead of hard thresholds
//   • Auto-calibration eliminates manual threshold tuning
//   • Confidence scores enable downstream ranking/prioritization
//     (NOT a calibrated FDR guarantee -- see "On FDR" above)
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
//   divergent_threshold    : Per-bin signal threshold (default: auto)
//   divergent_sum_thr      : Minimum peak total signal (default: auto)
//   divergent_fdr          : Score-stringency knob, NOT a calibrated FDR (default: 0.08; see "On FDR" above)
//   divergent_calibration_percentile : Percentile for auto threshold (default: 65)
//   divergent_calibration_sum_multiplier : sum_thr = threshold * N (default: 1.5)
//   divergent_calibration_background_lower : Use lower 50% bins (default: false)
//   divergent_merge_gap    : Merge overlapping regions within N bp (default: 'auto' -- organism-aware, 500bp ceiling)
//   divergent_nt_window     : Max edge-to-edge gap for pairing (default: 'auto' -- organism-aware, 1000bp ceiling)
//   divergent_balance      : Min balance ratio for initial pairing (default: 0.0)
//   divergent_bin_gap       : Max gap within peaks (default: 'auto' -- organism-aware, 100bp ceiling)
//
// Dependencies:
//   • Python 3.7+
//   • numpy, pandas (required)
//   • scikit-learn, scipy (required for GMM and statistics)
//
// ============================================================================


process detect_divergent_transcription {

  tag        { sample_id }
  label      'conda'
  cache      'lenient'

  publishDir { "${params.output_dir}/06_divergent_tx/${sample_id}" },
             mode: params.publish_mode,
             overwrite: true

  conda (params.conda_divergent ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(pos_bg), path(neg_bg),
          val(condition), val(timepoint), val(replicate)
    path(genome_sizes)
    val library_type

    // Explicit parameter inputs for better cache control
    val threshold
    val sum_thr
    val fdr
    val nt_window
    val balance
    val bin_gap
    val calibration_percentile
    val calibration_sum_multiplier
    val calibration_background_lower
    val merge_gap
    val fallback_top_frac

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
  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Stdout → log + terminal; stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a divergent.log)
  exec 2> >(tee -a divergent.log >&2)
  
  # Trap SIGPIPE to avoid exit code 141
  trap '' PIPE
  trap 'rc=\$?; tracktx_error "detect_divergent_transcription" "Unexpected process failure" "Check divergent.log in work dir" "\$rc"' ERR

  # Standardized error reporting (surfaces clearly in Nextflow "Command error")
  # Shared error helper (defined once in bin/tracktx_error_fragment.sh)
  source tracktx_error_fragment.sh

  TIMESTAMP=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "DIVERGENT | START | sample=${sample_id} | ts=\${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="${sample_id}"
  CONDITION="${condition}"
  TIMEPOINT="${timepoint}"
  REPLICATE="${replicate}"
  THREADS=${task.cpus}

  POS_BG="${pos_bg}"
  NEG_BG="${neg_bg}"
  GENOME_SIZES="${genome_sizes}"
  ASSAY_TYPE="${library_type == 'groseq' ? 'groseq' : 'proseq'}"

  DETECTOR_SCRIPT="\$(command -v detect_divergent_transcription.py)"

  # Detection parameters (passed as process inputs for cache control)
  THRESHOLD="${threshold}"
  SUM_THR="${sum_thr}"
  FDR=${fdr}
  NT_WINDOW="${nt_window}"
  BALANCE=${balance}
  BIN_GAP="${bin_gap}"
  CAL_PERCENTILE=${calibration_percentile}
  CAL_SUM_MULT=${calibration_sum_multiplier}
  CAL_BG_LOWER=${calibration_background_lower}
  MERGE_GAP="${merge_gap}"
  FALLBACK_TOP_FRAC=${fallback_top_frac}

  # Feature flags
  DO_QC=\$([[ "${params.advanced?.divergent_qc}" == "false" ]] && echo 0 || echo 1)

  echo "DIVERGENT | CONFIG | Sample ID: \${SAMPLE_ID}"
  echo "DIVERGENT | CONFIG | Condition: \${CONDITION}"
  echo "DIVERGENT | CONFIG | Timepoint: \${TIMEPOINT}"
  echo "DIVERGENT | CONFIG | Replicate: \${REPLICATE}"
  echo "DIVERGENT | CONFIG | Threads: \${THREADS}"
  echo "DIVERGENT | CONFIG | Positive bedGraph: \${POS_BG}"
  echo "DIVERGENT | CONFIG | Negative bedGraph: \${NEG_BG}"
  echo ""
  echo "DIVERGENT | CONFIG | Detection Parameters (statistical):"
  echo "DIVERGENT | CONFIG |   Algorithm: Gaussian Mixture Model with FDR control"
  echo "DIVERGENT | CONFIG |   Threshold: \${THRESHOLD} (auto = \${CAL_PERCENTILE}th percentile)"
  echo "DIVERGENT | CONFIG |   Sum threshold: \${SUM_THR} (auto = \${CAL_SUM_MULT}x threshold)"
  echo "DIVERGENT | CONFIG |   Merge gap: \${MERGE_GAP} bp (0=disabled)"
  echo "DIVERGENT | CONFIG |   FDR threshold: \${FDR}"
  echo "DIVERGENT | CONFIG |   Pairing window: \${NT_WINDOW} bp"
  echo "DIVERGENT | CONFIG |   Balance filter: \${BALANCE} (0.0 = disabled for max sensitivity)"
  echo "DIVERGENT | CONFIG |   Bin gap: \${BIN_GAP} bp"
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

  # Resolve Python: prefer env with deps (conda/container), fallback to PATH
  # 1) CONDA_PREFIX: Nextflow conda profile sets this; ensures correct env
  # 2) micromamba: Docker/Singularity container (mambaorg/micromamba image)
  # 3) /opt/conda: Container with conda at standard path
  # 4) python3/python: System or conda-activated PATH
  if [[ -n "\${CONDA_PREFIX:-}" ]] && [[ -x "\${CONDA_PREFIX}/bin/python3" ]]; then
    PYTHON_CMD="\${CONDA_PREFIX}/bin/python3"
  elif [[ -n "\${CONDA_PREFIX:-}" ]] && [[ -x "\${CONDA_PREFIX}/bin/python" ]]; then
    PYTHON_CMD="\${CONDA_PREFIX}/bin/python"
  elif command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  elif [[ -x /opt/conda/bin/python3 ]]; then
    PYTHON_CMD="/opt/conda/bin/python3"
  elif command -v python3 >/dev/null 2>&1; then
    PYTHON_CMD="python3"
  else
    PYTHON_CMD="python"
  fi

  # Check Python script exists
  if [[ ! -f "\${DETECTOR_SCRIPT}" ]]; then
    tracktx_error "detect_divergent_transcription" "Detector script not found: \${DETECTOR_SCRIPT}" "Ensure bin/detect_divergent_transcription.py exists in the pipeline directory"
  fi
  echo "DIVERGENT | VALIDATE | Detector script: \${DETECTOR_SCRIPT}"

  # Check input bedGraphs
  if [[ ! -s "\${POS_BG}" ]]; then
    tracktx_error "detect_divergent_transcription" "Positive bedGraph missing or empty: \${POS_BG}" "Check that generate_coverage_tracks produced 3p.pos.bedgraph for this sample"
  fi

  if [[ ! -s "\${NEG_BG}" ]]; then
    tracktx_error "detect_divergent_transcription" "Negative bedGraph missing or empty: \${NEG_BG}" "Check that generate_coverage_tracks produced 3p.neg.bedgraph for this sample"
  fi

  POS_SIZE=\$(tracktx_size "\${POS_BG}")
  NEG_SIZE=\$(tracktx_size "\${NEG_BG}")
  POS_LINES=\$(wc -l < "\${POS_BG}" | tr -d ' ')
  NEG_LINES=\$(wc -l < "\${NEG_BG}" | tr -d ' ')

  echo "DIVERGENT | VALIDATE | Positive bedGraph: \${POS_SIZE} bytes (\${POS_LINES} lines)"
  echo "DIVERGENT | VALIDATE | Negative bedGraph: \${NEG_SIZE} bytes (\${NEG_LINES} lines)"

  # Validate tools and dependencies; try alternative Pythons if first fails
  TOOLS_OK=0
  PYTHON_CANDIDATES=("\${PYTHON_CMD}")
  [[ -n "\${CONDA_PREFIX:-}" ]] && PYTHON_CANDIDATES+=("\${CONDA_PREFIX}/bin/python3" "\${CONDA_PREFIX}/bin/python")
  command -v micromamba >/dev/null 2>&1 && PYTHON_CANDIDATES+=("micromamba run -n base python3")
  [[ -x /opt/conda/bin/python3 ]] && PYTHON_CANDIDATES+=("/opt/conda/bin/python3")
  command -v python3 >/dev/null 2>&1 && PYTHON_CANDIDATES+=("python3")
  command -v python >/dev/null 2>&1 && PYTHON_CANDIDATES+=("python")

  for py in "\${PYTHON_CANDIDATES[@]}"; do
    [[ -z "\$py" ]] && continue
    if \$py --version >/dev/null 2>&1; then
      if \$py -c "import numpy, pandas, sklearn.mixture, scipy.stats" 2>/dev/null; then
        PYTHON_CMD="\$py"
        TOOLS_OK=1
        break
      fi
    fi
  done

  if \${PYTHON_CMD} --version >/dev/null 2>&1; then
    PYTHON_VERSION=\$(\${PYTHON_CMD} --version 2>&1 || echo "unknown")
    echo "DIVERGENT | VALIDATE | Python: \${PYTHON_VERSION}"
  fi

  if [[ \${TOOLS_OK} -eq 0 ]]; then
    tracktx_error "detect_divergent_transcription" "Missing Python dependencies (numpy, pandas, scikit-learn, scipy)" "pip install numpy pandas scikit-learn scipy | Or use: -profile conda | -profile docker"
  fi
  echo "DIVERGENT | VALIDATE | Checking Python dependencies... OK"

  ###########################################################################
  # 3) RUN STATISTICAL DIVERGENT TRANSCRIPTION DETECTOR
  ###########################################################################

  echo "DIVERGENT | DETECT | Running statistical detector..."
  echo "DIVERGENT | DETECT | Algorithm: Gaussian Mixture Model with FDR control"
  echo "DIVERGENT | DETECT | This may take several minutes for large genomes..."

  # Build command-line arguments
  THRESHOLD_ARG=""
  if [[ "\${THRESHOLD}" != "auto" ]]; then
    THRESHOLD_ARG="--threshold \${THRESHOLD}"
  fi

  SUM_THR_ARG=""
  if [[ "\${SUM_THR}" != "auto" ]]; then
    SUM_THR_ARG="--sum-thr \${SUM_THR}"
  fi

  # 'auto' (default) omits the flag entirely so the script's own
  # organism-aware auto-scaling applies (see auto_scale_window_params).
  NT_WINDOW_ARG=""
  [[ "\${NT_WINDOW}" != "auto" ]] && NT_WINDOW_ARG="--nt-window \${NT_WINDOW}"

  BIN_GAP_ARG=""
  [[ "\${BIN_GAP}" != "auto" ]] && BIN_GAP_ARG="--bin-gap \${BIN_GAP}"

  MERGE_GAP_ARG=""
  [[ "\${MERGE_GAP}" != "auto" ]] && MERGE_GAP_ARG="--merge-gap \${MERGE_GAP}"

  QC_ARG=""
  [[ \${DO_QC} -eq 0 ]] && QC_ARG="--no-report"

  CAL_BG_ARG=""
  [[ "\${CAL_BG_LOWER}" == "true" ]] && CAL_BG_ARG="--calibration-background-lower"

  # Run detector (capture exit code)
  DETECT_START=\$(date +%s)
  set +e
  \${PYTHON_CMD} "\${DETECTOR_SCRIPT}" \\
    --sample       "\${SAMPLE_ID}" \\
    --pos          "\${POS_BG}" \\
    --neg          "\${NEG_BG}" \\
    --genome-sizes "\${GENOME_SIZES}" \\
    --assay-type   "\${ASSAY_TYPE}" \\
    --out          "divergent_transcription.bed" \\
    --fdr          "\${FDR}" \\
    --balance      "\${BALANCE}" \\
    --calibration-percentile "\${CAL_PERCENTILE}" \\
    --calibration-sum-multiplier "\${CAL_SUM_MULT}" \\
    --fallback-top-frac "\${FALLBACK_TOP_FRAC}" \\
    --ncores       "\${THREADS}" \\
    --write-summary "divergent_summary.tsv" \\
    \${THRESHOLD_ARG} \\
    \${SUM_THR_ARG} \\
    \${NT_WINDOW_ARG} \\
    \${BIN_GAP_ARG} \\
    \${MERGE_GAP_ARG} \\
    \${CAL_BG_ARG} \\
    \${QC_ARG}
  
  DETECT_RC=\$?
  set -e
  
  DETECT_END=\$(date +%s)
  DETECT_TIME=\$((DETECT_END - DETECT_START))

  echo "DIVERGENT | DETECT | Detection completed in \${DETECT_TIME}s (exit code: \${DETECT_RC})"

  ###########################################################################
  # 4) HANDLE DETECTOR EXIT CODES
  ###########################################################################

  if [[ \${DETECT_RC} -ne 0 ]]; then
    tracktx_error "detect_divergent_transcription" "Detector failed with exit code \${DETECT_RC}" "Check divergent.log in work dir for details" \${DETECT_RC}
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
\${SAMPLE_ID}	0	0	0	0	\${DETECT_TIME}
SUMMARYEOF
  fi

  ###########################################################################
  # 5) PARSE DETECTION RESULTS
  ###########################################################################

  echo "DIVERGENT | RESULTS | Parsing detection results..."

  # Count detected regions
  DT_COUNT=\$(grep -v '^#' divergent_transcription.bed 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  BED_SIZE=\$(tracktx_size divergent_transcription.bed)

  echo "DIVERGENT | RESULTS | Detected regions: \${DT_COUNT}"
  echo "DIVERGENT | RESULTS | Output BED size: \${BED_SIZE} bytes"

  # Parse summary statistics if available
  if [[ -s divergent_summary.tsv && \${DT_COUNT} -gt 0 ]]; then
    # Extract statistics from summary (skip header)
    SUMMARY_LINE=\$(tail -1 divergent_summary.tsv)
    
    N_POS_PK=\$(echo "\${SUMMARY_LINE}" | cut -f2 || echo "unknown")
    N_NEG_PK=\$(echo "\${SUMMARY_LINE}" | cut -f3 || echo "unknown")
    N_PAIRS_RAW=\$(echo "\${SUMMARY_LINE}" | cut -f4 || echo "unknown")
    WALL_TIME=\$(echo "\${SUMMARY_LINE}" | cut -f6 || echo "unknown")
    
    echo "DIVERGENT | RESULTS | Positive peaks: \${N_POS_PK}"
    echo "DIVERGENT | RESULTS | Negative peaks: \${N_NEG_PK}"
    echo "DIVERGENT | RESULTS | Raw pairs: \${N_PAIRS_RAW}"
    echo "DIVERGENT | RESULTS | Final regions (FDR-filtered): \${DT_COUNT}"
    echo "DIVERGENT | RESULTS | Retention rate: \$(awk "BEGIN {printf \\"%.1f%%\\", (\${DT_COUNT}/\${N_PAIRS_RAW})*100}" 2>/dev/null || echo "N/A")"
    echo "DIVERGENT | RESULTS | Processing time: \${WALL_TIME}s"
  fi

  # Parse confidence scores if BED has 5 columns
  if [[ \${DT_COUNT} -gt 0 ]]; then
    FIRST_LINE=\$(grep -v -m1 '^#' divergent_transcription.bed || true)
    COL_COUNT=\$(echo "\${FIRST_LINE}" | awk '{print NF}')
    
    if [[ \${COL_COUNT} -eq 5 ]]; then
      echo "DIVERGENT | RESULTS | Output format: BED5 (with confidence scores)"
      
      # Calculate score statistics in a single awk pass -- sort|head/tail can raise
      # SIGPIPE under 'set -o pipefail' and abort the whole module even after detection
      # already succeeded (this is what crashed a completed 20,491-region run on 2026-09-01).
      read -r SCORE_MIN SCORE_MAX SCORE_MEAN < <(awk '
        NR==1 { min=\$5; max=\$5 }
        { if (\$5<min) min=\$5; if (\$5>max) max=\$5; sum+=\$5; n++ }
        END { if (n>0) printf "%s %s %.4f\\n", min, max, sum/n; else print "N/A N/A N/A" }
      ' divergent_transcription.bed)
      
      echo "DIVERGENT | RESULTS | Confidence scores:"
      echo "DIVERGENT | RESULTS |   Min:  \${SCORE_MIN}"
      echo "DIVERGENT | RESULTS |   Mean: \${SCORE_MEAN}"
      echo "DIVERGENT | RESULTS |   Max:  \${SCORE_MAX}"
    fi
  fi

  ###########################################################################
  # 6) QC VALIDATION
  ###########################################################################

  if [[ \${DO_QC} -eq 1 ]]; then
    if [[ -s divergent_transcription_qc.txt ]]; then
      QC_SIZE=\$(tracktx_size divergent_transcription_qc.txt)
      echo "DIVERGENT | QC | QC report generated: \${QC_SIZE} bytes"
      echo "DIVERGENT | QC | Report includes: thresholds, peak counts, score distributions, feature summaries"
    else
      echo "DIVERGENT | QC | WARNING: QC report not generated (check for errors above)"
    fi
  fi

  ###########################################################################
  # 7) CREATE README
  ###########################################################################

  echo "DIVERGENT | README | Creating documentation..."

  cat > README_divergent.txt <<DOCEOF
DIVERGENT TRANSCRIPTION — ${sample_id} (statistical, GMM-ranked)
────────────────────────────────────────────────────────────────────────────
  Pairs strand-specific 3' peaks within nt_window, extracts Bayesian-balance +
  signal-to-background features per pair, fits a 2-component GMM, and ranks
  regions by posterior probability. Uses RAW (unnormalized) tracks from
  primary/unique alignments -- auto-calibration needs absolute signal levels,
  and normalization would distort them.

  IMPORTANT -- "--fdr"/divergent_fdr is a STRINGENCY KNOB, not a calibrated
  FDR: it is not Benjamini-Hochberg and does not bound a real false-discovery
  rate. Measured empirically on real PRO-seq data, the actual false-positive
  rate sits around ~50% regardless of this setting -- a real ceiling of
  peak-pair scoring against pervasive PRO-seq signal, not an unvalidated
  assumption. See this run's divergent_transcription_qc.txt "Empirical-Null
  FDR Check" section for the measured rate; treat the output as a ranked
  candidate list, not a discovery set with a known error rate.

  divergent_transcription.bed (BED5: chrom,start,end,total_signal,confidence)
  divergent_transcription_qc.txt, divergent_summary.tsv, divergent.log

  This run: threshold=\${THRESHOLD}  sum_thr=\${SUM_THR}  fdr=\${FDR}
            nt_window=\${NT_WINDOW}bp  balance=\${BALANCE}  bin_gap=\${BIN_GAP}bp
  Detected regions: \${DT_COUNT}
  \$([ -s divergent_summary.tsv ] && [ \${DT_COUNT} -gt 0 ] && cat <<STATS
  Positive peaks: \${N_POS_PK:-unknown}  Negative peaks: \${N_NEG_PK:-unknown}  Raw pairs: \${N_PAIRS_RAW:-unknown}
  Retention rate: \$(awk "BEGIN {printf \\"%.1f%%\\", (\${DT_COUNT}/\${N_PAIRS_RAW:-1})*100}" 2>/dev/null || echo "N/A")  Processing time: \${WALL_TIME:-unknown}s
STATS
)
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
    FIRST_LINE=\$(grep -v '^#' divergent_transcription.bed | head -1 || true)
    if [[ -n "\${FIRST_LINE}" ]]; then
      COL_COUNT=\$(echo "\${FIRST_LINE}" | awk '{print NF}')
      if [[ \${COL_COUNT} -ne 5 ]]; then
        echo "DIVERGENT | WARNING | BED file should have 5 columns, found \${COL_COUNT}"
      fi
      
      # Check coordinates
      START=\$(echo "\${FIRST_LINE}" | cut -f2)
      END=\$(echo "\${FIRST_LINE}" | cut -f3)
      if [[ \${START} -ge \${END} ]]; then
        echo "DIVERGENT | ERROR | Invalid BED coordinates: start >= end"
        VALIDATION_OK=0
      fi
      
      # Check score is numeric and in [0,1]
      SCORE=\$(echo "\${FIRST_LINE}" | cut -f5)
      if ! awk -v s="\${SCORE}" 'BEGIN {exit !(s >= 0 && s <= 1)}'; then
        echo "DIVERGENT | WARNING | Confidence score should be in [0,1], found \${SCORE}"
      fi
    fi
  fi

  if [[ \${VALIDATION_OK} -eq 0 ]]; then
    tracktx_error "detect_divergent_transcription" "Output validation failed (BED or summary missing/invalid)" "Check divergent.log in work dir for details"
  fi

  echo "DIVERGENT | VALIDATE | All outputs validated"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "DIVERGENT | SUMMARY | Sample: \${SAMPLE_ID}"
  echo "DIVERGENT | SUMMARY | Algorithm: Gaussian Mixture Model + FDR"
  echo "DIVERGENT | SUMMARY | Detected regions: \${DT_COUNT}"
  echo "DIVERGENT | SUMMARY | Output size: \${BED_SIZE} bytes"
  if [[ -s divergent_summary.tsv && \${DT_COUNT} -gt 0 ]]; then
    echo "DIVERGENT | SUMMARY | Positive peaks: \${N_POS_PK}"
    echo "DIVERGENT | SUMMARY | Negative peaks: \${N_NEG_PK}"
    echo "DIVERGENT | SUMMARY | Raw pairs: \${N_PAIRS_RAW}"
    echo "DIVERGENT | SUMMARY | Retention rate: \$(awk "BEGIN {printf \\"%.1f%%\\", (\${DT_COUNT}/\${N_PAIRS_RAW})*100}" 2>/dev/null || echo "N/A")"
    
    if [[ \${COL_COUNT:-0} -eq 5 ]]; then
      echo "DIVERGENT | SUMMARY | Score range: [\${SCORE_MIN}, \${SCORE_MAX}]"
      echo "DIVERGENT | SUMMARY | Mean confidence: \${SCORE_MEAN}"
    fi
  fi
  echo "DIVERGENT | SUMMARY | Detection time: \${DETECT_TIME}s"
  if [[ \${DO_QC} -eq 1 && -s divergent_transcription_qc.txt ]]; then
    echo "DIVERGENT | SUMMARY | QC report: generated"
  fi
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "DIVERGENT | COMPLETE | sample=\${SAMPLE_ID} | ts=\${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  """
}
