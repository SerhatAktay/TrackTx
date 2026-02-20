// ============================================================================
// call_functional_regions.nf — Functional Region Assignment and Quantification
// ============================================================================
//
// Purpose:
//   Assigns 3' signal to functional genomic regions hierarchically
//
// Features:
//   • Signal-based region assignment (not just read counting)
//   • Hierarchical masking ensures mutually exclusive regions
//   • Active genes only (promoter intersects divergent transcription)
//   • Sequential masking prevents double-counting
//   • Comprehensive region definitions (promoter, gene body, CPS, etc.)
//
// Region Hierarchy (sequential masking):
//   1. Promoter (TSS ± window)
//   2. Gene Body (promoter to CPS)
//   3. CPS (Cleavage/Polyadenylation Site, TES ± 500bp)
//   4. Enhancers (from divergent transcription)
//   5. Termination Window (downstream of CPS)
//   6. Non-localized (remaining signal)
//
// Note: v10.0 removed "Short Genes" category - all genes processed uniformly
//
// Active Gene Definition:
//   Gene is "active" if promoter intersects divergent transcription regions
//
// Why RAW (Unnormalized) Tracks?
//   • Matches original bash script logic
//   • Absolute signal values for consistent thresholds
//   • Signal-based counting via bedtools map -o sum
//
// Inputs:
//   tuple(sample_id, divergent_bed,
//         pos3_cpm_bg, neg3_cpm_bg,  # Note: actually RAW tracks despite names
//         pos3_sicpm_bg, neg3_sicpm_bg,  # Placeholders
//         condition, timepoint, replicate)
//   path(gtf_file)              # Interface stability (unused)
//   path(functional_regions_py) # Python driver script
//   path(genes_tsv)             # Gene annotations
//   path(tss_bed)               # TSS positions (optional override)
//   path(tes_bed)               # TES positions (optional override)
//
// Outputs:
//   ${params.output_dir}/07_functional_regions/${sample_id}/
//     ├── functional_regions.bed              — BED9 with per-region signal
//     ├── functional_regions_summary.tsv      — Category totals
//     ├── README_functional_regions.txt       — Documentation
//     └── functional_regions.log              — Processing log
//
// Parameters (params.functional_regions.*):
//   prom_up             : Promoter upstream (default: 250 bp)
//   prom_down           : Promoter downstream (default: 250 bp)
//   div_inner           : Divergent inner bound (default: 250 bp)
//   div_outer           : Divergent outer bound (default: 750 bp)
//   tw_length           : Termination window length (default: 10,000 bp)
//   min_signal          : Minimum signal threshold (default: 0.0)
//   allow_unstranded    : Allow unstranded genes (default: true)
//   count_mode          : "signal" or "event" (default: "signal")
//   tss_active_pm       : TSS active window ±bp (default: 600 bp)
//   div_fallback_enable : Enable fallback for active genes (default: false)
//
// ============================================================================

nextflow.enable.dsl = 2

process call_functional_regions {

  tag        { sample_id }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/07_functional_regions/${sample_id}",
             mode: 'copy',
             overwrite: true

  conda (params.conda_fgr ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(divergent_bed),
          path(pos3_cpm_bg), path(neg3_cpm_bg),
          path(pos3_sicpm_bg), path(neg3_sicpm_bg),
          val(condition), val(timepoint), val(replicate)
    path gtf_file
    path functional_regions_py
    path genes_tsv
    path tss_bed
    path tes_bed

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path("functional_regions.bed"),
          path("functional_regions_summary.tsv"),
          val(condition), val(timepoint), val(replicate),
          emit: main

    path "README_functional_regions.txt", emit: readme
    path "functional_regions.log",        emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Redirect all output to log file
  exec > >(tee -a functional_regions.log) 2>&1

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "FUNCREGION | START | sample=!{sample_id} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="!{sample_id}"
  CONDITION="!{condition}"
  TIMEPOINT="!{timepoint}"
  REPLICATE="!{replicate}"
  THREADS=!{task.cpus}

  # Input files
  DIV_BED="!{divergent_bed}"
  POS_BG="!{pos3_cpm_bg}"    # Note: Actually RAW tracks despite name
  NEG_BG="!{neg3_cpm_bg}"    # Note: Actually RAW tracks despite name
  GENES_TSV="!{genes_tsv}"
  TSS_BED="!{tss_bed}"
  TES_BED="!{tes_bed}"
  FGR_SCRIPT="!{functional_regions_py}"

  # Region geometry parameters
  PROM_UP=!{params.functional_regions?.prom_up ?: 250}
  PROM_DOWN=!{params.functional_regions?.prom_down ?: 250}
  DIV_INNER=!{params.functional_regions?.div_inner ?: 250}
  DIV_OUTER=!{params.functional_regions?.div_outer ?: 750}
  TW_LENGTH=!{params.functional_regions?.tw_length ?: 10000}
  TSS_ACTIVE_PM=!{params.functional_regions?.tss_active_pm ?: 600}
  ACTIVE_SLOP=!{params.functional_regions?.active_slop ?: 0}

  # Signal parameters
  MIN_SIGNAL=!{params.functional_regions?.min_signal ?: 0.0}
  MIN_SIGNAL_MODE="!{params.functional_regions?.min_signal_mode ?: 'absolute'}"
  MIN_SIGNAL_QUANTILE=!{params.functional_regions?.min_signal_quantile ?: 0.90}
  COUNT_MODE="!{params.functional_regions?.count_mode ?: 'signal'}"

  # Divergent fallback parameters
  DIV_FALLBACK_ENABLE=$([[ "!{params.functional_regions?.div_fallback_enable}" == "true" ]] && echo 1 || echo 0)
  DIV_FALLBACK_THRESHOLD=!{params.functional_regions?.div_fallback_threshold ?: 0.30}
  DIV_FALLBACK_MAX_FRAC=!{params.functional_regions?.div_fallback_max_frac ?: 0.25}

  # Feature flags
  ALLOW_UNSTRANDED=$([[ "!{params.functional_regions?.allow_unstranded}" == "false" ]] && echo 0 || echo 1)

  echo "FUNCREGION | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "FUNCREGION | CONFIG | Condition: ${CONDITION}"
  echo "FUNCREGION | CONFIG | Timepoint: ${TIMEPOINT}"
  echo "FUNCREGION | CONFIG | Replicate: ${REPLICATE}"
  echo "FUNCREGION | CONFIG | Threads: ${THREADS}"
  echo ""
  echo "FUNCREGION | CONFIG | Input Files:"
  echo "FUNCREGION | CONFIG |   Divergent BED: $(basename ${DIV_BED})"
  echo "FUNCREGION | CONFIG |   Positive track: $(basename ${POS_BG}) [RAW]"
  echo "FUNCREGION | CONFIG |   Negative track: $(basename ${NEG_BG}) [RAW]"
  echo "FUNCREGION | CONFIG |   Genes: $(basename ${GENES_TSV})"
  echo "FUNCREGION | CONFIG |   TSS: $(basename ${TSS_BED})"
  echo "FUNCREGION | CONFIG |   TES: $(basename ${TES_BED})"
  echo ""
  echo "FUNCREGION | CONFIG | Region Geometry:"
  echo "FUNCREGION | CONFIG |   Promoter: TSS -${PROM_UP} to +${PROM_DOWN} bp"
  echo "FUNCREGION | CONFIG |   Divergent: TSS -${DIV_OUTER} to -${DIV_INNER} bp (opposite strand)"
  echo "FUNCREGION | CONFIG |   CPS: TES ±500 bp (fixed)"
  echo "FUNCREGION | CONFIG |   Termination window: ${TW_LENGTH} bp downstream of CPS"
  echo "FUNCREGION | CONFIG |   TSS active window: ±${TSS_ACTIVE_PM} bp"
  echo "FUNCREGION | CONFIG |   Active slop: ${ACTIVE_SLOP} bp"
  echo ""
  echo "FUNCREGION | CONFIG | Signal Parameters:"
  echo "FUNCREGION | CONFIG |   Count mode: ${COUNT_MODE}"
  echo "FUNCREGION | CONFIG |   Min signal: ${MIN_SIGNAL} (${MIN_SIGNAL_MODE})"
  if [[ "${MIN_SIGNAL_MODE}" == "quantile" ]]; then
    echo "FUNCREGION | CONFIG |   Min signal quantile: ${MIN_SIGNAL_QUANTILE}"
  fi
  echo ""
  echo "FUNCREGION | CONFIG | Feature Flags:"
  echo "FUNCREGION | CONFIG |   Allow unstranded: $([ ${ALLOW_UNSTRANDED} -eq 1 ] && echo "yes" || echo "no")"
  echo "FUNCREGION | CONFIG |   Divergent fallback: $([ ${DIV_FALLBACK_ENABLE} -eq 1 ] && echo "enabled" || echo "disabled")"
  if [[ ${DIV_FALLBACK_ENABLE} -eq 1 ]]; then
    echo "FUNCREGION | CONFIG |   Fallback threshold: ${DIV_FALLBACK_THRESHOLD}"
    echo "FUNCREGION | CONFIG |   Fallback max fraction: ${DIV_FALLBACK_MAX_FRAC}"
  fi

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "FUNCREGION | VALIDATE | Checking input files..."

  # Use micromamba run to ensure correct Python env when in container (Docker/Singularity)
  if command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  elif [[ -x /opt/conda/bin/python3 ]]; then
    PYTHON_CMD="/opt/conda/bin/python3"
  else
    PYTHON_CMD="python3"
  fi

  VALIDATION_OK=1

  # Check Python script
  if [[ ! -f "${FGR_SCRIPT}" ]]; then
    echo "FUNCREGION | ERROR | Python script not found: ${FGR_SCRIPT}"
    VALIDATION_OK=0
  else
    echo "FUNCREGION | VALIDATE | Python script: ${FGR_SCRIPT}"
  fi

  # Check required input files
  for FILE in "${DIV_BED}" "${POS_BG}" "${NEG_BG}" "${GENES_TSV}" "${TSS_BED}" "${TES_BED}"; do
    if [[ ! -s "${FILE}" ]]; then
      echo "FUNCREGION | ERROR | Missing or empty input: ${FILE}"
      VALIDATION_OK=0
    else
      FILE_SIZE=$(stat -c%s "${FILE}" 2>/dev/null || stat -f%z "${FILE}" 2>/dev/null || echo "unknown")
      echo "FUNCREGION | VALIDATE | $(basename ${FILE}): ${FILE_SIZE} bytes"
    fi
  done

  # Count input features
  DIV_COUNT=$(grep -v '^#' "${DIV_BED}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  GENE_COUNT=$(tail -n +2 "${GENES_TSV}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  TSS_COUNT=$(wc -l < "${TSS_BED}" 2>/dev/null | tr -d ' ' || echo 0)
  TES_COUNT=$(wc -l < "${TES_BED}" 2>/dev/null | tr -d ' ' || echo 0)

  echo "FUNCREGION | VALIDATE | Divergent regions: ${DIV_COUNT}"
  echo "FUNCREGION | VALIDATE | Genes: ${GENE_COUNT}"
  echo "FUNCREGION | VALIDATE | TSS sites: ${TSS_COUNT}"
  echo "FUNCREGION | VALIDATE | TES sites: ${TES_COUNT}"

  # Validate tools
  if ${PYTHON_CMD} --version >/dev/null 2>&1; then
    PYTHON_VERSION=$(${PYTHON_CMD} --version 2>&1 || echo "unknown")
    echo "FUNCREGION | VALIDATE | Python: ${PYTHON_VERSION}"
  else
    echo "FUNCREGION | ERROR | Python not found (tried: ${PYTHON_CMD})"
    VALIDATION_OK=0
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "FUNCREGION | ERROR | Validation failed"
    exit 1
  fi

  ###########################################################################
  # 3) RUN FUNCTIONAL REGION CALLER
  ###########################################################################

  echo "FUNCREGION | CALL | Running functional region caller..."
  echo "FUNCREGION | CALL | This may take several minutes..."

  # Build command with all parameters
  CALL_START=$(date +%s)
  
  EXTRA_ARGS=""
  [[ ${ALLOW_UNSTRANDED} -eq 1 ]] && EXTRA_ARGS="${EXTRA_ARGS} --allow-unstranded"
  [[ ${DIV_FALLBACK_ENABLE} -eq 1 ]] && EXTRA_ARGS="${EXTRA_ARGS} --div-fallback-enable"

  set +e
  ${PYTHON_CMD} "${FGR_SCRIPT}" \
    --sid "${SAMPLE_ID}" \
    --genes "${GENES_TSV}" \
    --divergent "${DIV_BED}" \
    --pos "${POS_BG}" \
    --neg "${NEG_BG}" \
    --tss "${TSS_BED}" \
    --tes "${TES_BED}" \
    --prom-up "${PROM_UP}" \
    --prom-down "${PROM_DOWN}" \
    --div-inner "${DIV_INNER}" \
    --div-outer "${DIV_OUTER}" \
    --tss-active-pm "${TSS_ACTIVE_PM}" \
    --tw-length "${TW_LENGTH}" \
    --min-signal "${MIN_SIGNAL}" \
    --min-signal-mode "${MIN_SIGNAL_MODE}" \
    --min-signal-quantile "${MIN_SIGNAL_QUANTILE}" \
    --div-fallback-threshold "${DIV_FALLBACK_THRESHOLD}" \
    --div-fallback-max-frac "${DIV_FALLBACK_MAX_FRAC}" \
    --active-slop "${ACTIVE_SLOP}" \
    --count-mode "${COUNT_MODE}" \
    --outdir "." \
    ${EXTRA_ARGS}
  
  CALL_RC=$?
  set -e
  
  CALL_END=$(date +%s)
  CALL_TIME=$((CALL_END - CALL_START))

  echo "FUNCREGION | CALL | Processing completed in ${CALL_TIME}s"

  # Handle failures gracefully
  if [[ ${CALL_RC} -ne 0 ]]; then
    echo "FUNCREGION | ERROR | Functional region caller failed with exit code ${CALL_RC}"
    echo "FUNCREGION | ERROR | Creating empty output files"
    
    : > functional_regions.bed
    
    cat > functional_regions_summary.tsv <<SUMMARYEOF
region  signal  region_count
promoter  0.0 0
gene_body 0.0 0
cps 0.0 0
enhancer  0.0 0
termination_window  0.0 0
non_localized 0.0 0
SUMMARYEOF
    
    exit 1
  fi

  # Ensure output files exist
  [[ -s functional_regions.bed ]] || : > functional_regions.bed
  
  if [[ ! -s functional_regions_summary.tsv ]]; then
    cat > functional_regions_summary.tsv <<SUMMARYEOF
region  signal  region_count
promoter  0.0 0
gene_body 0.0 0
cps 0.0 0
enhancer  0.0 0
termination_window  0.0 0
non_localized 0.0 0
SUMMARYEOF
  fi

  ###########################################################################
  # 4) PARSE RESULTS
  ###########################################################################

  echo "FUNCREGION | RESULTS | Parsing functional region results..."

  # Count regions
  REGION_COUNT=$(grep -v '^#' functional_regions.bed 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  BED_SIZE=$(stat -c%s functional_regions.bed 2>/dev/null || stat -f%z functional_regions.bed 2>/dev/null || echo "unknown")

  echo "FUNCREGION | RESULTS | Total regions: ${REGION_COUNT}"
  echo "FUNCREGION | RESULTS | Output BED size: ${BED_SIZE} bytes"

  # Parse summary statistics
  if [[ -s functional_regions_summary.tsv ]]; then
    echo "FUNCREGION | RESULTS | Signal distribution by category:"
    
    # Skip header and display each category
    tail -n +2 functional_regions_summary.tsv | while IFS=$'\t' read -r CATEGORY SIGNAL COUNT; do
      echo "FUNCREGION | RESULTS |   ${CATEGORY}: signal=${SIGNAL}, regions=${COUNT}"
    done
    
    # Calculate total signal
    TOTAL_SIGNAL=$(tail -n +2 functional_regions_summary.tsv | awk -F'\t' '{sum+=$2} END{print sum}')
    TOTAL_REGIONS=$(tail -n +2 functional_regions_summary.tsv | awk -F'\t' '{sum+=$3} END{print sum}')
    
    echo "FUNCREGION | RESULTS | Total signal: ${TOTAL_SIGNAL}"
    echo "FUNCREGION | RESULTS | Total regions: ${TOTAL_REGIONS}"
  fi

  ###########################################################################
  # 5) CREATE README
  ###########################################################################

  echo "FUNCREGION | README | Creating documentation..."

  cat > README_functional_regions.txt <<'DOCEOF'
================================================================================
FUNCTIONAL REGIONS — !{sample_id}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Assignment of 3' signal to functional genomic regions using hierarchical
  masking to ensure mutually exclusive categories.
  
  Each signal unit is assigned to exactly ONE functional category based on
  sequential priority: Promoter → Gene Body → CPS → Enhancers
  → Termination Window → Non-localized

METHOD
────────────────────────────────────────────────────────────────────────────
  1. Define active genes (promoter intersects divergent transcription)
  2. Define functional regions for active genes
  3. Sequentially mask and assign signal:
     a. Promoter regions
     b. Gene body regions
     c. CPS (Cleavage/Polyadenylation Site) regions
     d. Enhancer regions (from divergent transcription)
     e. Termination windows
     g. Remaining signal → Non-localized
  4. Quantify signal per region using bedtools map
  5. Summarize by category

ACTIVE GENE DEFINITION
────────────────────────────────────────────────────────────────────────────
  A gene is considered "active" if its promoter region (TSS ±${TSS_ACTIVE_PM} bp)
  intersects at least one divergent transcription region.
  
  Divergent fallback: $([ ${DIV_FALLBACK_ENABLE} -eq 1 ] && echo "Enabled" || echo "Disabled")
  $([ ${DIV_FALLBACK_ENABLE} -eq 1 ] && cat <<FALLBACK
  When enabled, genes with promoter signal ≥ ${DIV_FALLBACK_THRESHOLD} × max
  and representing ≤ ${DIV_FALLBACK_MAX_FRAC} of genes are marked active
  even without divergent transcription overlap.
FALLBACK
)

INPUT DATA
────────────────────────────────────────────────────────────────────────────
  Signal Tracks: RAW (unnormalized) 3' bedGraphs
    • Positive strand: $(basename ${POS_BG})
    • Negative strand: $(basename ${NEG_BG})
  
  Why RAW tracks?
    • Matches original bash script logic
    • Absolute signal values for consistent thresholds
    • Signal-based counting via bedtools map -o sum
  
  Annotations:
    • Divergent transcription: ${DIV_COUNT} regions
    • Genes: ${GENE_COUNT} annotations
    • TSS sites: ${TSS_COUNT}
    • TES sites: ${TES_COUNT}

SAMPLE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Sample:     !{sample_id}
  Condition:  !{condition}
  Timepoint:  !{timepoint}
  Replicate:  !{replicate}

REGION DEFINITIONS
────────────────────────────────────────────────────────────────────────────

1. Promoter:
   Location: TSS -${PROM_UP} to +${PROM_DOWN} bp
   Strand: Gene strand
   Purpose: Transcription initiation region

2. Divergent:
   Location: TSS -${DIV_OUTER} to -${DIV_INNER} bp
   Strand: OPPOSITE of gene strand
   Purpose: Divergent transcription detection
   Note: Used for active gene definition, not signal assignment

3. Gene Body:
   Location: End of promoter to start of CPS
   Strand: Gene strand
   Purpose: Elongation and productive transcription

4. CPS (Cleavage/Polyadenylation Site):
   Location: TES -500 to +500 bp (fixed)
   Strand: Gene strand
   Purpose: Transcription termination signal

5. Enhancers:
   Source: Divergent transcription regions not overlapping promoters
   Purpose: Putative enhancer elements

7. Termination Window:
   Location: ${TW_LENGTH} bp downstream of CPS end
   Strand: Gene strand
   Purpose: Readthrough and termination region

8. Non-localized:
   Definition: Signal not assigned to any category above
   Purpose: Intergenic, intronic, or unannotated transcription

SIGNAL QUANTIFICATION
────────────────────────────────────────────────────────────────────────────
  Count Mode: ${COUNT_MODE}
  
  signal mode (default):
    • Sum of all signal values overlapping region
    • Accounts for signal strength
    • bedtools map -o sum
  
  event mode (legacy):
    • Count of discrete 1bp events in region
    • Binary presence/absence
    • bedtools map -o count

  Minimum Signal: ${MIN_SIGNAL} (${MIN_SIGNAL_MODE})
  $([ "${MIN_SIGNAL_MODE}" == "quantile" ] && echo "  Quantile threshold: ${MIN_SIGNAL_QUANTILE}")

HIERARCHICAL MASKING
────────────────────────────────────────────────────────────────────────────
  Signal is assigned sequentially with masking to prevent double-counting:
  
  Step 1: Assign promoter signal → mask promoters
  Step 2: Assign gene body signal → mask gene bodies
  Step 3: Assign CPS signal → mask CPS regions
  Step 4: Assign enhancer signal → mask enhancers
  Step 5: Assign termination window signal → mask TW
  Step 7: Assign remaining signal → non-localized
  
  Result: Each signal unit counted exactly once

FILES
────────────────────────────────────────────────────────────────────────────
  functional_regions.bed          — BED9 format with per-region assignments
  functional_regions_summary.tsv  — Category totals (region, signal, count)
  README_functional_regions.txt   — This documentation
  functional_regions.log          — Processing log

OUTPUT FORMAT
────────────────────────────────────────────────────────────────────────────

BED9 Format (functional_regions.bed):
  Column 1: chromosome
  Column 2: start
  Column 3: end
  Column 4: region_name (e.g., "promoter", "gene_body")
  Column 5: signal (quantified value)
  Column 6: strand
  Column 7-9: BED9 fields (RGB color, thickStart, thickEnd)

Summary TSV (functional_regions_summary.tsv):
  Column 1: region      — Category name
  Column 2: signal      — Total signal in category
  Column 3: region_count— Number of regions in category

RESULTS SUMMARY
────────────────────────────────────────────────────────────────────────────
  Total regions: ${REGION_COUNT}
  Processing time: ${CALL_TIME}s
  
  Signal distribution:
$([ -s functional_regions_summary.tsv ] && tail -n +2 functional_regions_summary.tsv | \
  awk -F'\t' '{printf "    %-20s signal=%-12s regions=%s\n", $1":", $2, $3}' || echo "    (empty)")

QUALITY CONTROL
────────────────────────────────────────────────────────────────────────────

Expected Results:
  • Promoter signal: 10-30% of total (active transcription)
  • Gene body signal: 40-60% of total (elongation)
  • CPS signal: 5-15% of total (termination)
  • Non-localized: 10-30% of total (background)

Interpretation:
  • High promoter signal: Paused Pol-II, initiation activity
  • High gene body signal: Productive elongation
  • High non-localized: Intergenic transcription, enhancers
  • Low CPS signal: Efficient termination

Troubleshooting:
  • Zero regions: Check active gene definition, divergent transcription
  • All signal non-localized: Check gene annotations, strand matching
  • Low promoter signal: Check TSS coordinates, prom_up/prom_down
  • Unexpected distribution: Check min_signal threshold

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  These functional region assignments are used for:
  
  1. Pol-II Metrics Calculation:
     - Pausing Index (promoter/gene body ratio)
     - Termination Index (CPS/gene body ratio)
     - Region-specific coverage
  
  2. Quality Control:
     - Signal distribution patterns
     - Active gene counts
     - Coverage uniformity
  
  3. Differential Analysis:
     - Compare functional region changes
     - Region-specific differential expression
  
  4. Genome Annotation:
     - Validate gene boundaries
     - Identify novel elements

PARAMETER TUNING
────────────────────────────────────────────────────────────────────────────

More sensitive active gene detection:
  • Increase tss_active_pm (e.g., 1000)
  • Enable div_fallback_enable
  • Increase div_fallback_threshold (e.g., 0.5)

More specific active gene detection:
  • Decrease tss_active_pm (e.g., 300)
  • Disable div_fallback_enable
  • Increase min_signal

Adjust region sizes:
  • Promoter: prom_up, prom_down
  • Gene body: Automatic (promoter to CPS)
  • Termination: tw_length

PARAMETERS USED
────────────────────────────────────────────────────────────────────────────
  Promoter:             TSS -${PROM_UP} to +${PROM_DOWN} bp
  Divergent:            TSS -${DIV_OUTER} to -${DIV_INNER} bp
  CPS:                  TES ±500 bp (fixed)
  Termination window:   ${TW_LENGTH} bp
  TSS active window:    ±${TSS_ACTIVE_PM} bp
  Active slop:          ${ACTIVE_SLOP} bp
  Count mode:           ${COUNT_MODE}
  Min signal:           ${MIN_SIGNAL} (${MIN_SIGNAL_MODE})
  Allow unstranded:     $([ ${ALLOW_UNSTRANDED} -eq 1 ] && echo "yes" || echo "no")
  Divergent fallback:   $([ ${DIV_FALLBACK_ENABLE} -eq 1 ] && echo "enabled" || echo "disabled")

TECHNICAL NOTES
────────────────────────────────────────────────────────────────────────────
  • Uses bedtools for interval operations
  • Sequential masking ensures no double-counting
  • Strand-specific signal assignment
  • Active genes determined by divergent transcription overlap
  • All genes processed uniformly (v10.0: removed short gene special handling)
  • Python driver handles all region logic and masking

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Date: $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample: !{sample_id}
  Module: 10_call_functional_regions

================================================================================
DOCEOF

  echo "FUNCREGION | README | Documentation created"

  ###########################################################################
  # 6) VALIDATION
  ###########################################################################

  echo "FUNCREGION | VALIDATE | Validating outputs..."

  VALIDATION_OK=1

  # Check BED file
  if [[ ! -e functional_regions.bed ]]; then
    echo "FUNCREGION | ERROR | Output BED file missing"
    VALIDATION_OK=0
  fi

  # Check summary
  if [[ ! -s functional_regions_summary.tsv ]]; then
    echo "FUNCREGION | ERROR | Summary file missing or empty"
    VALIDATION_OK=0
  else
    # Validate summary format
    SUMMARY_LINES=$(wc -l < functional_regions_summary.tsv | tr -d ' ')
    if [[ ${SUMMARY_LINES} -lt 2 ]]; then
      echo "FUNCREGION | WARNING | Summary has fewer than expected lines"
    fi
  fi

  # Validate BED format if non-empty
  if [[ -s functional_regions.bed ]]; then
    FIRST_LINE=$(grep -v '^#' functional_regions.bed | head -1 || true)
    if [[ -n "${FIRST_LINE}" ]]; then
      COL_COUNT=$(echo "${FIRST_LINE}" | awk '{print NF}')
      if [[ ${COL_COUNT} -lt 6 ]]; then
        echo "FUNCREGION | WARNING | BED file has ${COL_COUNT} columns, expected at least 6"
      fi
    fi
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "FUNCREGION | ERROR | Validation failed"
    exit 1
  fi

  echo "FUNCREGION | VALIDATE | All outputs validated"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "FUNCREGION | SUMMARY | Sample: ${SAMPLE_ID}"
  echo "FUNCREGION | SUMMARY | Total regions: ${REGION_COUNT}"
  echo "FUNCREGION | SUMMARY | Output size: ${BED_SIZE} bytes"
  echo "FUNCREGION | SUMMARY | Processing time: ${CALL_TIME}s"
  
  if [[ -s functional_regions_summary.tsv ]]; then
    echo "FUNCREGION | SUMMARY | Signal by category:"
    tail -n +2 functional_regions_summary.tsv | while IFS=$'\t' read -r CAT SIG CNT; do
      echo "FUNCREGION | SUMMARY |   ${CAT}: ${SIG}"
    done
  fi
  
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "FUNCREGION | COMPLETE | sample=${SAMPLE_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}