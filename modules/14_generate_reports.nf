// ============================================================================
// generate_reports.nf — Per-Sample Comprehensive Report Generation
// ============================================================================
//
// Purpose:
//   Generates comprehensive per-sample reports integrating all analysis results
//
// Features:
//   • HTML report with interactive visualizations
//   • TSV summary for programmatic access
//   • JSON with structured data and schema versioning
//   • Optional plots page with embedded graphics
//   • Track file links for genome browser integration
//
// Report Contents:
//   1. Sample Metadata (condition, timepoint, replicate)
//   2. Quality Control Metrics (mapping, strand bias, coverage)
//   3. Divergent Transcription Summary (loci count, characteristics)
//   4. Functional Region Composition (promoter, body, CPS, etc.)
//   5. Pol-II Density Distribution (signal per region type)
//   6. Pausing Index Statistics (distribution, top genes)
//   7. Normalization Factors (CPM, siCPM)
//   8. Track File Links (for UCSC/IGV browsers)
//
// Track Links:
//   Supports both local file paths and HTTP(S) URLs
//   Types:
//     • Raw allMap tracks (unnormalized bedGraph)
//     • CPM normalized BigWig (main tracks)
//     • AllMap CPM normalized BigWig (including multimappers)
//
// Inputs:
//   tuple(sample_id,
//         divergent_bed, functional_summary, pol_density, pausing_index,
//         normalization_factors, dedup_stats, qc_json,
//         track_links (8 values: raw/normalized, pos/neg, main/allMap),
//         condition, timepoint, replicate)
//
// Outputs:
//   ${params.output_dir}/11_reports/samples/${sample_id}/
//     ├── ${sample_id}.report.html          — Main HTML report
//     ├── ${sample_id}.report.tsv           — Summary TSV
//     ├── ${sample_id}.report.json          — Structured JSON
//     ├── ${sample_id}.plots.html           — Plots page (optional)
//     ├── ${sample_id}.README_report.txt    — Documentation
//     └── ${sample_id}.report.log           — Processing log
//
//   Quick access symlink:
//     ${params.output_dir}/reports/${sample_id}.html → HTML report
//
// Parameters:
//   params.reports_plots : Enable plots page (default: 0/false)
//
// JSON Schema:
//   Versioned JSON output with stable schema for downstream parsing
//   Current version: 1.0
//
// ============================================================================

nextflow.enable.dsl = 2

process generate_reports {

  tag        { sample_id }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/11_reports/samples/${sample_id}",
             mode: 'copy',
             overwrite: true

  conda (params.conda_pol ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(div_bed),
          path(func_sum),
          path(pol_density),
          path(pausing_idx),
          path(norm_factors),
          path(dedup_stats),
          path(qc_json),
          val(allmap3p_pos_raw),
          val(allmap3p_neg_raw),
          val(pos3_cpm_bw),
          val(neg3_cpm_bw),
          val(allmap3p_pos_cpm_bw),
          val(allmap3p_neg_cpm_bw),
          val(condition),
          val(timepoint),
          val(replicate)

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    path("${sample_id}.report.html"),        emit: html
    path("${sample_id}.report.tsv"),         emit: tsv
    path("${sample_id}.report.json"),        emit: json
    path("${sample_id}.plots.html"),         optional: true, emit: plots
    path("${sample_id}.README_report.txt"),  emit: readme
    path("${sample_id}.report.log"),         emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  # NOTE: Using -e (not -u) because track links may be empty strings
  set -eo pipefail
  export LC_ALL=C

  # Redirect all output to log file
  exec > >(tee -a "!{sample_id}.report.log") 2>&1

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "REPORT | START | sample=!{sample_id} | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="!{sample_id}"
  CONDITION="!{condition}"
  TIMEPOINT="!{timepoint}"
  REPLICATE="!{replicate}"

  # Core input files
  DIV_BED="!{div_bed}"
  FUNC_SUM="!{func_sum}"
  POL_DENS="!{pol_density}"
  PAUSING_IDX="!{pausing_idx}"
  NORM_FACTORS="!{norm_factors}"
  DEDUP_STATS="!{dedup_stats}"
  QC_JSON="!{qc_json}"

  # Track links (may be empty strings)
  ALLMAP3P_POS_RAW="!{allmap3p_pos_raw}"
  ALLMAP3P_NEG_RAW="!{allmap3p_neg_raw}"
  POS3_CPM_BW="!{pos3_cpm_bw}"
  NEG3_CPM_BW="!{neg3_cpm_bw}"
  ALLMAP3P_POS_CPM_BW="!{allmap3p_pos_cpm_bw}"
  ALLMAP3P_NEG_CPM_BW="!{allmap3p_neg_cpm_bw}"

  # Renderer script
  RENDER_SCRIPT="!{projectDir}/bin/render_sample_report.py"

  # Parameters
  ENABLE_PLOTS=!{(params.reports_plots == null) ? 0 : (params.reports_plots as int)}

  # Output files
  OUT_HTML="${SAMPLE_ID}.report.html"
  OUT_TSV="${SAMPLE_ID}.report.tsv"
  OUT_JSON="${SAMPLE_ID}.report.json"
  OUT_PLOTS="${SAMPLE_ID}.plots.html"
  OUT_README="${SAMPLE_ID}.README_report.txt"

  echo "REPORT | CONFIG | Sample ID: ${SAMPLE_ID}"
  echo "REPORT | CONFIG | Condition: ${CONDITION}"
  echo "REPORT | CONFIG | Timepoint: ${TIMEPOINT}"
  echo "REPORT | CONFIG | Replicate: ${REPLICATE}"
  echo ""
  echo "REPORT | CONFIG | Renderer script: ${RENDER_SCRIPT}"
  echo "REPORT | CONFIG | Enable plots: $([ ${ENABLE_PLOTS} -eq 1 ] && echo "yes" || echo "no")"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "REPORT | VALIDATE | Checking required input files..."

  # Use micromamba run to ensure correct Python env when in container (Docker/Singularity)
  if command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  else
    PYTHON_CMD="python3"
  fi

  VALIDATION_OK=1

  # Check renderer script
  if [[ ! -e "${RENDER_SCRIPT}" ]]; then
    echo "REPORT | ERROR | Renderer script not found: ${RENDER_SCRIPT}"
    VALIDATION_OK=0
  else
    echo "REPORT | VALIDATE | Renderer script: ${RENDER_SCRIPT}"
  fi

  # Check core input files
  validate_file() {
    local label="$1"
    local file="$2"
    
    if [[ ! -s "${file}" ]]; then
      echo "REPORT | ERROR | ${label} missing or empty: ${file}"
      VALIDATION_OK=0
      return 1
    else
      FILE_SIZE=$(stat -c%s "${file}" 2>/dev/null || stat -f%z "${file}" 2>/dev/null || echo "unknown")
      FILE_LINES=$(wc -l < "${file}" 2>/dev/null | tr -d ' ' || echo 0)
      echo "REPORT | VALIDATE | ${label}: ${FILE_SIZE} bytes (${FILE_LINES} lines)"
      return 0
    fi
  }

  validate_file "Divergent bed" "${DIV_BED}"
  validate_file "Functional summary" "${FUNC_SUM}"
  validate_file "Pol-II density" "${POL_DENS}"
  validate_file "Pausing index" "${PAUSING_IDX}"
  validate_file "Normalization factors" "${NORM_FACTORS}"
  validate_file "QC JSON" "${QC_JSON}"

  # Dedup stats is optional
  if [[ -s "${DEDUP_STATS}" ]]; then
    DEDUP_SIZE=$(stat -c%s "${DEDUP_STATS}" 2>/dev/null || stat -f%z "${DEDUP_STATS}" 2>/dev/null || echo "unknown")
    echo "REPORT | VALIDATE | Dedup stats: ${DEDUP_SIZE} bytes"
  else
    echo "REPORT | VALIDATE | Dedup stats: not available"
  fi

  # Validate tools
  if ${PYTHON_CMD} --version >/dev/null 2>&1; then
    PYTHON_VERSION=$(${PYTHON_CMD} --version 2>&1 || echo "unknown")
    echo "REPORT | VALIDATE | Python: ${PYTHON_VERSION}"
  else
    echo "REPORT | ERROR | Python not found (tried: ${PYTHON_CMD})"
    VALIDATION_OK=0
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "REPORT | ERROR | Validation failed"
    exit 1
  fi

  ###########################################################################
  # 3) VALIDATE AND CATALOG TRACK LINKS
  ###########################################################################

  echo "REPORT | TRACKS | Validating track file links..."

  # Helper to validate track link (URL or file path)
  validate_track() {
    local label="$1"
    local link="$2"
    
    # Empty is OK (optional)
    if [[ -z "${link}" ]]; then
      echo "REPORT | TRACKS | ${label}: not provided"
      return 1
    fi
    
    # Check if it's a URL
    if [[ "${link}" =~ ^https?:// ]]; then
      echo "REPORT | TRACKS | ${label}: URL (${link})"
      return 0
    fi
    
    # Check if file exists
    if [[ -e "${link}" ]]; then
      LINK_SIZE=$(stat -c%s "${link}" 2>/dev/null || stat -f%z "${link}" 2>/dev/null || echo "unknown")
      echo "REPORT | TRACKS | ${label}: file (${LINK_SIZE} bytes)"
      return 0
    fi
    
    # Not found
    echo "REPORT | TRACKS | ${label}: not found (${link})"
    return 1
  }

  # Track availability flags
  HAVE_ALLMAP_POS_RAW=0
  HAVE_ALLMAP_NEG_RAW=0
  HAVE_POS_CPM_BW=0
  HAVE_NEG_CPM_BW=0
  HAVE_ALLMAP_POS_CPM_BW=0
  HAVE_ALLMAP_NEG_CPM_BW=0

  validate_track "AllMap 3p pos raw" "${ALLMAP3P_POS_RAW}" && HAVE_ALLMAP_POS_RAW=1
  validate_track "AllMap 3p neg raw" "${ALLMAP3P_NEG_RAW}" && HAVE_ALLMAP_NEG_RAW=1
  validate_track "3p pos CPM BigWig" "${POS3_CPM_BW}" && HAVE_POS_CPM_BW=1
  validate_track "3p neg CPM BigWig" "${NEG3_CPM_BW}" && HAVE_NEG_CPM_BW=1
  validate_track "AllMap 3p pos CPM BigWig" "${ALLMAP3P_POS_CPM_BW}" && HAVE_ALLMAP_POS_CPM_BW=1
  validate_track "AllMap 3p neg CPM BigWig" "${ALLMAP3P_NEG_CPM_BW}" && HAVE_ALLMAP_NEG_CPM_BW=1

  TOTAL_TRACKS=$((HAVE_ALLMAP_POS_RAW + HAVE_ALLMAP_NEG_RAW + \
                  HAVE_POS_CPM_BW + HAVE_NEG_CPM_BW + \
                  HAVE_ALLMAP_POS_CPM_BW + HAVE_ALLMAP_NEG_CPM_BW))
  
  echo "REPORT | TRACKS | Available tracks: ${TOTAL_TRACKS}/6"

  ###########################################################################
  # 4) PARSE INPUT DATA FOR SUMMARY
  ###########################################################################

  echo "REPORT | PARSE | Extracting summary statistics..."

  # Count divergent loci
  DIV_COUNT=$(grep -v '^#' "${DIV_BED}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  echo "REPORT | PARSE | Divergent loci: ${DIV_COUNT}"

  # Parse functional regions
  if [[ -s "${FUNC_SUM}" ]]; then
    FUNC_LINES=$(tail -n +2 "${FUNC_SUM}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
    echo "REPORT | PARSE | Functional regions: ${FUNC_LINES} categories"
  fi

  # Parse pausing index
  if [[ -s "${PAUSING_IDX}" ]]; then
    PAUSING_GENES=$(tail -n +2 "${PAUSING_IDX}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
    echo "REPORT | PARSE | Genes with pausing index: ${PAUSING_GENES}"
  fi

  # Parse QC JSON for key metrics
  if command -v jq >/dev/null 2>&1 && [[ -s "${QC_JSON}" ]]; then
    TOTAL_READS=$(jq -r '.total_reads_raw // 0' "${QC_JSON}" 2>/dev/null || echo 0)
    MAP_RATE=$(jq -r '.map_rate_percent // 0' "${QC_JSON}" 2>/dev/null || echo 0)
    DUP_RATE=$(jq -r '.duplicate_perc_of_total // 0' "${QC_JSON}" 2>/dev/null || echo 0)
    echo "REPORT | PARSE | Total reads: ${TOTAL_READS}"
    echo "REPORT | PARSE | Mapping rate: ${MAP_RATE}%"
    echo "REPORT | PARSE | Duplicate rate: ${DUP_RATE}%"
  else
    echo "REPORT | PARSE | jq not available, skipping QC metrics"
  fi

  ###########################################################################
  # 5) BUILD RENDERER COMMAND
  ###########################################################################

  echo "REPORT | BUILD | Building renderer command..."

  # Base arguments
  RENDERER_ARGS=(
    --sample "${SAMPLE_ID}"
    --condition "${CONDITION}"
    --timepoint "${TIMEPOINT}"
    --replicate "${REPLICATE}"
    --divergent-bed "${DIV_BED}"
    --functional-summary "${FUNC_SUM}"
    --pol-density "${POL_DENS}"
    --pausing-index "${PAUSING_IDX}"
    --norm-factors "${NORM_FACTORS}"
    --qc-json "${QC_JSON}"
    --out-html "${OUT_HTML}"
    --out-tsv "${OUT_TSV}"
    --out-json "${OUT_JSON}"
    --out-plots-html "${OUT_PLOTS}"
  )

  echo "REPORT | BUILD | Base arguments: ${#RENDERER_ARGS[@]}"

  # Helper to add track link if available
  add_track_link() {
    local flag="$1"
    local link="$2"
    local have_flag="$3"
    
    if [[ ${have_flag} -eq 1 ]]; then
      RENDERER_ARGS+=("${flag}" "${link}")
      echo "REPORT | BUILD | Added: ${flag}"
    fi
  }

  add_track_link "--allmap3p-pos-raw" "${ALLMAP3P_POS_RAW}" ${HAVE_ALLMAP_POS_RAW}
  add_track_link "--allmap3p-neg-raw" "${ALLMAP3P_NEG_RAW}" ${HAVE_ALLMAP_NEG_RAW}
  add_track_link "--pos3-cpm-bw" "${POS3_CPM_BW}" ${HAVE_POS_CPM_BW}
  add_track_link "--neg3-cpm-bw" "${NEG3_CPM_BW}" ${HAVE_NEG_CPM_BW}
  add_track_link "--allmap3p-pos-cpm-bw" "${ALLMAP3P_POS_CPM_BW}" ${HAVE_ALLMAP_POS_CPM_BW}
  add_track_link "--allmap3p-neg-cpm-bw" "${ALLMAP3P_NEG_CPM_BW}" ${HAVE_ALLMAP_NEG_CPM_BW}

  # Add plots flag if enabled
  if [[ ${ENABLE_PLOTS} -eq 1 ]]; then
    RENDERER_ARGS+=(--plots 1)
    echo "REPORT | BUILD | Plots enabled"
  fi

  echo "REPORT | BUILD | Total arguments: ${#RENDERER_ARGS[@]}"

  ###########################################################################
  # 6) RUN RENDERER
  ###########################################################################

  echo "REPORT | RENDER | Generating reports..."

  RENDER_START=$(date +%s)

  set +e
  ${PYTHON_CMD} "${RENDER_SCRIPT}" "${RENDERER_ARGS[@]}"
  RENDER_RC=$?
  set -e

  RENDER_END=$(date +%s)
  RENDER_TIME=$((RENDER_END - RENDER_START))

  echo "REPORT | RENDER | Rendering completed in ${RENDER_TIME}s"

  if [[ ${RENDER_RC} -ne 0 ]]; then
    echo "REPORT | ERROR | Renderer failed with exit code ${RENDER_RC}"
    exit 1
  fi

  ###########################################################################
  # 7) VALIDATE OUTPUTS
  ###########################################################################

  echo "REPORT | VALIDATE | Checking output files..."

  VALIDATION_OK=1

  # Check required outputs
  for OUTPUT in "${OUT_HTML}" "${OUT_TSV}" "${OUT_JSON}"; do
    if [[ ! -s "${OUTPUT}" ]]; then
      echo "REPORT | ERROR | Expected output missing or empty: ${OUTPUT}"
      VALIDATION_OK=0
    else
      OUTPUT_SIZE=$(stat -c%s "${OUTPUT}" 2>/dev/null || stat -f%z "${OUTPUT}" 2>/dev/null || echo "unknown")
      echo "REPORT | VALIDATE | $(basename ${OUTPUT}): ${OUTPUT_SIZE} bytes"
    fi
  done

  # Check plots file
  if [[ ${ENABLE_PLOTS} -eq 1 ]]; then
    if [[ -s "${OUT_PLOTS}" ]]; then
      PLOTS_SIZE=$(stat -c%s "${OUT_PLOTS}" 2>/dev/null || stat -f%z "${OUT_PLOTS}" 2>/dev/null || echo "unknown")
      echo "REPORT | VALIDATE | Plots HTML: ${PLOTS_SIZE} bytes"
    else
      echo "REPORT | WARNING | Plots enabled but file missing"
    fi
  else
    # Create placeholder if plots disabled
    if [[ ! -s "${OUT_PLOTS}" ]]; then
      cat > "${OUT_PLOTS}" <<'PLACEHOLDER'
<!DOCTYPE html>
<html>
<head>
  <meta charset="utf-8">
  <title>Plots Disabled</title>
  <style>
    body { font-family: sans-serif; padding: 2em; text-align: center; }
    .message { background: #f0f0f0; padding: 2em; border-radius: 8px; max-width: 600px; margin: 2em auto; }
  </style>
</head>
<body>
  <div class="message">
    <h1>Plots Page Disabled</h1>
    <p>To enable plots, set <code>params.reports_plots = 1</code> in your configuration.</p>
  </div>
</body>
</html>
PLACEHOLDER
      echo "REPORT | VALIDATE | Created plots placeholder"
    fi
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "REPORT | ERROR | Output validation failed"
    exit 1
  fi

  ###########################################################################
  # 8) CREATE README
  ###########################################################################

  echo "REPORT | README | Creating documentation..."

  cat > "${OUT_README}" <<'DOCEOF'
================================================================================
PER-SAMPLE REPORT — !{sample_id}
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Comprehensive report integrating all analysis results for a single sample.

SAMPLE INFORMATION
────────────────────────────────────────────────────────────────────────────
  Sample ID:    ${SAMPLE_ID}
  Condition:    ${CONDITION}
  Timepoint:    ${TIMEPOINT}
  Replicate:    ${REPLICATE}

REPORT FILES
────────────────────────────────────────────────────────────────────────────
  ${OUT_HTML}
    Interactive HTML report with:
      • Sample metadata and experimental design
      • Quality control summary
      • Divergent transcription statistics
      • Functional region composition
      • Pol-II density distribution
      • Pausing index summary
      • Normalization factors
      • Track file links for genome browsers
    
    View in web browser for best experience.
  
  ${OUT_TSV}
    Tab-separated summary table with key metrics:
      • One row per metric
      • Columns: metric_name, value, unit, category
      • Easy to parse programmatically
      • Compatible with R, Python, Excel
  
  ${OUT_JSON}
    Structured JSON with complete data:
      • Versioned schema (current: 1.0)
      • Nested structure by category
      • All numerical values and lists
      • Suitable for API integration
      • JSON Schema compliant
  
  ${OUT_PLOTS}
    $([ ${ENABLE_PLOTS} -eq 1 ] && cat <<PLOTINFO
Supplementary plots page with:
      • Functional region composition pie chart
      • Pausing index distribution histogram
      • QC metrics summary plots
      • Base64-encoded inline images
      • No external dependencies
PLOTINFO
|| echo "Placeholder (plots disabled)")
  
  ${OUT_README}
    This documentation file
  
  ${SAMPLE_ID}.report.log
    Processing log with timestamps

INPUT DATA FILES
────────────────────────────────────────────────────────────────────────────
  Core Inputs:
    • Divergent bed:          ${DIV_BED} (${DIV_COUNT} loci)
    • Functional summary:     ${FUNC_SUM}
    • Pol-II density:         ${POL_DENS}
    • Pausing index:          ${PAUSING_IDX} (${PAUSING_GENES} genes)
    • Normalization factors:  ${NORM_FACTORS}
    • QC JSON:                ${QC_JSON}
    • Dedup stats:            ${DEDUP_STATS}

TRACK FILE LINKS
────────────────────────────────────────────────────────────────────────────
  Raw Coverage (unnormalized bedGraph):
    • AllMap 3' pos:  $([ ${HAVE_ALLMAP_POS_RAW} -eq 1 ] && echo "${ALLMAP3P_POS_RAW}" || echo "Not available")
    • AllMap 3' neg:  $([ ${HAVE_ALLMAP_NEG_RAW} -eq 1 ] && echo "${ALLMAP3P_NEG_RAW}" || echo "Not available")
  
  Normalized BigWig (CPM):
    • 3' pos:         $([ ${HAVE_POS_CPM_BW} -eq 1 ] && echo "${POS3_CPM_BW}" || echo "Not available")
    • 3' neg:         $([ ${HAVE_NEG_CPM_BW} -eq 1 ] && echo "${NEG3_CPM_BW}" || echo "Not available")
  
  AllMap Normalized BigWig (CPM):
    • 3' pos:         $([ ${HAVE_ALLMAP_POS_CPM_BW} -eq 1 ] && echo "${ALLMAP3P_POS_CPM_BW}" || echo "Not available")
    • 3' neg:         $([ ${HAVE_ALLMAP_NEG_CPM_BW} -eq 1 ] && echo "${ALLMAP3P_NEG_CPM_BW}" || echo "Not available")
  
  Total tracks:       ${TOTAL_TRACKS}/6

USING THE REPORTS
────────────────────────────────────────────────────────────────────────────
  HTML Report:
    1. Open ${OUT_HTML} in web browser
    2. Navigate sections using table of contents
    3. Click track links to view in UCSC/IGV
    4. Review QC metrics for data quality
  
  TSV Summary:
    # In R
    data <- read.delim("${OUT_TSV}")
    
    # In Python
    import pandas as pd
    data = pd.read_csv("${OUT_TSV}", sep="\t")
  
  JSON Data:
    # In Python
    import json
    with open("${OUT_JSON}") as f:
        data = json.load(f)
    
    # Check schema version
    version = data.get("schema_version", "unknown")
    
    # Access metrics
    qc_metrics = data.get("qc", {})
    pausing = data.get("pausing", {})

TRACK LOADING IN GENOME BROWSERS
────────────────────────────────────────────────────────────────────────────
  UCSC Genome Browser:
    1. Open UCSC browser for your genome
    2. Click "add custom tracks"
    3. Paste track URLs from HTML report
    4. Click "submit"
  
  IGV (Integrative Genomics Viewer):
    1. Open IGV
    2. File → Load from URL
    3. Paste track URL or file path
    4. Or: File → Load from File (for local files)

REPORT CONTENTS DETAILS
────────────────────────────────────────────────────────────────────────────

1. Sample Metadata
   • Sample identification
   • Experimental conditions
   • Processing date
   • Pipeline version

2. Quality Control
   • Total reads: ${TOTAL_READS:-NA}
   • Mapping rate: ${MAP_RATE:-NA}%
   • Duplicate rate: ${DUP_RATE:-NA}%
   • Strand balance
   • Coverage depth
   • UMI deduplication (if applicable)

3. Divergent Transcription
   • Number of loci: ${DIV_COUNT}
   • Genomic distribution
   • Signal characteristics
   • Annotation overlap

4. Functional Regions
   • Promoter signal
   • Gene body signal
   • CPS (Cleavage/Polyadenylation Site) signal
   • Enhancer signal
   • Termination window signal
   • Non-localized signal

5. Pol-II Density
   • Signal per functional region type
   • Normalization method (CPM or siCPM)
   • Total signal quantification

6. Pausing Index
   • Per-gene pausing indices
   • Distribution statistics
   • Top paused genes
   • Length-normalized vs raw

7. Normalization
   • CPM factors
   • siCPM factors (if spike-in used)
   • Read counts per BAM type

QUALITY METRICS
────────────────────────────────────────────────────────────────────────────
  Check the HTML report for:
    ✓ Mapping rate >70%
    ✓ Duplicate rate <30%
    ✓ Strand balance 40-60%
    ✓ Mean depth >10×
    ✓ Expected pausing index distribution

  Red flags:
    ✗ Very low mapping rate
    ✗ Extreme strand bias
    ✗ Very high duplicate rate
    ✗ Insufficient coverage

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  1. Quality Assessment:
     - Review QC metrics
     - Identify failed samples
     - Compare across batch
  
  2. Data Export:
     - Load tracks in genome browser
     - Export TSV for further analysis
     - Parse JSON for custom scripts
  
  3. Publication:
     - Include QC metrics in methods
     - Reference track URLs in data availability
     - Share JSON for reproducibility

PROCESSING DETAILS
────────────────────────────────────────────────────────────────────────────
  Renderer:         ${RENDER_SCRIPT}
  Python version:   ${PYTHON_VERSION}
  Processing time:  ${RENDER_TIME}s
  Plots enabled:    $([ ${ENABLE_PLOTS} -eq 1 ] && echo "Yes" || echo "No")

TROUBLESHOOTING
────────────────────────────────────────────────────────────────────────────
  HTML not displaying correctly:
    → Use modern browser (Chrome, Firefox, Safari)
    → Enable JavaScript
    → Check file isn't corrupted
  
  Track links not working:
    → Verify files exist at specified paths
    → Check file permissions
    → Use absolute paths or URLs
    → Verify genome browser compatibility
  
  JSON parsing errors:
    → Check schema version compatibility
    → Validate JSON syntax
    → Check for special characters in strings

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Date: $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Sample: ${SAMPLE_ID}
  Module: 14_generate_reports

================================================================================
DOCEOF

  echo "REPORT | README | Documentation created"

  ###########################################################################
  # 9) CREATE QUICK ACCESS SYMLINK
  ###########################################################################

  echo "REPORT | SYMLINK | Creating quick access link..."

  SYMLINK_DIR="!{params.output_dir}/reports"
  mkdir -p "${SYMLINK_DIR}" || true

  # Create symlink for easy access
  ln -sf "$(pwd)/${OUT_HTML}" "${SYMLINK_DIR}/${SAMPLE_ID}.html" || \
    echo "REPORT | WARNING | Could not create symlink"

  if [[ -L "${SYMLINK_DIR}/${SAMPLE_ID}.html" ]]; then
    echo "REPORT | SYMLINK | Quick access: ${SYMLINK_DIR}/${SAMPLE_ID}.html"
  fi

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "REPORT | SUMMARY | Sample: ${SAMPLE_ID}"
  echo "REPORT | SUMMARY | Divergent loci: ${DIV_COUNT}"
  echo "REPORT | SUMMARY | Genes with PI: ${PAUSING_GENES}"
  echo "REPORT | SUMMARY | Available tracks: ${TOTAL_TRACKS}/6"
  echo "REPORT | SUMMARY | Output files: 5"
  echo "REPORT | SUMMARY | Processing time: ${RENDER_TIME}s"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "REPORT | COMPLETE | sample=${SAMPLE_ID} | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}