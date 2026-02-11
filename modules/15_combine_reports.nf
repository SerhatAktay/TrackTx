// ============================================================================
// combine_reports.nf — Cohort-Level Report Aggregation
// ============================================================================
//
// Purpose:
//   Aggregates per-sample reports into unified cohort-level summaries
//
// Features:
//   • Robust JSON intake with flexible naming (*.summary.json, *.report.json)
//   • Single-page application (SPA) HTML with interactive navigation
//   • Comprehensive TSV and JSON exports
//   • Optional functional region totals aggregation
//   • Landing page generation for easy navigation
//   • Pipeline metadata integration
//
// Report Contents:
//   1. Cohort Overview (samples, conditions, timepoints)
//   2. Quality Control Summary (across all samples)
//   3. Aggregate Statistics (divergent loci, pausing indices)
//   4. Per-Sample Metrics Table (interactive, sortable)
//   5. Functional Region Composition (cohort-wide)
//   6. Pipeline Execution Summary (runtime, resources)
//
// Input Discovery:
//   Accepts files and/or directories
//   Automatically finds: *.summary.json, *.report.json (case-insensitive)
//   Handles gzipped files: *.json.gz
//
// Inputs:
//   path('report_*.json') : Per-sample JSON reports (staged sequentially)
//
// Outputs:
//   ${params.output_dir}/11_reports/cohort/
//     ├── global_summary.html          — Interactive SPA report
//     ├── global_summary.tsv           — Cohort metrics table
//     ├── global_summary.json          — Structured cohort data
//     ├── global_region_totals.tsv     — Aggregated region counts (optional)
//     ├── README_cohort.txt            — Documentation
//     └── combine.log                  — Processing log
//
//   Quick access files:
//     ${params.output_dir}/
//       ├── index.html                 — Landing page with links
//       └── global_summary.html        — Copy of cohort report
//
// Output Formats:
//   HTML: Single-file SPA, no external dependencies, works offline
//   TSV:  Tab-separated, compatible with R/Python/Excel
//   JSON: Versioned schema, programmatic access
//
// ============================================================================

nextflow.enable.dsl = 2

process combine_reports {

  tag        { "cohort" }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/11_reports/cohort",
             mode: 'copy',
             overwrite: true

  conda (params.conda_pol2 ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  // Stage files with sequential names to avoid collisions (like metric_1, metric_2, etc.)
  input:
    path 'report_*.json'

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    path "global_summary.html",           emit: html
    path "global_summary.tsv",            emit: tsv
    path "global_summary.json",           emit: json
    path "global_region_totals.tsv",      optional: true, emit: regions
    path "README_cohort.txt",             emit: readme
    path "combine.log",                   emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Redirect all output to log file
  exec > >(tee -a combine.log) 2>&1

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "COHORT | START | report aggregation | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  COMBINER_SCRIPT="!{projectDir}/bin/combine_reports.py"
  
  # Pipeline metadata
  PIPELINE_VERSION="!{workflow.manifest.version ?: 'dev'}"
  RUN_NAME="!{workflow.runName ?: 'unnamed'}"
  DURATION="!{workflow.duration ?: 'unknown'}"
  PROFILE="!{workflow.profile ?: 'unknown'}"
  
  # Output files
  OUT_HTML="global_summary.html"
  OUT_TSV="global_summary.tsv"
  OUT_JSON="global_summary.json"
  OUT_REGIONS="global_region_totals.tsv"
  OUT_README="README_cohort.txt"

  echo "COHORT | CONFIG | Combiner script: ${COMBINER_SCRIPT}"
  echo "COHORT | CONFIG | Pipeline version: ${PIPELINE_VERSION}"
  echo "COHORT | CONFIG | Run name: ${RUN_NAME}"
  echo "COHORT | CONFIG | Duration: ${DURATION}"
  echo "COHORT | CONFIG | Profile: ${PROFILE}"

  ###########################################################################
  # 2) DISCOVER INPUT FILES
  ###########################################################################

  echo "COHORT | DISCOVER | Finding per-sample JSON files..."

  # Files are staged as report_1.json, report_2.json, etc. by Nextflow
  JSON_FILES=()
  for json_file in report_*.json; do
    if [[ -f "${json_file}" ]]; then
      JSON_FILES+=("${json_file}")
      echo "COHORT | DISCOVER | Found: ${json_file}"
    fi
  done

  TOTAL_JSON=${#JSON_FILES[@]}
  echo "COHORT | DISCOVER | JSON files found: ${TOTAL_JSON}"

  if [[ ${TOTAL_JSON} -eq 0 ]]; then
    echo "COHORT | ERROR | No JSON files found"
    echo "COHORT | ERROR | Expected files: *.summary.json, *.report.json, or *.json"
    exit 1
  fi

  ###########################################################################
  # 3) VALIDATE INPUT FILES
  ###########################################################################

  echo "COHORT | VALIDATE | Checking input files..."

  VALIDATION_OK=1
  MISSING_COUNT=0
  EMPTY_COUNT=0
  TOTAL_SIZE=0

  for JSON_FILE in "${JSON_FILES[@]}"; do
    if [[ ! -e "${JSON_FILE}" ]]; then
      echo "COHORT | ERROR | Missing file: ${JSON_FILE}"
      MISSING_COUNT=$((MISSING_COUNT + 1))
      VALIDATION_OK=0
      continue
    fi
    
    FILE_SIZE=$(stat -c%s "${JSON_FILE}" 2>/dev/null || stat -f%z "${JSON_FILE}" 2>/dev/null || echo 0)
    
    if [[ ${FILE_SIZE} -eq 0 ]]; then
      echo "COHORT | ERROR | Empty file: ${JSON_FILE}"
      EMPTY_COUNT=$((EMPTY_COUNT + 1))
      VALIDATION_OK=0
      continue
    fi
    
    TOTAL_SIZE=$((TOTAL_SIZE + FILE_SIZE))
    echo "COHORT | VALIDATE | $(basename ${JSON_FILE}): ${FILE_SIZE} bytes"
  done

  echo "COHORT | VALIDATE | Total input size: ${TOTAL_SIZE} bytes"

  if [[ ${MISSING_COUNT} -gt 0 ]]; then
    echo "COHORT | ERROR | Missing files: ${MISSING_COUNT}"
  fi

  if [[ ${EMPTY_COUNT} -gt 0 ]]; then
    echo "COHORT | ERROR | Empty files: ${EMPTY_COUNT}"
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "COHORT | ERROR | Validation failed"
    exit 1
  fi

  ###########################################################################
  # 4) VALIDATE TOOLS
  ###########################################################################

  echo "COHORT | VALIDATE | Checking required tools..."

  # Check combiner script
  if [[ ! -e "${COMBINER_SCRIPT}" ]]; then
    echo "COHORT | ERROR | Combiner script not found: ${COMBINER_SCRIPT}"
    exit 1
  else
    echo "COHORT | VALIDATE | Combiner script: ${COMBINER_SCRIPT}"
  fi

  # Check Python
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_VERSION=$(python3 --version 2>&1 || echo "unknown")
    echo "COHORT | VALIDATE | Python: ${PYTHON_VERSION}"
  else
    echo "COHORT | ERROR | python3 not found in PATH"
    exit 1
  fi

  ###########################################################################
  # 5) RUN REPORT COMBINER
  ###########################################################################

  echo "COHORT | COMBINE | Aggregating per-sample reports..."
  echo "COHORT | COMBINE | Processing ${TOTAL_JSON} samples..."

  COMBINE_START=$(date +%s)

  set +e
  python3 "${COMBINER_SCRIPT}" \
    --inputs "${JSON_FILES[@]}" \
    --out-tsv "${OUT_TSV}" \
    --out-json "${OUT_JSON}" \
    --out-html "${OUT_HTML}" \
    --out-regions "${OUT_REGIONS}" \
    --pipeline-version "${PIPELINE_VERSION}" \
    --run-name "${RUN_NAME}" \
    --duration "${DURATION}" \
    --profile "${PROFILE}"
  
  COMBINE_RC=$?
  set -e

  COMBINE_END=$(date +%s)
  COMBINE_TIME=$((COMBINE_END - COMBINE_START))

  echo "COHORT | COMBINE | Processing completed in ${COMBINE_TIME}s"

  if [[ ${COMBINE_RC} -ne 0 ]]; then
    echo "COHORT | ERROR | Combiner failed with exit code ${COMBINE_RC}"
    exit 1
  fi

  ###########################################################################
  # 6) VALIDATE OUTPUTS
  ###########################################################################

  echo "COHORT | VALIDATE | Checking output files..."

  VALIDATION_OK=1

  # Check required outputs
  for OUTPUT in "${OUT_HTML}" "${OUT_TSV}" "${OUT_JSON}"; do
    if [[ ! -s "${OUTPUT}" ]]; then
      echo "COHORT | ERROR | Missing or empty output: ${OUTPUT}"
      VALIDATION_OK=0
    else
      OUTPUT_SIZE=$(stat -c%s "${OUTPUT}" 2>/dev/null || stat -f%z "${OUTPUT}" 2>/dev/null || echo "unknown")
      echo "COHORT | VALIDATE | $(basename ${OUTPUT}): ${OUTPUT_SIZE} bytes"
    fi
  done

  # Check optional regions file
  if [[ -s "${OUT_REGIONS}" ]]; then
    REGIONS_SIZE=$(stat -c%s "${OUT_REGIONS}" 2>/dev/null || stat -f%z "${OUT_REGIONS}" 2>/dev/null || echo "unknown")
    REGIONS_LINES=$(wc -l < "${OUT_REGIONS}" | tr -d ' ')
    echo "COHORT | VALIDATE | Region totals: ${REGIONS_SIZE} bytes (${REGIONS_LINES} lines)"
  else
    echo "COHORT | VALIDATE | Region totals: not generated"
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "COHORT | ERROR | Output validation failed"
    exit 1
  fi

  ###########################################################################
  # 7) PARSE OUTPUT STATISTICS
  ###########################################################################

  echo "COHORT | STATS | Extracting summary statistics..."

  # Count samples in TSV
  if [[ -s "${OUT_TSV}" ]]; then
    TSV_LINES=$(wc -l < "${OUT_TSV}" | tr -d ' ')
    TSV_SAMPLES=$((TSV_LINES - 1))  # Exclude header
    echo "COHORT | STATS | TSV samples: ${TSV_SAMPLES}"
  fi

  # Parse JSON for key metrics (if jq available)
  if command -v jq >/dev/null 2>&1 && [[ -s "${OUT_JSON}" ]]; then
    SAMPLE_COUNT=$(jq -r '.samples | length' "${OUT_JSON}" 2>/dev/null || echo "NA")
    echo "COHORT | STATS | JSON samples: ${SAMPLE_COUNT}"
    
    # Try to get condition/timepoint info
    CONDITIONS=$(jq -r '.samples | map(.condition) | unique | length' "${OUT_JSON}" 2>/dev/null || echo "NA")
    echo "COHORT | STATS | Conditions: ${CONDITIONS}"
  else
    echo "COHORT | STATS | jq not available, skipping JSON parsing"
  fi

  ###########################################################################
  # 8) CREATE QUICK ACCESS COPIES
  ###########################################################################

  echo "COHORT | COPY | Creating quick access files..."

  # Copy cohort HTML to output root
  ROOT_HTML="!{params.output_dir}/global_summary.html"
  mkdir -p "!{params.output_dir}" || true
  cp -f "${OUT_HTML}" "${ROOT_HTML}" || \
    echo "COHORT | WARNING | Could not copy to ${ROOT_HTML}"
  
  if [[ -f "${ROOT_HTML}" ]]; then
    echo "COHORT | COPY | Cohort report: ${ROOT_HTML}"
  fi

  ###########################################################################
  # 9) CREATE LANDING PAGE
  ###########################################################################

  echo "COHORT | LANDING | Creating landing page..."

  ROOT_DIR="!{params.output_dir}"
  INDEX_HTML="${ROOT_DIR}/index.html"
  
  # Relative paths from root
  COHORT_REL="11_reports/cohort/global_summary.html"
  REPORTS_REL="reports/"
  NF_REPORT_REL="nf_report.html"
  NF_TIMELINE_REL="nf_timeline.html"
  NF_DAG_REL="nf_dag.png"
  NF_TRACE_REL="trace/trace.txt"

  cat > "${INDEX_HTML}" <<'LANDINGEOF'
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>TrackTx Pipeline Results</title>
  <style>
    * { margin: 0; padding: 0; box-sizing: border-box; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, 
                   "Helvetica Neue", Arial, sans-serif;
      line-height: 1.6;
      color: #333;
      background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
      min-height: 100vh;
      padding: 2rem;
    }
    .container {
      max-width: 900px;
      margin: 0 auto;
      background: white;
      border-radius: 12px;
      box-shadow: 0 20px 60px rgba(0,0,0,0.3);
      overflow: hidden;
    }
    .header {
      background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
      color: white;
      padding: 3rem 2rem;
      text-align: center;
    }
    .header h1 {
      font-size: 2.5rem;
      margin-bottom: 0.5rem;
      font-weight: 700;
    }
    .header p {
      font-size: 1.1rem;
      opacity: 0.9;
    }
    .content {
      padding: 2rem;
    }
    .section {
      margin-bottom: 2.5rem;
    }
    .section h2 {
      font-size: 1.5rem;
      color: #667eea;
      margin-bottom: 1rem;
      padding-bottom: 0.5rem;
      border-bottom: 2px solid #e0e0e0;
    }
    .link-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
      gap: 1rem;
      margin-top: 1rem;
    }
    .link-card {
      display: block;
      padding: 1.5rem;
      background: #f8f9fa;
      border: 2px solid #e0e0e0;
      border-radius: 8px;
      text-decoration: none;
      color: inherit;
      transition: all 0.3s ease;
    }
    .link-card:hover {
      border-color: #667eea;
      background: #f0f4ff;
      transform: translateY(-2px);
      box-shadow: 0 4px 12px rgba(102,126,234,0.15);
    }
    .link-card h3 {
      font-size: 1.1rem;
      color: #667eea;
      margin-bottom: 0.5rem;
    }
    .link-card p {
      font-size: 0.9rem;
      color: #666;
    }
    .meta {
      background: #f8f9fa;
      padding: 1.5rem;
      border-radius: 8px;
      font-size: 0.9rem;
      color: #666;
    }
    .meta-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 1rem;
    }
    .meta-item {
      display: flex;
      flex-direction: column;
    }
    .meta-label {
      font-weight: 600;
      color: #333;
      margin-bottom: 0.25rem;
    }
    .footer {
      text-align: center;
      padding: 1.5rem;
      color: #999;
      font-size: 0.85rem;
      border-top: 1px solid #e0e0e0;
    }
  </style>
</head>
<body>
  <div class="container">
    <div class="header">
      <h1>TrackTx Pipeline Results</h1>
      <p>PRO-seq Analysis Pipeline</p>
    </div>
    
    <div class="content">
      <div class="section">
        <h2>Analysis Reports</h2>
        <div class="link-grid">
          <a href="COHORT_REL" class="link-card">
            <h3>📊 Cohort Summary</h3>
            <p>Aggregated metrics across all samples</p>
          </a>
          <a href="REPORTS_REL" class="link-card">
            <h3>📄 Sample Reports</h3>
            <p>Individual per-sample analysis reports</p>
          </a>
        </div>
      </div>

      <div class="section">
        <h2>Pipeline Reports</h2>
        <div class="link-grid">
          <a href="NF_REPORT_REL" class="link-card">
            <h3>📈 Execution Report</h3>
            <p>Nextflow execution summary</p>
          </a>
          <a href="NF_TIMELINE_REL" class="link-card">
            <h3>⏱️ Timeline</h3>
            <p>Task execution timeline</p>
          </a>
          <a href="NF_DAG_REL" class="link-card">
            <h3>🔀 Workflow DAG</h3>
            <p>Pipeline execution graph</p>
          </a>
          <a href="NF_TRACE_REL" class="link-card">
            <h3>📝 Trace Log</h3>
            <p>Detailed execution trace</p>
          </a>
        </div>
      </div>

      <div class="section">
        <h2>Pipeline Information</h2>
        <div class="meta">
          <div class="meta-grid">
            <div class="meta-item">
              <span class="meta-label">Pipeline</span>
              <span>TrackTx PRO-seq</span>
            </div>
            <div class="meta-item">
              <span class="meta-label">Version</span>
              <span>PIPELINE_VERSION</span>
            </div>
            <div class="meta-item">
              <span class="meta-label">Run Name</span>
              <span>RUN_NAME</span>
            </div>
            <div class="meta-item">
              <span class="meta-label">Duration</span>
              <span>DURATION</span>
            </div>
            <div class="meta-item">
              <span class="meta-label">Profile</span>
              <span>PROFILE</span>
            </div>
            <div class="meta-item">
              <span class="meta-label">Completed</span>
              <span>TIMESTAMP</span>
            </div>
          </div>
        </div>
      </div>
    </div>

    <div class="footer">
      Generated by TrackTx Pipeline • TIMESTAMP
    </div>
  </div>
</body>
</html>
LANDINGEOF

  # Substitute variables
  sed -i.bak \
    -e "s|COHORT_REL|${COHORT_REL}|g" \
    -e "s|REPORTS_REL|${REPORTS_REL}|g" \
    -e "s|NF_REPORT_REL|${NF_REPORT_REL}|g" \
    -e "s|NF_TIMELINE_REL|${NF_TIMELINE_REL}|g" \
    -e "s|NF_DAG_REL|${NF_DAG_REL}|g" \
    -e "s|NF_TRACE_REL|${NF_TRACE_REL}|g" \
    -e "s|PIPELINE_VERSION|${PIPELINE_VERSION}|g" \
    -e "s|RUN_NAME|${RUN_NAME}|g" \
    -e "s|DURATION|${DURATION}|g" \
    -e "s|PROFILE|${PROFILE}|g" \
    -e "s|TIMESTAMP|${TIMESTAMP}|g" \
    "${INDEX_HTML}"
  
  rm -f "${INDEX_HTML}.bak"

  echo "COHORT | LANDING | Landing page created: ${INDEX_HTML}"

  ###########################################################################
  # 10) CREATE README
  ###########################################################################

  echo "COHORT | README | Creating documentation..."

  cat > "${OUT_README}" <<'DOCEOF'
================================================================================
COHORT-LEVEL REPORT AGGREGATION
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Aggregates per-sample reports into unified cohort-level summaries with
  interactive visualizations and comprehensive data exports.

PROCESSING SUMMARY
────────────────────────────────────────────────────────────────────────────
  Samples processed:    ${TOTAL_JSON}
  Processing time:      ${COMBINE_TIME}s
  Pipeline version:     ${PIPELINE_VERSION}
  Run name:             ${RUN_NAME}
  Duration:             ${DURATION}
  Profile:              ${PROFILE}

COHORT REPORT FILES
────────────────────────────────────────────────────────────────────────────

global_summary.html:
  Single-page application (SPA) with:
    • Cohort overview and statistics
    • Per-sample metrics table (interactive, sortable, filterable)
    • Quality control summary across all samples
    • Aggregate divergent transcription statistics
    • Functional region composition
    • Pausing index distributions
    • Pipeline execution metadata
  
  Features:
    - Works offline (no external dependencies)
    - All CSS and JavaScript embedded
    - Interactive tables with search/sort
    - Responsive design for mobile/desktop
    - Print-friendly formatting
  
  View in web browser for best experience.

global_summary.tsv:
  Tab-separated table with per-sample metrics:
    Columns vary based on available data:
      • Sample identification (sample_id, condition, timepoint, replicate)
      • QC metrics (total_reads, map_rate, duplicate_rate, etc.)
      • Divergent transcription (loci_count)
      • Functional regions (promoter_signal, body_signal, etc.)
      • Pausing indices (median_pi, mean_pi)
      • Normalization factors
  
  Compatible with:
    - R: read.delim("global_summary.tsv")
    - Python: pandas.read_csv("global_summary.tsv", sep="\\t")
    - Excel/LibreOffice Calc
    - Command-line tools (awk, cut, grep)

global_summary.json:
  Structured JSON with complete cohort data:
    Schema includes:
      • schema_version: JSON schema version
      • pipeline: Pipeline metadata
      • samples: Array of per-sample data
      • aggregate: Cohort-wide statistics
      • qc_summary: Aggregate QC metrics
  
  Use cases:
    - API integration
    - Custom analysis scripts
    - Database import
    - Reproducibility records

global_region_totals.tsv (optional):
  Aggregated functional region counts:
    Columns:
      • region: Region name (promoter, body, CPS, etc.)
      • total_signal: Sum across all samples
      • sample_count: Number of samples contributing
      • mean_signal: Average per sample
      • std_signal: Standard deviation
  
  Generated if functional region data available.

QUICK ACCESS FILES
────────────────────────────────────────────────────────────────────────────
  Created in output root for convenience:
  
  global_summary.html:
    Copy of cohort report for quick access
    Location: !{params.output_dir}/global_summary.html
  
  index.html:
    Landing page with links to all reports:
      • Cohort summary
      • Per-sample reports
      • Nextflow execution reports
      • Pipeline metadata
    
    Location: !{params.output_dir}/index.html
    
    Features:
      - Clean, modern design
      - Responsive layout
      - Quick navigation
      - Pipeline information display

INPUT FILES
────────────────────────────────────────────────────────────────────────────
  Expected per-sample JSON files:
    • *.summary.json
    • *.report.json
    • *.json
  
  Total files processed: ${TOTAL_JSON}
  
  JSON Discovery:
    - Case-insensitive matching
    - Accepts files or directories
    - Validates existence and size
    - Robust to missing fields

COHORT REPORT CONTENTS
────────────────────────────────────────────────────────────────────────────

1. Cohort Overview
   • Total samples analyzed
   • Experimental design summary
   • Conditions and timepoints
   • Replicates per group

2. Quality Control Summary
   • Aggregate mapping statistics
   • Duplicate rates across cohort
   • Strand balance distribution
   • Coverage depth summary
   • Sample outlier detection

3. Divergent Transcription
   • Total loci across cohort
   • Distribution per sample
   • Genomic characteristics
   • Summary statistics

4. Functional Region Composition
   • Promoter signal distribution
   • Gene body signal
   • CPS signal
   • Enhancer signal
   • Other regions
   • Comparative analysis

5. Pausing Index Analysis
   • Distribution across samples
   • Mean/median pausing indices
   • Top paused genes
   • Condition comparisons

6. Per-Sample Metrics Table
   • Interactive, sortable table
   • All key metrics per sample
   • Search and filter capabilities
   • Export functionality

7. Pipeline Execution
   • Runtime information
   • Resource utilization
   • Software versions
   • Configuration parameters

USING THE COHORT REPORT
────────────────────────────────────────────────────────────────────────────

HTML Report:
  1. Open global_summary.html or index.html in web browser
  2. Navigate using table of contents
  3. Use interactive table features:
     - Click column headers to sort
     - Type in search box to filter
     - Export data if needed
  4. Review QC summary for outliers
  5. Compare metrics across conditions

TSV Export:
  # Load in R
  cohort <- read.delim("global_summary.tsv")
  
  # Filter by condition
  treatment <- subset(cohort, condition == "treatment")
  
  # Calculate statistics
  mean(cohort$map_rate_percent)
  
  # Load in Python
  import pandas as pd
  cohort = pd.read_csv("global_summary.tsv", sep="\\t")
  
  # Group by condition
  cohort.groupby("condition").mean()

JSON Data:
  # Python
  import json
  with open("global_summary.json") as f:
      cohort = json.load(f)
  
  samples = cohort["samples"]
  aggregate = cohort["aggregate"]
  
  # R
  library(jsonlite)
  cohort <- fromJSON("global_summary.json")

QUALITY ASSESSMENT
────────────────────────────────────────────────────────────────────────────

Expected Cohort Metrics:
  ✓ All samples have mapping rate >70%
  ✓ Duplicate rates relatively consistent
  ✓ Strand balance similar across samples
  ✓ Coverage depths adequate (>10×)
  ✓ Replicates cluster by condition

Red Flags:
  ✗ Individual samples with very low mapping
  ✗ Extreme outliers in duplicate rates
  ✗ Inconsistent strand bias across samples
  ✗ Large variability within replicates
  ✗ Batch effects visible in QC metrics

Troubleshooting:
  - Outlier samples: Review per-sample reports
  - Batch effects: Check processing dates, reagent lots
  - Low quality: Consider excluding samples
  - Inconsistent metrics: Verify library prep protocol

DOWNSTREAM ANALYSIS
────────────────────────────────────────────────────────────────────────────
  
  1. Comparative Analysis:
     - Compare pausing indices across conditions
     - Identify condition-specific transcription patterns
     - Assess replicate consistency
  
  2. Quality Control:
     - Identify failed or low-quality samples
     - Detect batch effects
     - Plan additional sequencing if needed
  
  3. Data Integration:
     - Export TSV for statistical analysis (R, Python)
     - Load JSON for programmatic access
     - Share HTML for collaborators
  
  4. Publication:
     - QC summary for methods section
     - Cohort statistics for results
     - Data availability statement

TECHNICAL NOTES
────────────────────────────────────────────────────────────────────────────
  • Robust to missing fields in per-sample JSONs
  • Single-file HTML works offline
  • No external dependencies required
  • TSV uses standard tab-separated format
  • JSON follows versioned schema
  • Handles large cohorts (100+ samples)

FILE FORMATS
────────────────────────────────────────────────────────────────────────────
  HTML: UTF-8 encoded, HTML5 standard
  TSV:  UTF-8 encoded, tab-separated (\\t), newline-terminated (\\n)
  JSON: UTF-8 encoded, RFC 8259 compliant

PROCESSING DETAILS
────────────────────────────────────────────────────────────────────────────
  Combiner script:  ${COMBINER_SCRIPT}
  Python version:   ${PYTHON_VERSION}
  Processing time:  ${COMBINE_TIME}s
  Samples:          ${TOTAL_JSON}
  Output size:      $(stat -c%s global_summary.html 2>/dev/null || echo "unknown") bytes (HTML)

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Version:  ${PIPELINE_VERSION}
  Date:     $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Run:      ${RUN_NAME}
  Module:   15_combine_reports

================================================================================
DOCEOF

  echo "COHORT | README | Documentation created"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "COHORT | SUMMARY | Cohort Report Generation Complete"
  echo "COHORT | SUMMARY | Samples: ${TOTAL_JSON}"
  echo "COHORT | SUMMARY | Outputs: 5 files"
  echo "COHORT | SUMMARY | Processing time: ${COMBINE_TIME}s"
  echo "COHORT | SUMMARY | Landing page: ${INDEX_HTML}"
  echo "COHORT | SUMMARY | Cohort report: ${ROOT_HTML}"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "COHORT | COMPLETE | report aggregation | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}