// ============================================================================
// summarize_pol2_metrics.nf — Cohort-Level Pol-II Metrics Aggregation
// ============================================================================
//
// Purpose:
//   Aggregates per-sample Pol-II metrics into cohort-level summaries
//
// Features:
//   • Merges gene metrics from all samples into tidy format
//   • Optional differential contrasts (e.g., treatment vs control)
//   • Optional visualization (heatmaps, MA plots)
//   • Identifies top variable genes
//   • Calculates fold changes and statistics
//
// Input Manifest Format (samples.tsv):
//   sample_id  condition  timepoint  replicate  file
//   sample1    control    0h         1          /path/to/pol2_gene_metrics.tsv
//   sample2    treatment  24h        1          /path/to/pol2_gene_metrics.tsv
//
// Contrast Specification:
//   Format: "group_variable:numerator,denominator"
//   Examples:
//     • "condition:treatment,control"
//     • "timepoint:24h,0h"
//   Multiple contrasts supported
//
// Inputs:
//   path(samples_tsv) : Manifest of per-sample metric files
//
// Outputs:
//   ${params.output_dir}/09_pol2_aggregate/
//     ├── pol2_gene_metrics_merged.tsv     — Combined tidy table
//     ├── pol2_gene_metrics_contrasts.tsv  — Differential results (optional)
//     ├── plots/                           — Visualizations (optional)
//     │   ├── heatmap_*.png
//     │   └── ma_plot_*.png
//     ├── README_aggregate.txt             — Documentation
//     └── aggregate.log                    — Processing log
//
// Parameters (params.pol2.*):
//   top_n      : Top N variable genes for plots (default: 100)
//   plots      : Generate plots (default: true)
//   contrasts  : List of contrast specifications (default: [])
//
// Example Configuration:
//   params {
//     pol2 {
//       top_n = 100
//       plots = true
//       contrasts = [
//         "condition:treatment,control",
//         "timepoint:24h,0h"
//       ]
//     }
//   }
//
// ============================================================================

nextflow.enable.dsl = 2

process summarize_pol2_metrics {

  tag        'pol2-aggregate'
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/09_pol2_aggregate",
             mode: 'copy',
             overwrite: true

  conda (params.conda_pol2 ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  // NOTE: Each file is staged as "metric_N" where N is the index
  input:
    path(samples_tsv, stageAs: 'samples.tsv')
    path('metric_*')  // Will stage as metric_1, metric_2, etc.

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    path 'pol2_gene_metrics_merged.tsv',              emit: merged
    path 'pol2_gene_metrics_contrasts.tsv', optional: true, emit: contrasts
    path 'plots/**',                        optional: true, emit: plots
    path 'README_aggregate.txt',                      emit: readme
    path 'aggregate.log',                             emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # Redirect all output to log file
  exec > >(tee -a aggregate.log) 2>&1

  TIMESTAMP=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "AGGREGATE | START | cohort analysis | ts=${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLES_TSV_ORIG="samples.tsv"
  SAMPLES_TSV="samples_rewritten.tsv"
  AGGREGATOR_SCRIPT="!{projectDir}/bin/compare_pol2_metrics.py"

  # Parameters
  TOP_N=!{params.pol2?.top_n ?: 100}
  ENABLE_PLOTS=$([[ "!{params.pol2?.plots}" == "false" ]] && echo 0 || echo 1)

  echo "AGGREGATE | CONFIG | Samples manifest: ${SAMPLES_TSV}"
  echo "AGGREGATE | CONFIG | Aggregator script: ${AGGREGATOR_SCRIPT}"
  echo "AGGREGATE | CONFIG | Top N genes: ${TOP_N}"
  echo "AGGREGATE | CONFIG | Generate plots: $([ ${ENABLE_PLOTS} -eq 1 ] && echo "yes" || echo "no")"

  ###########################################################################
  # 1.5) CREATE SYMLINKS FROM STAGED FILES
  ###########################################################################

  echo "AGGREGATE | STAGE | Creating symlinks to staged files..."

  # Files are staged as metric_1, metric_2, etc. by Nextflow
  # TSV has the mapping: sample_id -> metric_N
  # Create symlinks with meaningful names: Sample_ID.pol2_gene_metrics.tsv -> metric_N

  tail -n +2 "${SAMPLES_TSV_ORIG}" | while IFS=$'\t' read -r SAMPLE_ID CONDITION TIMEPOINT REPLICATE STAGED_NAME; do
    TARGET_NAME="${SAMPLE_ID}.pol2_gene_metrics.tsv"
    
    # Nextflow stages single files as "metric_" instead of "metric_1"
    # Check both patterns
    ACTUAL_FILE=""
    if [[ -e "${STAGED_NAME}" ]]; then
      ACTUAL_FILE="${STAGED_NAME}"
    elif [[ "${STAGED_NAME}" == "metric_"* ]]; then
      # Try without number suffix for single-file case
      BASE_NAME="${STAGED_NAME%_*}"
      if [[ -e "${BASE_NAME}_" ]]; then
        ACTUAL_FILE="${BASE_NAME}_"
      fi
    fi
    
    if [[ -n "${ACTUAL_FILE}" && -e "${ACTUAL_FILE}" ]]; then
      ln -sf "${ACTUAL_FILE}" "${TARGET_NAME}"
      SIZE=$(stat -c%s "${ACTUAL_FILE}" 2>/dev/null || stat -f%z "${ACTUAL_FILE}" 2>/dev/null || echo "unknown")
      echo "AGGREGATE | STAGE | ✓ ${ACTUAL_FILE} -> ${TARGET_NAME} (${SIZE} bytes)"
    else
      echo "AGGREGATE | WARNING | ✗ ${STAGED_NAME} not found (tried ${STAGED_NAME} and metric_)"
      # List what files actually exist for debugging
      echo "AGGREGATE | DEBUG | Files in work dir: $(ls -1 metric* 2>/dev/null || echo 'none')"
    fi
  done

  # Rewrite TSV to use the meaningful names
  echo -e "sample_id\tcondition\ttimepoint\treplicate\tfile" > "${SAMPLES_TSV}"
  
  tail -n +2 "${SAMPLES_TSV_ORIG}" | while IFS=$'\t' read -r SAMPLE_ID CONDITION TIMEPOINT REPLICATE STAGED_NAME; do
    TARGET_NAME="${SAMPLE_ID}.pol2_gene_metrics.tsv"
    echo -e "${SAMPLE_ID}\t${CONDITION}\t${TIMEPOINT}\t${REPLICATE}\t${TARGET_NAME}"
  done >> "${SAMPLES_TSV}"

  echo "AGGREGATE | STAGE | Manifest rewritten with symlink names"

  # Parse contrasts from params
  cat > contrasts.txt <<'CONTRASTEOF'
!{((params.pol2?.contrasts ?: []) as List).join('\n')}
CONTRASTEOF

  # Count contrasts
  CONTRAST_COUNT=$( (grep -v '^$' contrasts.txt || true) | wc -l | tr -d ' ' )
  echo "AGGREGATE | CONFIG | Contrasts: ${CONTRAST_COUNT}"

  if [[ ${CONTRAST_COUNT} -gt 0 ]]; then
    echo "AGGREGATE | CONFIG | Contrast specifications:"
    cat contrasts.txt | grep -v '^$' | while read -r contrast; do
      echo "AGGREGATE | CONFIG |   ${contrast}"
    done
  fi

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "AGGREGATE | VALIDATE | Checking inputs..."

  VALIDATION_OK=1

  # Check Python script
  if [[ ! -f "${AGGREGATOR_SCRIPT}" ]]; then
    echo "AGGREGATE | ERROR | Aggregator script not found: ${AGGREGATOR_SCRIPT}"
    VALIDATION_OK=0
  else
    echo "AGGREGATE | VALIDATE | Aggregator script: ${AGGREGATOR_SCRIPT}"
  fi

  # Check samples manifest
  if [[ ! -s "${SAMPLES_TSV}" ]]; then
    echo "AGGREGATE | ERROR | Samples manifest missing or empty: ${SAMPLES_TSV}"
    VALIDATION_OK=0
  else
    SAMPLES_SIZE=$(stat -c%s "${SAMPLES_TSV}" 2>/dev/null || stat -f%z "${SAMPLES_TSV}" 2>/dev/null || echo "unknown")
    SAMPLES_LINES=$(wc -l < "${SAMPLES_TSV}" | tr -d ' ')
    echo "AGGREGATE | VALIDATE | Samples manifest: ${SAMPLES_SIZE} bytes (${SAMPLES_LINES} lines)"
  fi

  # Validate tools
  if command -v python3 >/dev/null 2>&1; then
    PYTHON_VERSION=$(python3 --version 2>&1 || echo "unknown")
    echo "AGGREGATE | VALIDATE | Python: ${PYTHON_VERSION}"
  else
    echo "AGGREGATE | ERROR | python3 not found in PATH"
    VALIDATION_OK=0
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "AGGREGATE | ERROR | Validation failed"
    exit 1
  fi

  ###########################################################################
  # 3) VALIDATE SAMPLES MANIFEST FORMAT
  ###########################################################################

  echo "AGGREGATE | VALIDATE | Checking manifest format..."

  # Expected header
  EXPECTED_HEADER="sample_id  condition timepoint replicate file"

  # Check header
  # Check header (normalize whitespace)
  ACTUAL_HEADER=$(head -1 "${SAMPLES_TSV}" | awk '{$1=$1};1')
  EXPECTED_HEADER_NORM=$(echo "${EXPECTED_HEADER}" | awk '{$1=$1};1')

  if [[ "${ACTUAL_HEADER}" != "${EXPECTED_HEADER_NORM}" ]]; then
    echo "AGGREGATE | ERROR | Invalid manifest header"
    echo "AGGREGATE | ERROR | Expected: ${EXPECTED_HEADER}"
    echo "AGGREGATE | ERROR | Got: ${ACTUAL_HEADER}"
    exit 1
  fi

  echo "AGGREGATE | VALIDATE | Manifest header: OK"

  # Validate each data row
  SAMPLE_COUNT=0
  INVALID_ROWS=0

  tail -n +2 "${SAMPLES_TSV}" | while IFS=$'\t' read -r SAMPLE_ID CONDITION TIMEPOINT REPLICATE FILE; do
    SAMPLE_COUNT=$((SAMPLE_COUNT + 1))
    
    # Check field count
    if [[ -z "${SAMPLE_ID}" || -z "${CONDITION}" || -z "${TIMEPOINT}" || -z "${REPLICATE}" || -z "${FILE}" ]]; then
      echo "AGGREGATE | WARNING | Row ${SAMPLE_COUNT}: incomplete fields"
      INVALID_ROWS=$((INVALID_ROWS + 1))
      continue
    fi
    
    # Check file exists
    if [[ ! -s "${FILE}" ]]; then
      echo "AGGREGATE | WARNING | Row ${SAMPLE_COUNT}: file missing or empty: ${FILE}"
      INVALID_ROWS=$((INVALID_ROWS + 1))
    fi
  done

  # Count samples (excluding header)
  SAMPLE_COUNT=$((SAMPLES_LINES - 1))
  echo "AGGREGATE | VALIDATE | Samples: ${SAMPLE_COUNT}"

  if [[ ${INVALID_ROWS} -gt 0 ]]; then
    echo "AGGREGATE | WARNING | ${INVALID_ROWS} rows with issues"
  fi

  ###########################################################################
  # 4) PREPARE OUTPUT DIRECTORIES
  ###########################################################################

  echo "AGGREGATE | SETUP | Creating output directories..."

  mkdir -p plots
  echo "AGGREGATE | SETUP | Plots directory created"

  ###########################################################################
  # 5) BUILD AGGREGATOR COMMAND
  ###########################################################################

  echo "AGGREGATE | BUILD | Building aggregator command..."

  # Base arguments
  AGGREGATOR_ARGS=(
    --samples-tsv "${SAMPLES_TSV}"
    --out-merged pol2_gene_metrics_merged.tsv
    --top-n "${TOP_N}"
  )

  # Add contrasts if specified
  if [[ ${CONTRAST_COUNT} -gt 0 ]]; then
    echo "AGGREGATE | BUILD | Adding ${CONTRAST_COUNT} contrasts..."
    
    # Read contrasts into array
    mapfile -t CONTRASTS < <(grep -v '^$' contrasts.txt)
    
    if [[ ${#CONTRASTS[@]} -gt 0 ]]; then
      AGGREGATOR_ARGS+=(
        --contrasts "${CONTRASTS[@]}"
        --out-contrasts pol2_gene_metrics_contrasts.tsv
      )
      echo "AGGREGATE | BUILD | Contrast output: pol2_gene_metrics_contrasts.tsv"
    fi
  fi

  # Add plots directory if enabled
  if [[ ${ENABLE_PLOTS} -eq 1 ]]; then
    AGGREGATOR_ARGS+=(--plots-dir plots)
    echo "AGGREGATE | BUILD | Plots will be generated in: plots/"
  fi

  # Display command (for debugging)
  echo "AGGREGATE | BUILD | Command arguments: ${#AGGREGATOR_ARGS[@]} args"

  ###########################################################################
  # 6) RUN AGGREGATION
  ###########################################################################

  echo "AGGREGATE | RUN | Running aggregation..."
  echo "AGGREGATE | RUN | This may take several minutes for large datasets..."

  AGG_START=$(date +%s)

  set +e
  python3 "${AGGREGATOR_SCRIPT}" "${AGGREGATOR_ARGS[@]}"
  AGG_RC=$?
  set -e

  AGG_END=$(date +%s)
  AGG_TIME=$((AGG_END - AGG_START))

  echo "AGGREGATE | RUN | Processing completed in ${AGG_TIME}s"

  # Handle failures
  if [[ ${AGG_RC} -ne 0 ]]; then
    echo "AGGREGATE | ERROR | Aggregation failed with exit code ${AGG_RC}"
    exit 1
  fi

  ###########################################################################
  # 7) VALIDATE OUTPUTS
  ###########################################################################

  echo "AGGREGATE | VALIDATE | Checking outputs..."

  VALIDATION_OK=1

  # Check merged table
  if [[ ! -s pol2_gene_metrics_merged.tsv ]]; then
    echo "AGGREGATE | ERROR | Merged table missing or empty"
    VALIDATION_OK=0
  else
    MERGED_SIZE=$(stat -c%s pol2_gene_metrics_merged.tsv 2>/dev/null || stat -f%z pol2_gene_metrics_merged.tsv 2>/dev/null || echo "unknown")
    MERGED_LINES=$(wc -l < pol2_gene_metrics_merged.tsv | tr -d ' ')
    MERGED_GENES=$((MERGED_LINES - 1))  # Exclude header
    echo "AGGREGATE | VALIDATE | Merged table: ${MERGED_SIZE} bytes (${MERGED_GENES} genes)"
  fi

  # Check contrasts table if expected
  if [[ ${CONTRAST_COUNT} -gt 0 ]]; then
    if [[ -s pol2_gene_metrics_contrasts.tsv ]]; then
      CONTRAST_SIZE=$(stat -c%s pol2_gene_metrics_contrasts.tsv 2>/dev/null || stat -f%z pol2_gene_metrics_contrasts.tsv 2>/dev/null || echo "unknown")
      CONTRAST_LINES=$(wc -l < pol2_gene_metrics_contrasts.tsv | tr -d ' ')
      echo "AGGREGATE | VALIDATE | Contrasts table: ${CONTRAST_SIZE} bytes (${CONTRAST_LINES} lines)"
    else
      echo "AGGREGATE | WARNING | Contrasts requested but output missing"
    fi
  fi

  # Check plots if enabled
  if [[ ${ENABLE_PLOTS} -eq 1 ]]; then
    if [[ -d plots ]]; then
      PLOT_COUNT=$(find plots -name "*.png" 2>/dev/null | wc -l | tr -d ' ')
      if [[ ${PLOT_COUNT} -gt 0 ]]; then
        echo "AGGREGATE | VALIDATE | Plots generated: ${PLOT_COUNT} PNG files"
      else
        echo "AGGREGATE | VALIDATE | No plots generated (empty results or insufficient data)"
        rmdir plots 2>/dev/null || true
      fi
    fi
  fi

  if [[ ${VALIDATION_OK} -eq 0 ]]; then
    echo "AGGREGATE | ERROR | Output validation failed"
    exit 1
  fi

  ###########################################################################
  # 8) PARSE RESULTS
  ###########################################################################

  echo "AGGREGATE | RESULTS | Parsing aggregation results..."

  # Initialize variables to avoid unbound errors
  SIG_COUNT="NA"
  PLOT_COUNT=0

  # Get sample counts from merged table
  if [[ -s pol2_gene_metrics_merged.tsv ]]; then
    # Count unique samples (columns after gene info)
    HEADER_LINE=$(head -1 pol2_gene_metrics_merged.tsv)
    TOTAL_COLS=$(echo "${HEADER_LINE}" | awk -F'\t' '{print NF}')
    
    echo "AGGREGATE | RESULTS | Merged table columns: ${TOTAL_COLS}"
    echo "AGGREGATE | RESULTS | Genes: ${MERGED_GENES}"
  fi

  # Parse contrast results if available
  if [[ -s pol2_gene_metrics_contrasts.tsv ]]; then
    CONTRAST_GENES=$(tail -n +2 pol2_gene_metrics_contrasts.tsv | wc -l | tr -d ' ')
    echo "AGGREGATE | RESULTS | Contrasts: ${CONTRAST_GENES} genes analyzed"
    
    # Count significant genes if possible (log2FC and padj columns)
    SIG_COUNT=$(tail -n +2 pol2_gene_metrics_contrasts.tsv | \
                awk -F'\t' 'NF>=6 && $5!="NA" && $6!="NA" && ($5>1 || $5<-1) && $6<0.05' | \
                wc -l | tr -d ' ' || echo "NA")
    
    if [[ "${SIG_COUNT}" != "NA" ]]; then
      echo "AGGREGATE | RESULTS | Significant genes (|log2FC|>1, padj<0.05): ${SIG_COUNT}"
    fi
  fi

  ###########################################################################
  # 9) CREATE README
  ###########################################################################

  echo "AGGREGATE | README | Creating documentation..."

  cat > README_aggregate.txt <<'DOCEOF'
================================================================================
POL-II METRICS AGGREGATION — COHORT ANALYSIS
================================================================================

OVERVIEW
────────────────────────────────────────────────────────────────────────────
  Cohort-level aggregation of per-sample Pol-II metrics into unified tables
  for comparative analysis across conditions, timepoints, and replicates.

PROCESSING SUMMARY
────────────────────────────────────────────────────────────────────────────
  Samples processed:    ${SAMPLE_COUNT}
  Genes analyzed:       ${MERGED_GENES}
  Contrasts performed:  ${CONTRAST_COUNT}
  Processing time:      ${AGG_TIME}s

AGGREGATION METHOD
────────────────────────────────────────────────────────────────────────────
  1. Load per-sample metric files from manifest
  2. Extract key columns: gene_id, gene_name, pi_len_norm, pi_raw, 
     body_cpm, tss_cpm, tss_density, body_density
  3. Merge into tidy long-format table
  4. Calculate summary statistics per gene
  5. Identify top variable genes (CV, variance)
  6. Optional: Perform differential contrasts
  7. Optional: Generate visualizations

METRICS AGGREGATED
────────────────────────────────────────────────────────────────────────────
  Primary Metrics:
    • Pausing Index (pi_len_norm, pi_raw)
    • TSS Coverage (tss_cpm, tss_density_per_bp)
    • Body Coverage (body_cpm, body_density_per_bp)
  
  Derived Statistics:
    • Mean across replicates
    • Standard deviation
    • Coefficient of variation (CV)
    • Min/Max values

FILES
────────────────────────────────────────────────────────────────────────────

pol2_gene_metrics_merged.tsv:
  Combined tidy table with all samples
  
  Format: One row per gene-sample combination
  Columns:
    • gene_id         — Gene identifier
    • gene_name       — Gene symbol
    • sample_id       — Sample identifier
    • condition       — Experimental condition
    • timepoint       — Time point
    • replicate       — Biological replicate
    • pi_len_norm     — Length-normalized pausing index
    • pi_raw          — Raw pausing index
    • tss_cpm         — TSS coverage (CPM)
    • body_cpm        — Body coverage (CPM)
    • tss_density     — TSS density (reads/bp)
    • body_density    — Body density (reads/bp)
  
  Lines: ${MERGED_LINES} (${MERGED_GENES} genes × samples)
  Size: ${MERGED_SIZE} bytes

pol2_gene_metrics_contrasts.tsv (optional):
  Differential analysis results
  
  Generated when: params.pol2.contrasts specified
  Contrasts analyzed: ${CONTRAST_COUNT}
  
  Format: One row per gene-contrast combination
  Columns:
    • gene_id         — Gene identifier
    • gene_name       — Gene symbol
    • contrast        — Contrast specification
    • mean_numerator  — Mean in numerator group
    • mean_denominator— Mean in denominator group
    • log2_fold_change— Log2(numerator/denominator)
    • pvalue          — Statistical p-value
    • padj            — Adjusted p-value (FDR)
  
  $([ -s pol2_gene_metrics_contrasts.tsv ] && cat <<STATS
  Lines: ${CONTRAST_LINES}
  Significant genes: ${SIG_COUNT} (|log2FC|>1, padj<0.05)
STATS
)

plots/ directory (optional):
  Visualization files
  
  Generated when: params.pol2.plots = true
  File format: PNG (Matplotlib)
  
  Types of plots:
    • Heatmaps: Top N variable genes per metric
    • MA plots: log2FC vs mean expression per contrast
    • Distribution plots: Metric distributions per group
  
  $([ -d plots ] && echo "  Plot count: ${PLOT_COUNT}" || echo "  No plots generated")

README_aggregate.txt:
  This documentation file

aggregate.log:
  Complete processing log with timestamps

CONTRAST SPECIFICATION
────────────────────────────────────────────────────────────────────────────
  Format: "grouping_variable:numerator,denominator"
  
  Examples:
    condition:treatment,control
      → Compare treatment vs control condition
    
    timepoint:24h,0h
      → Compare 24h vs 0h timepoint
  
  Multiple contrasts supported:
    params.pol2.contrasts = [
      "condition:KO,WT",
      "timepoint:late,early"
    ]

  Contrasts used:
$([ ${CONTRAST_COUNT} -gt 0 ] && cat contrasts.txt | grep -v '^$' | sed 's/^/    • /' || echo "    (none)")

TOP VARIABLE GENES
────────────────────────────────────────────────────────────────────────────
  Selection Method:
    • Calculate coefficient of variation (CV = SD/mean) per gene
    • Rank genes by CV
    • Select top N (default: ${TOP_N})
  
  Use cases:
    • Identify most dynamic genes
    • Focus on responsive genes
    • Reduce visualization complexity

STATISTICAL METHODS
────────────────────────────────────────────────────────────────────────────
  Differential Analysis:
    • Test: Two-sample t-test (default)
    • Multiple testing correction: Benjamini-Hochberg (FDR)
    • Significance threshold: padj < 0.05, |log2FC| > 1
  
  Missing Data:
    • Genes with zero coverage excluded from ratios
    • NA values propagated through calculations
    • Minimum sample size checked per contrast

QUALITY CONTROL
────────────────────────────────────────────────────────────────────────────

Expected Results:
  • All samples successfully merged
  • Gene count matches input annotations
  • Replicates cluster by condition
  • Contrasts show expected directionality

Red Flags:
  • Very few genes in merged table
  • High proportion of NA values
  • Replicates don't cluster
  • No significant genes in expected contrasts

Troubleshooting:
  • Missing samples: Check file paths in manifest
  • Low gene count: Check per-sample metric quality
  • NA values: Check for zero coverage samples
  • Clustering issues: Check sample labeling

DOWNSTREAM USAGE
────────────────────────────────────────────────────────────────────────────
  These aggregated metrics can be used for:
  
  1. Comparative Analysis:
     - Compare Pol-II dynamics across conditions
     - Identify condition-specific pausing
     - Track changes over time courses
  
  2. Integration:
     - Correlate with gene expression (RNA-seq)
     - Compare with ChIP-seq (NELF, DSIF, Pol-II)
     - Integrate with epigenetic marks
  
  3. Functional Analysis:
     - Gene set enrichment on paused genes
     - Pathway analysis of differentially paused genes
     - Transcription factor target analysis
  
  4. Visualization:
     - Custom heatmaps and clustering
     - Time series analysis
     - Principal component analysis (PCA)

PARAMETERS USED
────────────────────────────────────────────────────────────────────────────
  Top N genes:      ${TOP_N}
  Generate plots:   $([ ${ENABLE_PLOTS} -eq 1 ] && echo "Yes" || echo "No")
  Contrasts:        ${CONTRAST_COUNT}

TECHNICAL NOTES
────────────────────────────────────────────────────────────────────────────
  • Efficient column selection (usecols) for large files
  • TSV format for universal compatibility
  • Memory-efficient processing
  • Pandas-based data handling
  • Matplotlib visualizations
  • FDR correction via statsmodels

DATA FORMAT COMPATIBILITY
────────────────────────────────────────────────────────────────────────────
  Output tables can be imported into:
  • R (read.delim, readr::read_tsv)
  • Python (pandas.read_csv)
  • Excel/LibreOffice Calc
  • GraphPad Prism
  • DESeq2 (for further analysis)

GENERATED
────────────────────────────────────────────────────────────────────────────
  Pipeline: TrackTx PRO-seq
  Date: $(date -u +"%Y-%m-%d %H:%M:%S UTC")
  Module: 12_summarize_pol2_metrics
  Samples: ${SAMPLE_COUNT}

================================================================================
DOCEOF

  echo "AGGREGATE | README | Documentation created"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "AGGREGATE | SUMMARY | Cohort Analysis Complete"
  echo "AGGREGATE | SUMMARY | Samples: ${SAMPLE_COUNT}"
  echo "AGGREGATE | SUMMARY | Genes: ${MERGED_GENES}"
  echo "AGGREGATE | SUMMARY | Contrasts: ${CONTRAST_COUNT}"
  if [[ "${SIG_COUNT}" != "NA" && ${CONTRAST_COUNT} -gt 0 ]]; then
    echo "AGGREGATE | SUMMARY | Significant genes: ${SIG_COUNT}"
  fi
  if [[ ${ENABLE_PLOTS} -eq 1 && ${PLOT_COUNT} -gt 0 ]]; then
    echo "AGGREGATE | SUMMARY | Plots: ${PLOT_COUNT}"
  fi
  echo "AGGREGATE | SUMMARY | Processing time: ${AGG_TIME}s"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "AGGREGATE | COMPLETE | cohort analysis | ts=${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  '''
}