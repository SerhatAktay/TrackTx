// ============================================================================
// summarize_polymerase_metrics.nf — Cohort-Level Polymerase Metrics Aggregation
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
//   sample1    control    0h         1          /path/to/pol_gene_metrics.tsv
//   sample2    treatment  24h        1          /path/to/pol_gene_metrics.tsv
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
//   ${params.output_dir}/09_pol_aggregate/
//     ├── pol_gene_metrics_merged.tsv     — Combined tidy table
//     ├── pol_gene_metrics_contrasts.tsv  — Differential results (optional)
//     ├── plots/                           — Visualizations (optional)
//     │   ├── heatmap_*.png
//     │   └── ma_plot_*.png
//     ├── README_aggregate.txt             — Documentation
//     └── aggregate.log                    — Processing log
//
// Parameters (params.pol.*):
//   top_n      : Top N variable genes for plots (default: 100)
//   plots      : Generate plots (default: true)
//   contrasts  : List of contrast specifications (default: [])
//
// Example Configuration:
//   params {
//     pol {
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


process summarize_polymerase_metrics {

  tag        'pol-aggregate'
  label      'conda'
  cache      'lenient'

  publishDir "${params.output_dir}/09_pol_aggregate",
             mode: params.publish_mode,
             overwrite: true

  conda (params.conda_pol ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  // NOTE: Each file is staged as "metric_N" where N is the index
  input:
    path(samples_tsv, stageAs: 'samples.tsv')
    path('metric_*')  // Will stage as metric_1, metric_2, etc.

  // ── Outputs ───────────────────────────────────────────────────────────────
  output:
    path 'pol_gene_metrics_merged.tsv',              emit: merged
    path 'pol_gene_metrics_contrasts.tsv', optional: true, emit: contrasts
    path 'plots/**',                        optional: true, emit: plots
    path 'README_aggregate.txt',                      emit: readme
    path 'aggregate.log',                             emit: log

  // ── Main Script ───────────────────────────────────────────────────────────
  script:
  """
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C
  # Matplotlib font cache: use TMPDIR (fast local disk) so tasks don't stall on "building font cache"
  export MPLCONFIGDIR="\${TMPDIR:-/tmp}/matplotlib"

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a aggregate.log)
  exec 2> >(tee -a aggregate.log >&2)

  # Shared error helper (defined once in bin/tracktx_error_fragment.sh)
  source tracktx_error_fragment.sh
  trap 'rc=\$?; tracktx_error "summarize_polymerase_metrics" "Unexpected process failure" "Check aggregate.log in work dir" "\$rc"' ERR

  TIMESTAMP=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "AGGREGATE | START | cohort analysis | ts=\${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLES_TSV_ORIG="samples.tsv"
  SAMPLES_TSV="samples_rewritten.tsv"
  AGGREGATOR_SCRIPT="\$(command -v compare_pol_metrics.py)"
  THREADS=${task.cpus}

  # Parameters
  TOP_N=${params.pol?.top_n ?: 100}
  ENABLE_PLOTS=\$([[ "${params.pol?.plots}" == "false" ]] && echo 0 || echo 1)
  # Differential filtering (assay-agnostic; script has the same defaults):
  #   prior_count : log2FC shrinkage prior; min_expr : expression filter floor.
  PRIOR_COUNT=${params.pol?.prior_count ?: 1.0}
  MIN_EXPR=${params.pol?.min_expr ?: 1.0}

  echo "AGGREGATE | CONFIG | Samples manifest: \${SAMPLES_TSV}"
  echo "AGGREGATE | CONFIG | Aggregator script: \${AGGREGATOR_SCRIPT}"
  echo "AGGREGATE | CONFIG | Top N genes: \${TOP_N}"
  echo "AGGREGATE | CONFIG | Generate plots: \$([ \${ENABLE_PLOTS} -eq 1 ] && echo "yes" || echo "no")"

  ###########################################################################
  # 1.5) CREATE SYMLINKS FROM STAGED FILES
  ###########################################################################

  echo "AGGREGATE | STAGE | Creating symlinks to staged files..."

  # Files are staged as metric_1, metric_2, etc. by Nextflow
  # TSV has the mapping: sample_id -> metric_N
  # Create symlinks with meaningful names: Sample_ID.pol_gene_metrics.tsv -> metric_N

  tail -n +2 "\${SAMPLES_TSV_ORIG}" | while IFS=\$'\\t' read -r SAMPLE_ID CONDITION TIMEPOINT REPLICATE STAGED_NAME; do
    TARGET_NAME="\${SAMPLE_ID}.pol_gene_metrics.tsv"
    
    # Nextflow stages single files as "metric_" instead of "metric_1"
    # Check both patterns
    ACTUAL_FILE=""
    if [[ -e "\${STAGED_NAME}" ]]; then
      ACTUAL_FILE="\${STAGED_NAME}"
    elif [[ "\${STAGED_NAME}" == "metric_"* ]]; then
      # Try without number suffix for single-file case
      BASE_NAME="\${STAGED_NAME%_*}"
      if [[ -e "\${BASE_NAME}_" ]]; then
        ACTUAL_FILE="\${BASE_NAME}_"
      fi
    fi
    
    if [[ -n "\${ACTUAL_FILE}" && -e "\${ACTUAL_FILE}" ]]; then
      ln -sf "\${ACTUAL_FILE}" "\${TARGET_NAME}"
      SIZE=\$(tracktx_size "\${ACTUAL_FILE}")
      echo "AGGREGATE | STAGE | ✓ \${ACTUAL_FILE} -> \${TARGET_NAME} (\${SIZE} bytes)"
    else
      echo "AGGREGATE | WARNING | ✗ \${STAGED_NAME} not found (tried \${STAGED_NAME} and metric_)"
      # List what files actually exist for debugging
      echo "AGGREGATE | DEBUG | Files in work dir: \$(ls -1 metric* 2>/dev/null || echo 'none')"
    fi
  done

  # Rewrite TSV to use the meaningful names
  echo -e "sample_id\\tcondition\\ttimepoint\\treplicate\\tfile" > "\${SAMPLES_TSV}"
  
  tail -n +2 "\${SAMPLES_TSV_ORIG}" | while IFS=\$'\\t' read -r SAMPLE_ID CONDITION TIMEPOINT REPLICATE STAGED_NAME; do
    TARGET_NAME="\${SAMPLE_ID}.pol_gene_metrics.tsv"
    echo -e "\${SAMPLE_ID}\\t\${CONDITION}\\t\${TIMEPOINT}\\t\${REPLICATE}\\t\${TARGET_NAME}"
  done >> "\${SAMPLES_TSV}"

  echo "AGGREGATE | STAGE | Manifest rewritten with symlink names"

  # Parse contrasts from params.
  # Accepts three spec shapes and normalizes them to the "variable:num,denom"
  # format expected by compare_pol_metrics.py:
  #   • 4-element [cond1, tp1, cond2, tp2]  -> group:cond1_tp1,cond2_tp2  (treatment vs control)
  #   • 3-element [variable, num, denom]    -> variable:num,denom
  #   • already-formatted "variable:num,denom" string -> passed through
  #
  # NOTE: when replicates are merged (params.replicates.merge=true) every
  # condition collapses to a single pooled track (n=1), so a per-gene contrast
  # has zero residual degrees of freedom and every p-value is forced to 1.0.
  # Emitting such a table (and its MA plots) is statistically meaningless and
  # misleading, so contrasts are force-disabled in the merged case. The merged
  # summary table and descriptive heatmaps are still produced. For real
  # differential testing use the per-replicate handoff table written by
  # 11b_collect_pol_metrics_per_replicate (08b_pol_metrics_per_replicate/) as
  # input to DESeq2/edgeR downstream.
  cat > contrasts.txt <<'CONTRASTEOF'
${(params.replicates?.merge == true) ? '' : ((params.pol?.contrasts ?: []) as List).collect { spec ->
    if (spec instanceof List) {
        def s = spec as List
        if (s.size() == 4) { 'group:' + s[0] + '_' + s[1] + ',' + s[2] + '_' + s[3] }
        else if (s.size() == 3) { '' + s[0] + ':' + s[1] + ',' + s[2] }
        else { s.join(',') }
    } else { spec.toString() }
}.join('\n')}
CONTRASTEOF

  # Count contrasts
  CONTRAST_COUNT=\$( (grep -v '^\$' contrasts.txt || true) | wc -l | tr -d ' ' )
  echo "AGGREGATE | CONFIG | Contrasts: \${CONTRAST_COUNT}"
${(params.replicates?.merge == true) ? '  echo "AGGREGATE | CONFIG | Contrasts force-disabled: replicates were merged (n=1 per group, p-values not estimable). Use 08b_pol_metrics_per_replicate/ for differential testing."' : ''}

  if [[ \${CONTRAST_COUNT} -gt 0 ]]; then
    echo "AGGREGATE | CONFIG | Contrast specifications:"
    cat contrasts.txt | grep -v '^\$' | while read -r contrast; do
      echo "AGGREGATE | CONFIG |   \${contrast}"
    done
  fi

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "AGGREGATE | VALIDATE | Checking inputs..."

  # Shared resolver (bin/tracktx_error_fragment.sh): micromamba (container) ->
  # /opt/conda (container fallback) -> bare python3 (conda profile/local)
  tracktx_resolve_python

  # Check Python script
  if [[ ! -f "\${AGGREGATOR_SCRIPT}" ]]; then
    tracktx_error "summarize_polymerase_metrics" "Aggregator script not found: \${AGGREGATOR_SCRIPT}" "Ensure bin/compare_pol_metrics.py exists"
  fi
  echo "AGGREGATE | VALIDATE | Aggregator script: \${AGGREGATOR_SCRIPT}"

  # Check samples manifest
  if [[ ! -s "\${SAMPLES_TSV}" ]]; then
    tracktx_error "summarize_polymerase_metrics" "Samples manifest missing or empty: \${SAMPLES_TSV}" "Check samples manifest input"
  fi
  SAMPLES_SIZE=\$(tracktx_size "\${SAMPLES_TSV}")
  SAMPLES_LINES=\$(wc -l < "\${SAMPLES_TSV}" | tr -d ' ')
  echo "AGGREGATE | VALIDATE | Samples manifest: \${SAMPLES_SIZE} bytes (\${SAMPLES_LINES} lines)"

  # Validate tools
  if \${PYTHON_CMD} --version >/dev/null 2>&1; then
    PYTHON_VERSION=\$(\${PYTHON_CMD} --version 2>&1 || echo "unknown")
    echo "AGGREGATE | VALIDATE | Python: \${PYTHON_VERSION}"
  else
    tracktx_error "summarize_polymerase_metrics" "Python not found (tried: \${PYTHON_CMD})" "Use -profile docker"
  fi

  ###########################################################################
  # 3) VALIDATE SAMPLES MANIFEST FORMAT
  ###########################################################################

  echo "AGGREGATE | VALIDATE | Checking manifest format..."

  # Expected header
  EXPECTED_HEADER="sample_id  condition timepoint replicate file"

  # Check header
  # Check header (normalize whitespace)
  ACTUAL_HEADER=\$(head -1 "\${SAMPLES_TSV}" | awk '{\$1=\$1};1')
  EXPECTED_HEADER_NORM=\$(echo "\${EXPECTED_HEADER}" | awk '{\$1=\$1};1')

  if [[ "\${ACTUAL_HEADER}" != "\${EXPECTED_HEADER_NORM}" ]]; then
    tracktx_error "summarize_polymerase_metrics" "Invalid manifest header (expected: \${EXPECTED_HEADER}, got: \${ACTUAL_HEADER})" "Fix samples manifest format"
  fi

  echo "AGGREGATE | VALIDATE | Manifest header: OK"

  # Validate each data row (use process substitution so loop runs in main shell)
  SAMPLE_COUNT=0
  INVALID_ROWS=0

  while IFS=\$'\\t' read -r SAMPLE_ID CONDITION TIMEPOINT REPLICATE FILE; do
    SAMPLE_COUNT=\$((SAMPLE_COUNT + 1))

    # Check field count
    if [[ -z "\${SAMPLE_ID}" || -z "\${CONDITION}" || -z "\${TIMEPOINT}" || -z "\${REPLICATE}" || -z "\${FILE}" ]]; then
      echo "AGGREGATE | WARNING | Row \${SAMPLE_COUNT}: incomplete fields"
      INVALID_ROWS=\$((INVALID_ROWS + 1))
      continue
    fi

    # Check file exists
    if [[ ! -s "\${FILE}" ]]; then
      echo "AGGREGATE | WARNING | Row \${SAMPLE_COUNT}: file missing or empty: \${FILE}"
      INVALID_ROWS=\$((INVALID_ROWS + 1))
    fi
  done < <(tail -n +2 "\${SAMPLES_TSV}")

  echo "AGGREGATE | VALIDATE | Samples: \${SAMPLE_COUNT}"

  if [[ \${INVALID_ROWS} -gt 0 ]]; then
    echo "AGGREGATE | WARNING | \${INVALID_ROWS} rows with issues"
    if [[ \${INVALID_ROWS} -ge \${SAMPLE_COUNT} && \${SAMPLE_COUNT} -gt 0 ]]; then
      tracktx_error "summarize_polymerase_metrics" "All \${SAMPLE_COUNT} manifest rows are invalid (missing files or incomplete fields)" "Fix samples manifest and ensure metric files exist"
    fi
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
    --samples-tsv "\${SAMPLES_TSV}"
    --out-merged pol_gene_metrics_merged.tsv
    --top-n "\${TOP_N}"
    --prior-count "\${PRIOR_COUNT}"
    --min-expr "\${MIN_EXPR}"
    --threads "\${THREADS}"
  )

  # Add contrasts if specified
  if [[ \${CONTRAST_COUNT} -gt 0 ]]; then
    echo "AGGREGATE | BUILD | Adding \${CONTRAST_COUNT} contrasts..."
    
    # Read contrasts into array
    mapfile -t CONTRASTS < <(grep -v '^\$' contrasts.txt)
    
    if [[ \${#CONTRASTS[@]} -gt 0 ]]; then
      AGGREGATOR_ARGS+=(
        --contrasts "\${CONTRASTS[@]}"
        --out-contrasts pol_gene_metrics_contrasts.tsv
      )
      echo "AGGREGATE | BUILD | Contrast output: pol_gene_metrics_contrasts.tsv"
    fi
  fi

  # Add plots directory if enabled
  if [[ \${ENABLE_PLOTS} -eq 1 ]]; then
    AGGREGATOR_ARGS+=(--plots-dir plots)
    echo "AGGREGATE | BUILD | Plots will be generated in: plots/"
  fi

  # Display command (for debugging)
  echo "AGGREGATE | BUILD | Command arguments: \${#AGGREGATOR_ARGS[@]} args"

  ###########################################################################
  # 6) RUN AGGREGATION
  ###########################################################################

  echo "AGGREGATE | RUN | Running aggregation..."
  echo "AGGREGATE | RUN | This may take several minutes for large datasets..."

  AGG_START=\$(date +%s)

  set +e
  \${PYTHON_CMD} "\${AGGREGATOR_SCRIPT}" "\${AGGREGATOR_ARGS[@]}"
  AGG_RC=\$?
  set -e

  AGG_END=\$(date +%s)
  AGG_TIME=\$((AGG_END - AGG_START))

  echo "AGGREGATE | RUN | Processing completed in \${AGG_TIME}s"

  # Handle failures
  if [[ \${AGG_RC} -ne 0 ]]; then
    tracktx_error "summarize_polymerase_metrics" "Aggregation failed with exit code \${AGG_RC}" "Check aggregate.log in work dir" \${AGG_RC}
  fi

  ###########################################################################
  # 7) VALIDATE OUTPUTS
  ###########################################################################

  echo "AGGREGATE | VALIDATE | Checking outputs..."

  # Check merged table
  if [[ ! -s pol_gene_metrics_merged.tsv ]]; then
    tracktx_error "summarize_polymerase_metrics" "Merged table missing or empty" "Check aggregate.log in work dir"
  fi
  MERGED_SIZE=\$(tracktx_size pol_gene_metrics_merged.tsv)
  MERGED_LINES=\$(wc -l < pol_gene_metrics_merged.tsv | tr -d ' ')
  MERGED_GENES=\$((MERGED_LINES - 1))  # Exclude header
  echo "AGGREGATE | VALIDATE | Merged table: \${MERGED_SIZE} bytes (\${MERGED_GENES} genes)"

  # Check contrasts table if expected
  if [[ \${CONTRAST_COUNT} -gt 0 ]]; then
    if [[ -s pol_gene_metrics_contrasts.tsv ]]; then
      CONTRAST_SIZE=\$(tracktx_size pol_gene_metrics_contrasts.tsv)
      CONTRAST_LINES=\$(wc -l < pol_gene_metrics_contrasts.tsv | tr -d ' ')
      echo "AGGREGATE | VALIDATE | Contrasts table: \${CONTRAST_SIZE} bytes (\${CONTRAST_LINES} lines)"
    else
      echo "AGGREGATE | WARNING | Contrasts requested but output missing"
    fi
  fi

  # Check plots if enabled
  if [[ \${ENABLE_PLOTS} -eq 1 ]]; then
    if [[ -d plots ]]; then
      PLOT_COUNT=\$(find plots -name "*.png" 2>/dev/null | wc -l | tr -d ' ')
      if [[ \${PLOT_COUNT} -gt 0 ]]; then
        echo "AGGREGATE | VALIDATE | Plots generated: \${PLOT_COUNT} PNG files"
      else
        echo "AGGREGATE | VALIDATE | No plots generated (empty results or insufficient data)"
        rmdir plots 2>/dev/null || true
      fi
    fi
  fi

  ###########################################################################
  # 8) PARSE RESULTS
  ###########################################################################

  echo "AGGREGATE | RESULTS | Parsing aggregation results..."

  # Initialize variables to avoid unbound errors
  SIG_COUNT="NA"
  PLOT_COUNT=0

  # Get sample counts from merged table
  if [[ -s pol_gene_metrics_merged.tsv ]]; then
    # Count unique samples (columns after gene info)
    HEADER_LINE=\$(head -1 pol_gene_metrics_merged.tsv)
    TOTAL_COLS=\$(echo "\${HEADER_LINE}" | awk -F'\\t' '{print NF}')
    
    echo "AGGREGATE | RESULTS | Merged table columns: \${TOTAL_COLS}"
    echo "AGGREGATE | RESULTS | Genes: \${MERGED_GENES}"
  fi

  # Parse contrast results if available
  if [[ -s pol_gene_metrics_contrasts.tsv ]]; then
    CONTRAST_GENES=\$(tail -n +2 pol_gene_metrics_contrasts.tsv | wc -l | tr -d ' ')
    echo "AGGREGATE | RESULTS | Contrasts: \${CONTRAST_GENES} genes analyzed"
    
    # Count significant genes: log2FC col 9, padj col 11
    SIG_COUNT=\$(tail -n +2 pol_gene_metrics_contrasts.tsv | \\
                awk -F'\\t' 'NF>=11 && \$9!="" && \$9!="NA" && \$11!="" && \$11!="NA" && (\$9+0>1 || \$9+0<-1) && \$11+0<0.05 && \$11+0>0' | \\
                wc -l | tr -d ' ' || echo "NA")
    
    if [[ "\${SIG_COUNT}" != "NA" ]]; then
      echo "AGGREGATE | RESULTS | Significant genes (|log2FC|>1, padj<0.05): \${SIG_COUNT}"
    fi
  fi

  ###########################################################################
  # 9) CREATE README
  ###########################################################################

  echo "AGGREGATE | README | Creating documentation..."

  cat > README_aggregate.txt <<DOCEOF
POL-II METRICS AGGREGATION — COHORT ANALYSIS
────────────────────────────────────────────────────────────────────────────
  \${SAMPLE_COUNT} samples, \${MERGED_GENES} genes, \${CONTRAST_COUNT} contrasts, \${AGG_TIME}s.

  pol_gene_metrics_merged.tsv (\${MERGED_LINES} lines): one row per gene-sample --
    gene_id, gene_name, sample_id, condition, timepoint, replicate, pi_len_norm,
    pi_raw, tss_cpm, body_cpm, tss_density, body_density

  pol_gene_metrics_contrasts.tsv (if params.pol.contrasts set; FORCE-DISABLED
    when params.replicates.merge=true -- merged tracks are n=1/condition, so
    contrasts have zero residual degrees of freedom; use
    08b_pol_metrics_per_replicate/ for real differential testing instead):
    gene_id, gene_name, contrast, mean_numerator, mean_denominator,
    log2_fold_change, pvalue, padj (Benjamini-Hochberg)
    \$([ -s pol_gene_metrics_contrasts.tsv ] && cat <<STATS
    This run: \${CONTRAST_LINES} lines, \${SIG_COUNT} significant (|log2FC|>1, padj<0.05)
STATS
)
  Contrast spec: "variable:numerator,denominator" e.g. "condition:treatment,control".
  Used: \$([ \${CONTRAST_COUNT} -gt 0 ] && cat contrasts.txt | grep -v '^\$' | tr '\\n' ';' || echo "(none)")

  plots/ (if params.pol.plots=true): heatmaps of top \${TOP_N} CV-ranked genes,
    MA plots per contrast. \$([ -d plots ] && echo "\${PLOT_COUNT} generated" || echo "none generated this run").
DOCEOF

  echo "AGGREGATE | README | Documentation created"

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "AGGREGATE | SUMMARY | Cohort Analysis Complete"
  echo "AGGREGATE | SUMMARY | Samples: \${SAMPLE_COUNT}"
  echo "AGGREGATE | SUMMARY | Genes: \${MERGED_GENES}"
  echo "AGGREGATE | SUMMARY | Contrasts: \${CONTRAST_COUNT}"
  if [[ "\${SIG_COUNT}" != "NA" && \${CONTRAST_COUNT} -gt 0 ]]; then
    echo "AGGREGATE | SUMMARY | Significant genes: \${SIG_COUNT}"
  fi
  if [[ \${ENABLE_PLOTS} -eq 1 && \${PLOT_COUNT} -gt 0 ]]; then
    echo "AGGREGATE | SUMMARY | Plots: \${PLOT_COUNT}"
  fi
  echo "AGGREGATE | SUMMARY | Processing time: \${AGG_TIME}s"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "AGGREGATE | COMPLETE | cohort analysis | ts=\${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  """
}