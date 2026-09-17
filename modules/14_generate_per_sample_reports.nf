// ============================================================================
// generate_per_sample_reports.nf — Per-Sample Comprehensive Report Generation
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
//   6. Pausing Index Statistics (strand-specific counting; auto body-offset for any organism)
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


process generate_per_sample_reports {

  tag        { sample_id }
  label      'conda'
  cache      'lenient'

  publishDir { "${params.output_dir}/11_reports/samples/${sample_id}" },
             mode: params.publish_mode,
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
  script:
  """
  #!/usr/bin/env bash
  # NOTE: Using -e (not -u) because track links may be empty strings
  set -eo pipefail
  export LC_ALL=C
  # Matplotlib font cache: use TMPDIR so tasks don't stall on "building font cache"
  export MPLCONFIGDIR="\${TMPDIR:-/tmp}/matplotlib"

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a "${sample_id}.report.log")
  exec 2> >(tee -a "${sample_id}.report.log" >&2)

  # Shared error helper (defined once in bin/tracktx_error_fragment.sh)
  source tracktx_error_fragment.sh
  trap 'rc=\$?; tracktx_error "generate_per_sample_reports" "Unexpected process failure" "Check *.report.log in work dir" "\$rc"' ERR

  TIMESTAMP=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "REPORT | START | sample=${sample_id} | ts=\${TIMESTAMP}"
  echo "════════════════════════════════════════════════════════════════════════"

  ###########################################################################
  # 1) CONFIGURATION
  ###########################################################################

  SAMPLE_ID="${sample_id}"
  CONDITION="${condition}"
  TIMEPOINT="${timepoint}"
  REPLICATE="${replicate}"

  # Core input files
  DIV_BED="${div_bed}"
  FUNC_SUM="${func_sum}"
  POL_DENS="${pol_density}"
  PAUSING_IDX="${pausing_idx}"
  NORM_FACTORS="${norm_factors}"
  DEDUP_STATS="${dedup_stats}"
  QC_JSON="${qc_json}"

  # Track links (may be empty strings)
  ALLMAP3P_POS_RAW="${allmap3p_pos_raw}"
  ALLMAP3P_NEG_RAW="${allmap3p_neg_raw}"
  POS3_CPM_BW="${pos3_cpm_bw}"
  NEG3_CPM_BW="${neg3_cpm_bw}"
  ALLMAP3P_POS_CPM_BW="${allmap3p_pos_cpm_bw}"
  ALLMAP3P_NEG_CPM_BW="${allmap3p_neg_cpm_bw}"

  # Renderer script
  RENDER_SCRIPT="\$(command -v render_sample_report.py)"

  # Parameters
  ENABLE_PLOTS=${(params.reports_plots == null) ? 0 : (params.reports_plots as int)}

  # Output files
  OUT_HTML="\${SAMPLE_ID}.report.html"
  OUT_TSV="\${SAMPLE_ID}.report.tsv"
  OUT_JSON="\${SAMPLE_ID}.report.json"
  OUT_PLOTS="\${SAMPLE_ID}.plots.html"
  OUT_README="\${SAMPLE_ID}.README_report.txt"

  echo "REPORT | CONFIG | Sample ID: \${SAMPLE_ID}"
  echo "REPORT | CONFIG | Condition: \${CONDITION}"
  echo "REPORT | CONFIG | Timepoint: \${TIMEPOINT}"
  echo "REPORT | CONFIG | Replicate: \${REPLICATE}"
  echo ""
  echo "REPORT | CONFIG | Renderer script: \${RENDER_SCRIPT}"
  echo "REPORT | CONFIG | Enable plots: \$([ \${ENABLE_PLOTS} -eq 1 ] && echo "yes" || echo "no")"

  ###########################################################################
  # 2) VALIDATE INPUTS
  ###########################################################################

  echo "REPORT | VALIDATE | Checking required input files..."

  # Shared resolver (bin/tracktx_error_fragment.sh): micromamba (container) ->
  # /opt/conda (container fallback) -> bare python3 (conda profile/local)
  tracktx_resolve_python

  # Check renderer script
  if [[ ! -e "\${RENDER_SCRIPT}" ]]; then
    tracktx_error "generate_per_sample_reports" "Renderer script not found: \${RENDER_SCRIPT}" "Ensure bin/render_sample_report.py exists"
  fi
  echo "REPORT | VALIDATE | Renderer script: \${RENDER_SCRIPT}"

  # Check core input files
  validate_file() {
    local label="\$1"
    local file="\$2"
    local allow_empty="\${3:-0}"

    if [[ "\${allow_empty}" -eq 1 ]]; then
      # 0 bytes is a legitimate upstream result here (e.g. detect_divergent_transcription
      # can correctly call 0 regions for a low-signal replicate), so only a missing file
      # is an error -- an empty one is not.
      if [[ ! -e "\${file}" ]]; then
        tracktx_error "generate_per_sample_reports" "\${label} missing: \${file}" "Check upstream modules"
      fi
    elif [[ ! -s "\${file}" ]]; then
      tracktx_error "generate_per_sample_reports" "\${label} missing or empty: \${file}" "Check upstream modules"
    fi
    FILE_SIZE=\$(tracktx_size "\${file}")
    FILE_LINES=\$(wc -l < "\${file}" 2>/dev/null | tr -d ' ' || echo 0)
    echo "REPORT | VALIDATE | \${label}: \${FILE_SIZE} bytes (\${FILE_LINES} lines)"
  }

  validate_file "Divergent bed" "\${DIV_BED}" 1
  validate_file "Functional summary" "\${FUNC_SUM}"
  validate_file "Pol-II density" "\${POL_DENS}"
  validate_file "Pausing index" "\${PAUSING_IDX}"
  validate_file "Normalization factors" "\${NORM_FACTORS}"
  validate_file "QC JSON" "\${QC_JSON}"

  # Dedup stats is optional
  if [[ -s "\${DEDUP_STATS}" ]]; then
    DEDUP_SIZE=\$(tracktx_size "\${DEDUP_STATS}")
    echo "REPORT | VALIDATE | Dedup stats: \${DEDUP_SIZE} bytes"
  else
    echo "REPORT | VALIDATE | Dedup stats: not available"
  fi

  # Validate tools
  if \${PYTHON_CMD} --version >/dev/null 2>&1; then
    PYTHON_VERSION=\$(\${PYTHON_CMD} --version 2>&1 || echo "unknown")
    echo "REPORT | VALIDATE | Python: \${PYTHON_VERSION}"
  else
    tracktx_error "generate_per_sample_reports" "Python not found (tried: \${PYTHON_CMD})" "Use -profile docker"
  fi

  ###########################################################################
  # 3) VALIDATE AND CATALOG TRACK LINKS
  ###########################################################################

  echo "REPORT | TRACKS | Validating track file links..."

  # Helper to validate track link (URL or published file path).
  # IMPORTANT: a NON-EMPTY link is treated as available WITHOUT a filesystem
  # existence check. Availability is decided upstream in main.nf from the real
  # produced files on Nextflow channels; the link itself is the *published*
  # results path (e.g. /…/05_normalized_tracks/…/x.cpm.bw), which is NOT visible
  # from inside this task's sandbox (work dir / Docker mount). The previous
  # existence re-check (\${link} stat) therefore always failed and dropped every
  # link ("Available tracks: 0/6"), even though the files exist in results. Empty
  # string = not provided (main.nf emits '' when the track wasn't produced or its
  # publish toggle is off).
  validate_track() {
    local label="\$1"
    local link="\$2"

    if [[ -z "\${link}" ]]; then
      echo "REPORT | TRACKS | \${label}: not provided"
      return 1
    fi

    if [[ "\${link}" =~ ^https?:// ]]; then
      echo "REPORT | TRACKS | \${label}: URL (\${link})"
      return 0
    fi

    # Non-empty published path → trust upstream validation (do not stat sandbox).
    echo "REPORT | TRACKS | \${label}: published path (\${link})"
    return 0
  }

  # Track availability flags
  HAVE_ALLMAP_POS_RAW=0
  HAVE_ALLMAP_NEG_RAW=0
  HAVE_POS_CPM_BW=0
  HAVE_NEG_CPM_BW=0
  HAVE_ALLMAP_POS_CPM_BW=0
  HAVE_ALLMAP_NEG_CPM_BW=0

  validate_track "AllMap 3p pos raw" "\${ALLMAP3P_POS_RAW}" && HAVE_ALLMAP_POS_RAW=1
  validate_track "AllMap 3p neg raw" "\${ALLMAP3P_NEG_RAW}" && HAVE_ALLMAP_NEG_RAW=1
  validate_track "3p pos CPM BigWig" "\${POS3_CPM_BW}" && HAVE_POS_CPM_BW=1
  validate_track "3p neg CPM BigWig" "\${NEG3_CPM_BW}" && HAVE_NEG_CPM_BW=1
  validate_track "AllMap 3p pos CPM BigWig" "\${ALLMAP3P_POS_CPM_BW}" && HAVE_ALLMAP_POS_CPM_BW=1
  validate_track "AllMap 3p neg CPM BigWig" "\${ALLMAP3P_NEG_CPM_BW}" && HAVE_ALLMAP_NEG_CPM_BW=1

  TOTAL_TRACKS=\$((HAVE_ALLMAP_POS_RAW + HAVE_ALLMAP_NEG_RAW + \\
                  HAVE_POS_CPM_BW + HAVE_NEG_CPM_BW + \\
                  HAVE_ALLMAP_POS_CPM_BW + HAVE_ALLMAP_NEG_CPM_BW))
  
  echo "REPORT | TRACKS | Available tracks: \${TOTAL_TRACKS}/6"

  ###########################################################################
  # 4) PARSE INPUT DATA FOR SUMMARY
  ###########################################################################

  echo "REPORT | PARSE | Extracting summary statistics..."

  # Count divergent loci
  DIV_COUNT=\$(grep -v '^#' "\${DIV_BED}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  echo "REPORT | PARSE | Divergent loci: \${DIV_COUNT}"

  # Parse functional regions
  if [[ -s "\${FUNC_SUM}" ]]; then
    FUNC_LINES=\$(tail -n +2 "\${FUNC_SUM}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
    echo "REPORT | PARSE | Functional regions: \${FUNC_LINES} categories"
  fi

  # Parse pausing index
  if [[ -s "\${PAUSING_IDX}" ]]; then
    PAUSING_GENES=\$(tail -n +2 "\${PAUSING_IDX}" 2>/dev/null | wc -l | tr -d ' ' || echo 0)
    echo "REPORT | PARSE | Genes with pausing index: \${PAUSING_GENES}"
  fi

  # Parse QC JSON for key metrics
  if command -v jq >/dev/null 2>&1 && [[ -s "\${QC_JSON}" ]]; then
    TOTAL_READS=\$(jq -r '.total_reads_raw // 0' "\${QC_JSON}" 2>/dev/null || echo 0)
    # NOTE: module 13's qc_pol.json has no 'map_rate_percent' key (the genuine
    # overall alignment rate lives in 02_alignments/alignment_rates_summary.tsv,
    # not here). Use the unique/MAPQ-pass percent that IS present so this stops
    # silently reporting 0.
    MAP_RATE=\$(jq -r '.mapq_pass_percent // .unique_pass_percent // 0' "\${QC_JSON}" 2>/dev/null || echo 0)
    DUP_RATE=\$(jq -r '.duplicate_perc_of_total // 0' "\${QC_JSON}" 2>/dev/null || echo 0)
    echo "REPORT | PARSE | Total reads: \${TOTAL_READS}"
    echo "REPORT | PARSE | Unique/MAPQ-pass rate: \${MAP_RATE}%"
    echo "REPORT | PARSE | Duplicate rate: \${DUP_RATE}%"
  else
    echo "REPORT | PARSE | jq not available, skipping QC metrics"
  fi

  ###########################################################################
  # 5) BUILD RENDERER COMMAND
  ###########################################################################

  echo "REPORT | BUILD | Building renderer command..."

  # Base arguments
  RENDERER_ARGS=(
    --sample "\${SAMPLE_ID}"
    --condition "\${CONDITION}"
    --timepoint "\${TIMEPOINT}"
    --replicate "\${REPLICATE}"
    --divergent-bed "\${DIV_BED}"
    --functional-summary "\${FUNC_SUM}"
    --pol-density "\${POL_DENS}"
    --pausing-index "\${PAUSING_IDX}"
    --norm-factors "\${NORM_FACTORS}"
    --qc-json "\${QC_JSON}"
    --out-html "\${OUT_HTML}"
    --out-tsv "\${OUT_TSV}"
    --out-json "\${OUT_JSON}"
    --out-plots-html "\${OUT_PLOTS}"
    --pi-min-body-count "${params.pol?.pi_min_body_count ?: 10}"
  )

  echo "REPORT | BUILD | Base arguments: \${#RENDERER_ARGS[@]}"

  # Helper to add track link if available
  add_track_link() {
    local flag="\$1"
    local link="\$2"
    local have_flag="\$3"
    
    if [[ \${have_flag} -eq 1 ]]; then
      RENDERER_ARGS+=("\${flag}" "\${link}")
      echo "REPORT | BUILD | Added: \${flag}"
    fi
  }

  add_track_link "--allmap3p-pos-raw" "\${ALLMAP3P_POS_RAW}" \${HAVE_ALLMAP_POS_RAW}
  add_track_link "--allmap3p-neg-raw" "\${ALLMAP3P_NEG_RAW}" \${HAVE_ALLMAP_NEG_RAW}
  add_track_link "--pos3-cpm-bw" "\${POS3_CPM_BW}" \${HAVE_POS_CPM_BW}
  add_track_link "--neg3-cpm-bw" "\${NEG3_CPM_BW}" \${HAVE_NEG_CPM_BW}
  add_track_link "--allmap3p-pos-cpm-bw" "\${ALLMAP3P_POS_CPM_BW}" \${HAVE_ALLMAP_POS_CPM_BW}
  add_track_link "--allmap3p-neg-cpm-bw" "\${ALLMAP3P_NEG_CPM_BW}" \${HAVE_ALLMAP_NEG_CPM_BW}

  # Add plots flag if enabled
  if [[ \${ENABLE_PLOTS} -eq 1 ]]; then
    RENDERER_ARGS+=(--plots 1)
    echo "REPORT | BUILD | Plots enabled"
  fi

  echo "REPORT | BUILD | Total arguments: \${#RENDERER_ARGS[@]}"

  ###########################################################################
  # 6) RUN RENDERER
  ###########################################################################

  echo "REPORT | RENDER | Generating reports..."

  RENDER_START=\$(date +%s)

  set +e
  \${PYTHON_CMD} "\${RENDER_SCRIPT}" "\${RENDERER_ARGS[@]}"
  RENDER_RC=\$?
  set -e

  RENDER_END=\$(date +%s)
  RENDER_TIME=\$((RENDER_END - RENDER_START))

  echo "REPORT | RENDER | Rendering completed in \${RENDER_TIME}s"

  if [[ \${RENDER_RC} -ne 0 ]]; then
    tracktx_error "generate_per_sample_reports" "Renderer failed with exit code \${RENDER_RC}" "Check report.log in work dir" \${RENDER_RC}
  fi

  ###########################################################################
  # 7) VALIDATE OUTPUTS
  ###########################################################################

  echo "REPORT | VALIDATE | Checking output files..."

  # Check required outputs
  for OUTPUT in "\${OUT_HTML}" "\${OUT_TSV}" "\${OUT_JSON}"; do
    if [[ ! -s "\${OUTPUT}" ]]; then
      tracktx_error "generate_per_sample_reports" "Expected output missing or empty: \${OUTPUT}" "Check report.log in work dir"
    fi
    OUTPUT_SIZE=\$(tracktx_size "\${OUTPUT}")
    echo "REPORT | VALIDATE | \$(basename \${OUTPUT}): \${OUTPUT_SIZE} bytes"
  done

  # Check plots file
  if [[ \${ENABLE_PLOTS} -eq 1 ]]; then
    if [[ -s "\${OUT_PLOTS}" ]]; then
      PLOTS_SIZE=\$(tracktx_size "\${OUT_PLOTS}")
      echo "REPORT | VALIDATE | Plots HTML: \${PLOTS_SIZE} bytes"
    else
      echo "REPORT | WARNING | Plots enabled but file missing"
    fi
  else
    # Create placeholder if plots disabled
    if [[ ! -s "\${OUT_PLOTS}" ]]; then
      cat > "\${OUT_PLOTS}" <<'PLACEHOLDER'
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

  ###########################################################################
  # 8) CREATE README
  ###########################################################################

  echo "REPORT | README | Creating documentation..."

  cat > "\${OUT_README}" <<DOCEOF
PER-SAMPLE REPORT — ${sample_id}
────────────────────────────────────────────────────────────────────────────
  \${OUT_HTML}   — QC, divergent-TX, functional regions, Pol-II density,
                   pausing index, normalization factors, track links
  \${OUT_TSV}    — same metrics, one row each (metric_name, value, unit, category)
  \${OUT_JSON}   — same data, versioned schema (schema_version key)
  \${OUT_PLOTS}  — \$([ \${ENABLE_PLOTS} -eq 1 ] && echo "region pie chart + pausing-index histogram (inline base64 PNGs)" || echo "placeholder (plots disabled)")

  This sample: \${DIV_COUNT} divergent loci, \${PAUSING_GENES} genes with pausing
  index, \${TOTAL_TRACKS}/6 track links available. Note: total_reads/map_rate
  above are unique/MAPQ-pass stats -- the genuine overall alignment rate is in
  02_alignments/alignment_rates_summary.tsv, not in this report.
DOCEOF

  echo "REPORT | README | Documentation created"

  ###########################################################################
  # 9) CREATE QUICK ACCESS SYMLINK
  ###########################################################################

  echo "REPORT | SYMLINK | Creating quick access link..."

  SYMLINK_DIR="${params.output_dir}/reports"
  mkdir -p "\${SYMLINK_DIR}" || true

  # Create symlink for easy access
  ln -sf "\$(pwd)/\${OUT_HTML}" "\${SYMLINK_DIR}/\${SAMPLE_ID}.html" || \\
    echo "REPORT | WARNING | Could not create symlink"

  if [[ -L "\${SYMLINK_DIR}/\${SAMPLE_ID}.html" ]]; then
    echo "REPORT | SYMLINK | Quick access: \${SYMLINK_DIR}/\${SAMPLE_ID}.html"
  fi

  ###########################################################################
  # SUMMARY
  ###########################################################################

  echo "────────────────────────────────────────────────────────────────────────"
  echo "REPORT | SUMMARY | Sample: \${SAMPLE_ID}"
  echo "REPORT | SUMMARY | Divergent loci: \${DIV_COUNT}"
  echo "REPORT | SUMMARY | Genes with PI: \${PAUSING_GENES}"
  echo "REPORT | SUMMARY | Available tracks: \${TOTAL_TRACKS}/6"
  echo "REPORT | SUMMARY | Output files: 5"
  echo "REPORT | SUMMARY | Processing time: \${RENDER_TIME}s"
  echo "────────────────────────────────────────────────────────────────────────"

  TIMESTAMP_END=\$(date -u +"%Y-%m-%dT%H:%M:%SZ")
  echo "════════════════════════════════════════════════════════════════════════"
  echo "REPORT | COMPLETE | sample=\${SAMPLE_ID} | ts=\${TIMESTAMP_END}"
  echo "════════════════════════════════════════════════════════════════════════"
  """
}