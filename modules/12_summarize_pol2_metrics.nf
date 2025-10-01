// ============================================================================
// summarize_pol2_metrics.nf — merge per-sample gene metrics (+ optional contrasts/plots)
// ----------------------------------------------------------------------------
// Overview
//   • Input sheet: sample_id, condition, timepoint, replicate, file
//   • Loads minimal cols from each per-sample pol2_gene_metrics.tsv
//   • Produces merged tidy table; optional contrasts + plots
//
// Why it’s faster
//   • Reads only required columns (usecols)
//   • No HDF5; TSV only
//
// Outputs (publishDir):
//   ${params.output_dir}/09_pol2_aggregate/
//     ├── pol2_gene_metrics_merged.tsv
//     ├── pol2_gene_metrics_contrasts.tsv  (optional)
//     ├── plots/                           (optional)
//     ├── README_aggregate.txt
//     └── aggregate.log                    (unified step log)
// ============================================================================

nextflow.enable.dsl = 2

process summarize_pol2_metrics {

  // ── Meta / resources ─────────────────────────────────────────────────────
  tag        'pol2-aggregate'
  label      'conda'
  cache      'lenient'
  // Resource allocation handled dynamically by base.config

  // Ordered, discoverable output location
  publishDir "${params.output_dir}/09_pol2_aggregate", mode: 'copy', overwrite: true

  // Use the main tracktx conda environment
  conda (params.conda_pol2 ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    path(samples_tsv, stageAs: 'samples.tsv')  // fixed name in workdir

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    path 'pol2_gene_metrics_merged.tsv'
    path 'pol2_gene_metrics_contrasts.tsv', optional: true
    path 'plots', optional: true
    path 'README_aggregate.txt'
    path 'aggregate.log', emit: log

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # ── Unified log ──────────────────────────────────────────────────────────
  exec > >(tee -a aggregate.log) 2>&1
  trap 'echo "ERROR [aggregate] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  

  # ── Bindings / params ────────────────────────────────────────────────────
  SAMPLES_TSV='samples.tsv'
  # Scripts should be in PATH or conda environment
  TOPN='!{ params.pol2?.top_n ?: 100 }'
  PLOTS='!{ (params.pol2?.plots == null) ? true : (params.pol2.plots as boolean) }'

  echo "[aggregate] ▶ samples=${SAMPLES_TSV}  plots=${PLOTS}  top_n=${TOPN} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

  command -v python3 >/dev/null 2>&1 || { echo "[aggregate] ERROR: python3 not found" >&2; exit 1; }

  # contrasts (one per line) from params.pol2.contrasts (if provided)
  cat > contrasts.txt <<'TXT'
!{ ((params.pol2?.contrasts ?: []) as List).join('\n') }
TXT

  # ── Quick schema check for samples.tsv ───────────────────────────────────
  awk 'BEGIN{FS="\\t"}
       NR==1{
         want="sample_id\\tcondition\\ttimepoint\\treplicate\\tfile";
         if(\$0!=want){print "[aggregate] ERROR: bad header. want: " want " got: " \$0 >"/dev/stderr"; exit 2}
         next
       }
       { if(NF!=5){print "[aggregate] ERROR: row with NF!=5 :: " \$0 >"/dev/stderr"; exit 2} }' "${SAMPLES_TSV}"

  # ── Run aggregator ───────────────────────────────────────────────────────
  mkdir -p plots
  ARGS=( --samples-tsv "\${SAMPLES_TSV}"
         --out-merged  pol2_gene_metrics_merged.tsv
         --top-n       "\${TOPN}" )
  if [ -s contrasts.txt ]; then
    mapfile -t C < contrasts.txt || true
    (( \${#C[@]} )) && ARGS+=( --contrasts "\${C[@]}" --out-contrasts pol2_gene_metrics_contrasts.tsv )
  fi
  if [ "\${PLOTS}" = "true" ]; then ARGS+=( --plots-dir plots ); fi

  CMP_PY="!{projectDir}/bin/compare_pol2_metrics.py"
  python3 "${CMP_PY}" "\${ARGS[@]}"

  # ── Post checks / tidy ───────────────────────────────────────────────────
  [ -s pol2_gene_metrics_merged.tsv ] || { echo "[aggregate] ERROR: merged TSV missing" >&2; exit 5; }
  if [ -d plots ] && ! ls -1 plots/* >/dev/null 2>&1; then rmdir plots || true; fi
  if [ -s contrasts.txt ] && [ ! -s pol2_gene_metrics_contrasts.tsv ]; then rm -f pol2_gene_metrics_contrasts.tsv || true; fi

  # ── README ───────────────────────────────────────────────────────────────
  cat > README_aggregate.txt <<'TXT'
TrackTx — Pol II aggregate metrics

Aggregates
  • Per-sample pol2_gene_metrics.tsv files → merged tidy table.
  • Optional contrasts (e.g., condition:KI,WT or timepoint:90,0) with log2FC.
  • Optional plots: top-variable heatmaps (per metric x group) and MA plots.

Notes
  • Density over regions is computed upstream and is not aggregated here.
  • Plots are small PNGs (Matplotlib); enable via params.pol2.plots=true.
  • Required per-sample columns: gene_id, gene_name, pi_len_norm, pi_raw, body_cpm, tss_cpm.
TXT

  echo "INFO  [aggregate] ✔ done ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
