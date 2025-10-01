// ============================================================================
// detect_divergent_tx.nf — Divergent transcription from RAW allMap 3′ (±)
// ----------------------------------------------------------------------------
// Overview
//   • Inputs: unnormalized allMap 3′ bedGraphs (POS ≥0, NEG ≤0 mirrored upstream)
//   • Runs detector on sanitized, strand-correct signals (optionally |NEG|)
//   • Optional QC PDF, optional gzip; emits per-sample summary TSV + README
//
// Why it’s faster
//   • Lets the Python tool handle I/O + sanitization in one pass
//   • Skips resorting when upstream is already sorted (--assume-sorted)
//   • Multicore pairing per chromosome (ncores <= cpus)
//
// Inputs
//   tuple( sample_id, allmap3p_pos_bg, allmap3p_neg_bg, condition, timepoint, replicate )
//
// Outputs (publishDir):
//   ${params.output_dir}/06_divergent_tx/${sample_id}/
//     ├── divergent_transcription.bed
//     ├── divergent_transcription.bed.gz            (if gzip enabled)
//     ├── divergent_transcription.qc.pdf            (if qc enabled)
//     ├── divergent_summary.tsv
//     ├── README_divergent.txt
//     └── divergent.log
//
// Feature knobs (params.advanced.* with defaults):
//   divergent_nt_window   = 800
//   divergent_threshold   = 1.0
//   divergent_signal_sum  = 5.0
//   divergent_bin_gap     = 100
//   divergent_merge_gap   = 50
//   divergent_balance     = 0.0
//   divergent_overlap_bp  = 0
//   divergent_valley_thr  = "Inf"
//   divergent_max_width   = 0
//   divergent_smooth_sd   = 0
//   divergent_qc          = false
//   divergent_gzip        = false
//   divergent_abs_neg     = false      // pass |NEG| to detector
//   divergent_force_sort  = false      // if false → --assume-sorted
// ============================================================================

nextflow.enable.dsl = 2

process detect_divergent_tx {

  // ── Meta / resources ─────────────────────────────────────────────────────
  // Resource allocation handled dynamically by base.config
  tag      { sample_id }
  label    'conda'
  cache    'lenient'

  // Ordered, discoverable output location
  publishDir "${params.output_dir}/06_divergent_tx/${sample_id}", mode: 'copy', overwrite: true

  // Use the main tracktx conda environment
  conda (params.conda_divergent ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(allmap3p_pos_bg), path(allmap3p_neg_bg),
          val(condition), val(timepoint), val(replicate)

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    tuple val(sample_id),
          path("divergent_transcription.bed"),
          val(condition), val(timepoint), val(replicate)

    path "divergent_transcription.bed.gz", optional: true, emit: bed_gz
    path "divergent_transcription.qc.pdf", optional: true, emit: qc_pdf
    path "divergent_summary.tsv",                          emit: summary
    path "README_divergent.txt",                           emit: readme
    path "divergent.log",                                  emit: log

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C

  # ── Unified log ──────────────────────────────────────────────────────────
  exec > >(tee -a divergent.log) 2>&1
  trap 'echo "ERROR  [divergent] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  

  # ── Bindings / locals ────────────────────────────────────────────────────
  SID="!{sample_id}"
  POS_BG="!{allmap3p_pos_bg}"
  NEG_BG="!{allmap3p_neg_bg}"
  THREADS=!{task.cpus}
  DETPY="!{projectDir}/bin/detect_divergent_tx.py"

  # Tunables
  NT_WINDOW=!{ params.advanced?.divergent_nt_window   ?: 800 }
  THR=!{       params.advanced?.divergent_threshold   ?: 1.0 }
  SUMTHR=!{    params.advanced?.divergent_signal_sum  ?: 5.0 }
  BIN_GAP=!{   params.advanced?.divergent_bin_gap     ?: 100 }
  MERGE_GAP=!{ params.advanced?.divergent_merge_gap   ?: 50 }
  BAL=!{       params.advanced?.divergent_balance     ?: 0.0 }
  OVERLAP_BP=!{params.advanced?.divergent_overlap_bp  ?: 0 }
  VALLEY_THR="!{ params.advanced?.divergent_valley_thr ?: 'Inf' }"
  SMOOTH_SD=!{ params.advanced?.divergent_smooth_sd   ?: 0 }
  MAX_DT_WIDTH=!{ params.advanced?.divergent_max_width ?: 0 }
  ABS_NEG=$([ '!{ params.advanced?.divergent_abs_neg    in [true, "true"] }' = 'true' ] && echo 1 || echo 0)
  FORCE_SORT=$([ '!{ params.advanced?.divergent_force_sort in [true, "true"] }' = 'true' ] && echo 1 || echo 0)
  DO_QC=$([ '!{ params.advanced?.divergent_qc          in [true, "true"] }' = 'true' ] && echo 1 || echo 0)
  DO_GZIP=$([ '!{ params.advanced?.divergent_gzip       in [true, "true"] }' = 'true' ] && echo 1 || echo 0)

  echo "INFO  [divergent] ▶ sample=\${SID}  cpus=\${THREADS} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo "INFO  files: pos=\${POS_BG}  neg=\${NEG_BG}  abs_neg=\${ABS_NEG} force_sort=\${FORCE_SORT} qc=\${DO_QC} gzip=\${DO_GZIP}"
  echo "INFO  params: thr=\${THR} sum_thr=\${SUMTHR} win=\${NT_WINDOW} bin_gap=\${BIN_GAP} merge_gap=\${MERGE_GAP} overlap_bp=\${OVERLAP_BP} valley_thr=\${VALLEY_THR} bal=\${BAL} smooth_sd=\${SMOOTH_SD} max_width=\${MAX_DT_WIDTH}"

  command -v python3 >/dev/null 2>&1 || { echo "ERROR python3 not found" >&2; exit 1; }

  # ─────────────────────────────────────────────────────────────────────────
  # 1) Detector
  # ─────────────────────────────────────────────────────────────────────────
  ASSUME_SORT_ARG=""
  [[ \${FORCE_SORT} -eq 1 ]] || ASSUME_SORT_ARG="--assume-sorted"

  NEG_ABS_ARG=""
  [[ \${ABS_NEG} -eq 1 ]] && NEG_ABS_ARG="--neg-abs"

  QC_ARG=""
  [[ \${DO_QC} -eq 1 ]] && QC_ARG="--qc 1"

  set +e
  python3 "\${DETPY}" \\
    --sample       "\${SID}" \\
    --pos          "\${POS_BG}" \\
    --neg          "\${NEG_BG}" \\
    --nt-window    "\${NT_WINDOW}" \\
    --threshold    "\${THR}" \\
    --sum-thr      "\${SUMTHR}" \\
    --bin-gap      "\${BIN_GAP}" \\
    --merge-gap    "\${MERGE_GAP}" \\
    --balance      "\${BAL}" \\
    --overlap-bp   "\${OVERLAP_BP}" \\
    --valley-thr   "\${VALLEY_THR}" \\
    --smooth-sd    "\${SMOOTH_SD}" \\
    --max-dt-width "\${MAX_DT_WIDTH}" \\
    --ncores       "\${THREADS}" \\
    --out-bed      "divergent_transcription.bed" \\
    --write-summary "divergent_summary.tsv" \\
    \${ASSUME_SORT_ARG} \\
    \${NEG_ABS_ARG} \\
    \${QC_ARG} 2>&1 | tee -a divergent.log
  rc=\${PIPESTATUS[0]}
  set -e

  if [[ \${rc} -ne 0 ]]; then
    echo "WARN  detector exited with code \${rc}; emitting empty BED" | tee -a divergent.log
    : > divergent_transcription.bed
    : > divergent_summary.tsv
    rm -f divergent_transcription.qc.pdf || true
  fi

  # Ensure files exist
  [[ -s divergent_transcription.bed       ]] || : > divergent_transcription.bed
  [[ -s divergent_summary.tsv             ]] || echo -e "sample\tn_pos_pk\tn_neg_pk\tn_pairs_raw\tn_dt\twall_s" > divergent_summary.tsv

  # ─────────────────────────────────────────────────────────────────────────
  # 2) Optional compression
  # ─────────────────────────────────────────────────────────────────────────
  if [[ \${DO_GZIP} -eq 1 ]]; then
    if command -v bgzip >/dev/null 2>&1; then
      echo "INFO  compressing with bgzip"
      bgzip -@ "\${THREADS}" -fc divergent_transcription.bed > divergent_transcription.bed.gz
    elif command -v pigz >/dev/null 2>&1; then
      echo "INFO  compressing with pigz"
      pigz -p "\${THREADS}" -fc divergent_transcription.bed > divergent_transcription.bed.gz
    else
      echo "INFO  compressing with gzip"
      gzip -fc divergent_transcription.bed > divergent_transcription.bed.gz
    fi
  fi

  # ─────────────────────────────────────────────────────────────────────────
  # 3) README
  # ─────────────────────────────────────────────────────────────────────────
  cat > README_divergent.txt <<TXT
Divergent transcription (sample: \${SID})

Inputs (RAW allMap 3′):
  - POS bedGraph (≥0)
  - NEG bedGraph (≤0 mirrored upstream; abs_neg=\${ABS_NEG})
Detector:
  - Window=\${NT_WINDOW}  thr=\${THR}  sum_thr=\${SUMTHR}  bin_gap=\${BIN_GAP}  merge_gap=\${MERGE_GAP}
  - balance=\${BAL}  overlap_bp=\${OVERLAP_BP}  valley_thr=\${VALLEY_THR}
  - smooth_sd=\${SMOOTH_SD}  max_dt_width=\${MAX_DT_WIDTH}
  - ncores=\${THREADS}  assume_sorted=$([[ \${FORCE_SORT} -eq 1 ]] && echo false || echo true)

Files:
  - divergent_transcription.bed  (+ .bed.gz if enabled)
  - divergent_transcription.qc.pdf (if QC enabled)
  - divergent_summary.tsv (counts + wall time)
  - divergent.log

Notes:
  - If abs_neg=1, NEG magnitudes are used (|NEG|); otherwise NEG (≤0) is used.
  - Upstream sorting is trusted unless divergent_force_sort=true.
TXT

  echo "INFO  [divergent] ✔ done (\${SID}) ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  '''
}
