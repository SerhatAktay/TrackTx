// ============================================================================
// 14_generate_reports.nf — Per-sample HTML / TSV / JSON (+ optional track links)
// ----------------------------------------------------------------------------
// Overview
//   • Summarizes divergent loci, functional-region composition, Pol II density,
//     pausing metrics, normalization factors, and QC into one report.
//   • Emits a stable JSON (with schema_version), a concise TSV, and an HTML page.
//   • Optional plots page (base64 inline); disabled by default.
//   • Track links can be local paths or http(s) URLs (ignored if missing).
//
// Inputs
//   tuple(
//     sample_id,
//     divergent_transcription.bed,
//     functional_regions_summary.tsv,
//     pol2_density.tsv,
//     pausing_index.tsv,
//     normalization_factors.tsv,
//     dedup_stats        // not parsed; kept for interface parity
//     qc_pol2.json,
//     allmap3p_pos_raw, allmap3p_neg_raw,    // strings; may be ""
//     pos3_cpm_bw, neg3_cpm_bw,              // strings; may be ""
//     allmap3p_pos_cpm_bw, allmap3p_neg_cpm_bw, // strings; may be ""
//     condition, timepoint, replicate
//   )
//
// Outputs (publishDir):
//   ${params.output_dir}/11_reports/samples/${sample_id}/
//     ├── ${sample_id}.report.html
//     ├── ${sample_id}.report.tsv
//     ├── ${sample_id}.report.json
//     ├── ${sample_id}.plots.html                  (placeholder if disabled)
//     ├── ${sample_id}.report.log
//     └── ${sample_id}.README_report.txt
// ============================================================================

nextflow.enable.dsl = 2

process generate_reports {

  // ── Meta / resources ─────────────────────────────────────────────────────
  tag        { sample_id }
  label      'conda'
  cache      'lenient'
  // Resource allocation handled dynamically by base.config

  // Use the main tracktx conda environment
  conda 'envs/tracktx.yaml'

  // Ordered, discoverable output location
  publishDir "${params.output_dir}/11_reports/samples/${sample_id}", mode: 'copy', overwrite: true


  // ── Inputs ──────────────────────────────────────────────────────────────
  input:
    tuple val(sample_id),
          path(div_bed),
          path(func_sum),
          path(pol2_density),
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

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    path("${sample_id}.report.html"), emit: html
    path("${sample_id}.report.tsv"),  emit: tsv
    path("${sample_id}.report.json"), emit: json
    path("${sample_id}.plots.html"),  optional: true, emit: plots
    path("${sample_id}.README_report.txt")
    path("${sample_id}.report.log"), emit: log

  // ── Script ──────────────────────────────────────────────────────────────
  shell:
  '''
  #!/usr/bin/env bash
  # NOTE: no `-u` because optional links may be empty strings.
  set -eo pipefail
  export LC_ALL=C

  # ── Unified log ──────────────────────────────────────────────────────────
  exec > >(tee -a "!{sample_id}.report.log") 2>&1
  trap 'echo "ERROR [report] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  

  # ── Bindings / locals ────────────────────────────────────────────────────
  SID="!{sample_id}"
  COND="!{condition}"; TP="!{timepoint}"; REP="!{replicate}"

  BED="!{div_bed}"
  FSUM="!{func_sum}"
  DENS="!{pol2_density}"
  PAUS="!{pausing_idx}"
  NORM="!{norm_factors}"
  QCJ="!{qc_json}"

  RAW_AP3="!{allmap3p_pos_raw}"
  RAW_AN3="!{allmap3p_neg_raw}"
  BW_P3="!{pos3_cpm_bw}"
  BW_N3="!{neg3_cpm_bw}"
  BW_AP3="!{allmap3p_pos_cpm_bw}"
  BW_AN3="!{allmap3p_neg_cpm_bw}"

  # Container-safe bin script location
  RENDER="!{projectDir}/bin/render_sample_report.py"

  PLOTS_FLAG=!{ (params.reports_plots == null) ? 0 : (params.reports_plots as int) }

  OUT_HTML="${SID}.report.html"
  OUT_TSV="${SID}.report.tsv"
  OUT_JSON="${SID}.report.json"
  OUT_PLOTS="${SID}.plots.html"

  echo "INFO  [report] ▶ ${SID}  plots=${PLOTS_FLAG} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

  # ── Sanity (core) ────────────────────────────────────────────────────────
  for f in "${BED}" "${FSUM}" "${DENS}" "${PAUS}" "${NORM}" "${QCJ}"; do
    [[ -s "$f" ]] || { echo "ERROR missing/empty: $f" >&2; exit 2; }
  done
  command -v python3 >/dev/null 2>&1 || { echo "ERROR python3 not found" >&2; exit 2; }
  [[ -e "${RENDER}" ]] || { echo "ERROR script not found: ${RENDER}" >&2; exit 2; }

  # ── Build CLI (links: allow URLs or existing files) ──────────────────────
  ARGS=(
    --sample             "${SID}"
    --condition          "${COND}"
    --timepoint          "${TP}"
    --replicate          "${REP}"
    --divergent-bed      "${BED}"
    --functional-summary "${FSUM}"
    --pol2-density       "${DENS}"
    --pausing-index      "${PAUS}"
    --norm-factors       "${NORM}"
    --qc-json            "${QCJ}"
    --out-html           "${OUT_HTML}"
    --out-tsv            "${OUT_TSV}"
    --out-json           "${OUT_JSON}"
    --out-plots-html     "${OUT_PLOTS}"
  )

  add_link () {
    local flag="$1"; local s="$2"
    [[ -z "$s" ]] && { echo "DEBUG skip ${flag} (empty)"; return; }
    if [[ "$s" =~ ^https?:// ]]; then
      ARGS+=( "${flag}" "${s}" ); echo "DEBUG link URL ${flag}=${s}"
    elif [[ -e "$s" ]]; then
      ARGS+=( "${flag}" "${s}" ); echo "DEBUG link file ${flag}=${s}"
    else
      echo "DEBUG skip ${flag} (not URL, missing file: ${s})"
    fi
  }
  add_link --allmap3p-pos-raw     "${RAW_AP3}"
  add_link --allmap3p-neg-raw     "${RAW_AN3}"
  add_link --pos3-cpm-bw          "${BW_P3}"
  add_link --neg3-cpm-bw          "${BW_N3}"
  add_link --allmap3p-pos-cpm-bw  "${BW_AP3}"
  add_link --allmap3p-neg-cpm-bw  "${BW_AN3}"

  # Compact summary of which links were included
  inc=( )
  [[ -n "${RAW_AP3}" && ( "${RAW_AP3}" =~ ^https?:// || -e "${RAW_AP3}" ) ]] && inc+=(raw_ap3)
  [[ -n "${RAW_AN3}" && ( "${RAW_AN3}" =~ ^https?:// || -e "${RAW_AN3}" ) ]] && inc+=(raw_an3)
  [[ -n "${BW_P3}"  && ( "${BW_P3}"  =~ ^https?:// || -e "${BW_P3}"  ) ]] && inc+=(bw_pos3)
  [[ -n "${BW_N3}"  && ( "${BW_N3}"  =~ ^https?:// || -e "${BW_N3}"  ) ]] && inc+=(bw_neg3)
  [[ -n "${BW_AP3}" && ( "${BW_AP3}" =~ ^https?:// || -e "${BW_AP3}" ) ]] && inc+=(bw_allMap_pos3)
  [[ -n "${BW_AN3}" && ( "${BW_AN3}" =~ ^https?:// || -e "${BW_AN3}" ) ]] && inc+=(bw_allMap_neg3)
  echo "INFO  [report] links_included: ${inc[*]:-none}"

  [[ "${PLOTS_FLAG}" == "1" ]] && ARGS+=( --plots 1 )

  # ── Render ───────────────────────────────────────────────────────────────
  python3 "${RENDER}" "${ARGS[@]}"

  # Ensure expected outputs exist (renderer is defensive, but double-check)
  for o in "${OUT_HTML}" "${OUT_TSV}" "${OUT_JSON}"; do
    [[ -s "$o" ]] || { echo "ERROR expected output missing/empty: $o" >&2; exit 3; }
  done
  if [[ "${PLOTS_FLAG}" != "1" && ! -s "${OUT_PLOTS}" ]]; then
    printf "<!doctype html><meta charset='utf-8'><title>Plots disabled</title><p>Set --reports_plots 1 to enable.</p>\n" > "${OUT_PLOTS}"
  fi

  # ── README ───────────────────────────────────────────────────────────────
  cat > "${SID}.README_report.txt" <<TXT
TrackTx — Per-sample report

Sample    : ${SID}
Condition : ${COND}
Timepoint : ${TP}
Replicate : ${REP}

Outputs
  • ${OUT_HTML}
  • ${OUT_TSV}
  • ${OUT_JSON}
  • ${OUT_PLOTS}
  • ${SID}.report.log

Core inputs
  • divergent_bed              = ${BED}
  • functional_regions_summary = ${FSUM}
  • pol2_density               = ${DENS}
  • pausing_index              = ${PAUS}
  • normalization_factors      = ${NORM}
  • qc_pol2.json               = ${QCJ}

Optional track links
  • allMap3p_pos_raw.bg        = ${RAW_AP3:-NA}
  • allMap3p_neg_raw.bg        = ${RAW_AN3:-NA}
  • pos3_cpm.bw                = ${BW_P3:-NA}
  • neg3_cpm.bw                = ${BW_N3:-NA}
  • allMap3p_pos_cpm.bw        = ${BW_AP3:-NA}
  • allMap3p_neg_cpm.bw        = ${BW_AN3:-NA}
TXT

  echo "INFO  [report] ✔ ${SID} ts=$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  
  # Convenience: create a flat symlink for quick access under results/reports/
  SYM_DIR="!{params.output_dir}/reports"
  mkdir -p "${SYM_DIR}"
  ln -sf "$(pwd)/${OUT_HTML}" "${SYM_DIR}/${SID}.html"
  '''
}
