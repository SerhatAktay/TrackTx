// ============================================================================
// combine_reports.nf — Cohort aggregation of per-sample summaries (SPA HTML)
// ----------------------------------------------------------------------------
// Overview
//   • Accepts files and/or directories containing per-sample JSON summaries
//   • Robust intake (case-insensitive): *.summary.json / *.report.json / *.json / *.json.gz
//   • Produces a polished single-file interactive HTML + TSV + JSON
//   • Optional global functional-region totals table
//
// Inputs
//   path input_files   // one path or a list; may be files and/or directories
//
// Outputs (publishDir):
//   ${params.output_dir}/11_reports/cohort/
//     ├── global_summary.html
//     ├── global_summary.tsv
//     ├── global_summary.json
//     ├── global_region_totals.tsv           (optional)
//     ├── README.txt
//     └── combine.log
// ============================================================================

nextflow.enable.dsl = 2

process combine_reports {

  // ── Meta / resources ─────────────────────────────────────────────────────
  tag        { (input_files instanceof List) ? "n=${input_files.size()}" : input_files.getBaseName() }
  label      'conda'
  cache      'lenient'
  // Resource allocation handled dynamically by base.config
  publishDir "${params.output_dir}/11_reports/cohort", mode: 'copy', overwrite: true

  // Use the main tracktx conda environment
  conda (params.conda_pol2 ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ───────────────────────────────────────────────────────────────
  input:
    path input_files

  // ── Declared outputs ────────────────────────────────────────────────────
  output:
    path "global_summary.html"
    path "global_summary.tsv"
    path "global_summary.json"
    path "global_region_totals.tsv", optional: true
    path "README.txt"
    path "combine.log", emit: log

  shell:
  '''
  #!/usr/bin/env bash
  set -euo pipefail
  export LC_ALL=C
  # Unified log
  exec > >(tee -a combine.log) 2>&1
  trap 'echo "ERROR [combine] failed at $(date -u +"%Y-%m-%dT%H:%M:%SZ")" >&2' ERR
  

  echo "INFO  [combine] received paths: !{ (input_files instanceof List) ? input_files.size() : 1 }"

  # ------------------------------------------------------------------------
  # Build a bash array of candidate inputs (safe for spaces/newlines)
  # ------------------------------------------------------------------------
  mapfile -t INFILES < <(printf '%s\n' !{ (input_files instanceof List) ? input_files.collect{ '"' + it.toString().replace('"','\\"') + '"' }.join(' ') : '"' + input_files.toString().replace('"','\\"') + '"' })

  # ------------------------------------------------------------------------
  # Keep only *.summary.json or *.report.json (case-insensitive)
  # ------------------------------------------------------------------------
  JSONS=()
  for f in "${INFILES[@]}"; do
    lf="${f,,}"
    if [[ "$lf" == *.summary.json || "$lf" == *.report.json ]]; then
      JSONS+=( "$f" )
    fi
  done

  echo "INFO  [combine] json candidates: ${#JSONS[@]}"
  if [[ ${#JSONS[@]} -eq 0 ]]; then
    echo "ERROR no *.summary.json/*.report.json inputs. Upstream should collect per-sample JSONs." >&2
    exit 2
  fi

  # ------------------------------------------------------------------------
  # Validate inputs (exist + non-empty) and log sizes
  # ------------------------------------------------------------------------
  missing=0
  for f in "${JSONS[@]}"; do
    if [[ -e "$f" ]]; then
      sz=$(stat -c%s "$f" 2>/dev/null || stat -f%z "$f" || echo 0)
      echo "INFO  [combine] file: $f  size=$sz"
      [[ "$sz" -gt 0 ]] || { echo "ERROR empty file: $f" >&2; missing=1; }
    else
      echo "ERROR missing file: $f" >&2; missing=1
    fi
  done
  [[ $missing -eq 0 ]] || exit 2

  # ------------------------------------------------------------------------
  # Run the combiner
  # ------------------------------------------------------------------------
  command -v python3 >/dev/null 2>&1 || { echo "ERROR python3 not found" >&2; exit 3; }

  COMBINE_PY="!{projectDir}/bin/combine_reports.py"
  [[ -e "${COMBINE_PY}" ]] || { echo "ERROR combiner not found: ${COMBINE_PY}" >&2; exit 3; }

  python3 "${COMBINE_PY}" \
    --inputs "${JSONS[@]}" \
    --out-tsv   global_summary.tsv \
    --out-json  global_summary.json \
    --out-html  global_summary.html \
    --out-regions  global_region_totals.tsv \
    --pipeline-version "!{workflow.manifest.version ?: 'dev'}" \
    --run-name "!{workflow.runName ?: 'unnamed'}" \
    --duration "!{workflow.duration ?: 'unknown'}" \
    --profile "!{workflow.profile ?: 'unknown'}" || true

  # ------------------------------------------------------------------------
  # Validate outputs
  # ------------------------------------------------------------------------
  [[ -s global_summary.tsv  ]]  || { echo "ERROR cohort TSV missing"  >&2; exit 2; }
  [[ -s global_summary.json ]]  || { echo "ERROR cohort JSON missing" >&2; exit 2; }
  [[ -s global_summary.html ]]  || { echo "ERROR cohort HTML missing" >&2; exit 2; }

  # ------------------------------------------------------------------------
  # Convenience copy: place HTML at results root for quick access
  # ------------------------------------------------------------------------
  ROOT_HTML="!{params.output_dir}/global_summary.html"
  mkdir -p "!{params.output_dir}"
  cp -f global_summary.html "${ROOT_HTML}"
  echo "INFO  [combine] copied cohort HTML to ${ROOT_HTML}"

  # ------------------------------------------------------------------------
  # README
  # ------------------------------------------------------------------------
  cat > README.txt <<'TXT'
TrackTx — Cohort summary

Inputs
  • <sid>.summary.json or <sid>.report.json (one per sample; upstream should collect the list)

Outputs
  • global_summary.tsv
  • global_summary.json
  • global_summary.html
  • global_region_totals.tsv (optional; summed functional-region counts)

Notes
  - Robust to missing fields in per-sample JSON.
  - HTML is a single-file SPA (Summary + per-sample pages).
TXT

  echo "INFO  [combine] ✔ done"
  
  # ------------------------------------------------------------------------
  # Landing page: create results/index.html with quick links
  # ------------------------------------------------------------------------
  ROOT_DIR="!{params.output_dir}"
  INDEX_HTML="${ROOT_DIR}/index.html"
  COHORT_REL="11_reports/cohort/global_summary.html"
  NF_RPT_REL="nf_report.html"
  NF_TIM_REL="nf_timeline.html"
  NF_DAG_REL="nf_dag.png"
  mkdir -p "${ROOT_DIR}"
  cat > "${INDEX_HTML}" <<HTML
<!doctype html>
<meta charset="utf-8">
<title>TrackTx Results</title>
<h1>TrackTx Results</h1>
<ul>
  <li><a href="${COHORT_REL}">Cohort summary</a></li>
  <li><a href="reports/">Per-sample reports</a></li>
  <li><a href="${NF_RPT_REL}">Nextflow report</a></li>
  <li><a href="${NF_TIM_REL}">Nextflow timeline</a></li>
  <li><a href="${NF_DAG_REL}">Nextflow DAG</a></li>
  <li><a href="trace/trace.txt">Nextflow trace</a></li>
 </ul>
HTML
  echo "INFO  [combine] wrote landing page: ${INDEX_HTML}"
  '''
}
