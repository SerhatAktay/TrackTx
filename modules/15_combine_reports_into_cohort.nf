// ============================================================================
// combine_reports_into_cohort.nf — Cohort-Level Report Aggregation
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


process combine_reports_into_cohort {

  tag        { "cohort" }
  label      'conda'
  cache      'deep'

  publishDir "${params.output_dir}/11_reports/cohort",
             mode: params.publish_mode,
             overwrite: true

  conda (params.conda_pol ?: "${projectDir}/envs/tracktx.yaml")

  // ── Inputs ────────────────────────────────────────────────────────────────
  // Stage files with sequential names to avoid collisions (like metric_1, metric_2, etc.)
  input:
    path 'report_*.json'
    path concordance_tsv
    // Module 16 outputs — staged here so the landing page can embed/link them.
    // Files named "NO_FILE" are sentinel placeholders for optional outputs.
    path qc_multiqc_html
    path qc_igv_session
    path qc_runon_tsv
    path qc_pca_plot
    path qc_corr_heatmap

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

  # Stdout/stderr → log + terminal (kept separate for Nextflow "Command error")
  exec > >(tee -a combine.log)
  exec 2> >(tee -a combine.log >&2)

  tracktx_error() {
    local module="\$1" problem="\$2" fix="\$3" code="\${4:-1}"
    echo "" >&2
    echo "═══════════════════════════════════════════════════════════════════════" >&2
    echo "TRACKTX ERROR" >&2
    echo "═══════════════════════════════════════════════════════════════════════" >&2
    echo "Module:  \${module}" >&2
    echo "Problem: \${problem}" >&2
    echo "Fix:     \${fix}" >&2
    echo "═══════════════════════════════════════════════════════════════════════" >&2
    exit "\$code"
  }
  trap 'tracktx_error "combine_reports_into_cohort" "Unexpected process failure" "Check combine.log in work dir"' ERR

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
    tracktx_error "combine_reports_into_cohort" "No JSON files found (expected: *.summary.json, *.report.json, or *.json)" "Check generate_per_sample_reports module outputs"
  fi

  ###########################################################################
  # 3) VALIDATE INPUT FILES
  ###########################################################################

  echo "COHORT | VALIDATE | Checking input files..."

  TOTAL_SIZE=0

  for JSON_FILE in "${JSON_FILES[@]}"; do
    if [[ ! -e "${JSON_FILE}" ]]; then
      tracktx_error "combine_reports_into_cohort" "Missing file: ${JSON_FILE}" "Check generate_per_sample_reports module outputs"
    fi
    
    FILE_SIZE=$(stat -c%s "${JSON_FILE}" 2>/dev/null || stat -f%z "${JSON_FILE}" 2>/dev/null || echo 0)
    
    if [[ ${FILE_SIZE} -eq 0 ]]; then
      tracktx_error "combine_reports_into_cohort" "Empty file: ${JSON_FILE}" "Check generate_per_sample_reports module outputs"
    fi
    
    TOTAL_SIZE=$((TOTAL_SIZE + FILE_SIZE))
    echo "COHORT | VALIDATE | $(basename ${JSON_FILE}): ${FILE_SIZE} bytes"
  done

  echo "COHORT | VALIDATE | Total input size: ${TOTAL_SIZE} bytes"

  ###########################################################################
  # 4) VALIDATE TOOLS
  ###########################################################################

  echo "COHORT | VALIDATE | Checking required tools..."

  # Use micromamba run to ensure correct Python env when in container (Docker/Singularity)
  if command -v micromamba >/dev/null 2>&1; then
    PYTHON_CMD="micromamba run -n base python3"
  elif [[ -x /opt/conda/bin/python3 ]]; then
    PYTHON_CMD="/opt/conda/bin/python3"
  else
    PYTHON_CMD="python3"
  fi

  # Check combiner script
  if [[ ! -e "${COMBINER_SCRIPT}" ]]; then
    tracktx_error "combine_reports_into_cohort" "Combiner script not found: ${COMBINER_SCRIPT}" "Ensure bin/combine_reports.py exists"
  fi
  echo "COHORT | VALIDATE | Combiner script: ${COMBINER_SCRIPT}"

  # Check Python
  if ${PYTHON_CMD} --version >/dev/null 2>&1; then
    PYTHON_VERSION=$(${PYTHON_CMD} --version 2>&1 || echo "unknown")
    echo "COHORT | VALIDATE | Python: ${PYTHON_VERSION}"
  else
    tracktx_error "combine_reports_into_cohort" "Python not found (tried: ${PYTHON_CMD})" "Use -profile docker"
  fi

  ###########################################################################
  # 5) RUN REPORT COMBINER
  ###########################################################################

  echo "COHORT | COMBINE | Aggregating per-sample reports..."
  echo "COHORT | COMBINE | Processing ${TOTAL_JSON} samples..."

  COMBINE_START=$(date +%s)

  set +e
  ${PYTHON_CMD} "${COMBINER_SCRIPT}" \
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
    tracktx_error "combine_reports_into_cohort" "Combiner failed with exit code ${COMBINE_RC}" "Check combine.log in work dir" ${COMBINE_RC}
  fi

  ###########################################################################
  # 6) VALIDATE OUTPUTS
  ###########################################################################

  echo "COHORT | VALIDATE | Checking output files..."

  # Check required outputs
  for OUTPUT in "${OUT_HTML}" "${OUT_TSV}" "${OUT_JSON}"; do
    if [[ ! -s "${OUTPUT}" ]]; then
      tracktx_error "combine_reports_into_cohort" "Missing or empty output: ${OUTPUT}" "Check combine.log in work dir"
    fi
    OUTPUT_SIZE=$(stat -c%s "${OUTPUT}" 2>/dev/null || stat -f%z "${OUTPUT}" 2>/dev/null || echo "unknown")
    echo "COHORT | VALIDATE | $(basename ${OUTPUT}): ${OUTPUT_SIZE} bytes"
  done

  # Check optional regions file
  if [[ -s "${OUT_REGIONS}" ]]; then
    REGIONS_SIZE=$(stat -c%s "${OUT_REGIONS}" 2>/dev/null || stat -f%z "${OUT_REGIONS}" 2>/dev/null || echo "unknown")
    REGIONS_LINES=$(wc -l < "${OUT_REGIONS}" | tr -d ' ')
    echo "COHORT | VALIDATE | Region totals: ${REGIONS_SIZE} bytes (${REGIONS_LINES} lines)"
  else
    echo "COHORT | VALIDATE | Region totals: not generated"
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
  # 9) CREATE LANDING PAGE  (Python-generated; embeds run-on data + KPIs)
  ###########################################################################

  echo "COHORT | LANDING | Generating landing page..."

  ROOT_DIR="!{params.output_dir}"
  INDEX_HTML="${ROOT_DIR}/index.html"

  # Resolve staged module-16 file paths (sentinel "NO_FILE" = optional/absent)
  QC_MULTIQC="!{qc_multiqc_html}"
  QC_IGV="!{qc_igv_session}"
  QC_RUNON="!{qc_runon_tsv}"
  QC_PCA="!{qc_pca_plot}"
  QC_CORR="!{qc_corr_heatmap}"

  # Write the generator to a temp file (avoids nested-heredoc issues)
  GEN_PY="/tmp/gen_landing_$$.py"
  cat > "${GEN_PY}" <<'PYEOF'
#!/usr/bin/env python3
"""TrackTx landing-page generator — called from module 15 shell."""
import csv, json, os, sys
from datetime import datetime, timezone

ROOT_DIR  = os.environ["ROOT_DIR"]
INDEX_OUT = os.path.join(ROOT_DIR, "index.html")
VER       = os.environ.get("PIPELINE_VERSION", "dev")
RUN_NAME  = os.environ.get("RUN_NAME",  "unnamed")
DURATION  = os.environ.get("DURATION",  "unknown")
PROFILE   = os.environ.get("PROFILE",   "unknown")
TIMESTAMP = os.environ.get("TIMESTAMP", datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"))

QC_MULTIQC = os.environ.get("QC_MULTIQC", "")
QC_IGV     = os.environ.get("QC_IGV",     "")
QC_RUNON   = os.environ.get("QC_RUNON",   "")
QC_PCA     = os.environ.get("QC_PCA",     "")
QC_CORR    = os.environ.get("QC_CORR",    "")

def present(p):
    return bool(p) and p != "NO_FILE" and os.path.isfile(p) and os.path.getsize(p) > 0

# ── KPIs from cohort JSON ────────────────────────────────────────────────────
kpis = dict(n_samples="?", n_conditions="?", conditions=[], median_div="?", median_pi="?")
cj   = "global_summary.json"
if os.path.isfile(cj):
    try:
        data    = json.load(open(cj))
        samples = data.get("samples", [])
        kpis["n_samples"]    = len(samples)
        conds = sorted(set(s.get("condition","?") for s in samples))
        kpis["n_conditions"] = len(conds)
        kpis["conditions"]   = conds
        div_v = sorted(x for s in samples
                       for x in [s.get("divergent_regions", s.get("div_regions"))]
                       if x is not None)
        if div_v:
            kpis["median_div"] = f"{div_v[len(div_v)//2]:,}"
        pi_v = sorted(x for s in samples
                      for x in [s.get("pi_median_len_norm", s.get("pi_median"))]
                      if x is not None)
        if pi_v:
            kpis["median_pi"] = f"{pi_v[len(pi_v)//2]:.2f}"
    except Exception as e:
        print(f"WARNING: cohort JSON parse error: {e}", file=sys.stderr)

# ── Run-on efficiency table ──────────────────────────────────────────────────
def runon_html():
    if not present(QC_RUNON):
        return "<p class='unavail'>Run-on efficiency data not available.</p>"
    try:
        with open(QC_RUNON, newline="") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        if not rows:
            return "<p class='unavail'>Run-on efficiency file is empty.</p>"
        hdrs = list(rows[0].keys())

        def badge(v):
            il = v.lower()
            cls = ("badge-excellent" if "excellent" in il else
                   "badge-good"      if "good"      in il else
                   "badge-moderate"  if "moderate"  in il else
                   "badge-poor"      if "poor"      in il else "")
            return f"<span class='badge {cls}'>{v}</span>"

        th  = "".join(f"<th>{h.replace('_',' ').title()}</th>" for h in hdrs)
        trs = []
        for row in rows:
            cells = []
            for h in hdrs:
                v = row.get(h, "")
                if h == "interpretation":
                    cells.append(f"<td>{badge(v)}</td>")
                elif h in ("ratio_5p_3p","efficiency_5p_3p","run_on_ratio"):
                    try:    cells.append(f"<td><strong>{float(v):.3f}</strong></td>")
                    except: cells.append(f"<td>{v}</td>")
                else:
                    cells.append(f"<td>{v}</td>")
            trs.append("<tr>" + "".join(cells) + "</tr>")
        return (f"<div class='table-wrap'><table class='data-table'>"
                f"<thead><tr>{th}</tr></thead>"
                f"<tbody>{''.join(trs)}</tbody></table></div>")
    except Exception as e:
        return f"<p class='unavail'>Could not read run-on TSV: {e}</p>"

# ── Link card helper ─────────────────────────────────────────────────────────
def card(href, icon, title, desc, ok=True):
    if ok:
        return (f"<a href='{href}' class='link-card'>"
                f"<div class='card-icon'>{icon}</div>"
                f"<div class='card-body'><h3>{title}</h3><p>{desc}</p></div>"
                f"</a>")
    return (f"<div class='link-card unavailable'>"
            f"<div class='card-icon'>{icon}</div>"
            f"<div class='card-body'>"
            f"<h3>{title} <span class='card-badge'>not available</span></h3>"
            f"<p>{desc}</p></div></div>")

cond_pills = "".join(f"<span class='cond-pill'>{c}</span>"
                     for c in kpis["conditions"])
cond_row   = (f"<div class='conditions-row'>{cond_pills}</div>"
              if kpis["conditions"] else "")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>TrackTx &middot; {RUN_NAME}</title>
  <style>
    :root{{--bg:#0f1117;--surface:#1a1d27;--surface2:#232636;--border:#2e3248;
          --accent:#6c8ef5;--accent2:#4ecdc4;--text:#e8eaf0;--muted:#7b82a0;}}
    *{{margin:0;padding:0;box-sizing:border-box;}}
    body{{font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Inter",sans-serif;
          background:var(--bg);color:var(--text);line-height:1.6;min-height:100vh;}}
    /* Hero */
    .hero{{background:linear-gradient(135deg,#1a1d27 0%,#0f1117 60%);
           border-bottom:1px solid var(--border);padding:3rem 2rem 2.5rem;text-align:center;}}
    .hero-eyebrow{{font-size:.75rem;letter-spacing:.15em;text-transform:uppercase;
                   color:var(--accent);margin-bottom:.75rem;}}
    .hero h1{{font-size:2.4rem;font-weight:700;
              background:linear-gradient(90deg,var(--accent),var(--accent2));
              -webkit-background-clip:text;-webkit-text-fill-color:transparent;
              background-clip:text;margin-bottom:.5rem;}}
    .hero-meta{{color:var(--muted);font-size:.9rem;}}
    .hero-meta span{{margin:0 .5rem;}}
    /* KPI bar */
    .kpi-bar{{display:flex;justify-content:center;flex-wrap:wrap;
              background:var(--surface);border-bottom:1px solid var(--border);}}
    .kpi{{flex:1;min-width:140px;padding:1.25rem 1.5rem;text-align:center;
          border-right:1px solid var(--border);}}
    .kpi:last-child{{border-right:none;}}
    .kpi-value{{font-size:2rem;font-weight:700;color:var(--accent);
                line-height:1;margin-bottom:.3rem;}}
    .kpi-label{{font-size:.75rem;text-transform:uppercase;letter-spacing:.08em;
                color:var(--muted);}}
    /* Conditions */
    .conditions-row{{display:flex;flex-wrap:wrap;gap:.5rem;justify-content:center;
                     padding:1rem 2rem;background:var(--surface);
                     border-bottom:1px solid var(--border);}}
    .cond-pill{{background:var(--surface2);border:1px solid var(--border);
                border-radius:999px;padding:.2rem .85rem;font-size:.8rem;
                color:var(--accent2);}}
    /* Layout */
    .main{{max-width:1100px;margin:0 auto;padding:2.5rem 2rem;}}
    /* Section headers */
    .section-header{{display:flex;align-items:center;gap:.75rem;
                     margin:2.5rem 0 1.25rem;}}
    .section-header h2{{font-size:1rem;font-weight:600;letter-spacing:.08em;
                        text-transform:uppercase;color:var(--muted);white-space:nowrap;}}
    .section-header::after{{content:"";flex:1;height:1px;background:var(--border);}}
    /* Card grid */
    .card-grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(230px,1fr));gap:1rem;}}
    .link-card{{display:flex;align-items:flex-start;gap:1rem;padding:1.25rem;
                background:var(--surface);border:1px solid var(--border);
                border-radius:10px;text-decoration:none;color:var(--text);
                transition:border-color .2s,transform .2s,box-shadow .2s;}}
    .link-card:hover{{border-color:var(--accent);transform:translateY(-2px);
                      box-shadow:0 6px 20px rgba(108,142,245,.15);}}
    .link-card.unavailable{{opacity:.42;cursor:not-allowed;}}
    .card-icon{{font-size:1.6rem;line-height:1;flex-shrink:0;}}
    .card-body h3{{font-size:.95rem;font-weight:600;color:var(--accent);margin-bottom:.25rem;}}
    .link-card.unavailable .card-body h3{{color:var(--muted);}}
    .card-body p{{font-size:.82rem;color:var(--muted);line-height:1.4;}}
    .card-badge{{display:inline-block;font-size:.65rem;font-weight:600;
                 text-transform:uppercase;letter-spacing:.05em;
                 background:var(--surface2);color:var(--muted);
                 border:1px solid var(--border);border-radius:4px;
                 padding:.1rem .4rem;vertical-align:middle;margin-left:.4rem;}}
    /* Run-on section */
    .runon-section{{background:var(--surface);border:1px solid var(--border);
                    border-radius:10px;padding:1.5rem;}}
    .runon-header{{display:flex;align-items:center;gap:.75rem;margin-bottom:1.25rem;}}
    .runon-header h3{{font-size:1rem;font-weight:600;}}
    .runon-header p{{font-size:.82rem;color:var(--muted);margin-left:auto;}}
    /* Table */
    .table-wrap{{overflow-x:auto;}}
    .data-table{{width:100%;border-collapse:collapse;font-size:.875rem;}}
    .data-table th{{text-align:left;padding:.6rem 1rem;background:var(--surface2);
                    color:var(--muted);font-size:.75rem;font-weight:600;
                    letter-spacing:.06em;text-transform:uppercase;
                    border-bottom:1px solid var(--border);}}
    .data-table td{{padding:.7rem 1rem;border-bottom:1px solid var(--border);}}
    .data-table tr:last-child td{{border-bottom:none;}}
    .data-table tr:hover td{{background:var(--surface2);}}
    /* Badges */
    .badge{{display:inline-block;font-size:.75rem;font-weight:600;
            padding:.2rem .65rem;border-radius:999px;white-space:nowrap;}}
    .badge-excellent{{background:rgba(76,175,130,.18);color:#4caf82;}}
    .badge-good{{background:rgba(108,142,245,.18);color:#6c8ef5;}}
    .badge-moderate{{background:rgba(245,197,66,.18);color:#f5c542;}}
    .badge-poor{{background:rgba(224,90,90,.18);color:#e05a5a;}}
    .unavail{{color:var(--muted);font-size:.9rem;padding:.75rem 0;}}
    /* Info grid */
    .info-grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(190px,1fr));gap:1rem;}}
    .info-item{{background:var(--surface);border:1px solid var(--border);
                border-radius:8px;padding:1rem 1.25rem;}}
    .info-label{{font-size:.72rem;text-transform:uppercase;letter-spacing:.08em;
                 color:var(--muted);margin-bottom:.35rem;}}
    .info-value{{font-size:.95rem;font-weight:600;word-break:break-all;}}
    /* Footer */
    .footer{{text-align:center;padding:2rem;color:var(--muted);font-size:.8rem;
             border-top:1px solid var(--border);margin-top:3rem;}}
    code{{background:var(--surface2);padding:.1rem .35rem;border-radius:3px;
          font-size:.85em;color:var(--accent2);}}
  </style>
</head>
<body>
  <div class="hero">
    <div class="hero-eyebrow">TrackTx PRO-seq Analysis Pipeline</div>
    <h1>Run Results</h1>
    <div class="hero-meta">
      <span>&#128203; {RUN_NAME}</span>
      <span>&middot;</span>
      <span>v{VER}</span>
      <span>&middot;</span>
      <span>&#9200; {DURATION}</span>
      <span>&middot;</span>
      <span>&#128197; {TIMESTAMP}</span>
    </div>
  </div>

  <div class="kpi-bar">
    <div class="kpi">
      <div class="kpi-value">{kpis["n_samples"]}</div>
      <div class="kpi-label">Samples</div>
    </div>
    <div class="kpi">
      <div class="kpi-value">{kpis["n_conditions"]}</div>
      <div class="kpi-label">Conditions</div>
    </div>
    <div class="kpi">
      <div class="kpi-value">{kpis["median_div"]}</div>
      <div class="kpi-label">Median Divergent Sites</div>
    </div>
    <div class="kpi">
      <div class="kpi-value">{kpis["median_pi"]}</div>
      <div class="kpi-label">Median Pausing Index</div>
    </div>
  </div>

  {cond_row}

  <div class="main">

    <div class="section-header"><h2>Analysis Reports</h2></div>
    <div class="card-grid">
      {card("11_reports/cohort/global_summary.html","&#128202;","Cohort Summary",
            "Aggregated metrics, QC charts, and Pol&nbsp;II pausing across all samples")}
      {card("12_cohort_qc/multiqc/multiqc_report.html","&#128300;","MultiQC Report",
            "Alignment rates, trimming stats, and library QC for all samples",
            present(QC_MULTIQC))}
    </div>

    <div class="section-header"><h2>Signal QC</h2></div>
    <div class="card-grid">
      {card("12_cohort_qc/deeptools/pca_plot.pdf","&#128201;","PCA Plot",
            "Sample clustering from genome-wide 3&prime; BigWig signal (deepTools)",
            present(QC_PCA))}
      {card("12_cohort_qc/deeptools/correlation_heatmap.pdf","&#128279;","Correlation Heatmap",
            "Pearson correlation between all samples (deepTools)",
            present(QC_CORR))}
      {card("12_cohort_qc/igv_session.xml","&#129520;","IGV Session",
            "All tracks colour-coded by condition &mdash; drag into IGV to open",
            present(QC_IGV))}
    </div>

    <div class="section-header"><h2>Per-Sample Reports</h2></div>
    <div class="card-grid">
      {card("11_reports/","&#128196;","Sample Reports",
            "Individual HTML reports: divergent TX sites, Pol&nbsp;II metrics, region composition")}
    </div>

    <div class="section-header"><h2>Run-on Efficiency</h2></div>
    <div class="runon-section">
      <div class="runon-header">
        <span style="font-size:1.4rem">&#9879;</span>
        <h3>5&prime; / 3&prime; Signal Ratio</h3>
        <p>Measures nascent RNA synthesis quality over gene bodies &ge;&thinsp;10&thinsp;kb</p>
      </div>
      {runon_html()}
    </div>

    <div class="section-header"><h2>Pipeline Execution</h2></div>
    <p style="color:var(--muted);font-size:.82rem;margin-bottom:1rem;">
      Generated with <code>-with-report</code>, <code>-with-timeline</code>,
      <code>-with-dag</code>. Links may be unavailable if flags were not used.
    </p>
    <div class="card-grid">
      {card("nf_report.html","&#128200;","Execution Report","CPU, memory, and I/O usage per task")}
      {card("nf_timeline.html","&#9200;","Timeline","Task execution timeline")}
      {card("nf_dag.png","&#128256;","Workflow DAG","Pipeline execution graph")}
      {card("trace/trace.txt","&#128221;","Trace Log","Detailed per-process execution log")}
    </div>

    <div class="section-header"><h2>Pipeline Information</h2></div>
    <div class="info-grid">
      <div class="info-item">
        <div class="info-label">Pipeline</div>
        <div class="info-value">TrackTx PRO-seq</div>
      </div>
      <div class="info-item">
        <div class="info-label">Version</div>
        <div class="info-value">{VER}</div>
      </div>
      <div class="info-item">
        <div class="info-label">Run Name</div>
        <div class="info-value">{RUN_NAME}</div>
      </div>
      <div class="info-item">
        <div class="info-label">Duration</div>
        <div class="info-value">{DURATION}</div>
      </div>
      <div class="info-item">
        <div class="info-label">Profile</div>
        <div class="info-value">{PROFILE}</div>
      </div>
      <div class="info-item">
        <div class="info-label">Completed</div>
        <div class="info-value">{TIMESTAMP}</div>
      </div>
    </div>

  </div>

  <div class="footer">
    TrackTx PRO-seq Analysis Pipeline &nbsp;&middot;&nbsp; {TIMESTAMP}
  </div>
</body>
</html>"""

os.makedirs(ROOT_DIR, exist_ok=True)
with open(INDEX_OUT, "w") as f:
    f.write(html)
print(f"Landing page written: {INDEX_OUT}")
PYEOF

  # Run the generator, passing all values as environment variables
  ROOT_DIR="${ROOT_DIR}" \
  PIPELINE_VERSION="${PIPELINE_VERSION}" \
  RUN_NAME="${RUN_NAME}" \
  DURATION="${DURATION}" \
  PROFILE="${PROFILE}" \
  TIMESTAMP="${TIMESTAMP}" \
  QC_MULTIQC="${QC_MULTIQC}" \
  QC_IGV="${QC_IGV}" \
  QC_RUNON="${QC_RUNON}" \
  QC_PCA="${QC_PCA}" \
  QC_CORR="${QC_CORR}" \
  ${PYTHON_CMD} "${GEN_PY}" || \
    echo "COHORT | WARNING | Landing page generation failed (check Python errors above)"

  rm -f "${GEN_PY}"

  echo "COHORT | LANDING | Landing page: ${INDEX_HTML}"

  ###########################################################################
  # 10) CREATE README
  ###########################################################################

  echo "COHORT | README | Creating documentation..."

  cat > "${OUT_README}" <<DOCEOF
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
   • Small-cohort summary card (n ≤ 5) with QC and design verdicts

3. Divergent Transcription
   • Total loci across cohort
   • Distribution per sample (histograms for larger cohorts; dot/text summaries for n < 6)
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
   • Distribution across samples (histograms for larger cohorts; dot/text summaries for n < 6)
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
  5. For small cohorts (≤5 samples), focus on per-sample QC plots and the small-cohort summary card.
  6. Compare metrics across conditions

TSV Export:
  # Load in R
  cohort <- read.delim("global_summary.tsv")
  
  # Filter by condition
  treatment <- subset(cohort, condition == "treatment")
  
  # Calculate statistics
  mean(cohort\\$map_rate_percent)
  
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
  - Very small cohorts (≤5 samples): Treat cohort-level distributions as qualitative; rely on per-sample QC and the small-cohort summary card for decisions

DOWNSTREAM ANALYSIS
────────────────────────────────────────────────────────────────────────────
  
  1. Comparative Analysis:
     - Compare pausing indices and divergent transcription across conditions
     - Identify condition-specific transcription patterns
     - Assess replicate consistency (especially via per-condition QC summaries)
  
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
  Module:   15_combine_reports_into_cohort

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