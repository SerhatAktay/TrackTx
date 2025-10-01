#!/usr/bin/env python3
# =============================================================================
# render_sample_report.py — Per-sample HTML / TSV / JSON (+ optional plots)
# -----------------------------------------------------------------------------
# • Science-first summary: divergent loci, functional regions, Pol II density,
#   pausing metrics (favor length-normalized), normalization factors, QC.
# • Stable JSON schema (schema_version) + TSV (key name/value pairs).
# • Track links can be file paths or URLs; echoed verbatim in JSON/HTML.
# • Robust parsing (header-agnostic fallbacks), no hard failure on partial data.
# =============================================================================

from __future__ import annotations
import argparse, io, json, os, base64, re, sys, datetime
from typing import Optional, Dict, Any, List, Tuple
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ── CLI ────────────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--condition", required=True)
ap.add_argument("--timepoint", required=True)
ap.add_argument("--replicate", required=True)

ap.add_argument("--divergent-bed", required=True)
ap.add_argument("--functional-summary", required=True)
ap.add_argument("--pol2-density", required=True)
ap.add_argument("--pausing-index", required=True)
ap.add_argument("--norm-factors", required=True)
ap.add_argument("--qc-json", required=True)

ap.add_argument("--out-html", required=True)
ap.add_argument("--out-tsv", required=True)
ap.add_argument("--out-json", required=True)
ap.add_argument("--out-plots-html", required=True)
ap.add_argument("--plots", type=int, default=0)

# Optional track links (strings; may be empty; can be URLs or paths)
ap.add_argument("--pos3-cpm-bw", default=None)
ap.add_argument("--neg3-cpm-bw", default=None)
ap.add_argument("--allmap3p-pos-raw", default=None)
ap.add_argument("--allmap3p-neg-raw", default=None)
ap.add_argument("--allmap3p-pos-cpm-bw", default=None)
ap.add_argument("--allmap3p-neg-cpm-bw", default=None)

args = ap.parse_args()
print(f"[render_py] start ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)

SID, COND, TP, REP = args.sample, args.condition or "", args.timepoint or "", args.replicate or ""

# ── helpers ────────────────────────────────────────────────────────────────
def exists_nonempty(p: Optional[str]) -> bool:
    try: return bool(p) and os.path.exists(p) and os.path.getsize(p) > 0
    except Exception: return False

def read_table_any(path: Optional[str]) -> pd.DataFrame:
    if not exists_nonempty(path): return pd.DataFrame()
    for sep in ["\t", ",", ";", r"\s+"]:
        try:
            df = pd.read_csv(path, sep=sep, engine="python")
            if df.shape[0]:
                return df
        except Exception:
            pass
    return pd.DataFrame()

def read_bed6(path: Optional[str]) -> pd.DataFrame:
    df = read_table_any(path)
    if df.empty: return df
    while df.shape[1] < 6:
        df[df.shape[1]] = "."
    df = df.iloc[:, :6]
    df.columns = ["chr","start","end","name","score","strand"]
    return df

def to_num(x) -> float:
    try:
        v = float(x);  return v if np.isfinite(v) else np.nan
    except Exception:
        return np.nan

def median_num(series: pd.Series) -> float:
    x = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    return float(np.median(x)) if len(x) else float("nan")

def b64(fig) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight"); plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode("ascii")

def img_html(fig) -> str:
    return f"<img src='data:image/png;base64,{b64(fig)}' alt='plot'>"

def pick_region_reads_columns(df: pd.DataFrame) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    if df.empty: return None, None, None
    low = [str(c).strip().lower() for c in df.columns]
    reg = next((df.columns[i] for i,c in enumerate(low) if c in {"region","regions","feature","name","id"}), None)
    rea = next((df.columns[i] for i,c in enumerate(low) if c in {"reads","count","raw_count","n_reads","nread","n","signal"}), None)
    cnt = next((df.columns[i] for i,c in enumerate(low) if c in {"region_count","count","n_regions","num_regions"}), None)
    if reg is None:
        scores = [(1.0 - pd.to_numeric(df[c], errors="coerce").notna().mean(), c) for c in df.columns]
        reg = sorted(scores, reverse=True)[0][1]
    if rea is None:
        scores = [(pd.to_numeric(df[c], errors="coerce").notna().mean(), c) for c in df.columns]
        rea = sorted(scores, reverse=True)[0][1]
    if cnt is None and len(df.columns) >= 3:
        # If we have 3+ columns and no explicit count column, try the third column
        cnt = df.columns[2] if df.columns[2] != reg and df.columns[2] != rea else None
    if reg == rea: rea = None
    return reg, rea, cnt

def norm_factor_lookup(df: pd.DataFrame, key: str) -> Optional[float]:
    if df.empty: return None
    # Accept either 2-col key/value or labeled columns
    if df.shape[1] == 2:
        kv = {str(k).strip().lower(): str(v) for k,v in df.values.tolist()}
        for k,v in kv.items():
            if key in k: return to_num(v)
        return None
    low = {c.lower(): c for c in df.columns}
    for cand in ["cpm","cpm_factor","fac_cpm"]:
        if cand in low and key == "cpm": return to_num(df[low[cand]].iloc[0])
    for cand in ["sicpm","crpmsi","sicpm_factor","fac_sicpm"]:
        if cand in low and key == "sicpm": return to_num(df[low[cand]].iloc[0])
    return None

def nl(s: Optional[str]) -> str:
    return "" if s is None or str(s).strip()=="" else str(s)

# ── load inputs ────────────────────────────────────────────────────────────
bed_div  = read_bed6(args.divergent_bed)
fun_sum  = read_table_any(args.functional_summary)
density  = read_table_any(args.pol2_density)
pausing  = read_table_any(args.pausing_index)
norm_df  = read_table_any(args.norm_factors)

qc: Optional[Dict[str, Any]] = None
if exists_nonempty(args.qc_json):
    try:
        with open(args.qc_json, "r", encoding="utf-8") as fh:
            qc = json.load(fh)
    except Exception:
        qc = None

# ── normalize density columns ──────────────────────────────────────────────
if not density.empty:
    low = [c.lower() for c in density.columns]
    ren = {}
    for a,b in (("chrom","chr"),("contig","chr")):
        if a in low: ren[density.columns[low.index(a)]] = b
    for k in ("start","end","name","score","strand","cpm","raw_count","count","raw"):
        if k in low: ren[density.columns[low.index(k)]] = k
    density = density.rename(columns=ren)
    for k in ("start","end","raw_count","count","cpm"):
        if k in density.columns:
            density[k] = pd.to_numeric(density[k], errors="coerce")

# ── region counts (robust) ────────────────────────────────────────────────
region_counts_df = pd.DataFrame()
regions_list: List[Dict[str, Any]] = []
if not fun_sum.empty:
    f = fun_sum.copy()
    if f.shape[1] == 2 and not any(re.search(r"(region|name|feature|reads|count)", str(c), re.I) for c in f.columns):
        f.columns = ["region", "reads"]
    elif f.shape[1] == 3 and not any(re.search(r"(region|name|feature|reads|count)", str(c), re.I) for c in f.columns):
        f.columns = ["region", "reads", "region_count"]
    reg_col, rea_col, cnt_col = pick_region_reads_columns(f)
    if reg_col and rea_col:
        cols_to_use = [reg_col, rea_col]
        rename_dict = {reg_col:"region", rea_col:"reads"}
        if cnt_col:
            cols_to_use.append(cnt_col)
            rename_dict[cnt_col] = "region_count"
        
        f2 = f[cols_to_use].rename(columns=rename_dict).copy()
        f2["region"] = f2["region"].astype(str).str.strip()
        f2["reads"]  = pd.to_numeric(f2["reads"], errors="coerce").fillna(0.0)
        if "region_count" in f2.columns:
            f2["region_count"] = pd.to_numeric(f2["region_count"], errors="coerce").fillna(0.0)
        
        f2 = f2.groupby("region", as_index=False).agg({
            "reads": "sum",
            **({"region_count": "sum"} if "region_count" in f2.columns else {})
        }).sort_values("reads", ascending=False)
        region_counts_df = f2
        regions_list = []
        for _, row in f2.iterrows():
            region_dict = {"region": row["region"], "reads": float(row["reads"])}
            if "region_count" in row:
                region_dict["region_count"] = int(row["region_count"])
            regions_list.append(region_dict)

# ── summary stats ─────────────────────────────────────────────────────────
divergent_regions = int(len(bed_div)) if not bed_div.empty else 0
total_regions     = int(region_counts_df["region"].nunique()) if not region_counts_df.empty else 0
reads_total_func  = float(pd.to_numeric(region_counts_df["reads"], errors="coerce").fillna(0).sum()) if not region_counts_df.empty else 0.0

median_pausing = float("nan")
if not pausing.empty:
    low = [c.lower() for c in pausing.columns]
    for cand in ("pi_len_norm","pausing_index","pi_raw","pi"):
        if cand in low:
            median_pausing = median_num(pausing[pausing.columns[low.index(cand)]])
            break

median_density = float("nan")
if not density.empty:
    cpm_col = next((density.columns[i] for i,c in enumerate([c.lower() for c in density.columns]) if c=="cpm"), None)
    if cpm_col:
        median_density = median_num(density[cpm_col])

# Norm factors (accept synonyms)
cpm_factor    = norm_factor_lookup(norm_df, "cpm")
sicpm_factor  = norm_factor_lookup(norm_df, "sicpm")

# QC fields used for status + unlocalized fraction
input_reads       = qc.get("total_reads_raw") if isinstance(qc, dict) else None
dedup_reads       = qc.get("mapq_ge_reads_nodup") if isinstance(qc, dict) else None
duplicate_percent = qc.get("duplicate_perc_of_total") if isinstance(qc, dict) else None

# UMI deduplication fields
umi_enabled = qc.get("umi_deduplication_enabled", False) if isinstance(qc, dict) else False
umi_input_reads = qc.get("umi_input_reads") if isinstance(qc, dict) else None
umi_output_reads = qc.get("umi_output_reads") if isinstance(qc, dict) else None
umi_dedup_percent = qc.get("umi_deduplication_percent") if isinstance(qc, dict) else None

def safe_int(x): 
    try: 
        v=int(x) 
        return v if v>=0 else None
    except: 
        return None

input_reads = safe_int(input_reads)
dedup_reads = safe_int(dedup_reads)
duplicate_percent = None if duplicate_percent is None else float(duplicate_percent)

# Unlocalized fraction (1 − functional_reads / dedup_reads), if both present
unlocalized_fraction = None
if dedup_reads and dedup_reads > 0:
    frac = max(0.0, 1.0 - min(1.0, float(reads_total_func) / float(dedup_reads)))
    unlocalized_fraction = round(frac, 6)

# Status badge (very light heuristic)
def status_strip() -> str:
    inp = float(input_reads or 0)
    
    # Use UMI deduplication rate if available, otherwise fall back to PCR duplicate rate
    if umi_enabled and umi_dedup_percent is not None:
        dup = float(umi_dedup_percent)
    else:
        dup = float(duplicate_percent if duplicate_percent is not None else 100.0)
    
    ok_regions = (divergent_regions > 0) and (total_regions > 0)
    if inp >= 5e6 and dup < 15 and ok_regions: s, color = "PASS", "#10b981"
    elif inp >= 2e6 and dup < 30 and ok_regions: s, color = "WARN", "#f59e0b"
    else: s, color = "FAIL", "#ef4444"
    return f"<div class='status' style='--status:{color}'>{s}</div>"

# ── JSON (stable schema) ───────────────────────────────────────────────────
schema_version = "1.0.0"

row = dict(
    schema_version=schema_version,
    sample=SID, condition=COND, timepoint=TP, replicate=REP,
    metrics=dict(
        divergent_regions=int(divergent_regions),
        total_functional_regions=int(total_regions),
        reads_total_functional=float(reads_total_func),
        median_pausing_index=None if np.isnan(median_pausing) else float(median_pausing),
        median_functional_cpm=None if np.isnan(median_density) else float(median_density),
        unlocalized_fraction=unlocalized_fraction,
        cpm_factor=None if (cpm_factor is None or np.isnan(cpm_factor)) else float(cpm_factor),
        sicpm_factor=None if (sicpm_factor is None or np.isnan(sicpm_factor)) else float(sicpm_factor),
    ),
    qc=dict(
        total_reads_raw=input_reads,
        dedup_reads_mapq_ge=dedup_reads,
        duplicate_percent=duplicate_percent
    ),
    regions=regions_list,
    tracks=dict(
        pos3_cpm_bw=args.pos3_cpm_bw,
        neg3_cpm_bw=args.neg3_cpm_bw,
        allmap3p_pos_raw=args.allmap3p_pos_raw,
        allmap3p_neg_raw=args.allmap3p_neg_raw,
        allmap3p_pos_cpm_bw=args.allmap3p_pos_cpm_bw,
        allmap3p_neg_cpm_bw=args.allmap3p_neg_cpm_bw
    ),
    report_ok=True
)

# ── TSV (key name/value) ───────────────────────────────────────────────────
pairs = [
    ("schema_version", schema_version),
    ("sample", SID), ("condition", COND), ("timepoint", TP), ("replicate", REP),
    ("divergent_regions", row["metrics"]["divergent_regions"]),
    ("total_functional_regions", row["metrics"]["total_functional_regions"]),
    ("reads_total_functional", row["metrics"]["reads_total_functional"]),
    ("median_pausing_index", row["metrics"]["median_pausing_index"]),
    ("median_functional_cpm", row["metrics"]["median_functional_cpm"]),
    ("unlocalized_fraction", row["metrics"]["unlocalized_fraction"]),
    ("cpm_factor", row["metrics"]["cpm_factor"]),
    ("sicpm_factor", row["metrics"]["sicpm_factor"]),
    ("qc_total_reads_raw", row["qc"]["total_reads_raw"]),
    ("qc_dedup_reads_mapq_ge", row["qc"]["dedup_reads_mapq_ge"]),
    ("qc_duplicate_percent", row["qc"]["duplicate_percent"]),
]
pd.DataFrame(pairs, columns=["Metric","Value"]).to_csv(args.out_tsv, sep="\t", index=False)

with open(args.out_json, "w", encoding="utf-8") as fh:
    json.dump(row, fh, indent=2)

# ── HTML rendering (style consistent with other modules) ───────────────────
def html_table(df: pd.DataFrame, cols=None, max_rows=25) -> str:
    if df is None or df.empty: return "<p class='muted'>n/a</p>"
    if cols: df = df[[c for c in cols if c in df.columns]]
    return df.head(max_rows).to_html(index=False, escape=True, classes=["tbl"])

def kpi(label: str, val: Any, sub: str="") -> str:
    v = "" if val is None or (isinstance(val,float) and np.isnan(val)) else val
    return f"""
    <div class="kpi">
      <div class="kpi-label">{label}</div>
      <div class="kpi-value">{v}</div>
      <div class="kpi-sub">{sub}</div>
    </div>"""

def stacked_bar_regions(reg_df: pd.DataFrame) -> str:
    if reg_df is None or reg_df.empty: return "<div class='muted'>n/a</div>"
    fig = plt.figure(figsize=(7.5,1.8), dpi=160)
    vals = pd.to_numeric(reg_df["reads"], errors="coerce").fillna(0).values
    total = float(vals.sum()) or 1.0
    left = 0.0
    for v in vals:
        w = v/total
        plt.barh([0], [w], left=[left]); left += w
    plt.gca().set_yticks([]); plt.gca().set_xticks([]); plt.title("Functional-region composition (reads, normalized)")
    return img_html(fig)

def pie_regions(reg_df: pd.DataFrame) -> str:
    if reg_df is None or reg_df.empty: return "<div class='muted'>n/a</div>"
    vals = pd.to_numeric(reg_df["reads"], errors="coerce").fillna(0).values
    if vals.sum() <= 0: return "<div class='muted'>n/a</div>"
    fig = plt.figure(figsize=(4.6,4.6), dpi=160)
    plt.pie(vals, labels=reg_df["region"].astype(str).tolist(), autopct="%1.1f%%")
    plt.title("Functional regions (reads)")
    return img_html(fig)

def hist_cpm(density: pd.DataFrame) -> str:
    if density.empty or not any(c.lower()=="cpm" for c in density.columns): return "<div class='muted'>n/a</div>"
    cpm_col = [c for c in density.columns if c.lower()=="cpm"][0]
    x = pd.to_numeric(density[cpm_col], errors="coerce").replace([np.inf,-np.inf], np.nan).dropna().values
    if x.size == 0: return "<div class='muted'>n/a</div>"
    fig = plt.figure(figsize=(6.5,3.2), dpi=160); plt.hist(np.log10(x+1e-9), bins=50)
    plt.xlabel("log10(CPM)"); plt.ylabel("Regions"); plt.title("Functional-region CPM distribution")
    return img_html(fig)

def hist_pausing(pausing: pd.DataFrame) -> str:
    if pausing.empty: return "<div class='muted'>n/a</div>"
    low = [c.lower() for c in pausing.columns]
    col = next((pausing.columns[low.index(c)] for c in ("pi_len_norm","pausing_index","pi_raw","pi") if c in low), None)
    if not col: return "<div class='muted'>n/a</div>"
    vals = pd.to_numeric(pausing[col], errors="coerce").replace([np.inf,-np.inf], np.nan).dropna().values
    if vals.size == 0: return "<div class='muted'>n/a</div>"
    fig = plt.figure(figsize=(6.5,3.2), dpi=160); plt.hist(np.log2(vals+1e-9), bins=60)
    plt.xlabel(f"log2({col} + 1e-9)"); plt.ylabel("Genes"); plt.title(f"{col} distribution")
    return img_html(fig)

track_links = [(k,v) for k,v in [
    ("pos3_cpm.bw",            nl(args.pos3_cpm_bw)),
    ("neg3_cpm.bw",            nl(args.neg3_cpm_bw)),
    ("allMap3p_pos_raw.bg",    nl(args.allmap3p_pos_raw)),
    ("allMap3p_neg_raw.bg",    nl(args.allmap3p_neg_raw)),
    ("allMap3p_pos_cpm.bw",    nl(args.allmap3p_pos_cpm_bw)),
    ("allMap3p_neg_cpm.bw",    nl(args.allmap3p_neg_cpm_bw)),
] if v]

html = io.StringIO()
html.write(f"""<!doctype html><meta charset="utf-8">
<title>{SID} — TrackTx sample</title>
<style>
:root {{
  --bg:#fff; --fg:#111; --muted:#667085; --card:#f8fafc; --line:#ececec;
  --ok:#10b981; --warn:#f59e0b; --fail:#ef4444; --chip:#e5e7eb;
}}
@media (prefers-color-scheme: dark) {{
  :root {{ --bg:#0b0d10; --fg:#e5e7eb; --muted:#9aa4b2; --card:#11161c; --line:#1f2937; --chip:#374151; }}
}}
html,body{{background:var(--bg);color:var(--fg);font-family:system-ui,-apple-system,Segoe UI,Roboto,sans-serif;margin:0}}
.wrap{{max-width:1040px;margin:0 auto;padding:24px}}
h1{{margin:0 0 6px 0;font-size:22px}}
h2{{margin:20px 0 10px 0;font-size:16px}}
.mono{{font-family:ui-monospace,SFMono-Regular,Consolas,monospace}}
.muted{{color:var(--muted)}}
hr{{border:none;border-top:1px solid var(--line);margin:16px 0}}
.grid{{display:grid;gap:14px;grid-template-columns:repeat(auto-fit,minmax(260px,1fr))}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px}}
.kpi{{background:var(--bg);border:1px solid var(--line);border-radius:10px;padding:12px}}
.kpi-label{{font-size:12px;color:var(--muted)}}
.kpi-value{{font-size:18px;font-weight:600;margin-top:2px}}
.kpi-sub{{font-size:11px;color:var(--muted)}}
.tbl{{border-collapse:collapse;width:100%;font-size:13px}}
.tbl th,.tbl td{{border-bottom:1px solid var(--line);padding:6px 8px;text-align:left}}
.status{{display:inline-flex;align-items:center;gap:8px;font-weight:600;margin:8px 0}}
.status::before{{content:'';width:10px;height:10px;border-radius:50%;background:var(--status,#bbb)}}
.meta small{{color:var(--muted)}}
.dl a{{color:inherit;text-decoration:underline}}
</style>
<div class="wrap">
  <div class="meta">
    <h1>{SID}</h1>
    <div>{status_strip()}</div>
    <div class="muted">Condition: <b>{COND}</b> • Timepoint: <b>{TP}</b> • Replicate: <b>{REP}</b></div>
  </div>

  <h2>At-a-glance</h2>
  <div class="grid">
    {kpi("Total input reads", row["qc"]["total_reads_raw"], "qc_pol2.json")}
    {kpi("De-dup reads (MAPQ≥)", row["qc"]["dedup_reads_mapq_ge"], "qc_pol2.json")}
    {kpi("UMI Dedup %" if row["qc"].get("umi_deduplication_enabled", False) else "Duplicate %", 
         row["qc"].get("umi_deduplication_percent") if row["qc"].get("umi_deduplication_enabled", False) else row["qc"].get("duplicate_percent"))}
    {kpi("# divergent loci", row["metrics"]["divergent_regions"])}
    {kpi("Functional regions (distinct)", row["metrics"]["total_functional_regions"])}
    {kpi("Reads in functional regions", int(row["metrics"]["reads_total_functional"]))}
    {kpi("Unlocalized fraction", row["metrics"]["unlocalized_fraction"])}
    {kpi("Median pausing (len-norm)", None if np.isnan(median_pausing) else round(median_pausing,3))}
    {kpi("Median functional CPM", None if np.isnan(median_density) else round(median_density,3))}
    {kpi("CPM factor", row["metrics"]["cpm_factor"])}
    {kpi("siCPM factor", row["metrics"]["sicpm_factor"])}
  </div>

  <h2>Functional regions — composition</h2>
  <div class="grid">
    <div class="card">{stacked_bar_regions(region_counts_df)}</div>
    <div class="card">{pie_regions(region_counts_df)}</div>
  </div>

  <h2>Tables</h2>
  <div class="card">
    <h3 style="margin:0 0 8px">Reads per functional region</h3>
    {html_table(region_counts_df.rename(columns={"region":"Region","reads":"Reads","region_count":"Region Count"}))}
  </div>

  <div class="grid">
    <div class="card">
      <h3 style="margin:0 0 8px">Top functional regions by CPM</h3>
      {"<p class='muted'>n/a</p>" if density.empty else html_table(
         (density.sort_values([c for c in density.columns if c.lower()=="cpm"][0], ascending=False)
          if any(c.lower()=="cpm" for c in density.columns) else density),
         cols=["chr","start","end","name","cpm"] if any(c.lower()=="cpm" for c in density.columns) else None, max_rows=25)}
    </div>
    <div class="card">
      <h3 style="margin:0 0 8px">Pausing preview</h3>
      {"<p class='muted'>n/a</p>" if pausing.empty else html_table(pausing, max_rows=25)}
    </div>
  </div>

  <h2>Distributions</h2>
  <div class="grid">
    <div class="card">{hist_cpm(density)}</div>
    <div class="card">{hist_pausing(pausing)}</div>
  </div>

  <h2>Downloads</h2>
  <div class="card dl">
    <ul>
      <li><span class="mono">{os.path.basename(args.divergent_bed)}</span></li>
      <li><span class="mono">{os.path.basename(args.functional_summary)}</span></li>
      <li><span class="mono">{os.path.basename(args.pol2_density)}</span></li>
      <li><span class="mono">{os.path.basename(args.pausing_index)}</span></li>
      <li><span class="mono">{os.path.basename(args.norm_factors)}</span></li>
      <li><span class="mono">{os.path.basename(args.qc_json)}</span></li>
      <li><span class="mono">{os.path.basename(args.out_tsv)}</span></li>
      <li><span class="mono">{os.path.basename(args.out_json)}</span></li>
    </ul>
  </div>

  {("<h2>Track links</h2><div class='card'><ul>" +
     "".join(f"<li><span class='mono'>{k}</span>: {v}</li>" for k,v in track_links) +
     "</ul></div>") if track_links else ""}

  <hr>
  <div class="muted"><small>Generated by TrackTx — render_sample_report.py (schema {schema_version})</small></div>
</div>
""")
with open(args.out_html, "w", encoding="utf-8") as fh:
    fh.write(html.getvalue())

# ── Optional plots page ────────────────────────────────────────────────────
if int(args.plots) == 1:
    page = io.StringIO()
    page.write(f"<!doctype html><meta charset='utf-8'><title>{SID} plots</title><h1>{SID} plots</h1>")

    if not region_counts_df.empty and region_counts_df["reads"].sum() > 0:
        fig = plt.figure(figsize=(6,6), dpi=160)
        plt.pie(region_counts_df["reads"], labels=region_counts_df["region"], autopct=lambda p: f"{p:.1f}%")
        plt.title("Reads per functional region")
        page.write(f"<h3>Reads per functional region</h3>{img_html(fig)}")

    if not density.empty and any(c.lower()=="cpm" for c in density.columns):
        cpm_col = [c for c in density.columns if c.lower()=="cpm"][0]
        x = pd.to_numeric(density[cpm_col], errors="coerce").replace([np.inf,-np.inf], np.nan).dropna().values
        if x.size:
            fig = plt.figure(figsize=(6,3.4), dpi=160)
            plt.hist(np.log10(x + 1e-9), bins=50)
            plt.xlabel("log10(CPM)"); plt.ylabel("Count"); plt.title("Functional-region CPM")
            page.write(f"<h3>Functional-region CPM</h3>{img_html(fig)}")

    if not pausing.empty:
        low = [c.lower() for c in pausing.columns]
        col = next((pausing.columns[low.index(c)] for c in ("pi_len_norm","pausing_index","pi_raw","pi") if c in low), None)
        if col:
            vals = pd.to_numeric(pausing[col], errors="coerce").replace([np.inf,-np.inf], np.nan).dropna().values
            if vals.size:
                fig = plt.figure(figsize=(6,3.4), dpi=160)
                plt.hist(np.log2(vals + 1e-9), bins=60)
                plt.xlabel(f"log2({col}+1e-9)"); plt.ylabel("Genes"); plt.title(f"{col} distribution")
                page.write(f"<h3>{col} distribution</h3>{img_html(fig)}")

    with open(args.out_plots_html, "w", encoding="utf-8") as fh:
        fh.write(page.getvalue())
print(f"[render_py] done ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)
