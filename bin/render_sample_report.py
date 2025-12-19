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
    """Read a table file with automatic separator detection"""
    if not exists_nonempty(path): 
        return pd.DataFrame()
    
    # Try multiple separators in order of likelihood
    # \s+ is most flexible and handles multiple spaces/tabs
    for sep in [r"\s+", "\t", ",", ";"]:
        try:
            df = pd.read_csv(path, sep=sep, engine="python", skip_blank_lines=True)
            # Valid if we got rows AND columns
            if df.shape[0] > 0 and df.shape[1] > 0:
                return df
        except Exception:
            continue
    
    return pd.DataFrame()

def read_bed6(path: Optional[str]) -> pd.DataFrame:
    df = read_table_any(path)
    if df.empty: return df
    while df.shape[1] < 6:
        df[df.shape[1]] = "."
    df = df.iloc[:, :6]
    df.columns = ["chr","start","end","name","score","strand"]
    return df

def to_num(x) -> Optional[float]:
    """Convert to float, returning None for invalid/infinite values"""
    if x is None:
        return None
    try:
        v = float(x)
        # Only return valid, finite numbers
        if np.isfinite(v):
            return v
        return None
    except (ValueError, TypeError, AttributeError):
        return None

def median_num(series: pd.Series) -> float:
    x = pd.to_numeric(series, errors="coerce").replace([np.inf, -np.inf], np.nan).dropna()
    return float(np.median(x)) if len(x) > 0 else float("nan")

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
    """
    Look up normalization factor from DataFrame.
    
    Handles two formats:
    1. Two-column key-value format: method | factor
    2. Multi-column format with explicit column names
    """
    if df.empty or df.shape[0] == 0:
        return None
    
    key_lower = key.lower()
    
    # Check if columns themselves are the method names (single row format)
    # Example: CPM  siCPM
    #          0.01 0.02
    col_lower = {str(c).strip().lower(): c for c in df.columns}
    if key_lower in col_lower and df.shape[0] > 0:
        result = to_num(df[col_lower[key_lower]].iloc[0])
        if result is not None:
            return result
    
    # Standard two-column key-value format (most common)
    # Example: method  factor
    #          CPM     0.0106
    #          siCPM   0.0106
    if df.shape[1] >= 2:
        # Try iterating through rows
        for idx in range(len(df)):
            try:
                k = str(df.iloc[idx, 0]).strip().lower()
                v = str(df.iloc[idx, 1]).strip()
                
                # Exact match or substring match
                if key_lower == k or key_lower in k:
                    result = to_num(v)
                    if result is not None:
                        return result
            except (IndexError, ValueError, AttributeError, TypeError):
                continue
    
    # Format 3: Column-based with explicit names (e.g., 'cpm_factor' column)
    for candidate in ["cpm", "cpm_factor", "fac_cpm"] if key == "cpm" else ["sicpm", "crpmsi", "sicpm_factor", "fac_sicpm", "si_cpm"]:
        if candidate in col_lower and df.shape[0] > 0:
            result = to_num(df[col_lower[candidate]].iloc[0])
            if result is not None:
                return result
    
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
unlocalized_reads = 0.0  # Track unlocalized separately
# Pattern as string for pandas (pandas will compile it with the flags we provide)
unloc_pattern = r"non[-\s_]?localized|unlocalized|unlocalised"

if not fun_sum.empty:
    f = fun_sum.copy()
    if f.shape[1] == 2 and not any(re.search(r"(region|name|feature|reads|count)", str(c), re.I) for c in f.columns):
        f.columns = ["region", "reads"]
    elif f.shape[1] >= 3 and not any(re.search(r"(region|name|feature|reads|count)", str(c), re.I) for c in f.columns):
        # New summary may include length columns; assign conservative defaults
        base = ["region", "reads", "region_count"]
        extra = [f"col{i}" for i in range(3, f.shape[1])]
        f.columns = base + extra
    reg_col, rea_col, cnt_col = pick_region_reads_columns(f)
    if reg_col and rea_col:
        cols_to_use = [reg_col, rea_col]
        rename_dict = {reg_col:"region", rea_col:"reads"}
        if cnt_col:
            cols_to_use.append(cnt_col)
            rename_dict[cnt_col] = "region_count"
        
        # Optional length columns
        length_total_col = next((c for c in f.columns if str(c).strip().lower() in {"region_length_total_bp","length_total","total_length"}), None)
        length_median_col = next((c for c in f.columns if str(c).strip().lower() in {"region_length_median_bp","length_median","median_length"}), None)

        f2 = f[cols_to_use + ([length_total_col] if length_total_col else []) + ([length_median_col] if length_median_col else [])].rename(columns=rename_dict).copy()
        f2["region"] = f2["region"].astype(str).str.strip()
        f2["reads"]  = pd.to_numeric(f2["reads"], errors="coerce").fillna(0.0)
        if "region_count" in f2.columns:
            f2["region_count"] = pd.to_numeric(f2["region_count"], errors="coerce").fillna(0.0)
        if length_total_col:
            f2["region_length_total_bp"] = pd.to_numeric(f[length_total_col], errors="coerce").fillna(0.0)
        if length_median_col:
            f2["region_length_median_bp"] = pd.to_numeric(f[length_median_col], errors="coerce").fillna(0.0)
        
        agg_dict = {"reads": "sum"}
        if "region_count" in f2.columns: agg_dict["region_count"] = "sum"
        if "region_length_total_bp" in f2.columns: agg_dict["region_length_total_bp"] = "sum"
        if "region_length_median_bp" in f2.columns: agg_dict["region_length_median_bp"] = "median"
        f2 = f2.groupby("region", as_index=False).agg(agg_dict).sort_values("reads", ascending=False)
        
        # Separate unlocalized reads from functional regions
        for _, row in f2.iterrows():
            if re.search(unloc_pattern, str(row["region"]), re.I):
                unlocalized_reads += float(row["reads"])
        
        # Filter out unlocalized from region_counts_df (keep for regions_list output)
        region_counts_df = f2
        regions_list = []
        for _, row in f2.iterrows():
            region_dict = {"region": row["region"], "reads": float(row["reads"])}
            if "region_count" in row:
                region_dict["region_count"] = int(row["region_count"])
            if "region_length_total_bp" in row:
                region_dict["region_length_total_bp"] = float(row["region_length_total_bp"])
            if "region_length_median_bp" in row:
                region_dict["region_length_median_bp"] = float(row["region_length_median_bp"])
            regions_list.append(region_dict)

# ── summary stats ─────────────────────────────────────────────────────────
divergent_regions = int(len(bed_div)) if not bed_div.empty else 0
# Exclude non-localized regions from total count
if not region_counts_df.empty:
    mask = ~region_counts_df["region"].str.contains(unloc_pattern, case=False, na=False, regex=True)
    total_regions = int(region_counts_df[mask]["region"].nunique())
    reads_total_func = float(pd.to_numeric(region_counts_df[mask]["reads"], errors="coerce").fillna(0).sum())
else:
    total_regions = 0
    reads_total_func = 0.0

median_pausing = float("nan")
if not pausing.empty:
    low = [c.lower() for c in pausing.columns]
    for cand in ("pi_len_norm","pausing_index","pi_raw","pi"):
        if cand in low:
            median_pausing = median_num(pausing[pausing.columns[low.index(cand)]])
            break

median_density = float("nan")
density_source = None
density_reason = None
if not density.empty:
    low = [c.lower() for c in density.columns]
    # Prefer explicit CPM column; otherwise fall back to 'signal'
    cpm_col = next((density.columns[i] for i,c in enumerate(low) if c=="cpm"), None)
    sig_col = next((density.columns[i] for i,c in enumerate(low) if c=="signal"), None)
    use_col = cpm_col or sig_col
    if use_col:
        median_density = median_num(density[use_col])
        density_source = "cpm" if cpm_col is not None else "signal"
    else:
        density_reason = "pol2_density.tsv missing CPM/signal column"
else:
    density_reason = "pol2_density.tsv not found or empty"

# Norm factors (accept synonyms)
cpm_factor    = norm_factor_lookup(norm_df, "cpm")
sicpm_factor  = norm_factor_lookup(norm_df, "sicpm")

# QC fields used for status + unlocalized fraction
input_reads       = qc.get("total_reads_raw") if isinstance(qc, dict) else None
dedup_reads       = qc.get("mapq_ge_reads_nodup") if isinstance(qc, dict) else None
# Duplicate percentage: use explicit None checks to handle 0.0 correctly
# (Python's 'or' treats 0.0 as falsy, which would skip valid zero values)
duplicate_percent = None
if isinstance(qc, dict):
    # Try different field names
    for field in ["duplicate_perc_of_total", "duplicate_percent", "dup_percent"]:
        if field in qc and qc[field] is not None:
            duplicate_percent = qc[field]
            break

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

# Unlocalized fraction = unlocalized_reads / (functional_reads + unlocalized_reads)
unlocalized_fraction = None
total_assigned_reads = reads_total_func + unlocalized_reads
if total_assigned_reads > 0:
    unlocalized_fraction = round(unlocalized_reads / total_assigned_reads, 6)
elif dedup_reads and dedup_reads > 0:
    # Fallback: if we don't have unlocalized reads tracked, use old calculation
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
        density_source=density_source,
        density_reason=density_reason,
    ),
    qc=dict(
        total_reads_raw=input_reads,
        dedup_reads_mapq_ge=dedup_reads,
        duplicate_percent=duplicate_percent,
        umi_deduplication_enabled=umi_enabled,
        umi_input_reads=umi_input_reads,
        umi_output_reads=umi_output_reads,
        umi_deduplication_percent=umi_dedup_percent
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

def html_density_table(density: pd.DataFrame) -> str:
    """Render CPM table, properly sorted by CPM/signal value"""
    if density is None or density.empty: return "<p class='muted'>n/a</p>"
    
    # Find signal/CPM column (try multiple names)
    signal_candidates = ["signal", "cpm", "density", "coverage"]
    signal_col = None
    for cand in signal_candidates:
        matching = [c for c in density.columns if c.lower()==cand.lower()]
        if matching:
            signal_col = matching[0]
            break
    
    if signal_col:
        # Filter out zero values and sort by signal descending, take top 25
        filtered_df = density[pd.to_numeric(density[signal_col], errors='coerce') > 0]
        sorted_df = filtered_df.sort_values(signal_col, ascending=False).head(25)
        # Show relevant columns
        cols_to_show = ["chr","start","end","name",signal_col,"norm_method"]
        cols_to_show = [c for c in cols_to_show if c in sorted_df.columns]
        return html_table(sorted_df, cols=cols_to_show, max_rows=25)
    else:
        return html_table(density.head(25), max_rows=25)

def html_pausing_table(pausing: pd.DataFrame) -> str:
    """Render pausing table, properly sorted by PI value"""
    if pausing is None or pausing.empty: return "<p class='muted'>n/a</p>"
    
    # Find pausing index column (try multiple names)
    pi_candidates = ["pausing_index", "pi_len_norm", "pi_raw", "pi"]
    pi_col = None
    for cand in pi_candidates:
        matching = [c for c in pausing.columns if c.lower()==cand.lower()]
        if matching:
            pi_col = matching[0]
            break
    
    if pi_col:
        # Filter out nan/inf values and sort by PI descending, take top 25
        filtered_df = pausing[pd.to_numeric(pausing[pi_col], errors='coerce').notna()]
        filtered_df = filtered_df[~filtered_df[pi_col].isin([float('inf'), float('-inf')])]
        if not filtered_df.empty:
            sorted_df = filtered_df.sort_values(pi_col, ascending=False).head(25)
            return html_table(sorted_df, max_rows=25)
    
    # Fallback: just show first 25 if we can't find valid PI data
    return html_table(pausing.head(25), max_rows=25)

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
    # Filter out non-localized polymerase for cleaner visualization
    unloc_pattern = r"non[-\s_]?localized|unlocalized|unlocalised"
    mask = ~reg_df["region"].str.contains(unloc_pattern, case=False, na=False, regex=True)
    plot_df = reg_df[mask] if mask.any() else reg_df
    if plot_df.empty: return "<div class='muted'>n/a</div>"
    
    # Vertical bar chart instead of horizontal stacked  
    fig, ax = plt.subplots(figsize=(8,5), dpi=120)
    vals = pd.to_numeric(plot_df["reads"], errors="coerce").fillna(0).values
    labels = plot_df["region"].astype(str).tolist()
    total = float(vals.sum()) or 1.0
    
    # Define color map (canonical from functional_regions.py RGB values)
    color_map = {
        'promoter': '#f38400', 'activepromoter': '#f38400', 'pppol': '#f38400',
        'divergenttx': '#b23bd4', 'divtx': '#b23bd4', 'divergent': '#b23bd4', 'ppdiv': '#b23bd4',
        'enhancers': '#73d47a', 'enhancer': '#73d47a', 'enh': '#73d47a',
        'genebody': '#000000', 'gene body': '#000000', 'body': '#000000', 'gb': '#000000',
        'cps': '#67c8f9', 'cleavagepolyadenylation': '#67c8f9',
        'terminationwindow': '#ff3662', 'termination window': '#ff3662', 'termination': '#ff3662', 'tw': '#ff3662'
    }
    
    colors = []
    for lab in labels:
        key = str(lab).lower().replace(' ', '').replace('_', '').replace('-', '')
        colors.append(color_map.get(key, '#3b82f6'))
    
    # Create regular bar chart with percentages
    x_pos = np.arange(len(labels))
    percentages = 100 * vals / total
    
    bars = ax.bar(x_pos, vals, color=colors, edgecolor='white', linewidth=1.2, alpha=0.9)
    
    # Add percentage labels on top of bars
    for i, (bar, pct) in enumerate(zip(bars, percentages)):
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width()/2., height,
                   f'{pct:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=12)
    ax.set_ylabel('Number of reads', fontsize=14, fontweight='bold')
    ax.set_title("Functional Region Composition", fontsize=16, fontweight='bold', pad=25)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(axis='both', which='major', labelsize=11)
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: f'{int(x):,}'))
    
    plt.tight_layout()
    return img_html(fig)

def pie_regions(reg_df: pd.DataFrame) -> str:
    if reg_df is None or reg_df.empty: return "<div class='muted'>n/a</div>"
    # Filter out non-localized polymerase for cleaner visualization
    unloc_pattern = r"non[-\s_]?localized|unlocalized|unlocalised"
    mask = ~reg_df["region"].str.contains(unloc_pattern, case=False, na=False, regex=True)
    plot_df = reg_df[mask] if mask.any() else reg_df
    if plot_df.empty: return "<div class='muted'>n/a</div>"
    
    vals = pd.to_numeric(plot_df["reads"], errors="coerce").fillna(0).values
    labels = plot_df["region"].astype(str).tolist()
    if vals.sum() <= 0: return "<div class='muted'>n/a</div>"
    
    # Define color map (canonical from functional_regions.py RGB values)
    color_map = {
        'promoter': '#f38400', 'activepromoter': '#f38400', 'pppol': '#f38400',
        'divergenttx': '#b23bd4', 'divtx': '#b23bd4', 'divergent': '#b23bd4', 'ppdiv': '#b23bd4',
        'enhancers': '#73d47a', 'enhancer': '#73d47a', 'enh': '#73d47a',
        'genebody': '#000000', 'gene body': '#000000', 'body': '#000000', 'gb': '#000000',
        'cps': '#67c8f9', 'cleavagepolyadenylation': '#67c8f9',
        'terminationwindow': '#ff3662', 'termination window': '#ff3662', 'termination': '#ff3662', 'tw': '#ff3662'
    }
    
    colors = []
    for lab in labels:
        key = str(lab).lower().replace(' ', '').replace('_', '').replace('-', '')
        colors.append(color_map.get(key, '#3b82f6'))
    
    fig, ax = plt.subplots(figsize=(6,6), dpi=120)
    
    # NO labels on slices - only in legend to avoid overlap
    wedges, texts = ax.pie(vals, labels=None, colors=colors, startangle=90)
    
    ax.set_title("Functional Regions (Reads)", fontsize=14, fontweight='bold', pad=20)
    
    # Add comprehensive legend with percentages outside the pie
    total = float(vals.sum())
    legend_labels = [f"{lab}\n{v:,.0f} reads ({100*v/total:.1f}%)" 
                    for lab, v in zip(labels, vals)]
    ax.legend(wedges, legend_labels, 
              loc="center left", bbox_to_anchor=(1.05, 0.5),
              fontsize=10, frameon=True, fancybox=True, shadow=True)
    
    plt.tight_layout()
    return img_html(fig)

def hist_cpm(density: pd.DataFrame) -> str:
    if density.empty or not any(c.lower()=="cpm" for c in density.columns): return "<div class='muted'>n/a</div>"
    cpm_col = [c for c in density.columns if c.lower()=="cpm"][0]
    x = pd.to_numeric(density[cpm_col], errors="coerce").replace([np.inf,-np.inf], np.nan).dropna().values
    if x.size == 0: return "<div class='muted'>n/a</div>"
    
    fig, ax = plt.subplots(figsize=(6,4), dpi=110)
    n, bins, patches = ax.hist(np.log10(x+1e-9), bins=50, color='#3b82f6', 
                                alpha=0.7, edgecolor='black', linewidth=0.5)
    
    ax.set_xlabel("log₁₀(CPM)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Number of regions", fontsize=12, fontweight='bold')
    ax.set_title("Functional-region CPM distribution", fontsize=13, fontweight='bold', pad=15)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add summary statistics
    median_val = np.median(x)
    mean_val = np.mean(x)
    ax.text(0.98, 0.97, f'Median: {median_val:.1f}\nMean: {mean_val:.1f}\nN regions: {len(x):,}',
            transform=ax.transAxes, fontsize=9, verticalalignment='top',
            horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    return img_html(fig)

def hist_pausing(pausing: pd.DataFrame) -> str:
    if pausing.empty: return "<div class='muted'>n/a</div>"
    low = [c.lower() for c in pausing.columns]
    col = next((pausing.columns[low.index(c)] for c in ("pi_len_norm","pausing_index","pi_raw","pi") if c in low), None)
    if not col: return "<div class='muted'>n/a</div>"
    vals = pd.to_numeric(pausing[col], errors="coerce").replace([np.inf,-np.inf], np.nan).dropna().values
    if vals.size == 0: return "<div class='muted'>n/a</div>"
    
    fig, ax = plt.subplots(figsize=(6,4), dpi=110)
    n, bins, patches = ax.hist(np.log2(vals+1e-9), bins=60, color='#10b981', 
                                alpha=0.7, edgecolor='black', linewidth=0.5)
    
    ax.set_xlabel(f"log₂({col} + 10⁻⁹)", fontsize=12, fontweight='bold')
    ax.set_ylabel("Number of genes", fontsize=12, fontweight='bold')
    ax.set_title(f"{col} distribution", fontsize=13, fontweight='bold', pad=15)
    ax.tick_params(axis='both', which='major', labelsize=10)
    ax.grid(axis='y', alpha=0.3, linestyle='--', linewidth=0.5)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add summary statistics
    median_val = np.median(vals)
    mean_val = np.mean(vals)
    ax.text(0.98, 0.97, f'Median: {median_val:.2f}\nMean: {mean_val:.2f}\nN genes: {len(vals):,}',
            transform=ax.transAxes, fontsize=9, verticalalignment='top',
            horizontalalignment='right', bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    plt.tight_layout()
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
.grid{{display:grid;gap:16px;grid-template-columns:repeat(2,1fr)}}
@media (max-width:768px){{.grid{{grid-template-columns:1fr}}}}
.card{{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:12px}}
.kpi{{background:var(--bg);border:1px solid var(--line);border-radius:10px;padding:12px}}
.kpi-label{{font-size:12px;color:var(--muted)}}
.kpi-value{{font-size:18px;font-weight:600;margin-top:2px}}
.kpi-sub{{font-size:11px;color:var(--muted)}}
.tbl{{border-collapse:collapse;width:100%;font-size:0.9rem}}
.tbl th,.tbl td{{border-bottom:1px solid var(--line);padding:0.75rem;text-align:left}}
.tbl th{{background:var(--bg);font-weight:600;position:sticky;top:0;z-index:10}}
.tbl tbody tr:hover{{background:var(--chip)}}
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
    {kpi("Median density", (None if np.isnan(median_density) else round(median_density,3)), ("source: "+str(row['metrics'].get('density_source') or 'n/a') + (" — "+str(row['metrics'].get('density_reason')) if (np.isnan(median_density) and row['metrics'].get('density_reason')) else "")))}
    {kpi("CPM factor", row["metrics"]["cpm_factor"])}
    {kpi("siCPM factor", row["metrics"]["sicpm_factor"])}
  </div>

  <h2>Functional regions — composition</h2>
  <div style="display:grid;gap:16px;grid-template-columns:1fr 1fr;">
    <div class="card">{stacked_bar_regions(region_counts_df)}</div>
    <div class="card">{pie_regions(region_counts_df)}</div>
  </div>

  <h2>📋 Functional Region Summary</h2>
  <p class="muted" style="margin:-8px 0 16px 0;">Aggregate read counts and region statistics across all functional categories</p>
  <div class="card">
    <div style="overflow-x:auto;">
      {html_table(region_counts_df.rename(columns={
          "region":"Region",
          "reads":"Reads",
          "region_count":"Region Count",
          "region_length_total_bp":"Total Length (bp)",
          "region_length_median_bp":"Median Length (bp)"
        }))}
    </div>
  </div>

  <h2>📊 Top Regions & Genes</h2>
  <p class="muted" style="margin:-8px 0 16px 0;">Highest signal density regions and strongest pausing genes</p>
  <div style="display:grid;gap:16px;grid-template-columns:1fr 1fr;">
    <div class="card">
      <h3 style="margin:0 0 8px">Top Functional Regions by Signal</h3>
      <p class="muted" style="font-size:0.85em;margin-bottom:12px;">Regions with highest normalized signal (CPM). Only non-zero regions shown.</p>
      <div style="overflow-x:auto;">
        {html_density_table(density)}
      </div>
    </div>
    <div class="card">
      <h3 style="margin:0 0 8px">Top Genes by Pausing Index</h3>
      <p class="muted" style="font-size:0.85em;margin-bottom:12px;">Genes ranked by Pol II pausing (promoter-proximal vs gene body ratio).</p>
      <div style="overflow-x:auto;">
        {html_pausing_table(pausing)}
      </div>
    </div>
  </div>

  <h2>Distributions</h2>
  <div style="display:grid;gap:16px;grid-template-columns:1fr 1fr;">
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
    page.write(f"""<!doctype html><meta charset='utf-8'>
<title>{SID} plots</title>
<style>
body{{font-family:system-ui,sans-serif;max-width:1200px;margin:0 auto;padding:20px;background:#f9fafb}}
h1{{color:#1f2937;border-bottom:3px solid #3b82f6;padding-bottom:10px}}
h3{{color:#374151;margin-top:40px}}
.plot-container{{background:white;padding:20px;margin:20px 0;border-radius:8px;box-shadow:0 1px 3px rgba(0,0,0,0.1)}}
</style>
<h1>{SID} — Detailed Plots</h1>""")

    if not region_counts_df.empty and region_counts_df["reads"].sum() > 0:
        page.write("<div class='plot-container'>")
        page.write(pie_regions(region_counts_df))
        page.write("</div>")

    if not density.empty and any(c.lower()=="cpm" for c in density.columns):
        page.write("<div class='plot-container'>")
        page.write(hist_cpm(density))
        page.write("</div>")

    if not pausing.empty:
        page.write("<div class='plot-container'>")
        page.write(hist_pausing(pausing))
        page.write("</div>")

    with open(args.out_plots_html, "w", encoding="utf-8") as fh:
        fh.write(page.getvalue())
print(f"[render_py] done ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)
