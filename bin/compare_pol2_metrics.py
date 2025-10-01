#!/usr/bin/env python3
# =============================================================================
# compare_pol2_metrics.py — merge per-sample metrics, contrasts, plots (no HDF5)
# -----------------------------------------------------------------------------
# Inputs
#   --samples-tsv : sample_id,condition,timepoint,replicate,file
#   --contrasts   : items like "condition:KI,WT" or "timepoint:90,0" (optional)
# Outputs
#   --out-merged (tidy TSV, always)
#   --out-contrasts (optional TSV)
#   --plots-dir (optional PNGs)
# =============================================================================

from __future__ import annotations
import argparse, os, sys, math, datetime
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pandas.api.types import is_numeric_dtype

REQUIRED_SHEET_COLS = ["sample_id","condition","timepoint","replicate","file"]
NEEDED_METRICS      = ["gene_id","gene_name","pi_len_norm","pi_raw","body_cpm","tss_cpm"]

# ── CLI ────────────────────────────────────────────────────────────────────
ap = argparse.ArgumentParser()
ap.add_argument("--samples-tsv", required=True)
ap.add_argument("--out-merged",  required=True)
ap.add_argument("--contrasts", nargs="*", default=[])
ap.add_argument("--out-contrasts", default=None)
ap.add_argument("--plots-dir", default=None)
ap.add_argument("--top-n", type=int, default=100)
args = ap.parse_args()
print(f"[aggregate_py] start ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)

def eprint(*a): print(*a, file=sys.stderr, flush=True)

# ── Load sample sheet ──────────────────────────────────────────────────────
samples = pd.read_csv(args.samples_tsv, sep="\t", dtype=str)
for col in REQUIRED_SHEET_COLS:
    if col not in samples.columns:
        eprint(f"[aggregate] ERROR: samples.tsv missing column: {col}")
        sys.exit(2)

# drop duplicates conservatively
samples = samples.drop_duplicates(subset=REQUIRED_SHEET_COLS, keep="first")

missing = samples.loc[~samples["file"].apply(lambda p: os.path.exists(str(p)))]
if not missing.empty:
    eprint("[aggregate] WARN: some per-sample TSV paths do not exist; they will be skipped:")
    for _, r in missing.iterrows():
        eprint(f"  {r['sample_id']}: {r['file']}")
    samples = samples.loc[samples["file"].apply(lambda p: os.path.exists(str(p)))]

if samples.empty:
    eprint("[aggregate] ERROR: no valid per-sample TSV paths."); sys.exit(3)

# ── Load per-sample tables (usecols) ───────────────────────────────────────
rows = []
missing_cols_report = {}
for _, r in samples.iterrows():
    path = str(r["file"])
    try:
        df = pd.read_csv(path, sep="\t", usecols=lambda c: (c in NEEDED_METRICS), dtype="object")
    except Exception as e:
        eprint(f"[aggregate] WARN: failed to read {path}: {e}")
        continue

    present = set(df.columns)
    miss = [c for c in NEEDED_METRICS if c not in present]
    if miss:
        missing_cols_report[r["sample_id"]] = miss
        # fill missing metrics as NA to keep shape
        for c in miss: df[c] = np.nan

    df = df[["gene_id","gene_name","pi_len_norm","pi_raw","body_cpm","tss_cpm"]]
    df["sample_id"] = r["sample_id"]
    df["condition"] = r["condition"]
    df["timepoint"] = r["timepoint"]
    df["replicate"] = r["replicate"]
    rows.append(df)

if missing_cols_report:
    eprint("[aggregate] WARN: missing metric columns detected:")
    for sid, miss in missing_cols_report.items():
        eprint(f"  {sid}: missing {','.join(miss)}")

merged = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(
    columns=["gene_id","gene_name","pi_len_norm","pi_raw","body_cpm","tss_cpm",
             "sample_id","condition","timepoint","replicate"]
)

# coerce numeric metrics
for c in ["pi_len_norm","pi_raw","body_cpm","tss_cpm"]:
    if c in merged.columns and not is_numeric_dtype(merged[c]):
        merged[c] = pd.to_numeric(merged[c], errors="coerce")

# ── Save merged ────────────────────────────────────────────────────────────
merged.to_csv(args.out_merged, sep="\t", index=False)
eprint(f"[aggregate] merged rows={len(merged):,} samples={merged['sample_id'].nunique()} genes={merged['gene_id'].nunique()}")

# ── Contrasts ──────────────────────────────────────────────────────────────
def parse_contrast(s: str):
    # "condition:KI,WT"  → ("condition","KI","WT")
    k, v = s.split(":", 1)
    a, b = v.split(",", 1)
    return (k.strip(), a.strip(), b.strip())

def safe_log2fc(a: pd.Series, b: pd.Series) -> pd.Series:
    return np.log2((a.astype(float) + 1e-9) / (b.astype(float) + 1e-9))

def compute_contrasts(mrg: pd.DataFrame, cons):
    if mrg.empty: return None
    x = mrg.copy()
    for col in ["pi_len_norm","pi_raw","body_cpm","tss_cpm"]:
        if col in x.columns and not is_numeric_dtype(x[col]):
            x[col] = pd.to_numeric(x[col], errors="coerce")

    # replicate collapse: median per (gene, condition, timepoint)
    med = x.groupby(["gene_id","gene_name","condition","timepoint"], dropna=False)\
           .median(numeric_only=True).reset_index()

    outs = []
    for key, A, B in cons:
        if key not in {"condition","timepoint"}:
            eprint(f"[aggregate] WARN: unsupported contrast key: {key} (skip)")
            continue

        for metric in ["pi_len_norm","body_cpm"]:
            if key == "condition":
                a = med.query("condition == @A")[["gene_id","gene_name","timepoint",metric]].rename(columns={metric:"A"})
                b = med.query("condition == @B")[["gene_id","gene_name","timepoint",metric]].rename(columns={metric:"B"})
                jo = a.merge(b, on=["gene_id","gene_name","timepoint"], how="inner")
                if jo.empty:
                    eprint(f"[aggregate] WARN: no overlap for contrast {key}:{A},{B} on {metric}")
                    continue
                jo["group_by"] = "timepoint"; jo["level"] = jo["timepoint"]
            else:  # timepoint
                a = med.query("timepoint == @A")[["gene_id","gene_name","condition",metric]].rename(columns={metric:"A"})
                b = med.query("timepoint == @B")[["gene_id","gene_name","condition",metric]].rename(columns={metric:"B"})
                jo = a.merge(b, on=["gene_id","gene_name","condition"], how="inner")
                if jo.empty:
                    eprint(f"[aggregate] WARN: no overlap for contrast {key}:{A},{B} on {metric}")
                    continue
                jo["group_by"] = "condition"; jo["level"] = jo["condition"]

            jo = jo.dropna(subset=["A","B"])
            if jo.empty: 
                eprint(f"[aggregate] WARN: contrast {key}:{A},{B} on {metric} → all NA after join")
                continue

            jo["contrast"] = f"{key}:{A}_vs_{B}"
            jo["metric"]   = metric
            jo["log2FC"]   = safe_log2fc(jo["A"], jo["B"])
            outs.append(jo[["gene_id","gene_name","group_by","level","contrast","metric","A","B","log2FC"]])

    return None if not outs else pd.concat(outs, ignore_index=True)

contrasts_df = None
if args.contrasts and args.out_contrasts:
    # parse & validate specs
    specs = []
    for c in args.contrasts:
        try:
            specs.append(parse_contrast(c))
        except Exception as e:
            eprint(f"[aggregate] WARN: could not parse contrast '{c}': {e}")
    if specs:
        contrasts_df = compute_contrasts(merged, specs)
        if contrasts_df is not None and not contrasts_df.empty:
            contrasts_df.to_csv(args.out_contrasts, sep="\t", index=False)
            eprint(f"[aggregate] contrasts rows={len(contrasts_df):,}")
        else:
            eprint("[aggregate] INFO: no contrasts produced (empty overlap or bad specs).")

# ── Plots (optional) ───────────────────────────────────────────────────────
def ma_plot(df: pd.DataFrame, title: str, out_png: str):
    if df is None or df.empty: return
    x = 0.5 * np.log2((df["A"].astype(float) + 1e-9) * (df["B"].astype(float) + 1e-9))
    y = np.log2((df["A"].astype(float) + 1e-9) / (df["B"].astype(float) + 1e-9))
    plt.figure(figsize=(6,5), dpi=130)
    plt.scatter(x, y, s=6, alpha=0.6)
    plt.axhline(0, linewidth=1)
    plt.xlabel("mean expression (log2)")
    plt.ylabel("log2 fold-change")
    plt.title(title)
    plt.tight_layout(); plt.savefig(out_png); plt.close()

def heatmap_top_var(df: pd.DataFrame, value_col: str, group_col: str, out_png: str, top_n: int = 100):
    if df is None or df.empty: return
    pv = df.pivot_table(index=["gene_id","gene_name"], columns=group_col, values=value_col, aggfunc="median")
    pv = pv.replace([np.inf, -np.inf], np.nan).dropna()
    if pv.empty: return
    # guard for very wide plots
    cols = list(pv.columns)
    if len(cols) == 0: return
    var = pv.var(axis=1, numeric_only=True).sort_values(ascending=False)
    top = pv.loc[var.index[:top_n]]
    if top.empty: return
    zs  = (top - top.mean(axis=1).values.reshape(-1,1)) / (top.std(axis=1, ddof=0).values.reshape(-1,1) + 1e-9)
    plt.figure(figsize=(max(6, zs.shape[1]*0.25), max(6, zs.shape[0]*0.06)), dpi=130)
    plt.imshow(zs.values, aspect='auto', interpolation='nearest')
    plt.xticks(range(zs.shape[1]), zs.columns, rotation=90)
    plt.yticks([])
    plt.colorbar(fraction=0.046, pad=0.04)
    plt.title(f"Top {min(top_n, zs.shape[0])} variable genes — {value_col} (z-score)")
    plt.tight_layout(); plt.savefig(out_png); plt.close()

if args.plots_dir:
    os.makedirs(args.plots_dir, exist_ok=True)
    if not merged.empty:
        for metric in ["pi_len_norm","body_cpm"]:
            for group in ["condition","timepoint"]:
                try:
                    heatmap_top_var(merged, metric, group,
                                    os.path.join(args.plots_dir, f"heatmap_{metric}_by_{group}.png"),
                                    top_n=max(10, min(args.top_n, 500)))
                except Exception as e:
                    eprint(f"[aggregate] WARN: heatmap failed: {e}")
    if contrasts_df is not None and not contrasts_df.empty:
        for (contrast, metric, level), sub in contrasts_df.groupby(["contrast","metric","level"]):
            try:
                ma_plot(sub, title=f"{contrast} @ {level} — {metric}",
                        out_png=os.path.join(args.plots_dir, f"MA_{metric}_{contrast}_@{level}.png"))
            except Exception as e:
                eprint(f"[aggregate] WARN: MA plot failed: {e}")
print(f"[aggregate_py] done ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)
