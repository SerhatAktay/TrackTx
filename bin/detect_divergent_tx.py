#!/usr/bin/env python3
# =============================================================================
# detect_divergent_tx.py — fast divergent transcription caller (v2.5)
# -----------------------------------------------------------------------------
# Inputs (bedGraph, 4-col: chr start end value; plain or .gz)
#   --pos : positive-strand 3′ signal (treated as |sig|)
#   --neg : negative-strand 3′ signal (treated as -|sig| by default; or |sig| with --neg-abs)
#
# Steps
#   1) Read & sanitize bedGraphs: enforce 4 cols, start<end, enforce signs
#      (skip sorting if --assume-sorted)
#   2) Optional smoothing (Gaussian SD in bp)
#   3) Call blocks on + and − via per-bin threshold + merge (bin_gap)
#   4) Pair +/− blocks per chromosome within ±nt_window (vectorized filters)
#   5) Optional “balance”, “overlap_bp”, and “valley_thr” constraints
#   6) Merge paired windows within merge_gap; write BED4: chr start end total
#   7) Optional QC PDF; optional summary TSV with counts/timings
# =============================================================================

from __future__ import annotations
import argparse, sys, os, time, math, pathlib, datetime
from typing import Dict, List, Tuple
import numpy as np
import pandas as pd
from multiprocessing import Pool, cpu_count

# keep hidden threading modest
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# ─────────────── CLI ───────────────
ap = argparse.ArgumentParser()
ap.add_argument("--sample", required=True)
ap.add_argument("--pos", required=True)
ap.add_argument("--neg", required=True)
ap.add_argument("--nt-window", type=int, default=800)
ap.add_argument("--threshold", type=float, default=1.0)
ap.add_argument("--sum-thr",  type=float, default=5.0)
ap.add_argument("--bin-gap",  type=int, default=100)
ap.add_argument("--merge-gap",type=int, default=50)
ap.add_argument("--balance",  type=float, default=0.0)

# Preferred:
ap.add_argument("--overlap-bp", type=int, default=0)
# Back-compat alias:
ap.add_argument("--valley-gap", type=int, default=None)

ap.add_argument("--valley-thr",  default="Inf")
ap.add_argument("--smooth-sd", type=int, default=0)
ap.add_argument("--max-dt-width", type=int, default=0)  # 0 = off
ap.add_argument("--ncores", type=int, default=1)
ap.add_argument("--qc", type=int, default=0)
ap.add_argument("--out-bed", required=True)

# New UX flags
ap.add_argument("--assume-sorted", action="store_true",
                help="Trust input sortedness; skip pandas sort")
ap.add_argument("--neg-abs", action="store_true",
                help="Treat --neg as magnitudes (|NEG|) instead of -|NEG|")
ap.add_argument("--write-summary", default=None,
                help="Path to a TSV with counts/timings (sample, n_pos_pk, n_neg_pk, n_pairs_raw, n_dt, wall_s)")

# Advanced
ap.add_argument("--max-neg-per-pos", type=int, default=64)
ap.add_argument("--report-every", type=int, default=20000)
ap.add_argument("--quiet", action="store_true")
args = ap.parse_args()
sys.stderr.write(f"[divergent_py] start ts={datetime.datetime.utcnow().isoformat()}Z\n")

# ───────── normalize / clamp args ─────────
def clamp_int(x, lo, hi=None):
    x = int(x); x = max(lo, x)
    return x if hi is None else min(hi, x)

def clamp_float(x, lo, hi=None):
    x = float(x); x = max(lo, x)
    return x if hi is None else min(hi, x)

THR       = clamp_float(args.threshold, 0.0)
SUMTHR    = clamp_float(args.sum_thr,   0.0)
BIN_GAP   = clamp_int(args.bin_gap,     0)
MERGE_GAP = clamp_int(args.merge_gap,   0)
NTWIN     = clamp_int(args.nt_window,   1)
BAL       = clamp_float(args.balance,   0.0, 1.0)
SMOOTH_SD = clamp_int(args.smooth_sd,   0)
NCORES    = clamp_int(args.ncores,      1)
MAXNEG    = clamp_int(args.max_neg_per_pos, 0)
REVERY    = clamp_int(args.report_every,    0)
MAXW      = clamp_int(args.max_dt_width,    0)

OVERLAP_BP = clamp_int(args.overlap_bp if args.overlap_bp is not None else 0, 0)
if args.valley_gap is not None:
    OVERLAP_BP = max(OVERLAP_BP, clamp_int(args.valley_gap, 0))

try:
    VALLEY_THR = float(args.valley_thr) if str(args.valley_thr).lower() != "inf" else float("inf")
    VALLEY_THR = clamp_float(VALLEY_THR, 0.0)
except Exception:
    VALLEY_THR = float("inf")

def log(msg: str):
    if not args.quiet:
        sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n"); sys.stderr.flush()

# ───────── I/O & sanitize ─────────
def read_bg_sane(path: str, role: str, assume_sorted: bool) -> pd.DataFrame:
    # Prefer pyarrow engine if present (faster/lower RAM), else C
    engine = "c"
    try:
        import pyarrow  # type: ignore
        engine = "pyarrow"
    except Exception:
        pass

    df = pd.read_csv(
        path,
        sep="\t", header=None, comment="#",
        names=["chr","start","end","sig"],
        usecols=[0,1,2,3],
        dtype={"chr":"category","start":"int64","end":"int64","sig":"float64"},
        compression="infer",
        engine=engine,
    )
    if df.empty:
        return df
    df = df.loc[(df["end"] > df["start"]) & df["chr"].notna()].copy()
    if role == "pos":
        df["sig"] = df["sig"].abs()
    elif role == "neg":
        df["sig"] = (df["sig"].abs() if args.neg_abs else -df["sig"].abs())
    else:
        raise ValueError("role must be 'pos' or 'neg'")
    if not assume_sorted:
        df.sort_values(["chr","start"], inplace=True, kind="mergesort")
    df.reset_index(drop=True, inplace=True)
    return df

def gaussian_smooth_inplace(df: pd.DataFrame, sd: int):
    if sd <= 0 or df.empty:
        return
    try:
        from scipy.ndimage import gaussian_filter1d as g1d  # type: ignore
        smoother = lambda v: g1d(v, sd, mode="constant")
    except Exception:
        half = max(1, int(5*sd))
        xs = np.arange(-half, half+1)
        ker = np.exp(-(xs*xs)/(2.0*sd*sd)); ker /= ker.sum()
        smoother = lambda v: np.convolve(v, ker, mode="same")
    for _, idx in df.groupby("chr", sort=False, observed=False).groups.items():
        v = df.loc[idx, "sig"].to_numpy()
        df.loc[idx, "sig"] = smoother(v)

# ───────── block calling (vectorized) ─────────
def call_blocks(df: pd.DataFrame, thr: float, sum_thr: float, bin_gap: int, anchor_mode: str = "adjacent") -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["chr","start","end","height","signal"])
    v = np.abs(df["sig"].to_numpy())
    keep = v >= thr
    if not keep.any():
        return pd.DataFrame(columns=["chr","start","end","height","signal"])
    chr_a = df["chr"].to_numpy()[keep]
    s     = df["start"].to_numpy()[keep]
    e     = df["end"].to_numpy()[keep]
    v     = v[keep]

    # Grouping strategies:
    # - adjacent (legacy): compare to previous interval end
    # - left:     left-anchored grouping by the first bin start within <= bin_gap of anchor
    # - right:    right-anchored grouping by the first bin end (process in reverse)
    if anchor_mode not in ("adjacent","left","right"):
        anchor_mode = "adjacent"

    if anchor_mode == "adjacent":
        grp = np.concatenate([[True], (chr_a[1:] != chr_a[:-1]) | (s[1:] > e[:-1] + bin_gap)]).cumsum()
    else:
        grp = np.empty(len(s), dtype=np.int64)
        gid = 0
        i = 0
        n = len(s)
        # Process per chromosome
        # Build indices per chromosome to avoid crossing boundaries
        # chr_a is categorical/strings; iterate spans of equal chr
        start_idx = 0
        while start_idx < n:
            cur_chr = chr_a[start_idx]
            end_idx = start_idx
            while end_idx < n and chr_a[end_idx] == cur_chr:
                end_idx += 1
            # Slice for this chromosome
            idx_range = np.arange(start_idx, end_idx)
            if anchor_mode == "left":
                anchor = s[idx_range[0]]
                for j in idx_range:
                    if s[j] > anchor + bin_gap:
                        gid += 1
                        anchor = s[j]
                    grp[j] = gid
                gid += 1
            else:  # right-anchored
                anchor = e[idx_range[-1]]
                for j in idx_range[::-1]:
                    if e[j] < anchor - bin_gap:
                        gid += 1
                        anchor = e[j]
                    grp[j] = gid
                gid += 1
            start_idx = end_idx

    s_ser = pd.Series(s); e_ser = pd.Series(e); v_ser = pd.Series(v); chr_ser = pd.Series(chr_a)
    out = pd.DataFrame({
        "chr":   chr_ser.groupby(grp).first(),
        "start": s_ser.groupby(grp).min(),
        "end":   e_ser.groupby(grp).max(),
        "height": v_ser.groupby(grp).max(),
        "signal": v_ser.groupby(grp).sum(),
    }).reset_index(drop=True)
    return out.loc[out["signal"] >= sum_thr].reset_index(drop=True)

def build_chr_index(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    if df.empty: return {}
    return {str(c): sub.sort_values("start").reset_index(drop=True)
            for c, sub in df.groupby("chr", sort=False, observed=False)}

def valley_ok(idx_df: pd.DataFrame, a: int, b: int, thr: float) -> bool:
    if not math.isfinite(thr) or idx_df is None or idx_df.empty:
        return True
    if a > b: a, b = b, a
    starts = idx_df["start"].to_numpy()
    ends   = idx_df["end"].to_numpy()
    sig    = np.abs(idx_df["sig"].to_numpy())
    i = np.searchsorted(ends, a, side="right")
    n = len(starts)
    while i < n and starts[i] <= b:
        if (ends[i] > a) and (starts[i] < b) and (sig[i] > thr):
            return False
        i += 1
    return True

# ───────── pairing per chromosome ─────────
def pair_on_chr_fast(task) -> pd.DataFrame:
    (chrom, p_df, n_df, win, bal, overlap_bp, vthr, pidx, nidx,
     max_neg, report_every) = task
    if p_df.empty or n_df.empty:
        return pd.DataFrame(columns=["chr","start","end","total","dist","balv"])
    ps, pe = p_df["start"].to_numpy(), p_df["end"].to_numpy()
    ph, pg = p_df["height"].to_numpy(), p_df["signal"].to_numpy()
    ns, ne = n_df["start"].to_numpy(), n_df["end"].to_numpy()
    nh, ng = n_df["height"].to_numpy(), n_df["signal"].to_numpy()
    out_rows: List[Tuple[str,int,int,float,int,float]] = []
    for i in range(len(ps)):
        if report_every and (i % report_every == 0):
            log(f"{chrom}: +peaks {i}/{len(ps)}")
        left  = np.searchsorted(ns, ps[i] - win, side="left")
        right = np.searchsorted(ns, ps[i] + win, side="right")
        if right <= left: continue
        if max_neg and (right - left) > max_neg:
            mid = (left + right) // 2; half = max_neg // 2
            left = max(left, mid - half); right = left + max_neg
        ns_sl, ne_sl = ns[left:right], ne[left:right]
        nh_sl, ng_sl = nh[left:right], ng[left:right]
        keep = np.abs(ns_sl - ps[i]) <= win
        if not keep.any(): continue
        ns_sl, ne_sl = ns_sl[keep], ne_sl[keep]
        nh_sl, ng_sl = nh_sl[keep], ng_sl[keep]
        bal_mask = np.ones_like(nh_sl, dtype=bool)
        balv_all = None
        if bal > 0.0:
            balv = np.minimum(ph[i], nh_sl) / np.maximum(ph[i], nh_sl)
            bal_mask = balv >= bal
            if not bal_mask.any(): continue
            ns_sl, ne_sl = ns_sl[bal_mask], ne_sl[bal_mask]
            ng_sl = ng_sl[bal_mask]; balv_all = balv[bal_mask]
        else:
            balv_all = np.minimum(ph[i], nh_sl) / np.maximum(ph[i], nh_sl)
        if overlap_bp > 0:
            overlap = (np.minimum(pe[i], ne_sl) - np.maximum(ps[i], ns_sl))
            ov_mask = overlap >= overlap_bp
            if not ov_mask.any(): continue
            ns_sl, ne_sl = ns_sl[ov_mask], ne_sl[ov_mask]
            ng_sl = ng_sl[ov_mask]; balv_all = balv_all[ov_mask]
            overlap = overlap[ov_mask]
        else:
            overlap = (np.minimum(pe[i], ne_sl) - np.maximum(ps[i], ns_sl))
        if math.isfinite(vthr):
            keep_mask = np.ones(len(ns_sl), dtype=bool)
            for j in range(len(ns_sl)):
                a_end = int(min(pe[i], ne_sl[j])); b_sta = int(max(ps[i], ns_sl[j]))
                if not (valley_ok(pidx, a_end, b_sta, vthr) and valley_ok(nidx, a_end, b_sta, vthr)):
                    keep_mask[j] = False
            if not keep_mask.any(): continue
            ns_sl, ne_sl = ns_sl[keep_mask], ne_sl[keep_mask]
            ng_sl = ng_sl[keep_mask]; balv_all = balv_all[keep_mask]
            overlap = overlap[keep_mask]
        start = np.minimum(ps[i], ns_sl)
        end   = np.maximum(pe[i], ne_sl)
        total = pg[i] + ng_sl
        dist  = ((ns_sl + ne_sl)//2) - ((ps[i] + pe[i])//2)
        for s,e,t,d,bv in zip(start.tolist(), end.tolist(), total.tolist(), dist.tolist(), balv_all.tolist()):
            out_rows.append((chrom, int(s), int(e), float(t), int(abs(d)), float(bv)))
    if not out_rows:
        return pd.DataFrame(columns=["chr","start","end","total","dist","balv"])
    return pd.DataFrame(out_rows, columns=["chr","start","end","total","dist","balv"])

def merge_windows(dt: pd.DataFrame, gap: int) -> pd.DataFrame:
    if dt is None or dt.empty:
        return pd.DataFrame(columns=["chr","start","end","total"])
    dt = dt.sort_values(["chr","start","end"]).reset_index(drop=True)
    rows: List[Tuple[str,int,int,float]] = []
    cur_chr, cur_s, cur_e, cur_t = None, None, None, 0.0
    for r in dt.itertuples(index=False):
        if (cur_chr is None) or (r.chr != cur_chr) or (r.start > cur_e + gap):
            if cur_chr is not None:
                rows.append((cur_chr, cur_s, cur_e, cur_t))
            cur_chr, cur_s, cur_e, cur_t = r.chr, int(r.start), int(r.end), float(r.total)
        else:
            cur_e = max(cur_e, int(r.end))
            cur_t += float(r.total)
    if cur_chr is not None:
        rows.append((cur_chr, cur_s, cur_e, cur_t))
    return pd.DataFrame(rows, columns=["chr","start","end","total"])

# ───────── main ─────────
start_t = time.time()
def tsec(): return f"{time.time()-start_t:.1f}s"

log(f"reading: {args.sample}")
plus  = read_bg_sane(args.pos, "pos", assume_sorted=args.assume_sorted)
minus = read_bg_sane(args.neg, "neg", assume_sorted=args.assume_sorted)

if plus.empty or minus.empty:
    pathlib.Path(args.out_bed).write_text("")
    if args.write-summary:
        pathlib.Path(args.write_summary).write_text(
            "sample\tn_pos_pk\tn_neg_pk\tn_pairs_raw\tn_dt\twall_s\n"
            f"{args.sample}\t0\t0\t0\t0\t{time.time()-start_t:.3f}\n"
        )
    log(f"✓ {args.sample}: DT=0 (empty input) → {args.out_bed}  (wall={tsec()})")
    sys.exit(0)

if SMOOTH_SD > 0:
    log("smoothing…")
    gaussian_smooth_inplace(plus,  SMOOTH_SD)
    gaussian_smooth_inplace(minus, SMOOTH_SD)

log("calling peaks…")
pk_pos = call_blocks(plus,  THR, SUMTHR, BIN_GAP, anchor_mode="left")
pk_neg = call_blocks(minus, THR, SUMTHR, BIN_GAP, anchor_mode="right")
log(f"peaks: +{len(pk_pos)}  -{len(pk_neg)}")

chroms = sorted(set(map(str, pk_pos["chr"])) & set(map(str, pk_neg["chr"])))
pos_idx = build_chr_index(plus)
neg_idx = build_chr_index(minus)

work = [(c,
         pk_pos[pk_pos["chr"].astype(str)==c].reset_index(drop=True),
         pk_neg[pk_neg["chr"].astype(str)==c].reset_index(drop=True),
         NTWIN, BAL, OVERLAP_BP, VALLEY_THR,
         pos_idx.get(c, pd.DataFrame(columns=plus.columns)),
         neg_idx.get(c, pd.DataFrame(columns=minus.columns)),
         max(0, MAXNEG),
         max(0, REVERY))
        for c in chroms]

pieces: List[pd.DataFrame] = []
if NCORES > 1 and len(work) > 1:
    n = min(NCORES, max(1, cpu_count()-1), len(work))
    chunksize = max(1, math.ceil(len(work) / (n * 4)))
    log(f"pairing (ncores={n}, chroms={len(work)}, chunksize={chunksize})…")
    with Pool(processes=n) as pool:
        for res in pool.imap_unordered(pair_on_chr_fast, work, chunksize=chunksize):
            if res is not None and not res.empty:
                pieces.append(res)
            log(f"…chrom done (pairs={0 if res is None or res.empty else len(res)})")
else:
    log("pairing (single core)…")
    for t in work:
        res = pair_on_chr_fast(t)
        if res is not None and not res.empty:
            pieces.append(res)
        log(f"…{t[0]} done (pairs={0 if res is None or res.empty else len(res)})")

paired = (pd.concat(pieces, ignore_index=True)
          if pieces else pd.DataFrame(columns=["chr","start","end","total","dist","balv"]))
log(f"raw pairs: {len(paired)}")

merged = merge_windows(paired[["chr","start","end","total"]], MERGE_GAP)

# optional max width filter
if MAXW > 0 and not merged.empty:
    merged = merged[(merged["end"] - merged["start"]) <= MAXW].reset_index(drop=True)

# deterministic sort before write
merged.sort_values(["chr","start","end"], inplace=True, kind="mergesort")

out = pathlib.Path(args.out_bed)
merged.to_csv(out, sep="\t", header=False, index=False)
log(f"✓ {args.sample}: DT={len(merged)} → {out}  (wall={tsec()})")

# ───────── summary TSV (optional) ─────────
if args.write_summary:
    wall = f"{time.time()-start_t:.3f}"
    with open(args.write_summary, "w") as fh:
        fh.write("sample\tn_pos_pk\tn_neg_pk\tn_pairs_raw\tn_dt\twall_s\n")
        fh.write(f"{args.sample}\t{len(pk_pos)}\t{len(pk_neg)}\t{len(paired)}\t{len(merged)}\t{wall}\n")

# ───────── QC (optional) ─────────
if int(args.qc) == 1:
    try:
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        qc_pdf = out.with_suffix(".qc.pdf")
        with PdfPages(qc_pdf) as pdf:
            def _hist(arr, title, bins=100, logy=True):
                plt.figure()
                a = np.asarray(arr)
                if a.size == 0: a = np.array([0.0])
                plt.hist(a, bins=bins)
                if logy: plt.yscale("log")
                plt.title(f"{args.sample} — {title}")
                pdf.savefig(); plt.close()
            _hist(pk_pos["signal"], "pos peak sum")
            _hist(pk_neg["signal"], "neg peak sum")
            widths = (merged["end"] - merged["start"]) if not merged.empty else np.array([])
            _hist(widths, "DT widths (bp)")
            if not paired.empty and "dist" in paired:
                _hist(paired["dist"], "pair center distance |Δ| (bp)")
            if not paired.empty and "balv" in paired:
                _hist(paired["balv"], "balance min/max", bins=50, logy=False)
        log(f"✓ QC {qc_pdf}")
    except Exception as e:
        log(f"WARN QC failed: {e}")
sys.stderr.write(f"[divergent_py] done ts={datetime.datetime.utcnow().isoformat()}Z\n")
