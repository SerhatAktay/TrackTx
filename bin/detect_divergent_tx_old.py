#!/usr/bin/env python3
# =============================================================================
# detect_divergent_tx.py — Divergent Transcription Detector
# =============================================================================
# Detects divergent transcription from strand-specific 3' PRO-seq bedGraphs
#
# Inputs: bedGraph files (4-col: chr start end signal; plain or .gz)
#   --pos : positive-strand 3' signal
#   --neg : negative-strand 3' signal (should already be negative values)
#
# Algorithm:
#   1. Load and sanitize bedGraphs
#   2. Optional Gaussian smoothing
#   3. Call peak blocks on each strand independently
#   4. Pair peaks using edge-to-edge distance (overlap-aware)
#   5. Apply balance, overlap, and valley filters
#   6. Merge paired windows within merge_gap
#   7. Output BED4: chr start end total_signal
#
# Exit Codes:
#   0 = Success
#   1 = Empty input (expected, creates empty outputs)
#   2 = Invalid parameters
#   3 = Missing dependencies
#   4 = File I/O error
#   5 = Processing error
# =============================================================================

from __future__ import annotations
import argparse
import sys
import os
import time
import math
import pathlib
import datetime
import re
from typing import Dict, List, Tuple, Optional
import warnings

# Suppress specific warnings
warnings.filterwarnings('ignore', category=FutureWarning)
warnings.filterwarnings('ignore', category=RuntimeWarning)

# Check dependencies early
MISSING_DEPS = []
try:
    import numpy as np
except ImportError:
    MISSING_DEPS.append("numpy")
    np = None

try:
    import pandas as pd
except ImportError:
    MISSING_DEPS.append("pandas")
    pd = None

if MISSING_DEPS:
    sys.stderr.write(f"ERROR: Missing required dependencies: {', '.join(MISSING_DEPS)}\n")
    sys.stderr.write("Install with: pip install numpy pandas\n")
    sys.exit(3)

# Control threading for better resource management
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def parse_args() -> argparse.Namespace:
    """Parse and validate command-line arguments."""
    ap = argparse.ArgumentParser(
        description="Detect divergent transcription from 3' PRO-seq bedGraphs",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Required arguments
    ap.add_argument("--sample", required=True,
                    help="Sample identifier")
    ap.add_argument("--pos", required=True,
                    help="Positive strand bedGraph (.gz supported)")
    ap.add_argument("--neg", required=True,
                    help="Negative strand bedGraph (.gz supported)")
    ap.add_argument("--out-bed", required=True,
                    help="Output BED file path")
    
    # Detection parameters
    ap.add_argument("--nt-window", type=int, default=500,
                    help="Max edge-to-edge gap for pairing (default: 500)")
    ap.add_argument("--threshold", type=float, default=1.0,
                    help="Per-bin signal threshold (default: 1.0)")
    ap.add_argument("--sum-thr", type=float, default=5.0,
                    help="Minimum peak total signal (default: 5.0)")
    ap.add_argument("--bin-gap", type=int, default=100,
                    help="Max gap to merge bins into peaks (default: 100)")
    ap.add_argument("--merge-gap", type=int, default=50,
                    help="Max gap to merge paired windows (default: 50)")
    ap.add_argument("--balance", type=float, default=0.0,
                    help="Min balance ratio min/max (default: 0.0, off)")
    ap.add_argument("--overlap-bp", type=int, default=0,
                    help="Min overlap requirement for overlapping peaks (default: 0)")
    ap.add_argument("--valley-thr", default="Inf",
                    help="Max signal in gap between peaks (default: Inf, off)")
    ap.add_argument("--smooth-sd", type=int, default=0,
                    help="Gaussian smoothing SD in bp (default: 0, off)")
    ap.add_argument("--max-dt-width", type=int, default=0,
                    help="Max divergent region width (default: 0, off)")
    
    # Processing options
    ap.add_argument("--ncores", type=int, default=1,
                    help="Number of CPU cores (default: 1)")
    ap.add_argument("--assume-sorted", action="store_true",
                    help="Trust input sortedness, skip sorting")
    ap.add_argument("--neg-abs", action="store_true",
                    help="Treat --neg as magnitudes |NEG| instead of -|NEG|")
    
    # Output options
    ap.add_argument("--qc", type=int, default=0,
                    help="Generate QC PDF (1=yes, 0=no)")
    ap.add_argument("--write-summary", default=None,
                    help="Path to summary TSV")
    ap.add_argument("--quiet", action="store_true",
                    help="Suppress progress messages")
    
    # Advanced options
    ap.add_argument("--max-neg-per-pos", type=int, default=64,
                    help="Max negative peaks to check per positive (default: 64)")
    ap.add_argument("--report-every", type=int, default=20000,
                    help="Progress report interval (default: 20000, 0=off)")
    
    return ap.parse_args()

# =============================================================================
# UTILITIES
# =============================================================================

class Timer:
    """Simple timer for tracking wall time."""
    def __init__(self):
        self.start_time = time.time()
    
    def elapsed(self) -> str:
        return f"{time.time() - self.start_time:.1f}s"

def log(msg: str, quiet: bool = False):
    """Thread-safe logging to stderr."""
    if not quiet:
        timestamp = time.strftime('%H:%M:%S')
        sys.stderr.write(f"[{timestamp}] {msg}\n")
        sys.stderr.flush()

def natural_sort_key(s: str) -> List:
    """Natural sorting key for chromosome names (chr1, chr2, ..., chr10)."""
    return [int(c) if c.isdigit() else c.lower() 
            for c in re.split(r'([0-9]+)', str(s))]

def clamp_int(x: int, lo: int, hi: Optional[int] = None) -> int:
    """Clamp integer to range [lo, hi]."""
    x = int(x)
    x = max(lo, x)
    return x if hi is None else min(hi, x)

def clamp_float(x: float, lo: float, hi: Optional[float] = None) -> float:
    """Clamp float to range [lo, hi]."""
    x = float(x)
    x = max(lo, x)
    return x if hi is None else min(hi, x)

def validate_parameters(args: argparse.Namespace) -> Tuple[bool, List[str]]:
    """Validate parameters and return (is_valid, warnings)."""
    warnings_list = []
    is_valid = True
    
    # Check file existence
    if not os.path.exists(args.pos):
        log(f"ERROR: Positive bedGraph not found: {args.pos}")
        is_valid = False
    if not os.path.exists(args.neg):
        log(f"ERROR: Negative bedGraph not found: {args.neg}")
        is_valid = False
    
    # Validate numeric parameters
    if args.threshold < 0:
        warnings_list.append(f"threshold < 0 ({args.threshold}), clamped to 0")
    if args.sum_thr < 0:
        warnings_list.append(f"sum_thr < 0 ({args.sum_thr}), clamped to 0")
    if args.nt_window < 1:
        warnings_list.append(f"nt_window < 1 ({args.nt_window}), clamped to 1")
    if args.balance < 0 or args.balance > 1:
        warnings_list.append(f"balance out of [0,1] ({args.balance}), clamped")
    
    return is_valid, warnings_list

# =============================================================================
# DATA I/O
# =============================================================================

def read_bedgraph(path: str, role: str, assume_sorted: bool, 
                  neg_abs: bool, quiet: bool = False) -> pd.DataFrame:
    """
    Read and sanitize bedGraph file.
    
    Args:
        path: Path to bedGraph file (.gz supported)
        role: "pos" or "neg" for strand
        assume_sorted: Skip sorting if True
        neg_abs: If True, treat neg strand as |NEG| instead of -|NEG|
        quiet: Suppress logging
    
    Returns:
        DataFrame with columns: chr, start, end, sig
    """
    log(f"Reading {role} strand: {os.path.basename(path)}", quiet)
    
    # Try pyarrow engine for speed, fallback to C
    engine = "c"
    try:
        import pyarrow
        engine = "pyarrow"
    except ImportError:
        pass
    
    try:
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            comment="#",
            names=["chr", "start", "end", "sig"],
            usecols=[0, 1, 2, 3],
            dtype={
                "chr": "category",
                "start": "int64",
                "end": "int64",
                "sig": "float64"
            },
            compression="infer",
            engine=engine,
        )
    except Exception as e:
        log(f"ERROR: Failed to read {path}: {e}")
        sys.exit(4)
    
    if df.empty:
        log(f"WARNING: Empty bedGraph: {path}", quiet)
        return df
    
    # Sanitize
    df = df.loc[(df["end"] > df["start"]) & df["chr"].notna()].copy()
    
    # Handle signs
    if role == "pos":
        df["sig"] = df["sig"].abs()
    elif role == "neg":
        if neg_abs:
            df["sig"] = df["sig"].abs()
        else:
            df["sig"] = -df["sig"].abs()
    else:
        raise ValueError(f"Invalid role: {role}")
    
    # Sort if needed
    if not assume_sorted and not df.empty:
        df.sort_values(["chr", "start"], inplace=True, kind="mergesort")
    
    df.reset_index(drop=True, inplace=True)
    log(f"  → {len(df):,} bins loaded", quiet)
    
    return df

def gaussian_smooth(df: pd.DataFrame, sd: int, quiet: bool = False):
    """Apply Gaussian smoothing to signal column in-place."""
    if sd <= 0 or df.empty:
        return
    
    log(f"Applying Gaussian smoothing (SD={sd}bp)...", quiet)
    
    # Try scipy first
    try:
        from scipy.ndimage import gaussian_filter1d
        smoother = lambda v: gaussian_filter1d(v, sd, mode="constant")
    except ImportError:
        # Fallback to numpy convolution
        half = max(1, int(5 * sd))
        xs = np.arange(-half, half + 1)
        kernel = np.exp(-(xs * xs) / (2.0 * sd * sd))
        kernel /= kernel.sum()
        smoother = lambda v: np.convolve(v, kernel, mode="same")
    
    # Smooth per chromosome
    for chrom, idx in df.groupby("chr", sort=False, observed=False).groups.items():
        v = df.loc[idx, "sig"].to_numpy()
        df.loc[idx, "sig"] = smoother(v)

# =============================================================================
# PEAK CALLING
# =============================================================================

def call_peaks(df: pd.DataFrame, threshold: float, sum_thr: float, 
               bin_gap: int, quiet: bool = False) -> pd.DataFrame:
    """
    Call peak blocks from bedGraph.
    
    Args:
        df: Input bedGraph DataFrame
        threshold: Per-bin minimum signal
        sum_thr: Minimum total signal for peak
        bin_gap: Max gap to merge bins
        quiet: Suppress logging
    
    Returns:
        DataFrame with columns: chr, start, end, height, signal
    """
    if df.empty:
        return pd.DataFrame(columns=["chr", "start", "end", "height", "signal"])
    
    # Filter by threshold
    v = np.abs(df["sig"].to_numpy())
    keep = v >= threshold
    
    if not keep.any():
        log("  → 0 peaks (no bins above threshold)", quiet)
        return pd.DataFrame(columns=["chr", "start", "end", "height", "signal"])
    
    # Extract passing bins
    chr_a = df["chr"].to_numpy()[keep]
    start = df["start"].to_numpy()[keep]
    end = df["end"].to_numpy()[keep]
    signal = v[keep]
    
    # Group bins into blocks
    # New block starts on: chromosome change OR gap > bin_gap
    grp = np.concatenate([
        [True],  # First bin starts new group
        (chr_a[1:] != chr_a[:-1]) | (start[1:] > end[:-1] + bin_gap)
    ]).cumsum()
    
    # Aggregate blocks
    chr_ser = pd.Series(chr_a)
    start_ser = pd.Series(start)
    end_ser = pd.Series(end)
    signal_ser = pd.Series(signal)
    
    peaks = pd.DataFrame({
        "chr": chr_ser.groupby(grp).first(),
        "start": start_ser.groupby(grp).min(),
        "end": end_ser.groupby(grp).max(),
        "height": signal_ser.groupby(grp).max(),
        "signal": signal_ser.groupby(grp).sum(),
    }).reset_index(drop=True)
    
    # Filter by total signal
    peaks = peaks.loc[peaks["signal"] >= sum_thr].reset_index(drop=True)
    
    log(f"  → {len(peaks):,} peaks called", quiet)
    return peaks

# =============================================================================
# PAIRING LOGIC
# =============================================================================

def check_valley(idx_df: pd.DataFrame, valley_start: int, valley_end: int, 
                 threshold: float) -> bool:
    """
    Check if valley region has signal below threshold.
    
    Args:
        idx_df: Chromosome bedGraph DataFrame
        valley_start: Valley start position
        valley_end: Valley end position
        threshold: Max allowed signal
    
    Returns:
        True if valley is clear (all signal <= threshold)
    """
    if not math.isfinite(threshold) or idx_df is None or idx_df.empty:
        return True
    
    if valley_start >= valley_end:
        # No valley (peaks overlap)
        return True
    
    # Find overlapping bins
    starts = idx_df["start"].to_numpy()
    ends = idx_df["end"].to_numpy()
    sig = np.abs(idx_df["sig"].to_numpy())
    
    # Binary search for potential overlaps
    i = np.searchsorted(ends, valley_start, side="right")
    n = len(starts)
    
    # Check each overlapping bin
    while i < n and starts[i] < valley_end:
        if ends[i] > valley_start:  # Overlaps valley
            if sig[i] > threshold:
                return False  # Valley has signal above threshold
        i += 1
    
    return True

def pair_peaks_on_chromosome(
    chrom: str,
    pos_peaks: pd.DataFrame,
    neg_peaks: pd.DataFrame,
    nt_window: int,
    balance: float,
    overlap_bp: int,
    valley_thr: float,
    pos_idx: pd.DataFrame,
    neg_idx: pd.DataFrame,
    max_neg_per_pos: int,
    report_every: int,
    quiet: bool = False
) -> pd.DataFrame:
    """
    Pair positive and negative peaks on a chromosome.
    
    Uses edge-to-edge distance with overlap awareness:
    - Overlapping peaks are always considered (gap=0)
    - Non-overlapping peaks paired if gap <= nt_window
    
    Args:
        chrom: Chromosome name
        pos_peaks: Positive strand peaks
        neg_peaks: Negative strand peaks
        nt_window: Max edge-to-edge gap
        balance: Min balance ratio
        overlap_bp: Min overlap requirement (if overlapping)
        valley_thr: Max signal in gap
        pos_idx: Positive strand bedGraph (for valley check)
        neg_idx: Negative strand bedGraph (for valley check)
        max_neg_per_pos: Max neg peaks to check per pos peak
        report_every: Progress report interval
        quiet: Suppress logging
    
    Returns:
        DataFrame with paired windows
    """
    if pos_peaks.empty or neg_peaks.empty:
        return pd.DataFrame(columns=["chr", "start", "end", "total", 
                                     "overlap", "gap", "balance"])
    
    # Extract arrays
    ps = pos_peaks["start"].to_numpy()
    pe = pos_peaks["end"].to_numpy()
    ph = pos_peaks["height"].to_numpy()
    pg = pos_peaks["signal"].to_numpy()
    
    ns = neg_peaks["start"].to_numpy()
    ne = neg_peaks["end"].to_numpy()
    nh = neg_peaks["height"].to_numpy()
    ng = neg_peaks["signal"].to_numpy()
    
    pairs: List[Tuple[str, int, int, float, int, int, float]] = []
    
    # For each positive peak
    for i in range(len(ps)):
        if report_every > 0 and i > 0 and i % report_every == 0:
            log(f"  {chrom}: processed {i:,}/{len(ps):,} + peaks", quiet)
        
        # Find candidate negative peaks
        # Use a wide search window (nt_window + max peak width)
        search_margin = nt_window + 1000  # Assume peaks < 1kb wide
        left = np.searchsorted(ns, ps[i] - search_margin, side="left")
        right = np.searchsorted(ns, pe[i] + search_margin, side="right")
        
        if right <= left:
            continue
        
        # Limit candidates if too many
        if max_neg_per_pos > 0 and (right - left) > max_neg_per_pos:
            mid = (left + right) // 2
            half = max_neg_per_pos // 2
            left = max(left, mid - half)
            right = left + max_neg_per_pos
        
        # Get candidate slices
        ns_slice = ns[left:right]
        ne_slice = ne[left:right]
        nh_slice = nh[left:right]
        ng_slice = ng[left:right]
        
        # Calculate overlaps and gaps (VECTORIZED)
        overlap = np.maximum(0, np.minimum(pe[i], ne_slice) - np.maximum(ps[i], ns_slice))
        gap = np.maximum(0, np.maximum(ps[i], ns_slice) - np.minimum(pe[i], ne_slice))
        
        # Distance filter: overlap > 0 OR gap <= nt_window
        dist_mask = (overlap > 0) | (gap <= nt_window)
        if not dist_mask.any():
            continue
        
        # Apply distance filter
        overlap = overlap[dist_mask]
        gap = gap[dist_mask]
        ns_slice = ns_slice[dist_mask]
        ne_slice = ne_slice[dist_mask]
        nh_slice = nh_slice[dist_mask]
        ng_slice = ng_slice[dist_mask]
        
        # Balance filter
        if balance > 0.0:
            balance_ratio = np.minimum(ph[i], nh_slice) / np.maximum(ph[i], nh_slice)
            bal_mask = balance_ratio >= balance
            if not bal_mask.any():
                continue
            
            overlap = overlap[bal_mask]
            gap = gap[bal_mask]
            ns_slice = ns_slice[bal_mask]
            ne_slice = ne_slice[bal_mask]
            ng_slice = ng_slice[bal_mask]
            balance_ratio = balance_ratio[bal_mask]
        else:
            balance_ratio = np.minimum(ph[i], nh_slice) / np.maximum(ph[i], nh_slice)
        
        # Overlap requirement filter
        if overlap_bp > 0:
            ov_mask = overlap >= overlap_bp
            if not ov_mask.any():
                continue
            
            overlap = overlap[ov_mask]
            gap = gap[ov_mask]
            ns_slice = ns_slice[ov_mask]
            ne_slice = ne_slice[ov_mask]
            ng_slice = ng_slice[ov_mask]
            balance_ratio = balance_ratio[ov_mask]
        
        # Valley filter (if enabled)
        if math.isfinite(valley_thr):
            keep_mask = np.ones(len(ns_slice), dtype=bool)
            for j in range(len(ns_slice)):
                if gap[j] > 0:  # Only check if there's a gap
                    valley_start = int(min(pe[i], ne_slice[j]))
                    valley_end = int(max(ps[i], ns_slice[j]))
                    
                    if not (check_valley(pos_idx, valley_start, valley_end, valley_thr) and
                            check_valley(neg_idx, valley_start, valley_end, valley_thr)):
                        keep_mask[j] = False
            
            if not keep_mask.any():
                continue
            
            overlap = overlap[keep_mask]
            gap = gap[keep_mask]
            ns_slice = ns_slice[keep_mask]
            ne_slice = ne_slice[keep_mask]
            ng_slice = ng_slice[keep_mask]
            balance_ratio = balance_ratio[keep_mask]
        
        # Create paired windows
        start = np.minimum(ps[i], ns_slice)
        end = np.maximum(pe[i], ne_slice)
        total = pg[i] + ng_slice
        
        # Store pairs
        for s, e, t, ov, g, bal in zip(start, end, total, overlap, gap, balance_ratio):
            pairs.append((chrom, int(s), int(e), float(t), int(ov), int(g), float(bal)))
    
    if not pairs:
        return pd.DataFrame(columns=["chr", "start", "end", "total", 
                                     "overlap", "gap", "balance"])
    
    return pd.DataFrame(
        pairs,
        columns=["chr", "start", "end", "total", "overlap", "gap", "balance"]
    )

# =============================================================================
# MERGING
# =============================================================================

def merge_windows(df: pd.DataFrame, merge_gap: int, quiet: bool = False) -> pd.DataFrame:
    """
    Merge overlapping or nearby paired windows.
    
    Args:
        df: Paired windows DataFrame
        merge_gap: Max gap to merge windows
        quiet: Suppress logging
    
    Returns:
        Merged DataFrame
    """
    if df is None or df.empty:
        return pd.DataFrame(columns=["chr", "start", "end", "total"])
    
    log(f"Merging windows (gap={merge_gap})...", quiet)
    
    # Sort by position
    df = df.sort_values(["chr", "start", "end"]).reset_index(drop=True)
    
    merged_rows: List[Tuple[str, int, int, float]] = []
    
    # Track current merged window
    cur_chr: Optional[str] = None
    cur_start: int = 0
    cur_end: int = 0
    cur_total: float = 0.0
    
    for row in df.itertuples(index=False):
        # Start new window if different chromosome or gap too large
        if (cur_chr is None or 
            row.chr != cur_chr or 
            row.start > cur_end + merge_gap):
            
            # Save previous window
            if cur_chr is not None:
                merged_rows.append((cur_chr, cur_start, cur_end, cur_total))
            
            # Start new window
            cur_chr = row.chr
            cur_start = int(row.start)
            cur_end = int(row.end)
            cur_total = float(row.total)
        else:
            # Extend current window
            cur_end = max(cur_end, int(row.end))
            cur_total += float(row.total)
    
    # Save final window
    if cur_chr is not None:
        merged_rows.append((cur_chr, cur_start, cur_end, cur_total))
    
    result = pd.DataFrame(
        merged_rows,
        columns=["chr", "start", "end", "total"]
    )
    
    log(f"  → {len(result):,} merged regions", quiet)
    return result

# =============================================================================
# QC PLOTTING
# =============================================================================

def generate_qc_pdf(
    out_bed: pathlib.Path,
    sample: str,
    pos_peaks: pd.DataFrame,
    neg_peaks: pd.DataFrame,
    paired: pd.DataFrame,
    merged: pd.DataFrame,
    quiet: bool = False
):
    """Generate QC PDF with diagnostic plots."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
    except ImportError:
        log("WARNING: matplotlib not available, skipping QC PDF", quiet)
        return
    
    log("Generating QC PDF...", quiet)
    qc_pdf = out_bed.with_suffix(".qc.pdf")
    
    try:
        with PdfPages(qc_pdf) as pdf:
            # Helper for histograms
            def plot_hist(data, title, xlabel, bins=100, log_y=True):
                plt.figure(figsize=(8, 5))
                arr = np.asarray(data)
                if arr.size == 0:
                    arr = np.array([0.0])
                plt.hist(arr, bins=bins, edgecolor='black', linewidth=0.5)
                if log_y:
                    plt.yscale('log')
                plt.xlabel(xlabel)
                plt.ylabel('Count (log scale)' if log_y else 'Count')
                plt.title(f"{sample} — {title}")
                plt.grid(True, alpha=0.3)
                plt.tight_layout()
                pdf.savefig()
                plt.close()
            
            # Plot 1: Positive peak signals
            if not pos_peaks.empty:
                plot_hist(pos_peaks["signal"], 
                         "Positive Strand Peak Signals",
                         "Total Signal")
            
            # Plot 2: Negative peak signals
            if not neg_peaks.empty:
                plot_hist(neg_peaks["signal"],
                         "Negative Strand Peak Signals",
                         "Total Signal")
            
            # Plot 3: Divergent region widths
            if not merged.empty:
                widths = merged["end"] - merged["start"]
                plot_hist(widths,
                         "Divergent Region Widths",
                         "Width (bp)")
            
            # Plot 4: Pair overlaps
            if not paired.empty and "overlap" in paired:
                plot_hist(paired["overlap"],
                         "Pair Overlap Sizes",
                         "Overlap (bp)")
            
            # Plot 5: Pair gaps (for non-overlapping)
            if not paired.empty and "gap" in paired:
                gaps_only = paired.loc[paired["gap"] > 0, "gap"]
                if len(gaps_only) > 0:
                    plot_hist(gaps_only,
                             "Pair Edge-to-Edge Gaps (non-overlapping)",
                             "Gap (bp)")
            
            # Plot 6: Balance ratios
            if not paired.empty and "balance" in paired:
                plot_hist(paired["balance"],
                         "Peak Balance Ratios",
                         "Balance (min/max)",
                         bins=50,
                         log_y=False)
        
        log(f"  → QC PDF: {qc_pdf}", quiet)
    
    except Exception as e:
        log(f"WARNING: QC PDF generation failed: {e}", quiet)

# =============================================================================
# MAIN WORKFLOW
# =============================================================================

def main():
    """Main workflow."""
    timer = Timer()
    
    # Parse arguments
    args = parse_args()
    
    log(f"╔════════════════════════════════════════════════════════════════╗")
    log(f"║ Divergent Transcription Detector                              ║")
    log(f"║ Sample: {args.sample:<53} ║")
    log(f"╚════════════════════════════════════════════════════════════════╝")
    log(f"Start: {datetime.datetime.utcnow().isoformat()}Z")
    
    # Validate parameters
    is_valid, warnings_list = validate_parameters(args)
    if not is_valid:
        log("ERROR: Invalid parameters")
        sys.exit(2)
    
    for warning in warnings_list:
        log(f"WARNING: {warning}")
    
    # Clamp parameters
    THR = clamp_float(args.threshold, 0.0)
    SUMTHR = clamp_float(args.sum_thr, 0.0)
    BIN_GAP = clamp_int(args.bin_gap, 0)
    MERGE_GAP = clamp_int(args.merge_gap, 0)
    NTWIN = clamp_int(args.nt_window, 1)
    BAL = clamp_float(args.balance, 0.0, 1.0)
    SMOOTH_SD = clamp_int(args.smooth_sd, 0)
    OVERLAP_BP = clamp_int(args.overlap_bp, 0)
    MAXW = clamp_int(args.max_dt_width, 0)
    
    # Parse valley threshold
    try:
        VALLEY_THR = (float("inf") if str(args.valley_thr).lower() == "inf" 
                     else clamp_float(float(args.valley_thr), 0.0))
    except:
        VALLEY_THR = float("inf")
    
    # Processing options
    NCORES = clamp_int(args.ncores, 1)
    MAXNEG = clamp_int(args.max_neg_per_pos, 0)
    REVERY = clamp_int(args.report_every, 0)
    
    # Log configuration
    log(f"Configuration:")
    log(f"  Pairing window: {NTWIN} bp (edge-to-edge)")
    log(f"  Signal threshold: {THR}")
    log(f"  Sum threshold: {SUMTHR}")
    log(f"  Bin gap: {BIN_GAP} bp")
    log(f"  Merge gap: {MERGE_GAP} bp")
    log(f"  Balance: {BAL}" + (" (off)" if BAL == 0.0 else ""))
    log(f"  Overlap requirement: {OVERLAP_BP} bp" + (" (off)" if OVERLAP_BP == 0 else ""))
    log(f"  Valley threshold: {VALLEY_THR}" + (" (off)" if not math.isfinite(VALLEY_THR) else ""))
    log(f"  Smooth SD: {SMOOTH_SD} bp" + (" (off)" if SMOOTH_SD == 0 else ""))
    log(f"  Max width: {MAXW} bp" + (" (off)" if MAXW == 0 else ""))
    log(f"  Cores: {NCORES}")
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1: Load Data
    # ─────────────────────────────────────────────────────────────────────────
    log(f"\n[1/7] Loading bedGraphs...")
    
    pos_df = read_bedgraph(args.pos, "pos", args.assume_sorted, 
                           args.neg_abs, args.quiet)
    neg_df = read_bedgraph(args.neg, "neg", args.assume_sorted,
                           args.neg_abs, args.quiet)
    
    # Handle empty input
    if pos_df.empty or neg_df.empty:
        log("\nWARNING: Empty input detected")
        
        # Create empty outputs
        pathlib.Path(args.out_bed).write_text("")
        
        if args.write_summary:
            with open(args.write_summary, "w") as f:
                f.write("sample\tn_pos_pk\tn_neg_pk\tn_pairs_raw\tn_dt\twall_s\n")
                f.write(f"{args.sample}\t0\t0\t0\t0\t{timer.elapsed().rstrip('s')}\n")
        
        log(f"✓ Complete: 0 divergent regions (wall={timer.elapsed()})")
        sys.exit(1)  # Exit code 1 = empty input (expected)
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2: Optional Smoothing
    # ─────────────────────────────────────────────────────────────────────────
    if SMOOTH_SD > 0:
        log(f"\n[2/7] Applying smoothing...")
        gaussian_smooth(pos_df, SMOOTH_SD, args.quiet)
        gaussian_smooth(neg_df, SMOOTH_SD, args.quiet)
    else:
        log(f"\n[2/7] Smoothing: skipped")
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3: Call Peaks
    # ─────────────────────────────────────────────────────────────────────────
    log(f"\n[3/7] Calling peaks...")
    log(f"  Positive strand:")
    pos_peaks = call_peaks(pos_df, THR, SUMTHR, BIN_GAP, args.quiet)
    log(f"  Negative strand:")
    neg_peaks = call_peaks(neg_df, THR, SUMTHR, BIN_GAP, args.quiet)
    
    # Build chromosome indices for valley checking
    def build_chr_index(df: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        if df.empty:
            return {}
        return {
            str(c): sub.sort_values("start").reset_index(drop=True)
            for c, sub in df.groupby("chr", sort=False, observed=False)
        }
    
    pos_idx = build_chr_index(pos_df)
    neg_idx = build_chr_index(neg_df)
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 4: Pair Peaks
    # ─────────────────────────────────────────────────────────────────────────
    log(f"\n[4/7] Pairing peaks...")
    
    # Get chromosomes present in both strands
    chroms = sorted(
        set(map(str, pos_peaks["chr"])) & set(map(str, neg_peaks["chr"])),
        key=natural_sort_key
    )
    
    log(f"  Chromosomes: {len(chroms)}")
    
    # Prepare tasks for parallel processing
    tasks = []
    for chrom in chroms:
        pos_sub = pos_peaks[pos_peaks["chr"].astype(str) == chrom].reset_index(drop=True)
        neg_sub = neg_peaks[neg_peaks["chr"].astype(str) == chrom].reset_index(drop=True)
        
        tasks.append((
            chrom, pos_sub, neg_sub,
            NTWIN, BAL, OVERLAP_BP, VALLEY_THR,
            pos_idx.get(chrom, pd.DataFrame(columns=pos_df.columns)),
            neg_idx.get(chrom, pd.DataFrame(columns=neg_df.columns)),
            MAXNEG, REVERY, args.quiet
        ))
    
    # Process chromosomes
    pieces: List[pd.DataFrame] = []
    
    if NCORES > 1 and len(tasks) > 1:
        from multiprocessing import Pool, cpu_count
        n_proc = min(NCORES, max(1, cpu_count() - 1), len(tasks))
        chunksize = max(1, math.ceil(len(tasks) / (n_proc * 4)))
        
        log(f"  Using {n_proc} cores (chunksize={chunksize})")
        
        with Pool(processes=n_proc) as pool:
            for result in pool.starmap(pair_peaks_on_chromosome, tasks, chunksize=chunksize):
                if result is not None and not result.empty:
                    pieces.append(result)
    else:
        log(f"  Single-threaded processing")
        for task in tasks:
            result = pair_peaks_on_chromosome(*task)
            if result is not None and not result.empty:
                pieces.append(result)
    
    # Combine results
    if pieces:
        paired = pd.concat(pieces, ignore_index=True)
    else:
        paired = pd.DataFrame(columns=["chr", "start", "end", "total", 
                                       "overlap", "gap", "balance"])
    
    log(f"  → {len(paired):,} raw pairs")
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 5: Merge Windows
    # ─────────────────────────────────────────────────────────────────────────
    log(f"\n[5/7] Merging windows...")
    
    merged = merge_windows(
        paired[["chr", "start", "end", "total"]] if not paired.empty else paired,
        MERGE_GAP,
        args.quiet
    )
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 6: Apply Width Filter
    # ─────────────────────────────────────────────────────────────────────────
    if MAXW > 0 and not merged.empty:
        log(f"\n[6/7] Applying width filter (max={MAXW}bp)...")
        n_before = len(merged)
        merged = merged[(merged["end"] - merged["start"]) <= MAXW].reset_index(drop=True)
        n_after = len(merged)
        log(f"  → Retained {n_after:,}/{n_before:,} regions")
    else:
        log(f"\n[6/7] Width filter: skipped")
    
    # ─────────────────────────────────────────────────────────────────────────
    # STEP 7: Write Outputs
    # ─────────────────────────────────────────────────────────────────────────
    log(f"\n[7/7] Writing outputs...")
    
    # Sort deterministically
    if not merged.empty:
        merged.sort_values(["chr", "start", "end"], inplace=True, kind="mergesort")
    
    # Write BED
    out_path = pathlib.Path(args.out_bed)
    merged.to_csv(out_path, sep="\t", header=False, index=False)
    log(f"  → BED: {out_path} ({len(merged):,} regions)")
    
    # Write summary
    if args.write_summary:
        with open(args.write_summary, "w") as f:
            f.write("sample\tn_pos_pk\tn_neg_pk\tn_pairs_raw\tn_dt\twall_s\n")
            f.write(f"{args.sample}\t{len(pos_peaks)}\t{len(neg_peaks)}\t"
                   f"{len(paired)}\t{len(merged)}\t{timer.elapsed().rstrip('s')}\n")
        log(f"  → Summary: {args.write_summary}")
    
    # Generate QC PDF
    if args.qc == 1:
        generate_qc_pdf(out_path, args.sample, pos_peaks, neg_peaks, 
                       paired, merged, args.quiet)
    
    # ─────────────────────────────────────────────────────────────────────────
    # COMPLETE
    # ─────────────────────────────────────────────────────────────────────────
    log(f"\n╔════════════════════════════════════════════════════════════════╗")
    log(f"║ ✓ COMPLETE                                                     ║")
    log(f"║   Sample: {args.sample:<52} ║")
    log(f"║   Divergent regions: {len(merged):<40} ║")
    log(f"║   Wall time: {timer.elapsed():<47} ║")
    log(f"╚════════════════════════════════════════════════════════════════╝")
    log(f"End: {datetime.datetime.utcnow().isoformat()}Z")
    
    sys.exit(0)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("\nInterrupted by user\n")
        sys.exit(130)
    except Exception as e:
        sys.stderr.write(f"\nFATAL ERROR: {e}\n")
        import traceback
        traceback.print_exc()
        sys.exit(5)