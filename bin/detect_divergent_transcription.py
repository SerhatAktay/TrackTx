#!/usr/bin/env python3
"""
detect_divergent_transcription.py

Statistical detection of divergent transcription from PRO-seq data.

Uses relaxed initial thresholds followed by Gaussian Mixture Model scoring
and FDR-controlled filtering to identify bidirectional transcription sites.

Output BED format: chr, start, end, total_signal, confidence_score

Author: Serhat Aktay (adapted from TrackTx pipeline)
Date: December 2025
Version: 1.0

Dependencies: pandas, numpy, scikit-learn, scipy

Usage:
    python detect_divergent_transcription.py \
        --pos sample_pos.bedgraph \
        --neg sample_neg.bedgraph \
        --out divergent.bed \
        --threshold 1.0 \
        --sum-thr 5.0 \
        --ncores 4 \
        --fdr 0.05
"""

import argparse
import sys
import os
import time
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')

__version__ = "1.0"


def parse_args():
    """Parse command-line arguments."""
    ap = argparse.ArgumentParser(
        description="Detect divergent transcription from PRO-seq bedGraphs",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic usage with auto-calibration
  python detect_divergent_transcription.py --sample MySample --pos pos.bedgraph --neg neg.bedgraph --out divergent.bed
  
  # Relaxed thresholds for maximum sensitivity
  python detect_divergent_transcription.py --sample MySample --pos pos.bg --neg neg.bg --out div.bed \\
      --threshold 1.0 --sum-thr 5.0 --fdr 0.05
  
  # Conservative mode
  python detect_divergent_transcription.py --sample MySample --pos pos.bg --neg neg.bg --out div.bed \\
      --threshold 2.0 --sum-thr 15.0 --fdr 0.01

Output:
  BED file with 5 columns: chr, start, end, total_signal, confidence_score
  - confidence_score: Posterior probability from mixture model (0-1)
  - Higher scores = higher confidence of true divergent transcription
        """
    )
    
    # Required arguments
    ap.add_argument("--sample", required=True,
                    help="Sample identifier")
    ap.add_argument("--pos", required=True, 
                    help="Positive strand bedGraph file (.gz supported)")
    ap.add_argument("--neg", required=True,
                    help="Negative strand bedGraph file (.gz supported)")
    ap.add_argument("--out", required=True,
                    help="Output BED file")
    
    # Detection parameters
    ap.add_argument("--threshold", type=float, default=None,
                    help="Per-bin signal threshold (default: auto-calibrate to 95th percentile)")
    ap.add_argument("--sum-thr", type=float, default=None,
                    help="Minimum peak total signal (default: 10x threshold)")
    ap.add_argument("--fdr", type=float, default=0.05,
                    help="False discovery rate threshold (default: 0.05)")
    
    # Pairing parameters
    ap.add_argument("--nt-window", type=int, default=1000,
                    help="Maximum edge-to-edge gap for pairing peaks (default: 1000bp)")
    ap.add_argument("--bin-gap", type=int, default=100,
                    help="Maximum gap within peaks (default: 100bp)")
    ap.add_argument("--balance", type=float, default=0.0,
                    help="Minimum strand balance for initial pairing (default: 0.0, filtering done statistically)")
    
    # Performance
    ap.add_argument("--ncores", type=int, default=1,
                    help="Number of CPU cores (default: 1, currently unused)")
    
    # Output options
    ap.add_argument("--report", default=None,
                    help="QC report output file (default: <out>_qc.txt)")
    ap.add_argument("--no-report", action="store_true",
                    help="Disable QC report generation")
    ap.add_argument("--write-summary", default=None,
                    help="Path to summary TSV file (for pipeline integration)")
    ap.add_argument("--quiet", action="store_true",
                    help="Suppress progress messages")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    
    return ap.parse_args()


def log(msg: str, quiet: bool = False):
    """Print timestamped log message to stderr."""
    if not quiet:
        timestamp = time.strftime('%H:%M:%S')
        sys.stderr.write(f"[{timestamp}] {msg}\n")
        sys.stderr.flush()


def natural_sort_key(s: str) -> List:
    """Natural sorting key for chromosome names (chr1, chr2, ..., chr10, ...)."""
    import re
    return [int(c) if c.isdigit() else c.lower() 
            for c in re.split(r'([0-9]+)', str(s))]


class FastBedGraph:
    """
    Chromosome-indexed bedGraph for fast region queries.
    
    Organizes bedGraph by chromosome for O(1) chromosome lookup
    followed by fast pandas boolean indexing.
    """
    
    def __init__(self, df: pd.DataFrame):
        """Build chromosome index from bedGraph dataframe."""
        self.by_chrom = {}
        for chrom, group in df.groupby('chr', sort=False):
            group = group.sort_values('start').reset_index(drop=True)
            self.by_chrom[str(chrom)] = group
    
    def query_sum(self, chrom: str, start: int, end: int) -> float:
        """
        Query total signal in genomic region.
        
        Args:
            chrom: Chromosome name
            start: Region start (0-based)
            end: Region end (exclusive)
            
        Returns:
            Sum of absolute signal values in region
        """
        if chrom not in self.by_chrom:
            return 0.0
        
        df = self.by_chrom[chrom]
        mask = (df['end'] > start) & (df['start'] < end)
        
        if not mask.any():
            return 0.0
        
        return float(df.loc[mask, 'sig'].abs().sum())


def read_bedgraph(path: str, role: str, quiet: bool = False) -> pd.DataFrame:
    """
    Read and process bedGraph file.
    
    Args:
        path: Path to bedGraph file (.gz supported)
        role: 'pos' or 'neg' for strand
        quiet: Suppress logging
        
    Returns:
        DataFrame with columns: chr, start, end, sig
    """
    log(f"Loading {role} strand: {os.path.basename(path)}", quiet)
    
    try:
        df = pd.read_csv(
            path, sep="\t", header=None, comment="#",
            names=["chr", "start", "end", "sig"],
            dtype={"chr": str, "start": int, "end": int, "sig": float},
            compression="infer"
        )
    except Exception as e:
        log(f"ERROR: Failed to read {path}: {e}")
        sys.exit(1)
    
    # Validate and clean
    df = df.loc[(df["end"] > df["start"]) & df["chr"].notna()].copy()
    
    if df.empty:
        log(f"ERROR: Empty bedGraph file: {path}")
        sys.exit(1)
    
    # Ensure correct sign convention
    if role == "pos":
        df["sig"] = df["sig"].abs()
    else:  # neg
        df["sig"] = -df["sig"].abs()
    
    log(f"  Loaded {len(df):,} bins", quiet)
    return df


def auto_calibrate(pos_df: pd.DataFrame, neg_df: pd.DataFrame, 
                   quiet: bool = False) -> Dict[str, float]:
    """
    Auto-calibrate detection thresholds from background distribution.
    
    Samples 100K random bins to estimate background, then sets:
    - threshold = 95th percentile of background
    - sum_thr = 10x threshold
    
    Args:
        pos_df: Positive strand bedGraph
        neg_df: Negative strand bedGraph
        quiet: Suppress logging
        
    Returns:
        Dict with keys: threshold, sum_thr, bg_mean, bg_std, bg_p95
    """
    log("Auto-calibrating detection thresholds...", quiet)
    
    sample_size = min(100000, len(pos_df) + len(neg_df))
    
    all_signals = np.concatenate([
        pos_df['sig'].abs().sample(
            min(sample_size//2, len(pos_df)), 
            random_state=42
        ).values,
        neg_df['sig'].abs().sample(
            min(sample_size//2, len(neg_df)), 
            random_state=42
        ).values
    ])
    
    bg_mean = float(np.mean(all_signals))
    bg_std = float(np.std(all_signals))
    bg_p95 = float(np.percentile(all_signals, 95))
    
    threshold = bg_p95
    sum_thr = bg_p95 * 10
    
    log(f"  Background mean: {bg_mean:.2f}, std: {bg_std:.2f}", quiet)
    log(f"  Background 95th percentile: {bg_p95:.2f}", quiet)
    log(f"  Recommended threshold: {threshold:.2f}", quiet)
    log(f"  Recommended sum_thr: {sum_thr:.2f}", quiet)
    
    return {
        'threshold': threshold,
        'sum_thr': sum_thr,
        'bg_mean': bg_mean,
        'bg_std': bg_std,
        'bg_p95': bg_p95
    }


def call_peaks(df: pd.DataFrame, threshold: float, sum_thr: float,
               bin_gap: int) -> pd.DataFrame:
    """
    Call peak blocks from bedGraph.
    
    Groups consecutive bins above threshold into peaks, then filters
    by total signal.
    
    Args:
        df: bedGraph dataframe
        threshold: Per-bin threshold
        sum_thr: Minimum peak total signal
        bin_gap: Maximum gap within peak
        
    Returns:
        DataFrame with columns: chr, start, end, height, signal
    """
    if df.empty:
        return pd.DataFrame(columns=["chr", "start", "end", "height", "signal"])
    
    # Filter bins by threshold
    v = np.abs(df["sig"].to_numpy())
    keep = v >= threshold
    
    if not keep.any():
        return pd.DataFrame(columns=["chr", "start", "end", "height", "signal"])
    
    # Extract passing bins
    chr_a = df["chr"].to_numpy()[keep]
    start = df["start"].to_numpy()[keep]
    end = df["end"].to_numpy()[keep]
    signal = v[keep]
    
    # Group into contiguous blocks
    grp = np.concatenate([
        [True],
        (chr_a[1:] != chr_a[:-1]) | (start[1:] > end[:-1] + bin_gap)
    ]).cumsum()
    
    # Aggregate blocks
    peaks = pd.DataFrame({
        "chr": pd.Series(chr_a).groupby(grp).first(),
        "start": pd.Series(start).groupby(grp).min(),
        "end": pd.Series(end).groupby(grp).max(),
        "height": pd.Series(signal).groupby(grp).max(),
        "signal": pd.Series(signal).groupby(grp).sum(),
    }).reset_index(drop=True)
    
    # Filter by total signal
    peaks = peaks.loc[peaks["signal"] >= sum_thr].reset_index(drop=True)
    
    return peaks


def pair_peaks_on_chromosome(
    chrom: str,
    pos_peaks: pd.DataFrame,
    neg_peaks: pd.DataFrame,
    nt_window: int,
    balance: float
) -> pd.DataFrame:
    """
    Pair positive and negative peaks on a single chromosome.
    
    Finds all positive-negative peak pairs that are:
    - Within nt_window bp of each other (edge-to-edge or overlapping)
    - Meet minimum balance requirement (if balance > 0)
    
    Args:
        chrom: Chromosome name
        pos_peaks: Positive peaks on this chromosome
        neg_peaks: Negative peaks on this chromosome  
        nt_window: Maximum edge-to-edge distance
        balance: Minimum balance ratio (0.0 = no filter)
        
    Returns:
        DataFrame with columns: chr, start, end, total
    """
    if pos_peaks.empty or neg_peaks.empty:
        return pd.DataFrame(columns=["chr", "start", "end", "total"])
    
    ps = pos_peaks["start"].to_numpy()
    pe = pos_peaks["end"].to_numpy()
    pg = pos_peaks["signal"].to_numpy()
    ph = pos_peaks["height"].to_numpy()
    
    ns = neg_peaks["start"].to_numpy()
    ne = neg_peaks["end"].to_numpy()
    ng = neg_peaks["signal"].to_numpy()
    nh = neg_peaks["height"].to_numpy()
    
    pairs = []
    
    # For each positive peak, find compatible negative peaks
    for i in range(len(ps)):
        # Binary search for candidates
        search_margin = nt_window + 1000
        left = np.searchsorted(ns, ps[i] - search_margin, side="left")
        right = np.searchsorted(ns, pe[i] + search_margin, side="right")
        
        if right <= left:
            continue
        
        ns_slice = ns[left:right]
        ne_slice = ne[left:right]
        ng_slice = ng[left:right]
        nh_slice = nh[left:right]
        
        # Calculate overlaps and gaps
        overlap = np.maximum(0, np.minimum(pe[i], ne_slice) - np.maximum(ps[i], ns_slice))
        gap = np.maximum(0, np.maximum(ps[i], ns_slice) - np.minimum(pe[i], ne_slice))
        
        # Distance filter
        dist_mask = (overlap > 0) | (gap <= nt_window)
        if not dist_mask.any():
            continue
        
        ns_slice = ns_slice[dist_mask]
        ne_slice = ne_slice[dist_mask]
        ng_slice = ng_slice[dist_mask]
        nh_slice = nh_slice[dist_mask]
        
        # Optional balance filter
        if balance > 0.0:
            balance_ratio = np.minimum(ph[i], nh_slice) / np.maximum(ph[i], nh_slice)
            bal_mask = balance_ratio >= balance
            if not bal_mask.any():
                continue
            ns_slice = ns_slice[bal_mask]
            ne_slice = ne_slice[bal_mask]
            ng_slice = ng_slice[bal_mask]
        
        # Create paired regions
        start = np.minimum(ps[i], ns_slice)
        end = np.maximum(pe[i], ne_slice)
        total = pg[i] + ng_slice
        
        for s, e, t in zip(start, end, total):
            pairs.append((chrom, int(s), int(e), float(t)))
    
    if not pairs:
        return pd.DataFrame(columns=["chr", "start", "end", "total"])
    
    return pd.DataFrame(pairs, columns=["chr", "start", "end", "total"])


def extract_features(
    paired_df: pd.DataFrame,
    pos_idx: FastBedGraph,
    neg_idx: FastBedGraph,
    quiet: bool = False
) -> pd.DataFrame:
    """
    Extract features for statistical scoring.
    
    For each paired region, computes:
    - Total signal and strand-specific sums
    - Bayesian balance score (Beta-Binomial model)
    - Local background and signal-to-background ratio
    - Region width and signal density
    
    Args:
        paired_df: Paired regions
        pos_idx: Indexed positive bedGraph
        neg_idx: Indexed negative bedGraph
        quiet: Suppress logging
        
    Returns:
        DataFrame with feature columns
    """
    log(f"Extracting features for {len(paired_df):,} regions...", quiet)
    
    try:
        from scipy.stats import beta as beta_dist
        has_scipy = True
    except ImportError:
        has_scipy = False
        log("  WARNING: scipy not available, using simple balance", quiet)
    
    features_list = []
    
    for idx, row in paired_df.iterrows():
        if (idx + 1) % 5000 == 0:
            log(f"  Progress: {idx+1:,} / {len(paired_df):,}", quiet)
        
        chrom = str(row['chr'])
        start = int(row['start'])
        end = int(row['end'])
        total_from_pairing = float(row['total'])
        
        # Query strand-specific signals
        pos_sum = max(pos_idx.query_sum(chrom, start, end), 0.01)
        neg_sum = max(neg_idx.query_sum(chrom, start, end), 0.01)
        
        total_signal = max(total_from_pairing, pos_sum + neg_sum)
        
        # Bayesian balance score
        if has_scipy and total_signal > 1:
            try:
                # Beta(2,2) prior, posterior mean balance probability
                alpha_post = 2 + pos_sum
                beta_post = 2 + neg_sum
                balance_bayesian = float(
                    beta_dist.cdf(0.7, alpha_post, beta_post) - 
                    beta_dist.cdf(0.3, alpha_post, beta_post)
                )
            except:
                balance_bayesian = min(pos_sum, neg_sum) / max(pos_sum, neg_sum)
        else:
            balance_bayesian = min(pos_sum, neg_sum) / max(pos_sum, neg_sum)
        
        # Local background (±5kb excluding peak)
        window = 5000
        local_start = max(0, start - window)
        local_end = end + window
        
        pos_bg_sum = (pos_idx.query_sum(chrom, local_start, start) + 
                      pos_idx.query_sum(chrom, end, local_end))
        neg_bg_sum = (neg_idx.query_sum(chrom, local_start, start) + 
                      neg_idx.query_sum(chrom, end, local_end))
        
        bg_length = 2 * window
        local_bg = (pos_bg_sum + neg_bg_sum) / max(1, bg_length)
        local_bg = max(local_bg, 0.1)
        
        # Derived features
        width = end - start
        signal_to_bg = total_signal / max(0.1, local_bg * width)
        
        features = {
            'total_signal': total_signal,
            'log_total': np.log1p(total_signal),
            'balance_bayesian': balance_bayesian,
            'width': width,
            'signal_density': total_signal / max(1, width),
            'local_bg': local_bg,
            'signal_to_bg': signal_to_bg,
            'log_snr': np.log1p(signal_to_bg)
        }
        
        features_list.append(features)
    
    return pd.DataFrame(features_list)


def score_with_mixture_model(
    features_df: pd.DataFrame,
    fdr_threshold: float,
    quiet: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Score regions using Gaussian Mixture Model and apply FDR control.
    
    Fits 2-component GMM on key features to identify true divergent TX
    vs noise/artifacts. Computes posterior probabilities and applies
    FDR-controlled filtering.
    
    Args:
        features_df: Feature matrix
        fdr_threshold: FDR cutoff (e.g., 0.05 = 5%)
        quiet: Suppress logging
        
    Returns:
        Tuple of (passing_mask, scores) where:
        - passing_mask: Boolean array of regions passing FDR
        - scores: Posterior probability for each region (0-1)
    """
    log("Fitting statistical model...", quiet)
    
    try:
        from sklearn.mixture import GaussianMixture
    except ImportError:
        log("ERROR: scikit-learn required for statistical scoring")
        log("Install with: pip install scikit-learn scipy")
        sys.exit(1)
    
    # Select features for model
    feature_cols = ['log_total', 'balance_bayesian', 'log_snr']
    
    X = features_df[feature_cols].copy()
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median())
    
    # Fit 2-component Gaussian mixture
    gmm = GaussianMixture(
        n_components=2,
        covariance_type='full',
        random_state=42,
        max_iter=200,
        n_init=10
    )
    gmm.fit(X.values)
    
    # Identify positive component (higher mean signal)
    means = gmm.means_[:, 0]  # log_total means
    pos_component = int(np.argmax(means))
    
    log(f"  Positive component: {pos_component} (mean log_total={gmm.means_[pos_component][0]:.2f})", quiet)
    log(f"  Negative component: {1-pos_component} (mean log_total={gmm.means_[1-pos_component][0]:.2f})", quiet)
    
    # Compute posterior probabilities
    probs = gmm.predict_proba(X.values)
    scores = probs[:, pos_component]
    
    # FDR control via Benjamini-Hochberg-like procedure
    sorted_idx = np.argsort(-scores)
    sorted_scores = scores[sorted_idx]
    
    cumulative_fp = np.cumsum(1 - sorted_scores)
    cumulative_calls = np.arange(1, len(sorted_scores) + 1)
    fdr = cumulative_fp / cumulative_calls
    
    passing = fdr <= fdr_threshold
    
    if not passing.any():
        log(f"  WARNING: No regions pass FDR={fdr_threshold}", quiet)
        log(f"  Taking top 10% by score as fallback", quiet)
        n_pass = max(100, len(sorted_scores) // 10)
    else:
        n_pass = np.where(passing)[0][-1] + 1
    
    mask = np.array([False] * len(features_df))
    mask[sorted_idx[:n_pass]] = True
    
    log(f"  → {n_pass:,} / {len(features_df):,} regions pass FDR (retention: {100*n_pass/len(features_df):.1f}%)", quiet)
    
    return mask, scores


def generate_qc_report(
    report_path: str,
    stats: Dict,
    features_df: pd.DataFrame,
    passing_mask: np.ndarray,
    scores: np.ndarray,
    calibration: Dict,
    actual_threshold: float,
    actual_sum_thr: float
):
    """Generate QC report with pipeline statistics."""
    with open(report_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("DIVERGENT TRANSCRIPTION DETECTION - QC REPORT\n")
        f.write("="*70 + "\n\n")
        
        f.write("Threshold Parameters\n")
        f.write("-"*70 + "\n")
        f.write(f"  Threshold used:           {actual_threshold:.3f}\n")
        f.write(f"  Sum threshold used:       {actual_sum_thr:.3f}\n")
        f.write(f"  (Auto-calibrated bg p95:  {calibration['bg_p95']:.3f})\n")
        f.write(f"  (Auto-calibrated values:  {calibration['threshold']:.3f}, {calibration['sum_thr']:.3f})\n\n")
        
        f.write("Peak Calling\n")
        f.write("-"*70 + "\n")
        f.write(f"  Positive peaks:           {stats['pos_peaks']:,}\n")
        f.write(f"  Negative peaks:           {stats['neg_peaks']:,}\n")
        f.write(f"  Strand ratio (pos/neg):   {stats['pos_peaks']/max(1,stats['neg_peaks']):.2f}\n\n")
        
        f.write("Pairing\n")
        f.write("-"*70 + "\n")
        f.write(f"  Candidate pairs:          {stats['paired']:,}\n\n")
        
        f.write("Statistical Filtering\n")
        f.write("-"*70 + "\n")
        f.write(f"  FDR threshold:            {stats['fdr']:.4f}\n")
        f.write(f"  Regions passing:          {np.sum(passing_mask):,}\n")
        f.write(f"  Retention rate:           {100*np.sum(passing_mask)/len(passing_mask):.1f}%\n\n")
        
        f.write("Score Distribution\n")
        f.write("-"*70 + "\n")
        f.write(f"  Min:                      {scores.min():.4f}\n")
        f.write(f"  25th percentile:          {np.percentile(scores, 25):.4f}\n")
        f.write(f"  Median:                   {np.median(scores):.4f}\n")
        f.write(f"  75th percentile:          {np.percentile(scores, 75):.4f}\n")
        f.write(f"  Max:                      {scores.max():.4f}\n\n")
        
        f.write("Passing Regions - Feature Summary\n")
        f.write("-"*70 + "\n")
        passing_features = features_df[passing_mask]
        
        f.write(f"  Total signal:\n")
        f.write(f"    Mean:                   {passing_features['total_signal'].mean():.0f}\n")
        f.write(f"    Median:                 {passing_features['total_signal'].median():.0f}\n")
        
        f.write(f"  Balance (Bayesian):\n")
        f.write(f"    Mean:                   {passing_features['balance_bayesian'].mean():.3f}\n")
        f.write(f"    Median:                 {passing_features['balance_bayesian'].median():.3f}\n")
        
        f.write(f"  Width (bp):\n")
        f.write(f"    Mean:                   {passing_features['width'].mean():.0f}\n")
        f.write(f"    Median:                 {passing_features['width'].median():.0f}\n")
        
        f.write(f"  Signal-to-background:\n")
        f.write(f"    Mean:                   {passing_features['signal_to_bg'].mean():.2f}\n")
        f.write(f"    Median:                 {passing_features['signal_to_bg'].median():.2f}\n")
        
        f.write("\n" + "="*70 + "\n")
        f.write("Score Interpretation:\n")
        f.write("  Score = P(true divergent TX | features) from Gaussian mixture model\n")
        f.write("  Higher scores = higher confidence\n")
        f.write("  All regions collectively pass FDR control\n")
        f.write("="*70 + "\n")


def main():
    """Main pipeline."""
    args = parse_args()
    start_time = time.time()
    
    print("="*70)
    print("DIVERGENT TRANSCRIPTION DETECTOR v" + __version__)
    print("="*70)
    print(f"Sample:           {args.sample}")
    print(f"Positive strand:  {os.path.basename(args.pos)}")
    print(f"Negative strand:  {os.path.basename(args.neg)}")
    print(f"Output:           {args.out}")
    print(f"FDR threshold:    {args.fdr}")
    print()
    
    # Check dependencies
    try:
        from sklearn.mixture import GaussianMixture
    except ImportError:
        print("ERROR: scikit-learn required")
        print("Install with: pip install scikit-learn scipy")
        sys.exit(1)
    
    # Load bedGraphs
    log("[1/7] Loading bedGraphs...")
    pos_df = read_bedgraph(args.pos, "pos", args.quiet)
    neg_df = read_bedgraph(args.neg, "neg", args.quiet)
    
    # Build chromosome indices
    log("[2/7] Indexing chromosomes...")
    pos_idx = FastBedGraph(pos_df)
    neg_idx = FastBedGraph(neg_df)
    log(f"  Indexed {len(pos_idx.by_chrom)} chromosomes", args.quiet)
    
    # Auto-calibrate thresholds
    log("[3/7] Calibrating thresholds...")
    calibration = auto_calibrate(pos_df, neg_df, args.quiet)
    
    # Use user-specified or auto-calibrated values
    if args.threshold is not None:
        threshold = args.threshold
        log(f"  Using user-specified threshold: {threshold:.2f}", args.quiet)
    else:
        threshold = calibration['threshold']
        log(f"  Using auto-calibrated threshold: {threshold:.2f}", args.quiet)
    
    if args.sum_thr is not None:
        sum_thr = args.sum_thr
        log(f"  Using user-specified sum_thr: {sum_thr:.2f}", args.quiet)
    else:
        sum_thr = calibration['sum_thr']
        log(f"  Using auto-calibrated sum_thr: {sum_thr:.2f}", args.quiet)
    
    # Call peaks on both strands
    log("[4/7] Calling peaks...")
    pos_peaks = call_peaks(pos_df, threshold, sum_thr, args.bin_gap)
    neg_peaks = call_peaks(neg_df, threshold, sum_thr, args.bin_gap)
    log(f"  Positive: {len(pos_peaks):,} peaks", args.quiet)
    log(f"  Negative: {len(neg_peaks):,} peaks", args.quiet)
    
    # Pair peaks chromosome-by-chromosome
    log("[5/7] Pairing peaks...")
    chroms = sorted(
        set(pos_peaks['chr'].unique()) & set(neg_peaks['chr'].unique()),
        key=natural_sort_key
    )
    log(f"  Processing {len(chroms)} chromosomes", args.quiet)
    
    pieces = []
    for chrom in chroms:
        result = pair_peaks_on_chromosome(
            chrom,
            pos_peaks[pos_peaks['chr'] == chrom].reset_index(drop=True),
            neg_peaks[neg_peaks['chr'] == chrom].reset_index(drop=True),
            args.nt_window,
            args.balance
        )
        if not result.empty:
            pieces.append(result)
    
    if not pieces:
        log("ERROR: No paired peaks found")
        sys.exit(1)
    
    paired = pd.concat(pieces, ignore_index=True)
    log(f"  → {len(paired):,} candidate pairs", args.quiet)
    
    # Extract features
    log("[6/7] Extracting features...")
    features_df = extract_features(paired, pos_idx, neg_idx, args.quiet)
    
    # Statistical scoring and FDR filtering
    log("[7/7] Statistical scoring...")
    passing_mask, scores = score_with_mixture_model(
        features_df, args.fdr, args.quiet
    )
    
    # Output results
    final = paired[passing_mask].copy()
    final['score'] = scores[passing_mask]
    final = final.sort_values('score', ascending=False).reset_index(drop=True)
    
    final[['chr', 'start', 'end', 'total', 'score']].to_csv(
        args.out, sep='\t', header=False, index=False, float_format='%.4f'
    )
    
    log(f"\nOutput: {args.out} ({len(final):,} regions)", args.quiet)
    
    # Generate QC report
    if not args.no_report:
        report_path = args.report if args.report else args.out.replace('.bed', '_qc.txt')
        
        stats = {
            'pos_peaks': len(pos_peaks),
            'neg_peaks': len(neg_peaks),
            'paired': len(paired),
            'fdr': args.fdr
        }
        
        generate_qc_report(
            report_path, stats, features_df, passing_mask, 
            scores, calibration, threshold, sum_thr
        )
        log(f"QC report: {report_path}", args.quiet)
    
    # Generate summary TSV (for pipeline integration)
    if args.write_summary:
        elapsed = time.time() - start_time
        with open(args.write_summary, 'w') as f:
            f.write("sample\tn_pos_pk\tn_neg_pk\tn_pairs_raw\tn_dt\twall_s\n")
            f.write(f"{args.sample}\t{len(pos_peaks)}\t{len(neg_peaks)}\t"
                   f"{len(paired)}\t{len(final)}\t{elapsed:.1f}\n")
        log(f"Summary: {args.write_summary}", args.quiet)
    
    elapsed = time.time() - start_time
    
    print("\n" + "="*70)
    print("COMPLETE")
    print(f"  Sample:             {args.sample}")
    print(f"  Divergent regions:  {len(final):,}")
    print(f"  Mean score:         {final['score'].mean():.3f}")
    print(f"  Score range:        [{final['score'].min():.4f}, {final['score'].max():.4f}]")
    print(f"  Elapsed time:       {elapsed/60:.1f} min")
    print("="*70)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n\nInterrupted by user")
        sys.exit(130)
    except Exception as e:
        print(f"\n\nERROR: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
