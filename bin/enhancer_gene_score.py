#!/usr/bin/env python3
"""
enhancer_gene_score.py

Assign a continuous enhancer-vs-gene score in [0,1] to each divergent
transcription site, based on positional features relative to annotated genes.

Inputs
------
- divergent_transcription.bed  (BED5: chr, start, end, total_signal, dt_score)
- genes.tsv                    (from gtf_to_catalog.py)
- tss.bed                      (1bp TSS sites; BED6)

Outputs
-------
1) TSV with features and scores (one row per divergent site)
2) BED-like file with enhancer_score appended

The score has the interpretation:
  0   → perfectly gene-like (classic promoter / gene-proximal)
  1   → perfectly enhancer-like (distal / intergenic)

Notes
-----
- This script is intentionally self-contained and light-weight:
  it computes interpretable features, a simple rule-based score,
  and then (optionally) trains a small gradient-boosted tree on
  pseudo-labels derived from the rule for a smoother score.
"""

from __future__ import annotations

import argparse
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class GeneRecord:
    chrom: str
    start: int
    end: int
    tss: int


@dataclass
class GeneIndex:
    """Per-chromosome index for fast positional queries."""

    starts: np.ndarray   # gene starts (sorted)
    ends: np.ndarray     # gene ends (aligned with starts)
    tss: np.ndarray      # TSS positions (subset or same order as starts)


def build_gene_indices(genes_tsv: str) -> Dict[str, GeneIndex]:
    """
    Build per-chromosome indices from genes.tsv.

    genes.tsv columns (at minimum):
      gene_id, gene_name, chr, strand, start, end, tss, tes, ...
    """
    df = pd.read_csv(genes_tsv, sep="\t")

    # Normalize column names
    cols = {c.lower(): c for c in df.columns}
    required = ["chr", "start", "end", "tss"]
    missing = [c for c in required if c not in cols]
    if missing:
        raise ValueError(f"genes.tsv missing required columns: {missing}")

    chr_col = cols["chr"]
    start_col = cols["start"]
    end_col = cols["end"]
    tss_col = cols["tss"]

    gene_indices: Dict[str, GeneIndex] = {}

    for chrom, sub in df.groupby(chr_col, sort=False):
        starts = sub[start_col].astype(int).to_numpy()
        ends = sub[end_col].astype(int).to_numpy()
        tss = sub[tss_col].astype(int).to_numpy()

        order = np.argsort(starts, kind="mergesort")
        starts_sorted = starts[order]
        ends_sorted = ends[order]
        tss_sorted = tss[order]

        gene_indices[str(chrom)] = GeneIndex(
            starts=starts_sorted,
            ends=ends_sorted,
            tss=tss_sorted,
        )

    return gene_indices


def build_tss_index(tss_bed: str) -> Dict[str, np.ndarray]:
    """
    Build per-chromosome sorted arrays of TSS positions from BED.

    Assumes 0-based BED coordinates; we use the start column as position.
    """
    pos_by_chrom: Dict[str, List[int]] = {}
    with open(tss_bed) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 3:
                continue
            chrom = f[0]
            try:
                start = int(f[1])
            except ValueError:
                continue
            pos_by_chrom.setdefault(chrom, []).append(start)

    idx: Dict[str, np.ndarray] = {}
    for chrom, arr in pos_by_chrom.items():
        a = np.array(sorted(arr), dtype=int)
        idx[chrom] = a
    return idx


# ---------------------------------------------------------------------------
# Feature computation
# ---------------------------------------------------------------------------

def nearest_tss_distance(
    chrom: str,
    center: int,
    tss_index: Dict[str, np.ndarray],
) -> float:
    arr = tss_index.get(chrom)
    if arr is None or arr.size == 0:
        return float("inf")
    # Find nearest via binary search
    i = np.searchsorted(arr, center)
    candidates = []
    if i > 0:
        candidates.append(abs(center - int(arr[i - 1])))
    if i < arr.size:
        candidates.append(abs(center - int(arr[i])))
    return float(min(candidates)) if candidates else float("inf")


def gene_body_features_for_point(
    chrom: str,
    center: int,
    genes: Dict[str, GeneIndex],
) -> Tuple[bool, float, Optional[int], Optional[int], Optional[float]]:
    """
    Compute gene-related features for a single genomic point:

    Returns:
      - overlaps_gene: bool, whether point lies within any gene body
      - dist_to_gene_body: float, 0 if overlapping, else distance to nearest body
      - upstream_gene_end: int or None
      - downstream_gene_start: int or None
      - intergenic_rel_pos: float in [0,1] if between two genes, else None
    """
    gi = genes.get(chrom)
    if gi is None or gi.starts.size == 0:
        return False, float("inf"), None, None, None

    starts = gi.starts
    ends = gi.ends

    # Locate rightmost gene with start <= center
    i = int(np.searchsorted(starts, center, side="right") - 1)

    overlaps = False
    dist_body = float("inf")
    upstream_end: Optional[int] = None
    downstream_start: Optional[int] = None
    rel_pos: Optional[float] = None

    # Check overlap with gene i (max start <= center)
    if i >= 0:
        g_start = int(starts[i])
        g_end = int(ends[i])
        if g_start <= center <= g_end:
            overlaps = True
            dist_body = 0.0
        upstream_end = g_end

    # Candidate downstream gene: i+1 (smallest start > center)
    j = i + 1
    if j < starts.size:
        downstream_start = int(starts[j])

    # If not overlapping, estimate distance to nearest gene body
    if not overlaps:
        dists = []
        if upstream_end is not None:
            dists.append(max(0, center - upstream_end))
        if downstream_start is not None:
            dists.append(max(0, downstream_start - center))
        if dists:
            dist_body = float(min(dists))

    # Intergenic relative position: between upstream_end and downstream_start
    if upstream_end is not None and downstream_start is not None:
        if upstream_end < center < downstream_start:
            span = downstream_start - upstream_end
            if span > 0:
                rel_pos = float((center - upstream_end) / span)

    return overlaps, dist_body, upstream_end, downstream_start, rel_pos


def compute_features(
    divergent_bed: str,
    genes_tsv: str,
    tss_bed: str,
    prom_up: int,
    prom_down: int,
    tss_active_pm: int,
) -> pd.DataFrame:
    """
    Compute positional features for each divergent site.

    The divergent BED is expected to have 5 columns:
      chr, start, end, total_signal, dt_score
    """
    gene_idx = build_gene_indices(genes_tsv)
    tss_idx = build_tss_index(tss_bed)

    records = []
    with open(divergent_bed) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            chrom = f[0]
            try:
                start = int(f[1])
                end = int(f[2])
                total_signal = float(f[3])
                dt_score = float(f[4])
            except ValueError:
                continue
            if end <= start:
                continue
            center = (start + end) // 2
            length = end - start

            # Distances to TSS and promoter windows
            dist_tss = nearest_tss_distance(chrom, center, tss_idx)

            # Active gene detection window: TSS ± tss_active_pm
            # Approximate by checking against per-chromosome TSS list
            arr_tss = tss_idx.get(chrom)
            in_active_window = False
            if arr_tss is not None and arr_tss.size > 0:
                i = np.searchsorted(arr_tss, center)
                candidates = []
                if i > 0:
                    candidates.append(arr_tss[i - 1])
                if i < arr_tss.size:
                    candidates.append(arr_tss[i])
                for pos in candidates:
                    if abs(int(pos) - center) <= tss_active_pm:
                        in_active_window = True
                        break

            # Promoter geometry (TSS - prom_up .. TSS + prom_down)
            in_promoter_window = False
            if arr_tss is not None and arr_tss.size > 0:
                # Same nearest-TSS candidates as above
                i = np.searchsorted(arr_tss, center)
                candidates = []
                if i > 0:
                    candidates.append(arr_tss[i - 1])
                if i < arr_tss.size:
                    candidates.append(arr_tss[i])
                for pos in candidates:
                    prom_start = int(pos) - prom_up
                    prom_end = int(pos) + prom_down
                    if prom_start <= center <= prom_end:
                        in_promoter_window = True
                        break

            # Gene body / intergenic features
            (
                overlaps_gene,
                dist_gene_body,
                upstream_end,
                downstream_start,
                rel_pos,
            ) = gene_body_features_for_point(chrom, center, gene_idx)

            # Local TSS density (±100 kb)
            local_tss_density = 0
            window = 100_000
            if arr_tss is not None and arr_tss.size > 0:
                lo = max(0, center - window)
                hi = center + window
                left = np.searchsorted(arr_tss, lo, side="left")
                right = np.searchsorted(arr_tss, hi, side="right")
                local_tss_density = int(max(0, right - left))

            records.append(
                {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "center": center,
                    "length": length,
                    "total_signal": total_signal,
                    "dt_score": dt_score,
                    "dist_to_nearest_tss": dist_tss,
                    "in_active_tss_window": int(in_active_window),
                    "in_promoter_window": int(in_promoter_window),
                    "overlaps_gene_body": int(overlaps_gene),
                    "dist_to_gene_body": dist_gene_body,
                    "upstream_gene_end": upstream_end if upstream_end is not None else -1,
                    "downstream_gene_start": (
                        downstream_start if downstream_start is not None else -1
                    ),
                    "intergenic_rel_pos": rel_pos if rel_pos is not None else np.nan,
                    "local_tss_density_100kb": local_tss_density,
                }
            )

    if not records:
        return pd.DataFrame(
            columns=[
                "chrom",
                "start",
                "end",
                "center",
                "length",
                "total_signal",
                "dt_score",
                "dist_to_nearest_tss",
                "in_active_tss_window",
                "in_promoter_window",
                "overlaps_gene_body",
                "dist_to_gene_body",
                "upstream_gene_end",
                "downstream_gene_start",
                "intergenic_rel_pos",
                "local_tss_density_100kb",
            ]
        )

    df = pd.DataFrame.from_records(records)
    return df


# ---------------------------------------------------------------------------
# Rule-based score
# ---------------------------------------------------------------------------

def rule_based_enhancer_score(df: pd.DataFrame) -> pd.Series:
    """
    Compute a heuristic enhancer_score in [0,1] encoding simple priors:

    - Strongly gene-like near TSS/promoters and inside gene bodies
    - More enhancer-like when intergenic and far from TSS,
      especially mid-way between neighboring genes
    """
    s = np.full(len(df), 0.5, dtype=float)

    dist_tss = df["dist_to_nearest_tss"].to_numpy()
    in_prom = df["in_promoter_window"].to_numpy().astype(bool)
    in_active = df["in_active_tss_window"].to_numpy().astype(bool)
    overlaps_gene = df["overlaps_gene_body"].to_numpy().astype(bool)
    dist_gene = df["dist_to_gene_body"].to_numpy()
    rel_pos = df["intergenic_rel_pos"].to_numpy()
    local_tss = df["local_tss_density_100kb"].to_numpy()

    # 1) Strongly down-weight promoter-proximal sites
    mask_prom = in_prom
    s[mask_prom] -= 0.4

    # Extra penalty if also inside gene body (canonical promoters)
    mask_prom_gene = mask_prom & overlaps_gene
    s[mask_prom_gene] -= 0.1

    # 2) Other gene-body sites (intronic / gene-proximal)
    #    Slightly gene-like, especially if close to TSS
    mask_gene_not_prom = overlaps_gene & ~mask_prom
    close_gene = (mask_gene_not_prom) & (dist_tss < 2000)
    far_gene = (mask_gene_not_prom) & (dist_tss >= 2000)
    s[close_gene] -= 0.2
    s[far_gene] -= 0.05

    # 3) Intergenic sites (no gene-body overlap)
    intergenic = ~overlaps_gene
    s[intergenic] += 0.15

    # 3a) Intergenic and mid-way between flanking genes
    mid_intergenic = intergenic & np.isfinite(rel_pos) & (rel_pos >= 0.3) & (rel_pos <= 0.7)
    s[mid_intergenic] += 0.2

    # 3b) Far from any TSS (>5 kb) → more enhancer-like
    far_from_tss = (dist_tss > 5000)
    s[intergenic & far_from_tss] += 0.1

    # 4) Very high local TSS density (promoter clusters)
    # Slightly pull towards gene-like, as many nearby promoters
    high_tss_density = local_tss >= 20
    s[high_tss_density] -= 0.05

    # 5) Sites in active TSS window but not overlapping gene body:
    # likely promoter-like but perhaps unannotated; gentle penalty.
    in_active_only = in_active & ~overlaps_gene
    s[in_active_only] -= 0.15

    # Clamp to [0,1]
    s = np.clip(s, 0.0, 1.0)
    return pd.Series(s, index=df.index, name="enhancer_score_rule")


# ---------------------------------------------------------------------------
# Pseudo-labels and ML model
# ---------------------------------------------------------------------------

def create_pseudo_labels(
    df: pd.DataFrame,
    rule_col: str = "enhancer_score_rule",
    high_thr: float = 0.8,
    low_thr: float = 0.2,
) -> pd.Series:
    """
    Create weak labels from the rule-based score:

    - label 1 (enhancer-like): rule_score >= high_thr and intergenic
    - label 0 (gene-like):     rule_score <= low_thr and promoter/gene-body
    - NaN: ambiguous samples (not used for training)
    """
    rule = df[rule_col].to_numpy()
    overlaps_gene = df["overlaps_gene_body"].to_numpy().astype(bool)
    in_prom = df["in_promoter_window"].to_numpy().astype(bool)

    labels = np.full(len(df), np.nan, dtype=float)

    intergenic = ~overlaps_gene

    pos_mask = intergenic & (rule >= high_thr)
    neg_mask = (overlaps_gene | in_prom) & (rule <= low_thr)

    labels[pos_mask] = 1.0
    labels[neg_mask] = 0.0

    return pd.Series(labels, index=df.index, name="pseudo_label")


def train_ml_model_and_score(
    df: pd.DataFrame,
    label_col: str = "pseudo_label",
    max_depth: int = 3,
    n_estimators: int = 200,
    random_state: int = 42,
) -> Optional[pd.Series]:
    """
    Train a small gradient-boosted tree on pseudo-labels and
    return predicted enhancer probabilities for all rows.
    """
    try:
        from sklearn.ensemble import GradientBoostingClassifier
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics import roc_auc_score
    except ImportError:
        # Fallback: no ML dependencies available
        return None

    y = df[label_col].to_numpy()
    mask = np.isfinite(y)
    if mask.sum() < 50:
        # Not enough pseudo-labeled points to fit a reasonable model
        return None

    # Feature set (exclude coordinates and directly derived labels/scores)
    feature_cols = [
        "length",
        "total_signal",
        "dt_score",
        "dist_to_nearest_tss",
        "in_active_tss_window",
        "in_promoter_window",
        "overlaps_gene_body",
        "dist_to_gene_body",
        "intergenic_rel_pos",
        "local_tss_density_100kb",
    ]

    X = df[feature_cols].copy()
    # Replace NaNs in intergenic_rel_pos (and any other) with column medians
    X = X.replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median(numeric_only=True))

    X_train = X.loc[mask]
    y_train = y[mask]

    clf = GradientBoostingClassifier(
        max_depth=max_depth,
        n_estimators=n_estimators,
        random_state=random_state,
    )

    # Optional: simple cross-validation for sanity; errors are non-fatal
    try:
        if len(np.unique(y_train)) > 1 and len(y_train) >= 100:
            skf = StratifiedKFold(n_splits=3, shuffle=True, random_state=random_state)
            aucs: List[float] = []
            for train_idx, val_idx in skf.split(X_train, y_train):
                clf_cv = GradientBoostingClassifier(
                    max_depth=max_depth,
                    n_estimators=n_estimators,
                    random_state=random_state,
                )
                clf_cv.fit(X_train.iloc[train_idx], y_train[train_idx])
                proba_val = clf_cv.predict_proba(X_train.iloc[val_idx])[:, 1]
                aucs.append(roc_auc_score(y_train[val_idx], proba_val))
            mean_auc = float(np.mean(aucs))
            print(
                f"[enhancer_gene_score] Weak-label CV AUC (3-fold, n={len(y_train)}): "
                f"{mean_auc:.3f}",
                file=sys.stderr,
            )
    except Exception as e:
        print(
            f"[enhancer_gene_score] WARNING: CV evaluation failed: {e}",
            file=sys.stderr,
        )

    # Fit final model on all pseudo-labeled data
    clf.fit(X_train, y_train)
    proba = clf.predict_proba(X)[:, 1]
    proba = np.clip(proba, 0.0, 1.0)

    return pd.Series(proba, index=df.index, name="enhancer_score_ml")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description=(
            "Compute enhancer-vs-gene scores for divergent transcription sites "
            "using positional features and weakly supervised ML."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    ap.add_argument(
        "--divergent",
        required=True,
        help="Input divergent_transcription.bed (BED5)",
    )
    ap.add_argument(
        "--genes",
        required=True,
        help="genes.tsv from gtf_to_catalog.py",
    )
    ap.add_argument(
        "--tss",
        required=True,
        help="tss.bed (1bp sites, BED6)",
    )
    ap.add_argument(
        "--prom-up",
        type=int,
        default=250,
        help="Promoter upstream window (bp) from TSS",
    )
    ap.add_argument(
        "--prom-down",
        type=int,
        default=250,
        help="Promoter downstream window (bp) from TSS",
    )
    ap.add_argument(
        "--tss-active-pm",
        type=int,
        default=500,
        help="Half-window for active TSS detection (bp)",
    )
    ap.add_argument(
        "--no-ml",
        action="store_true",
        help="Disable ML refinement; output rule-based score only",
    )
    ap.add_argument(
        "--out-tsv",
        required=True,
        help="Output TSV with features and scores",
    )
    ap.add_argument(
        "--out-bed",
        required=True,
        help=(
            "Output BED-like file with enhancer_score appended "
            "(chr, start, end, total_signal, dt_score, enhancer_score)"
        ),
    )

    return ap.parse_args()


def main() -> None:
    args = parse_args()

    div_path = Path(args.divergent)
    if not div_path.exists():
        raise SystemExit(f"Input divergent BED not found: {div_path}")

    print(
        f"[enhancer_gene_score] Computing features from {div_path.name} ...",
        file=sys.stderr,
    )
    df = compute_features(
        divergent_bed=str(div_path),
        genes_tsv=args.genes,
        tss_bed=args.tss,
        prom_up=args.prom_up,
        prom_down=args.prom_down,
        tss_active_pm=args.tss_active_pm,
    )

    if df.empty:
        print(
            "[enhancer_gene_score] No divergent sites found; writing empty outputs.",
            file=sys.stderr,
        )
        # Still write empty files with headers for robustness
        Path(args.out_tsv).write_text(
            "chrom\tstart\tend\tcenter\tlength\t"
            "total_signal\tdt_score\t"
            "dist_to_nearest_tss\tin_active_tss_window\tin_promoter_window\t"
            "overlaps_gene_body\tdist_to_gene_body\t"
            "upstream_gene_end\tdownstream_gene_start\t"
            "intergenic_rel_pos\tlocal_tss_density_100kb\t"
            "enhancer_score_rule\tenhancer_score\n",
            encoding="utf-8",
        )
        Path(args.out_bed).write_text("", encoding="utf-8")
        return

    # Rule-based score
    df["enhancer_score_rule"] = rule_based_enhancer_score(df)

    # Pseudo-labels
    df["pseudo_label"] = create_pseudo_labels(df)

    # ML refinement
    enhancer_score = df["enhancer_score_rule"].copy()
    if not args.no_ml:
        ml_scores = train_ml_model_and_score(df)
        if ml_scores is not None:
            enhancer_score = ml_scores
            print(
                "[enhancer_gene_score] Using ML-refined enhancer_score.",
                file=sys.stderr,
            )
        else:
            print(
                "[enhancer_gene_score] ML refinement unavailable; "
                "falling back to rule-based enhancer_score.",
                file=sys.stderr,
            )
    else:
        print(
            "[enhancer_gene_score] ML refinement disabled (--no-ml); "
            "using rule-based enhancer_score.",
            file=sys.stderr,
        )

    df["enhancer_score"] = enhancer_score

    # Write TSV
    out_tsv_path = Path(args.out_tsv)
    out_tsv_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_tsv_path, sep="\t", index=False)

    # Write BED-like file: chr, start, end, total_signal, dt_score, enhancer_score
    out_bed_path = Path(args.out_bed)
    with open(div_path) as src, open(out_bed_path, "w") as dst:
        idx = 0
        for ln in src:
            if not ln.strip() or ln.startswith(("#", "track", "browser")):
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            if idx >= len(df):
                break
            score_val = float(df.loc[df.index[idx], "enhancer_score"])
            dst.write(
                f"{f[0]}\t{f[1]}\t{f[2]}\t{f[3]}\t{f[4]}\t{score_val:.4f}\n"
            )
            idx += 1

    print(
        f"[enhancer_gene_score] Wrote features to {out_tsv_path} "
        f"and scored BED to {out_bed_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

