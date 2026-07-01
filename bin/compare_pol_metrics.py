#!/usr/bin/env python3
# =============================================================================
# compare_pol_metrics.py — Pol-II Metrics Aggregation and Comparison
# =============================================================================
#
# Purpose:
#   Aggregates per-sample Pol-II metrics into cohort-level summaries with
#   optional differential contrasts and visualizations.
#
# Features:
#   • Merges per-sample metrics into tidy long-format table
#   • Efficient column selection (usecols) for large files
#   • Replicate collapsing via median aggregation
#   • Differential contrasts (condition or timepoint comparisons)
#   • Log2 fold-change calculations
#   • Heatmaps for top variable genes
#   • MA plots for differential analysis
#
# Inputs:
#   • Samples manifest (TSV): sample_id, condition, timepoint, replicate, file
#   • Per-sample metrics files: gene-level pausing and expression data
#
# Outputs:
#   • Merged tidy table (TSV): All samples combined
#   • Contrasts table (TSV): Differential analysis results (optional)
#   • Plots directory (PNG): Heatmaps and MA plots (optional)
#
# Contrast Specification:
#   Format: "variable:numerator,denominator"
#   Examples:
#     • "condition:treatment,control"
#     • "timepoint:24h,0h"
#
# =============================================================================

from __future__ import annotations
import argparse
import datetime
import math
import os
import sys
from pathlib import Path
from typing import List, Tuple, Optional

import numpy as np
import pandas as pd
from pandas.api.types import is_numeric_dtype
from scipy import stats

# Matplotlib backend must be set before importing pyplot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# CONSTANTS
# =============================================================================

VERSION = "2.1.0"
LOG_PREFIX = "[COMPARE]"

# Required columns in samples manifest
REQUIRED_MANIFEST_COLS = ["sample_id", "condition", "timepoint", "replicate", "file"]

# Required metric columns in per-sample files
REQUIRED_METRIC_COLS = ["gene_id", "gene_name", "pi_len_norm", "pi_raw", "body_cpm", "tss_cpm"]

# Metrics available for analysis
AVAILABLE_METRICS = ["pi_len_norm", "pi_raw", "body_cpm", "tss_cpm"]

# --- Differential contrast defaults (assay-agnostic; no biology hardcoded) ---
# DEFAULT_PRIOR_COUNT: additive prior ("prior count", edgeR-style) used in the
#   log2 fold-change so that genes with near-zero signal in one group cannot
#   produce extreme fold-changes. log2((num + prior)/(denom + prior)). In the
#   units of the contrasted metric (CPM for *_cpm; index units for pausing).
# DEFAULT_MIN_EXPR: independent expression filter. A gene is flagged
#   passes_filter=True when max(numerator, denominator) >= this value, i.e. it is
#   detectable in at least one of the two groups (keeps genuine on/off switches,
#   excludes genes that are ~0 in both). Reported as a flag — rows are NOT dropped.
DEFAULT_PRIOR_COUNT = 1.0
DEFAULT_MIN_EXPR = 1.0

# =============================================================================
# LOGGING UTILITIES
# =============================================================================

def log(section: str, message: str):
    """Consistent logging format"""
    timestamp = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    print(f"{LOG_PREFIX} {section} | {message} | ts={timestamp}", flush=True)

def log_info(message: str):
    """Log informational message"""
    print(f"{LOG_PREFIX} INFO | {message}", flush=True)

def log_error(message: str):
    """Log error message"""
    print(f"{LOG_PREFIX} ERROR | {message}", file=sys.stderr, flush=True)

def log_warning(message: str):
    """Log warning message"""
    print(f"{LOG_PREFIX} WARNING | {message}", flush=True)

# =============================================================================
# MANIFEST LOADING AND VALIDATION
# =============================================================================

def load_samples_manifest(manifest_path: str) -> pd.DataFrame:
    """
    Load and validate samples manifest
    
    Args:
        manifest_path: Path to samples TSV file
        
    Returns:
        DataFrame with validated manifest
        
    Raises:
        SystemExit on validation failure
    """
    log("MANIFEST", f"Loading: {manifest_path}")
    
    try:
        samples = pd.read_csv(manifest_path, sep="\t", dtype=str)
    except Exception as e:
        log_error(f"Failed to read manifest: {e}")
        raise SystemExit(2)
    
    # Validate required columns
    missing_cols = [col for col in REQUIRED_MANIFEST_COLS if col not in samples.columns]
    if missing_cols:
        log_error(f"Missing required columns: {', '.join(missing_cols)}")
        raise SystemExit(2)
    
    log("MANIFEST", f"Loaded {len(samples)} entries")
    
    # Remove duplicates
    original_count = len(samples)
    samples = samples.drop_duplicates(subset=REQUIRED_MANIFEST_COLS, keep="first")
    
    if len(samples) < original_count:
        log_warning(f"Removed {original_count - len(samples)} duplicate entries")
    
    # Validate file paths exist
    samples["file_exists"] = samples["file"].apply(
        lambda p: os.path.exists(str(p)) if pd.notna(p) else False
    )
    
    missing_files = samples[~samples["file_exists"]]
    if not missing_files.empty:
        log_warning(f"Found {len(missing_files)} missing files:")
        for _, row in missing_files.head(5).iterrows():
            log_warning(f"  {row['sample_id']}: {row['file']}")
        if len(missing_files) > 5:
            log_warning(f"  ... and {len(missing_files) - 5} more")
        
        samples = samples[samples["file_exists"]]
    
    samples = samples.drop("file_exists", axis=1)
    
    if samples.empty:
        log_error("No valid file paths in manifest")
        raise SystemExit(3)
    
    log("MANIFEST", f"Validated {len(samples)} samples with existing files")
    return samples

# =============================================================================
# PER-SAMPLE METRICS LOADING
# =============================================================================

def load_sample_metrics(
    sample_row: pd.Series,
    required_cols: List[str]
) -> Optional[pd.DataFrame]:
    """
    Load metrics for a single sample
    
    Args:
        sample_row: Row from samples manifest
        required_cols: Required metric columns
        
    Returns:
        DataFrame with sample metrics or None on failure
    """
    filepath = str(sample_row["file"])
    
    try:
        # Efficient loading: only read needed columns
        df = pd.read_csv(
            filepath,
            sep="\t",
            usecols=lambda c: c in required_cols,
            dtype="object"
        )
        
        # Check for missing columns
        present_cols = set(df.columns)
        missing_cols = [col for col in required_cols if col not in present_cols]
        
        if missing_cols:
            log_warning(f"Sample {sample_row['sample_id']}: missing columns {', '.join(missing_cols)}")
            # Add missing columns as NaN
            for col in missing_cols:
                df[col] = np.nan
        
        # Ensure column order
        df = df[required_cols]
        
        # Add sample metadata
        df["sample_id"] = sample_row["sample_id"]
        df["condition"] = sample_row["condition"]
        df["timepoint"] = sample_row["timepoint"]
        df["replicate"] = sample_row["replicate"]
        
        return df
        
    except Exception as e:
        log_error(f"Failed to load {filepath}: {e}")
        return None

def load_all_samples(samples_df: pd.DataFrame) -> pd.DataFrame:
    """
    Load metrics from all samples
    
    Args:
        samples_df: Samples manifest DataFrame
        
    Returns:
        Merged DataFrame with all samples
    """
    log("LOAD", f"Loading metrics from {len(samples_df)} samples...")
    
    sample_dfs = []
    failed_samples = []
    
    for i, (_, sample_row) in enumerate(samples_df.iterrows(), 1):
        if i % 10 == 0:
            log_info(f"Loaded {i}/{len(samples_df)} samples...")
        
        df = load_sample_metrics(sample_row, REQUIRED_METRIC_COLS)
        
        if df is not None:
            sample_dfs.append(df)
        else:
            failed_samples.append(sample_row["sample_id"])
    
    if failed_samples:
        log_warning(f"Failed to load {len(failed_samples)} samples")
    
    if not sample_dfs:
        log_error("No samples loaded successfully")
        return pd.DataFrame(
            columns=REQUIRED_METRIC_COLS + ["sample_id", "condition", "timepoint", "replicate"]
        )
    
    # Concatenate all samples
    merged = pd.concat(sample_dfs, ignore_index=True)
    
    log("LOAD", f"Loaded {len(merged):,} gene-sample combinations")
    log("LOAD", f"Unique samples: {merged['sample_id'].nunique()}")
    log("LOAD", f"Unique genes: {merged['gene_id'].nunique()}")
    
    return merged

def coerce_numeric_columns(df: pd.DataFrame, numeric_cols: List[str]) -> pd.DataFrame:
    """
    Convert columns to numeric type, handling errors
    
    Args:
        df: DataFrame to process
        numeric_cols: Columns to convert
        
    Returns:
        DataFrame with numeric columns
    """
    log("CONVERT", "Converting metric columns to numeric...")
    
    for col in numeric_cols:
        if col in df.columns and not is_numeric_dtype(df[col]):
            df[col] = pd.to_numeric(df[col], errors="coerce")
            
            nan_count = df[col].isna().sum()
            if nan_count > 0:
                log_warning(f"Column {col}: {nan_count:,} values could not be converted")
    
    return df

# =============================================================================
# CONTRAST PARSING AND COMPUTATION
# =============================================================================

def parse_contrast_spec(contrast_str: str) -> Tuple[str, str, str]:
    """
    Parse contrast specification string
    
    Format: "variable:numerator,denominator"
    Example: "condition:treatment,control"
    
    Args:
        contrast_str: Contrast specification
        
    Returns:
        Tuple of (variable, numerator, denominator)
        
    Raises:
        ValueError on invalid format
    """
    if ":" not in contrast_str:
        raise ValueError("Contrast must contain ':'")
    
    variable, comparison = contrast_str.split(":", 1)
    
    if "," not in comparison:
        raise ValueError("Comparison must contain ','")
    
    numerator, denominator = comparison.split(",", 1)
    
    return (variable.strip(), numerator.strip(), denominator.strip())

def safe_log2_fold_change(numerator: pd.Series, denominator: pd.Series) -> pd.Series:
    """
    Calculate log2 fold change with pseudocount
    
    Args:
        numerator: Numerator values
        denominator: Denominator values
        
    Returns:
        Log2 fold change values
    """
    return safe_log2_fold_change_prior(numerator, denominator, DEFAULT_PRIOR_COUNT)


def safe_log2_fold_change_prior(numerator, denominator, prior_count: float):
    """
    Log2 fold change with an additive prior count (edgeR-style shrinkage).

    Using a prior on the scale of the data (rather than a vanishing 1e-9
    pseudocount) prevents genes with near-zero signal in one group from
    producing extreme, meaningless fold-changes: a gene going 0 -> 0.7 with
    prior=1 yields log2(1.7/1.0)=0.77 instead of ~30. Assay-agnostic.
    """
    num = np.asarray(numerator, dtype=float)
    den = np.asarray(denominator, dtype=float)
    return np.log2((num + prior_count) / (den + prior_count))


def benjamini_hochberg(pvalues: np.ndarray) -> np.ndarray:
    """
    Benjamini-Hochberg FDR correction.
    
    Args:
        pvalues: Array of p-values (may contain NaN)
        
    Returns:
        Adjusted p-values (padj), same shape as input
    """
    p = np.asarray(pvalues, dtype=float)
    padj = np.full_like(p, np.nan)
    valid = np.isfinite(p)
    if not np.any(valid):
        return padj

    # Rank ONLY the finite p-values. NaNs (e.g. genes where a group had no signal,
    # or merged-replicate contrasts where a Mann-Whitney U test is undefined for
    # n=1 vs n=1) are excluded from the test count and left as NaN in the output.
    # Operating on the full array while counting only finite values broadcasts
    # mismatched shapes — the bug this replaces.
    idx_valid = np.flatnonzero(valid)
    pv = p[idx_valid]
    n = pv.size

    order = np.argsort(pv)
    pv_sorted = pv[order]
    ranks = np.arange(1, n + 1, dtype=float)
    padj_sorted = np.minimum(1.0, pv_sorted * n / ranks)

    # Enforce BH monotonicity as a step-up: sweep from the largest p-value down so
    # each adjusted value is no greater than the next-larger rank's value.
    for i in range(n - 2, -1, -1):
        padj_sorted[i] = min(padj_sorted[i], padj_sorted[i + 1])

    padj_valid = np.empty(n, dtype=float)
    padj_valid[order] = padj_sorted
    padj[idx_valid] = padj_valid
    return padj


def _compute_contrast_with_stats(
    merged_df: pd.DataFrame,
    variable: str,
    numerator: str,
    denominator: str,
    metric: str,
    group_by: str,
    level_col: str,
    prior_count: float = DEFAULT_PRIOR_COUNT,
    min_expr: float = DEFAULT_MIN_EXPR,
) -> Optional[pd.DataFrame]:
    """
    Compute contrast for one metric using replicate-level data with Mann-Whitney U test.
    """
    num_df = merged_df[merged_df[variable] == numerator].copy()
    denom_df = merged_df[merged_df[variable] == denominator].copy()

    if num_df.empty or denom_df.empty:
        return None

    # Aggregate per (gene_id, gene_name, level): collect values from all replicates.
    # When level_col is None (direct group-vs-group, e.g. an unpaired
    # treatment-vs-control contrast), there is no shared level to pair on, so we
    # compare the two groups directly per gene.
    paired = level_col is not None
    id_cols = ["gene_id", "gene_name", level_col] if paired else ["gene_id", "gene_name"]
    num_agg = (
        num_df.groupby(id_cols)[metric]
        .apply(lambda x: x.dropna().tolist())
        .reset_index()
        .rename(columns={metric: "num_vals"})
    )
    denom_agg = (
        denom_df.groupby(id_cols)[metric]
        .apply(lambda x: x.dropna().tolist())
        .reset_index()
        .rename(columns={metric: "denom_vals"})
    )

    joined = num_agg.merge(denom_agg, on=id_cols, how="inner")
    if joined.empty:
        return None

    # Compute median, log2FC, and pvalue per gene
    def row_stats(row):
        nv = row["num_vals"]
        dv = row["denom_vals"]
        if not nv or not dv:
            return pd.Series({"numerator": np.nan, "denominator": np.nan, "log2FC": np.nan, "pvalue": np.nan})
        med_n = np.median(nv)
        med_d = np.median(dv)
        # Prior-count shrinkage: stops near-zero groups from exploding log2FC.
        log2fc = float(np.log2((med_n + prior_count) / (med_d + prior_count)))
        try:
            _, pval = stats.mannwhitneyu(nv, dv, alternative="two-sided")
        except Exception:
            pval = np.nan
        return pd.Series({"numerator": med_n, "denominator": med_d, "log2FC": log2fc, "pvalue": pval})

    stats_df = joined.apply(row_stats, axis=1)
    result = joined[id_cols].copy()
    result["numerator"] = stats_df["numerator"]
    result["denominator"] = stats_df["denominator"]
    result["log2FC"] = stats_df["log2FC"]
    result["pvalue"] = stats_df["pvalue"]
    result["group_by"] = group_by
    result["level"] = result[level_col] if paired else "all"
    result["contrast"] = f"{variable}:{numerator}_vs_{denominator}"
    result["metric"] = metric
    # Independent expression filter (reported as a flag — rows are NOT dropped):
    #   expr_mean    : mean of the two group values (DESeq2 baseMean analogue)
    #   passes_filter: detectable in >=1 group (max >= min_expr). Keeps genuine
    #                  on/off switches; excludes genes that are ~0 in both groups,
    #                  which are the source of pseudocount-driven extreme log2FC.
    num_v = result["numerator"].astype(float)
    den_v = result["denominator"].astype(float)
    result["expr_mean"] = 0.5 * (num_v + den_v)
    result["passes_filter"] = (np.maximum(num_v, den_v) >= float(min_expr))
    return result[[
        "gene_id", "gene_name", "group_by", "level",
        "contrast", "metric", "numerator", "denominator", "log2FC", "pvalue",
        "expr_mean", "passes_filter"
    ]]


def compute_single_contrast(
    merged_df: pd.DataFrame,
    variable: str,
    numerator: str,
    denominator: str,
    metrics: List[str],
    prior_count: float = DEFAULT_PRIOR_COUNT,
    min_expr: float = DEFAULT_MIN_EXPR,
) -> Optional[pd.DataFrame]:
    """
    Compute a single contrast with replicate-level statistical testing.

    Uses Mann-Whitney U test for significance. log2FC from group medians.
    padj is applied in compute_all_contrasts (Benjamini-Hochberg per contrast+metric).

    Args:
        merged_df: Replicate-level metrics (sample_id, condition, timepoint, etc.)
        variable: Variable to contrast (condition or timepoint)
        numerator: Numerator group name
        denominator: Denominator group name
        metrics: Metrics to compute contrasts for

    Returns:
        DataFrame with contrast results including pvalue, or None
    """
    if variable not in ["condition", "timepoint", "group"]:
        log_warning(f"Unsupported contrast variable: {variable}")
        return None

    results = []
    for metric in metrics:
        if metric not in merged_df.columns:
            continue
        if variable == "condition":
            df = _compute_contrast_with_stats(
                merged_df, variable, numerator, denominator, metric,
                group_by="timepoint", level_col="timepoint",
                prior_count=prior_count, min_expr=min_expr,
            )
        elif variable == "timepoint":
            df = _compute_contrast_with_stats(
                merged_df, variable, numerator, denominator, metric,
                group_by="condition", level_col="condition",
                prior_count=prior_count, min_expr=min_expr,
            )
        else:  # group — direct, unpaired group-vs-group (e.g. treatment vs control)
            df = _compute_contrast_with_stats(
                merged_df, variable, numerator, denominator, metric,
                group_by="group", level_col=None,
                prior_count=prior_count, min_expr=min_expr,
            )
        if df is not None and not df.empty:
            results.append(df)

    if not results:
        return None
    return pd.concat(results, ignore_index=True)

def compute_all_contrasts(
    merged_df: pd.DataFrame,
    contrast_specs: List[Tuple[str, str, str]],
    metrics: List[str],
    prior_count: float = DEFAULT_PRIOR_COUNT,
    min_expr: float = DEFAULT_MIN_EXPR,
) -> Optional[pd.DataFrame]:
    """
    Compute all specified contrasts with replicate-level statistical testing.

    Uses Mann-Whitney U for pvalue, Benjamini-Hochberg for padj (per contrast+metric).

    Args:
        merged_df: Replicate-level merged metrics
        contrast_specs: List of (variable, numerator, denominator) tuples
        metrics: Metrics to contrast (e.g. pi_len_norm, pi_raw, body_cpm, tss_cpm)

    Returns:
        DataFrame with log2FC, pvalue, padj or None
    """
    if merged_df.empty:
        log_warning("Empty DataFrame, cannot compute contrasts")
        return None

    log("CONTRAST", f"Computing {len(contrast_specs)} contrasts...")
    log("CONTRAST", f"Metrics: {', '.join(metrics)}")

    working_df = merged_df.copy()
    working_df = coerce_numeric_columns(working_df, metrics)

    # Synthetic composite key for direct group-vs-group contrasts (condition+timepoint),
    # e.g. "celastrol_40" vs "no_treatment_0". Mirrors how the pipeline names groups.
    if "condition" in working_df.columns and "timepoint" in working_df.columns:
        working_df["group"] = (
            working_df["condition"].astype(str) + "_" + working_df["timepoint"].astype(str)
        )

    all_results = []
    for i, (variable, numerator, denominator) in enumerate(contrast_specs, 1):
        log_info(f"Computing contrast {i}/{len(contrast_specs)}: {variable}:{numerator} vs {denominator}")

        result = compute_single_contrast(
            working_df,
            variable,
            numerator,
            denominator,
            metrics,
            prior_count=prior_count,
            min_expr=min_expr,
        )

        if result is not None and not result.empty:
            all_results.append(result)

    if not all_results:
        log_warning("No contrasts produced")
        return None

    contrasts_df = pd.concat(all_results, ignore_index=True)

    # Benjamini-Hochberg correction per (contrast, metric)
    log("CONTRAST", "Applying Benjamini-Hochberg FDR correction...")
    padj_arr = np.full(len(contrasts_df), np.nan, dtype=float)
    for (_, _), grp in contrasts_df.groupby(["contrast", "metric"]):
        idx = grp.index
        padj_arr[idx] = benjamini_hochberg(grp["pvalue"].values)
    contrasts_df["padj"] = padj_arr

    # Stable column order. Keep log2FC/pvalue/padj at positions 9/10/11 (downstream
    # log-parsing in module 12 reads those by index); append the new filter columns.
    col_order = [
        "gene_id", "gene_name", "group_by", "level", "contrast", "metric",
        "numerator", "denominator", "log2FC", "pvalue", "padj",
        "expr_mean", "passes_filter",
    ]
    contrasts_df = contrasts_df[[c for c in col_order if c in contrasts_df.columns]]

    n_pass = int(contrasts_df["passes_filter"].sum()) if "passes_filter" in contrasts_df else 0
    log("CONTRAST", f"Generated {len(contrasts_df):,} contrast results (log2FC, pvalue, padj)")
    log("CONTRAST", f"Passing expression filter (max>={min_expr}): {n_pass:,} of {len(contrasts_df):,}")
    return contrasts_df

# =============================================================================
# VISUALIZATION
# =============================================================================

def create_ma_plot(
    contrast_df: pd.DataFrame,
    title: str,
    output_path: str
):
    """
    Create MA plot (mean vs log2 fold change)
    
    Args:
        contrast_df: DataFrame with contrast results
        title: Plot title
        output_path: Output PNG path
    """
    if contrast_df is None or contrast_df.empty:
        return

    # Restrict to genes passing the expression filter so the plot is not dominated
    # by undetectable genes whose fold-changes are pseudocount artifacts. Falls
    # back to all rows if the column is absent (older inputs).
    if "passes_filter" in contrast_df.columns:
        contrast_df = contrast_df[contrast_df["passes_filter"] == True]  # noqa: E712
    if contrast_df.empty:
        return

    # Mean expression (log2). Use the precomputed expr_mean when available.
    if "expr_mean" in contrast_df.columns:
        mean_expr = np.log2(contrast_df["expr_mean"].astype(float) + 1.0)
    else:
        mean_expr = 0.5 * np.log2(
            (contrast_df["numerator"].astype(float) + 1.0) *
            (contrast_df["denominator"].astype(float) + 1.0)
        )

    log2fc = contrast_df["log2FC"]

    # Create plot
    plt.figure(figsize=(6, 5), dpi=130)
    plt.scatter(mean_expr, log2fc, s=6, alpha=0.6, edgecolors="none")
    plt.axhline(0, color="red", linewidth=1, linestyle="--", alpha=0.7)
    plt.xlabel("Mean Expression (log2)", fontsize=10)
    plt.ylabel("Log2 Fold Change", fontsize=10)
    plt.title(title, fontsize=11, weight="bold")
    plt.grid(True, alpha=0.3, linestyle=":", linewidth=0.5)
    plt.tight_layout()
    plt.savefig(output_path, dpi=130, bbox_inches="tight")
    plt.close()
    
    log_info(f"Created MA plot: {output_path}")

def create_heatmap(
    merged_df: pd.DataFrame,
    value_col: str,
    group_col: str,
    output_path: str,
    top_n: int = 100
):
    """
    Create heatmap of top variable genes
    
    Args:
        merged_df: Merged metrics DataFrame
        value_col: Column to use for values
        group_col: Column to group by
        output_path: Output PNG path
        top_n: Number of top variable genes
    """
    if merged_df is None or merged_df.empty:
        return
    
    log_info(f"Creating heatmap: {value_col} by {group_col}...")
    
    # Pivot table: genes × groups
    pivot = merged_df.pivot_table(
        index=["gene_id", "gene_name"],
        columns=group_col,
        values=value_col,
        aggfunc="median"
    )
    
    # Remove infinite values and NAs
    pivot = pivot.replace([np.inf, -np.inf], np.nan).dropna()

    # Sort columns numerically when the group labels are numeric-like (e.g.
    # timepoint "0","10","160","20","40","60"), which otherwise sort as text
    # and scramble the chronological order in the plot (0,10,160,20,40,60).
    # Falls back to a plain string sort for genuinely categorical labels
    # (e.g. condition names) that aren't all numeric.
    try:
        sorted_cols = sorted(pivot.columns, key=lambda x: float(x))
    except (TypeError, ValueError):
        sorted_cols = sorted(pivot.columns, key=lambda x: str(x))
    pivot = pivot[sorted_cols]

    if pivot.empty:
        log_warning(f"Empty pivot table for {value_col} by {group_col}")
        return
    
    if len(pivot.columns) == 0:
        log_warning(f"No columns in pivot for {value_col} by {group_col}")
        return
    
    # Select top variable genes
    variance = pivot.var(axis=1, numeric_only=True).sort_values(ascending=False)
    top_genes = pivot.loc[variance.index[:top_n]]
    
    if top_genes.empty:
        log_warning(f"No genes remaining after variance filtering")
        return
    
    # Calculate z-scores
    z_scores = (
        top_genes.subtract(top_genes.mean(axis=1), axis=0)
        .divide(top_genes.std(axis=1, ddof=0) + 1e-9, axis=0)
    )
    
    # Create figure
    fig_width = max(6, z_scores.shape[1] * 0.25)
    fig_height = max(6, z_scores.shape[0] * 0.06)
    
    plt.figure(figsize=(fig_width, fig_height), dpi=130)
    
    # Plot heatmap
    im = plt.imshow(
        z_scores.values,
        aspect="auto",
        interpolation="nearest",
        cmap="RdBu_r",
        vmin=-3,
        vmax=3
    )
    
    # Configure axes
    plt.xticks(
        range(z_scores.shape[1]),
        z_scores.columns,
        rotation=90,
        fontsize=8
    )
    plt.yticks([])
    plt.ylabel(f"{len(z_scores)} genes", fontsize=9)
    
    # Add colorbar
    cbar = plt.colorbar(im, fraction=0.046, pad=0.04)
    cbar.set_label("Z-score", fontsize=9)
    
    # Title
    actual_n = min(top_n, z_scores.shape[0])
    plt.title(
        f"Top {actual_n} Variable Genes — {value_col} (by {group_col})",
        fontsize=11,
        weight="bold"
    )
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=130, bbox_inches="tight")
    plt.close()
    
    log_info(f"Created heatmap: {output_path}")

def generate_all_plots(
    merged_df: pd.DataFrame,
    contrasts_df: Optional[pd.DataFrame],
    plots_dir: str,
    top_n: int
):
    """
    Generate all plots
    
    Args:
        merged_df: Merged metrics DataFrame
        contrasts_df: Contrasts DataFrame (optional)
        plots_dir: Output directory for plots
        top_n: Number of top genes for heatmaps
    """
    log("PLOTS", f"Generating plots in: {plots_dir}")
    
    os.makedirs(plots_dir, exist_ok=True)
    
    # Heatmaps for merged data
    if not merged_df.empty:
        metrics = ["pi_len_norm", "pi_raw", "body_cpm", "tss_cpm"]
        metrics = [m for m in metrics if m in merged_df.columns]
        groups = ["condition", "timepoint"]
        
        for metric in metrics:
            for group in groups:
                try:
                    output_path = os.path.join(
                        plots_dir,
                        f"heatmap_{metric}_by_{group}.png"
                    )
                    
                    create_heatmap(
                        merged_df,
                        metric,
                        group,
                        output_path,
                        top_n=max(10, min(top_n, 500))
                    )
                    
                except Exception as e:
                    log_error(f"Heatmap failed ({metric} by {group}): {e}")
    
    # MA plots for contrasts
    if contrasts_df is not None and not contrasts_df.empty:
        for (contrast, metric, level), subset in contrasts_df.groupby(
            ["contrast", "metric", "level"]
        ):
            try:
                safe_level = str(level).replace("/", "_").replace(" ", "_")
                safe_contrast = str(contrast).replace("/", "_").replace(" ", "_")
                
                output_path = os.path.join(
                    plots_dir,
                    f"MA_{metric}_{safe_contrast}_at_{safe_level}.png"
                )
                
                title = f"{contrast} @ {level} — {metric}"
                
                create_ma_plot(subset, title, output_path)
                
            except Exception as e:
                log_error(f"MA plot failed ({contrast}, {metric}, {level}): {e}")
    
    log("PLOTS", "Plot generation complete")

# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main execution function"""
    
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Aggregate and compare Pol-II metrics across samples",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "--samples-tsv",
        required=True,
        help="Samples manifest (TSV): sample_id, condition, timepoint, replicate, file"
    )
    parser.add_argument(
        "--out-merged",
        required=True,
        help="Output merged TSV path"
    )
    parser.add_argument(
        "--contrasts",
        nargs="*",
        default=[],
        help="Contrast specifications (e.g., 'condition:treatment,control')"
    )
    parser.add_argument(
        "--out-contrasts",
        default=None,
        help="Output contrasts TSV path (optional)"
    )
    parser.add_argument(
        "--plots-dir",
        default=None,
        help="Output directory for plots (optional)"
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=100,
        help="Number of top variable genes for heatmaps [default: 100]"
    )
    parser.add_argument(
        "--metrics",
        default=",".join(AVAILABLE_METRICS),
        help=f"Comma-separated metrics for contrasts [default: {','.join(AVAILABLE_METRICS)}]"
    )
    parser.add_argument(
        "--prior-count",
        type=float,
        default=DEFAULT_PRIOR_COUNT,
        help=(
            "Additive prior count for log2 fold-change shrinkage (edgeR-style); "
            f"prevents extreme FC from near-zero groups [default: {DEFAULT_PRIOR_COUNT}]"
        ),
    )
    parser.add_argument(
        "--min-expr",
        type=float,
        default=DEFAULT_MIN_EXPR,
        help=(
            "Independent expression filter: gene flagged passes_filter when "
            f"max(numerator, denominator) >= this value [default: {DEFAULT_MIN_EXPR}]"
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}"
    )
    
    args = parser.parse_args()
    
    # Start
    log("START", f"compare_pol_metrics.py v{VERSION}")
    log("START", f"Samples manifest: {args.samples_tsv}")
    log("START", f"Output merged: {args.out_merged}")
    
    # Load samples manifest
    log("═" * 70, "")
    samples_df = load_samples_manifest(args.samples_tsv)
    
    # Load per-sample metrics
    log("═" * 70, "")
    merged_df = load_all_samples(samples_df)
    
    # Convert to numeric
    merged_df = coerce_numeric_columns(merged_df, AVAILABLE_METRICS)
    
    # Save merged table
    log("═" * 70, "")
    log("OUTPUT", f"Writing merged table: {args.out_merged}")
    merged_df.to_csv(args.out_merged, sep="\t", index=False)
    
    merged_size = os.path.getsize(args.out_merged)
    log("OUTPUT", f"Written: {merged_size:,} bytes")
    log("OUTPUT", f"Rows: {len(merged_df):,}")
    log("OUTPUT", f"Samples: {merged_df['sample_id'].nunique()}")
    log("OUTPUT", f"Genes: {merged_df['gene_id'].nunique()}")
    
    # Compute contrasts if specified
    contrasts_df = None
    if args.contrasts and args.out_contrasts:
        log("═" * 70, "")
        
        # Parse contrast specifications
        contrast_specs = []
        for contrast_str in args.contrasts:
            try:
                spec = parse_contrast_spec(contrast_str)
                contrast_specs.append(spec)
                log_info(f"Parsed contrast: {spec[0]}:{spec[1]} vs {spec[2]}")
            except Exception as e:
                log_warning(f"Could not parse contrast '{contrast_str}': {e}")
        
        if contrast_specs:
            # Parse metrics (filter to available only)
            contrast_metrics = [
                m.strip() for m in args.metrics.split(",")
                if m.strip() in AVAILABLE_METRICS
            ]
            if not contrast_metrics:
                contrast_metrics = AVAILABLE_METRICS
            # Compute contrasts
            contrasts_df = compute_all_contrasts(
                merged_df,
                contrast_specs,
                contrast_metrics,
                prior_count=args.prior_count,
                min_expr=args.min_expr,
            )
            
            if contrasts_df is not None and not contrasts_df.empty:
                # Save contrasts table
                log("OUTPUT", f"Writing contrasts table: {args.out_contrasts}")
                contrasts_df.to_csv(args.out_contrasts, sep="\t", index=False)
                
                contrasts_size = os.path.getsize(args.out_contrasts)
                log("OUTPUT", f"Written: {contrasts_size:,} bytes")
                log("OUTPUT", f"Rows: {len(contrasts_df):,}")
            else:
                log_warning("No contrasts generated (empty overlap or invalid specs)")
    
    # Generate plots if requested
    if args.plots_dir:
        log("═" * 70, "")
        generate_all_plots(
            merged_df,
            contrasts_df,
            args.plots_dir,
            args.top_n
        )
    
    # Success
    log("═" * 70, "")
    log("COMPLETE", "Processing complete")
    log("COMPLETE", f"Samples: {merged_df['sample_id'].nunique()}")
    log("COMPLETE", f"Genes: {merged_df['gene_id'].nunique()}")
    if contrasts_df is not None:
        log("COMPLETE", f"Contrasts: {len(contrasts_df):,} results")
    
    return 0

# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        log_error("Interrupted by user")
        sys.exit(130)
    except Exception as e:
        log_error(f"Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)