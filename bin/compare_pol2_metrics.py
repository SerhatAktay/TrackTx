#!/usr/bin/env python3
# =============================================================================
# compare_pol2_metrics.py — Pol-II Metrics Aggregation and Comparison
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

# Matplotlib backend must be set before importing pyplot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =============================================================================
# CONSTANTS
# =============================================================================

VERSION = "2.0.0"
LOG_PREFIX = "[COMPARE]"

# Required columns in samples manifest
REQUIRED_MANIFEST_COLS = ["sample_id", "condition", "timepoint", "replicate", "file"]

# Required metric columns in per-sample files
REQUIRED_METRIC_COLS = ["gene_id", "gene_name", "pi_len_norm", "pi_raw", "body_cpm", "tss_cpm"]

# Metrics available for analysis
AVAILABLE_METRICS = ["pi_len_norm", "pi_raw", "body_cpm", "tss_cpm"]

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
    return np.log2((numerator.astype(float) + 1e-9) / (denominator.astype(float) + 1e-9))

def compute_single_contrast(
    medians_df: pd.DataFrame,
    variable: str,
    numerator: str,
    denominator: str,
    metrics: List[str]
) -> Optional[pd.DataFrame]:
    """
    Compute a single contrast
    
    Args:
        medians_df: DataFrame with median values per group
        variable: Variable to contrast (condition or timepoint)
        numerator: Numerator group name
        denominator: Denominator group name
        metrics: Metrics to compute contrasts for
        
    Returns:
        DataFrame with contrast results or None
    """
    if variable not in ["condition", "timepoint"]:
        log_warning(f"Unsupported contrast variable: {variable}")
        return None
    
    results = []
    
    for metric in metrics:
        if metric not in medians_df.columns:
            continue
        
        if variable == "condition":
            # Compare conditions across timepoints
            numerator_df = medians_df.query("condition == @numerator")[
                ["gene_id", "gene_name", "timepoint", metric]
            ].rename(columns={metric: "numerator"})
            
            denominator_df = medians_df.query("condition == @denominator")[
                ["gene_id", "gene_name", "timepoint", metric]
            ].rename(columns={metric: "denominator"})
            
            joined = numerator_df.merge(
                denominator_df,
                on=["gene_id", "gene_name", "timepoint"],
                how="inner"
            )
            
            if joined.empty:
                log_warning(f"No overlap for {variable}:{numerator} vs {denominator} on {metric}")
                continue
            
            joined["group_by"] = "timepoint"
            joined["level"] = joined["timepoint"]
            
        else:  # timepoint
            # Compare timepoints across conditions
            numerator_df = medians_df.query("timepoint == @numerator")[
                ["gene_id", "gene_name", "condition", metric]
            ].rename(columns={metric: "numerator"})
            
            denominator_df = medians_df.query("timepoint == @denominator")[
                ["gene_id", "gene_name", "condition", metric]
            ].rename(columns={metric: "denominator"})
            
            joined = numerator_df.merge(
                denominator_df,
                on=["gene_id", "gene_name", "condition"],
                how="inner"
            )
            
            if joined.empty:
                log_warning(f"No overlap for {variable}:{numerator} vs {denominator} on {metric}")
                continue
            
            joined["group_by"] = "condition"
            joined["level"] = joined["condition"]
        
        # Remove rows with missing values
        joined = joined.dropna(subset=["numerator", "denominator"])
        
        if joined.empty:
            log_warning(f"All NA after join for {variable}:{numerator} vs {denominator} on {metric}")
            continue
        
        # Calculate fold change
        joined["contrast"] = f"{variable}:{numerator}_vs_{denominator}"
        joined["metric"] = metric
        joined["log2FC"] = safe_log2_fold_change(joined["numerator"], joined["denominator"])
        
        # Select output columns
        result = joined[[
            "gene_id", "gene_name", "group_by", "level",
            "contrast", "metric", "numerator", "denominator", "log2FC"
        ]]
        
        results.append(result)
    
    if not results:
        return None
    
    return pd.concat(results, ignore_index=True)

def compute_all_contrasts(
    merged_df: pd.DataFrame,
    contrast_specs: List[Tuple[str, str, str]],
    metrics: List[str]
) -> Optional[pd.DataFrame]:
    """
    Compute all specified contrasts
    
    Args:
        merged_df: Merged metrics DataFrame
        contrast_specs: List of contrast tuples
        metrics: Metrics to analyze
        
    Returns:
        DataFrame with all contrasts or None
    """
    if merged_df.empty:
        log_warning("Empty DataFrame, cannot compute contrasts")
        return None
    
    log("CONTRAST", f"Computing {len(contrast_specs)} contrasts...")
    
    # Ensure numeric columns
    working_df = merged_df.copy()
    working_df = coerce_numeric_columns(working_df, metrics)
    
    # Collapse replicates: median per (gene, condition, timepoint)
    log("CONTRAST", "Collapsing replicates (median)...")
    medians = working_df.groupby(
        ["gene_id", "gene_name", "condition", "timepoint"],
        dropna=False
    ).median(numeric_only=True).reset_index()
    
    log("CONTRAST", f"Median aggregation: {len(medians):,} unique combinations")
    
    # Compute each contrast
    all_results = []
    
    for i, (variable, numerator, denominator) in enumerate(contrast_specs, 1):
        log_info(f"Computing contrast {i}/{len(contrast_specs)}: {variable}:{numerator} vs {denominator}")
        
        result = compute_single_contrast(
            medians,
            variable,
            numerator,
            denominator,
            metrics
        )
        
        if result is not None:
            all_results.append(result)
    
    if not all_results:
        log_warning("No contrasts produced")
        return None
    
    contrasts_df = pd.concat(all_results, ignore_index=True)
    
    log("CONTRAST", f"Generated {len(contrasts_df):,} contrast results")
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
    
    # Calculate mean expression (log2)
    mean_expr = 0.5 * np.log2(
        (contrast_df["numerator"].astype(float) + 1e-9) *
        (contrast_df["denominator"].astype(float) + 1e-9)
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
        metrics = ["pi_len_norm", "body_cpm"]
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
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}"
    )
    
    args = parser.parse_args()
    
    # Start
    log("START", f"compare_pol2_metrics.py v{VERSION}")
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
            # Compute contrasts
            contrasts_df = compute_all_contrasts(
                merged_df,
                contrast_specs,
                ["pi_len_norm", "body_cpm"]
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