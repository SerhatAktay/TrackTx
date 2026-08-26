#!/usr/bin/env python3
# =============================================================================
# combine_reports.py — TrackTx Cohort Report Aggregation
# =============================================================================
#
# Purpose:
#   Aggregates per-sample JSON reports into unified cohort-level summaries
#   with interactive single-page application (SPA) HTML dashboard.
#
# Features:
#   • Robust JSON intake (*.summary.json, *.report.json, *.json.gz)
#   • Cohort-wide quality control assessment
#   • Interactive visualizations (no external dependencies)
#   • Functional region aggregation
#   • Export capabilities (TSV, JSON, CSV)
#   • Offline-capable HTML dashboard
#
# Inputs:
#   Per-sample JSON reports (files or directories)
#   Supported formats:
#     - *.summary.json
#     - *.report.json  
#     - *.json
#     - *.json.gz (gzipped)
#
# Outputs:
#   • global_summary.html: Interactive SPA dashboard
#   • global_summary.tsv: Cohort metrics table
#   • global_summary.json: Structured cohort data
#   • global_region_totals.tsv: Aggregated region counts (optional)
#
# Design:
#   • Single-file HTML (embedded CSS/JS, no CDNs)
#   • Works offline
#   • Responsive design
#   • Export functionality
#
# =============================================================================

from __future__ import annotations
import argparse
import datetime
import gzip
import io
import json
import math
import os
import re
import shlex
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

# =============================================================================
# CONSTANTS
# =============================================================================

VERSION = "2.0.0"
LOG_PREFIX = "[COMBINE]"

# Region name pattern for unlocalized reads
UNLOCALIZED_PATTERN = re.compile(
    r"non[-\s_]?localized|unlocalized|unlocalised",
    re.IGNORECASE
)

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
# FILE DISCOVERY AND LOADING
# =============================================================================

def discover_json_files(input_paths: List[str]) -> List[str]:
    """
    Discover JSON files from input paths (files or directories)
    
    Args:
        input_paths: List of file paths or directory paths
        
    Returns:
        Sorted list of discovered JSON file paths
    """
    log("DISCOVER", f"Searching {len(input_paths)} input path(s)...")
    
    discovered = []
    
    for path_str in input_paths:
        if not path_str:
            continue
        
        path = Path(path_str)
        
        # Handle directories
        if path.is_dir():
            log_info(f"Scanning directory: {path}")
            
            for root, _, files in os.walk(path):
                for filename in files:
                    if is_json_file(filename):
                        full_path = os.path.join(root, filename)
                        discovered.append(full_path)
                        
        # Handle individual files
        elif path.is_file():
            if is_json_file(path.name):
                discovered.append(str(path))
            else:
                log_warning(f"Skipping non-JSON file: {path}")
        else:
            log_warning(f"Path not found: {path}")
    
    # Remove duplicates and sort
    unique_files = sorted(set(discovered))
    
    log("DISCOVER", f"Found {len(unique_files)} JSON files")
    return unique_files

def is_json_file(filename: str) -> bool:
    """
    Check if filename matches JSON patterns
    
    Args:
        filename: File name to check
        
    Returns:
        True if matches JSON patterns
    """
    lowercase = filename.lower()
    return lowercase.endswith((
        ".summary.json",
        ".report.json",
        ".json",
        ".json.gz"
    ))

def read_json_file(filepath: str) -> Optional[Dict[str, Any]]:
    """
    Read JSON file (handles gzip compression)
    
    Args:
        filepath: Path to JSON file
        
    Returns:
        Parsed JSON as dictionary, or None on failure
    """
    if not os.path.exists(filepath):
        log_warning(f"File not found: {filepath}")
        return None
    
    try:
        # Handle gzipped files
        if filepath.lower().endswith(".gz"):
            with gzip.open(filepath, "rb") as f:
                content = f.read().decode("utf-8", errors="replace")
                return json.loads(content)
        
        # Handle regular files
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            return json.load(f)
            
    except json.JSONDecodeError as e:
        log_error(f"JSON parse error in {filepath}: {e}")
        return None
    except Exception as e:
        log_error(f"Failed to read {filepath}: {e}")
        return None

# =============================================================================
# DATA VALIDATION AND SANITIZATION
# =============================================================================

def is_valid_number(value: Any) -> bool:
    """
    Check if value is a valid finite number
    
    Args:
        value: Value to check
        
    Returns:
        True if valid finite number
    """
    if not isinstance(value, (int, float)):
        return False
    
    if isinstance(value, float):
        return not (math.isnan(value) or math.isinf(value))
    
    return True

def to_float(value: Any) -> Optional[float]:
    """
    Convert value to float, handling NaN/Inf
    
    Args:
        value: Value to convert
        
    Returns:
        Float value or None if invalid
    """
    try:
        num = float(value)
        if math.isnan(num) or math.isinf(num):
            return None
        return num
    except (ValueError, TypeError):
        return None

def sanitize_data(obj: Any) -> Any:
    """
    Recursively sanitize data structure (remove NaN/Inf)
    
    Args:
        obj: Object to sanitize
        
    Returns:
        Sanitized object
    """
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
    
    if isinstance(obj, dict):
        return {k: sanitize_data(v) for k, v in obj.items()}
    
    if isinstance(obj, list):
        return [sanitize_data(v) for v in obj]
    
    return obj

# =============================================================================
# SAMPLE DATA NORMALIZATION
# =============================================================================

def _resolve_input_reads(
    qc: Optional[Dict[str, Any]],
    json_data: Dict[str, Any],
    reads_total_functional: Optional[float],
) -> Dict[str, Any]:
    """
    Resolve input_reads with fallback chain; track source for UI labeling.

    Returns dict with input_reads and input_reads_source.
    Source: "qc" | "fallback_dedup" | "fallback_functional"
    """
    qc_val = (qc.get("total_reads_raw") if qc else None) or json_data.get("input_reads")
    dedup_val = (qc.get("dedup_reads_mapq_ge") if qc else None) or json_data.get("dedup_reads")
    func_val = int(reads_total_functional) if reads_total_functional else None

    if qc_val is not None:
        return {"input_reads": qc_val, "input_reads_source": "qc"}
    if dedup_val is not None:
        return {"input_reads": dedup_val, "input_reads_source": "fallback_dedup"}
    if func_val is not None:
        return {"input_reads": func_val, "input_reads_source": "fallback_functional"}
    return {"input_reads": None, "input_reads_source": None}


def normalize_sample_data(
    json_data: Dict[str, Any],
    fallback_name: str
) -> Optional[Dict[str, Any]]:
    """
    Normalize per-sample JSON into consistent format
    
    Handles multiple schema versions and missing fields
    
    Args:
        json_data: Raw JSON data from sample report
        fallback_name: Fallback sample name if not in JSON
        
    Returns:
        Normalized sample dictionary or None on failure
    """
    # Extract sample identification
    sample_id = str(
        json_data.get("sample") or 
        json_data.get("sample_id") or 
        fallback_name or 
        "NA"
    )
    
    condition = json_data.get("condition", "NA")
    timepoint = json_data.get("timepoint", "NA")
    replicate = json_data.get("replicate", "NA")
    
    # Extract regions data
    regions = json_data.get("regions") or []
    
    # Parse functional regions
    func_totals = {}
    region_counts = {}
    region_len_totals = {}
    region_len_medians = {}
    
    # Every region entry that fails to parse is tracked here (with why) so a
    # malformed report doesn't just quietly under-count a sample's regions
    # with no trace -- mirrors the aggregated-warning pattern already used in
    # calculate_pol_metrics.py's count_reads_pysam().
    skipped_regions: list = []
    for region_item in regions:
        try:
            region_name = str(region_item.get("region", ""))
            reads = float(region_item.get("reads", 0) or 0.0)
            count = int(region_item.get("region_count", 0) or 0)
            len_total = float(region_item.get("region_length_total_bp", 0) or 0.0)
            len_median = float(region_item.get("region_length_median_bp", 0) or 0.0)
            
            if region_name:
                func_totals[region_name] = func_totals.get(region_name, 0.0) + reads
                region_counts[region_name] = region_counts.get(region_name, 0) + count
                region_len_totals[region_name] = region_len_totals.get(region_name, 0.0) + len_total
                
                # Handle multiple entries: average medians
                if region_name in region_len_medians and region_len_medians[region_name] != 0:
                    region_len_medians[region_name] = (region_len_medians[region_name] + len_median) / 2.0
                else:
                    region_len_medians[region_name] = len_median
        except Exception as e:
            skipped_regions.append((str(region_item.get("region", "?")), str(e)))
            continue

    if skipped_regions:
        examples = ", ".join(f"{name!r} ({reason})" for name, reason in skipped_regions[:5])
        log_warning(
            f"{sample_id}: skipped {len(skipped_regions)} malformed region "
            f"entr{'y' if len(skipped_regions) == 1 else 'ies'} while combining reports "
            f"(first {min(5, len(skipped_regions))}: {examples})"
        )
    
    # Separate localized from unlocalized reads
    unlocalized_reads = sum(
        reads for name, reads in func_totals.items()
        if UNLOCALIZED_PATTERN.search(name or "")
    )
    
    localized_reads = sum(
        reads for name, reads in func_totals.items()
        if not UNLOCALIZED_PATTERN.search(name or "")
    )
    
    reads_total_functional = localized_reads
    
    # Calculate unlocalized fraction
    unloc_frac = json_data.get("unlocalized_fraction")
    if unloc_frac is None:
        total_with_unloc = localized_reads + unlocalized_reads
        if total_with_unloc > 0:
            unloc_frac = unlocalized_reads / total_with_unloc
    
    unloc_frac = to_float(unloc_frac)
    
    # Extract nested metrics and QC
    metrics = json_data.get("metrics") or {}
    qc = json_data.get("qc") or {}
    
    # Determine duplicate percentage (prefer UMI if available)
    # Use explicit None checks to handle 0.0 values correctly
    dup_percent = None
    
    # Priority 1: UMI deduplication if enabled
    if qc.get("umi_deduplication_enabled", False):
        dup_percent = qc.get("umi_deduplication_percent")
    
    # Priority 2: Try QC dict fields
    if dup_percent is None:
        for field in ["duplicate_percent", "duplicate_perc_of_total", "dup_percent"]:
            if field in qc and qc[field] is not None:
                dup_percent = qc[field]
                break
    
    # Priority 3: Try top-level JSON fields
    if dup_percent is None:
        for field in ["duplicate_percent", "duplicate_perc_of_total"]:
            if field in json_data and json_data[field] is not None:
                dup_percent = json_data[field]
                break
    
    # Build normalized record
    record = {
        "sample_id": sample_id,
        "condition": condition,
        "timepoint": timepoint,
        "replicate": replicate,
        "divergent_regions": metrics.get("divergent_regions"),
        "total_regions": metrics.get("total_functional_regions") or metrics.get("total_regions"),
        "reads_total_functional": reads_total_functional,
        "median_pausing_index": (
            metrics.get("median_pausing_index") if metrics.get("median_pausing_index") is not None 
            else metrics.get("median_pi")
        ),
        "median_density": (
            metrics.get("median_functional_cpm") if metrics.get("median_functional_cpm") is not None 
            else metrics.get("median_density")
        ),
        "cpm_factor": metrics.get("cpm_factor"),
        "crpmsi_factor": (
            metrics.get("sicpm_factor") if metrics.get("sicpm_factor") is not None 
            else metrics.get("crpmsi_factor")
        ),
        "density_source": metrics.get("density_source"),
        "density_reason": metrics.get("density_reason"),
        # Read depth: prefer QC total_reads_raw; fallback to dedup_reads or reads_total_functional
        **(_resolve_input_reads(qc, json_data, reads_total_functional)),
        # Uniquely-mapped reads (NH==1 / MAPQ>=threshold). 'unique_reads_nh1' is
        # the correctly named field; 'dedup_reads' is kept as a back-compat alias
        # (it is NOT PCR-deduplicated unless umi_deduplication_enabled is true).
        "unique_reads_nh1": (
            (qc.get("unique_reads_nh1") if qc else None)
            or (qc.get("dedup_reads_mapq_ge") if qc else None)
            or json_data.get("dedup_reads")
        ),
        "dedup_reads": (
            (qc.get("dedup_reads_mapq_ge") if qc else None) or json_data.get("dedup_reads")
        ),
        "duplicate_percent": dup_percent,
        "umi_deduplication_enabled": qc.get("umi_deduplication_enabled", False),
        "umi_deduplication_percent": qc.get("umi_deduplication_percent"),
        "multimapper_percent": qc.get("multimapper_percent"),
        "uniqueness_method": qc.get("uniqueness_method"),
        "multimap_k": qc.get("multimap_k"),
        "unlocalized_fraction": unloc_frac,
        "func_totals": func_totals,
        "region_counts": region_counts,
        "region_length_totals": region_len_totals,
        "region_length_medians": region_len_medians,
        "regions": regions
    }
    
    return record

# =============================================================================
# REGION AGGREGATION
# =============================================================================

def collect_region_keys(samples: List[Dict[str, Any]]) -> List[str]:
    """
    Collect all unique region names from samples
    
    Args:
        samples: List of normalized sample dicts
        
    Returns:
        Sorted list of region names
    """
    region_set = set()
    
    for sample in samples:
        func_totals = sample.get("func_totals") or {}
        region_set.update(func_totals.keys())
    
    return sorted(region_set)

def aggregate_region_totals(samples: List[Dict[str, Any]]) -> Dict[str, float]:
    """
    Sum region read counts across all samples
    
    Args:
        samples: List of normalized sample dicts
        
    Returns:
        Dictionary mapping region name to total reads
    """
    totals = {}
    
    for sample in samples:
        func_totals = sample.get("func_totals") or {}
        for region_name, reads in func_totals.items():
            try:
                totals[region_name] = totals.get(region_name, 0.0) + float(reads or 0.0)
            except (ValueError, TypeError):
                continue
    
    return totals


def compute_aggregate_summary(rows: List[Dict[str, Any]]) -> Dict[str, Optional[float]]:
    """
    Compute simple cohort-level summary statistics from per-sample rows.
    """
    def collect_numeric(key: str) -> List[float]:
        vals: List[float] = []
        for r in rows:
            v = r.get(key)
            try:
                num = float(v)
            except (TypeError, ValueError):
                continue
            if math.isnan(num) or math.isinf(num):
                continue
            vals.append(num)
        return vals

    def agg_min(vals: List[float]) -> Optional[float]:
        return float(min(vals)) if vals else None

    def agg_max(vals: List[float]) -> Optional[float]:
        return float(max(vals)) if vals else None

    def agg_sum(vals: List[float]) -> Optional[float]:
        return float(sum(vals)) if vals else None

    def agg_mean(vals: List[float]) -> Optional[float]:
        return float(sum(vals) / len(vals)) if vals else None

    def agg_median(vals: List[float]) -> Optional[float]:
        return float(statistics.median(vals)) if vals else None

    div_vals = collect_numeric("divergent_regions")
    func_vals = collect_numeric("reads_total_functional")
    depth_vals = collect_numeric("input_reads")
    unloc_vals = collect_numeric("unlocalized_fraction")

    return {
        "divergent_regions_total": int(agg_sum(div_vals)) if div_vals else None,
        "divergent_regions_min_per_sample": agg_min(div_vals),
        "divergent_regions_max_per_sample": agg_max(div_vals),
        "divergent_regions_mean_per_sample": agg_mean(div_vals),
        "divergent_regions_median_per_sample": agg_median(div_vals),
        "reads_total_functional_total": agg_sum(func_vals),
        "reads_total_functional_mean": agg_mean(func_vals),
        "reads_total_functional_median": agg_median(func_vals),
        "input_reads_mean": agg_mean(depth_vals),
        "input_reads_median": agg_median(depth_vals),
        "unlocalized_fraction_median": agg_median(unloc_vals),
    }

# =============================================================================
# OUTPUT GENERATION
# =============================================================================

def build_cohort_dataframe(
    samples: List[Dict[str, Any]],
    region_keys: List[str]
) -> pd.DataFrame:
    """
    Build pandas DataFrame with cohort data
    
    Args:
        samples: List of normalized sample dicts
        region_keys: List of all region names
        
    Returns:
        DataFrame with cohort metrics
    """
    log("DATAFRAME", "Building cohort table...")
    
    # Core columns
    core_cols = [
        "sample_id", "condition", "timepoint", "replicate",
        "input_reads", "input_reads_source", "unique_reads_nh1", "dedup_reads", "duplicate_percent",
        "multimapper_percent", "uniqueness_method",
        "divergent_regions", "total_regions", "reads_total_functional",
        "median_pausing_index", "median_density",
        "cpm_factor", "crpmsi_factor", "unlocalized_fraction"
    ]
    
    rows = []
    for sample in samples:
        # Start with core columns
        row = {col: sample.get(col) for col in core_cols}
        
        # Add optional diagnostic fields
        row["density_source"] = sample.get("density_source")
        row["density_reason"] = sample.get("density_reason")
        
        # Add per-region columns
        for region_name in region_keys:
            func_totals = sample.get("func_totals") or {}
            region_counts = sample.get("region_counts") or {}
            region_len_totals = sample.get("region_length_totals") or {}
            region_len_medians = sample.get("region_length_medians") or {}
            
            row[f"func_{region_name}"] = func_totals.get(region_name, 0)
            row[f"count_{region_name}"] = region_counts.get(region_name, 0)
            row[f"len_total_{region_name}"] = region_len_totals.get(region_name, 0)
            row[f"len_median_{region_name}"] = region_len_medians.get(region_name, 0)
        
        rows.append(row)
    
    df = pd.DataFrame(rows)
    log("DATAFRAME", f"Created table with {len(df)} rows × {len(df.columns)} columns")
    
    return df

def write_tsv_output(df: pd.DataFrame, output_path: str):
    """
    Write DataFrame to TSV file
    
    Args:
        df: DataFrame to write
        output_path: Output file path
    """
    log("TSV", f"Writing: {output_path}")
    df.to_csv(output_path, sep="\t", index=False)
    
    file_size = os.path.getsize(output_path)
    log("TSV", f"Written: {file_size:,} bytes")

def write_json_output(
    samples: List[Dict[str, Any]],
    rows: List[Dict[str, Any]],
    region_totals: Dict[str, float],
    region_keys: List[str],
    skipped_files: List[str],
    aggregate: Dict[str, Optional[float]],
    output_path: str
):
    """
    Write cohort JSON output
    
    Args:
        samples: List of sample dicts
        rows: List of row dicts
        region_totals: Aggregated region totals
        region_keys: List of region names
        skipped_files: List of skipped file paths
        output_path: Output file path
    """
    log("JSON", f"Writing: {output_path}")
    
    # Build column list
    core_cols = [
        "sample_id", "condition", "timepoint", "replicate",
        "input_reads", "input_reads_source", "unique_reads_nh1", "dedup_reads", "duplicate_percent",
        "multimapper_percent", "uniqueness_method",
        "divergent_regions", "total_regions", "reads_total_functional",
        "median_pausing_index", "median_density",
        "cpm_factor", "crpmsi_factor", "unlocalized_fraction"
    ]
    
    all_columns = (
        core_cols +
        [f"func_{k}" for k in region_keys] +
        [f"count_{k}" for k in region_keys] +
        [f"len_total_{k}" for k in region_keys] +
        [f"len_median_{k}" for k in region_keys]
    )
    
    # Build payload
    payload = {
        "n_samples": len(samples),
        "rows": rows,
        "samples": samples,
        "region_totals": region_totals,
        "region_keys": region_keys,
        "columns": all_columns,
        "skipped_inputs": skipped_files,
        "aggregate": aggregate,
    }
    
    # Sanitize and write
    sanitized = sanitize_data(payload)
    
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(sanitized, f, indent=2)
    
    file_size = os.path.getsize(output_path)
    log("JSON", f"Written: {file_size:,} bytes")

def write_region_totals_tsv(
    region_totals: Dict[str, float],
    region_keys: List[str],
    output_path: str
):
    """
    Write optional region totals TSV
    
    Args:
        region_totals: Dictionary of region totals
        region_keys: Sorted list of region names
        output_path: Output file path
    """
    if not region_totals:
        log_warning("No region totals to write")
        return
    
    log("REGIONS", f"Writing: {output_path}")
    
    with open(output_path, "w", encoding="utf-8") as f:
        f.write("region\treads\n")
        for region_name in region_keys:
            reads = region_totals.get(region_name, 0)
            f.write(f"{region_name}\t{reads}\n")
    
    file_size = os.path.getsize(output_path)
    log("REGIONS", f"Written: {len(region_keys)} regions ({file_size:,} bytes)")

# =============================================================================
# HTML GENERATION
# =============================================================================

def generate_html_report(
    data_json: str,
    css: str,
    js: str,
    args: argparse.Namespace,
    run_command: str,
    output_path: str
):
    """Generate single-file interactive HTML cohort report (embedded CSS/JS)."""
    import html as _html
    log("HTML", "Generating interactive dashboard...")
    timestamp = datetime.datetime.now().strftime("%B %d, %Y at %H:%M")

    template = r'''<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>TrackTx Cohort Report</title>
%%CSS%%
</head>
<body>
<div id="tt"></div>
<div class="container">

  <header class="page-header">
    <div class="eyebrow">TrackTx · PRO-seq · Cohort report</div>
    <h1>%%RUN_NAME%%</h1>
    <div class="meta">
      <span>Profile <b>%%PROFILE%%</b></span><span class="dot">·</span>
      <span id="hdr-samples">— samples</span><span class="dot">·</span>
      <span>Duration %%DURATION%%</span><span class="dot">·</span>
      <span>Generated %%TIMESTAMP%%</span>
    </div>
    <div class="kpi-strip" id="kpi-strip"></div>
  </header>

  <div class="toolbar">
    <div class="toolbar-inner">
      <nav class="nav" id="nav">
        <a href="#overview">Overview</a>
        <a href="#trends">Trends</a>
        <a href="#qc">Quality</a>
        <a href="#divergent">Divergent TX</a>
        <a href="#pausing">Pausing</a>
        <a href="#regions">Functional regions</a>
        <a href="#normalization">Normalization</a>
        <a href="#samples">Samples</a>
      </nav>
      <div class="controls">
        <div class="control">
          <label>Order</label>
          <span class="seg" id="order-seg">
            <button data-order="timepoint" class="on">Timepoint</button>
            <button data-order="name">Name</button>
          </span>
        </div>
        <label class="control toggle"><input type="checkbox" id="color-by-cond"> Color by condition</label>
      </div>
    </div>
  </div>

  <div class="hl-banner" id="hl-banner">
    <span>Highlighting <b id="hl-name"></b> across all charts.</span>
    <button class="btn" id="hl-clear">Clear</button>
  </div>

  <!-- OVERVIEW -->
  <section class="section" id="overview">
    <div class="head">
      <h2>Overview</h2>
      <div class="sub">Experimental design and per-condition summary</div>
    </div>
    <details class="note-box">
      <summary>About this report</summary>
      <div class="body">
        Aggregates metrics from every per-sample report in this run. Use it to compare conditions, spot outliers,
        and read transcriptional trends across the cohort. Click any bar, point, or table row to highlight that
        sample everywhere; click again to clear. For gene-level statistics, open the individual sample reports.
      </div>
    </details>
    <div class="grid cols-2" style="align-items:start;">
      <div class="card">
        <h3>Design matrix</h3>
        <div class="h3sub">Samples per condition × timepoint</div>
        <div style="overflow-x:auto;"><table class="design-table" id="design-table"></table></div>
      </div>
      <div>
        <div class="grid cols-2" id="condition-cards"></div>
      </div>
    </div>
  </section>

  <!-- TRENDS -->
  <section class="section" id="trends">
    <div class="head">
      <h2>Trends across timepoints</h2>
      <div class="sub">How each metric moves along the ordered timepoint series, split by condition</div>
    </div>
    <div class="card">
      <div class="table-tools" style="margin-bottom:0.8rem;">
        <div class="control">
          <label style="font-size:0.8rem;color:var(--muted);font-weight:600;">Metric</label>
          <select id="trend-metric"></select>
        </div>
        <div class="spacer"></div>
        <div class="legend" id="trend-legend"></div>
      </div>
      <div class="chart" id="chart-trend"></div>
      <p class="h3sub" id="trend-note" style="margin-top:0.6rem;"></p>
    </div>
  </section>

  <!-- QC -->
  <section class="section" id="qc">
    <div class="head">
      <h2>Quality control</h2>
      <div class="sub">Sequencing depth, duplication, multimapping and unlocalized signal per sample</div>
    </div>
    <details class="note-box">
      <summary>How to read these</summary>
      <div class="body">
        <ul>
          <li><strong>Read depth</strong> — PRO-seq usually wants 5–20M+ usable reads. Bars are colored by a rough PASS/WARN/FAIL on depth.</li>
          <li><strong>Duplication / UMI-dedup</strong> — &lt;15% great, 15–30% acceptable, &gt;30% suggests over-amplification or low complexity.</li>
          <li><strong>Multimapper %</strong> — high values flag repetitive genomes; those reads stay in allMap tracks but are excluded from quantification.</li>
          <li><strong>Unlocalized %</strong> — reads outside annotated features; &lt;20% typical, higher may mean rRNA or annotation gaps.</li>
        </ul>
      </div>
    </details>
    <div class="grid cols-4" id="qc-charts">
      <div class="card"><h3>Read depth</h3><div class="h3sub">million reads</div><div class="chart" id="chart-depth"></div></div>
      <div class="card"><h3 id="dup-title">Duplication</h3><div class="h3sub">percent</div><div class="chart" id="chart-dup"></div></div>
      <div class="card"><h3>Multimapper %</h3><div class="h3sub">1 − unique/mapped</div><div class="chart" id="chart-mm"></div></div>
      <div class="card"><h3>Unlocalized %</h3><div class="h3sub">outside features</div><div class="chart" id="chart-unloc"></div></div>
    </div>
  </section>

  <!-- DIVERGENT -->
  <section class="section" id="divergent">
    <div class="head">
      <h2>Divergent transcription</h2>
      <div class="sub">High-confidence bidirectional transcription loci (GMM + FDR)</div>
    </div>
    <details class="note-box">
      <summary>What this means</summary>
      <div class="body">
        Divergent transcription marks active promoters and enhancers where Pol II initiates in both directions.
        Counts are per sample (not de-duplicated across samples). Condition-level shifts can reflect
        stimulus-dependent enhancer activity.
      </div>
    </details>
    <div class="chips" id="div-chips"></div>
    <div class="grid cols-2">
      <div class="card"><h3>Divergent loci per sample</h3><div class="chart" id="chart-div"></div></div>
      <div class="card"><h3>Mean per condition</h3><div class="chart" id="chart-div-cond"></div></div>
    </div>
  </section>

  <!-- PAUSING -->
  <section class="section" id="pausing">
    <div class="head">
      <h2>Pol II pausing</h2>
      <div class="sub">Length-normalized pausing index (promoter density ÷ gene-body density)</div>
    </div>
    <details class="note-box">
      <summary>Interpreting the pausing index</summary>
      <div class="body">
        PI &gt; 1.5 indicates strong promoter-proximal pausing; ≈1 is balanced; &lt;1 favors productive elongation.
        PI should be roughly independent of sequencing depth — a strong depth correlation in the scatter would suggest technical bias.
      </div>
    </details>
    <div class="chips" id="pi-chips"></div>
    <div class="grid cols-2">
      <div class="card"><h3>Median PI per sample</h3><div class="chart" id="chart-pi"></div></div>
      <div class="card"><h3>PI vs read depth</h3><div class="h3sub">bias check</div><div class="chart" id="chart-pi-depth"></div></div>
    </div>
  </section>

  <!-- FUNCTIONAL REGIONS -->
  <section class="section" id="regions">
    <div class="head">
      <h2>Functional region composition</h2>
      <div class="sub">Where Pol II signal lands across genomic features</div>
    </div>
    <div class="grid cols-2" style="align-items:start;">
      <div class="card">
        <h3>Composition per sample</h3>
        <div class="h3sub">share of localized functional signal (%)</div>
        <div class="chart" id="chart-composition"></div>
        <div class="legend" id="composition-legend"></div>
      </div>
      <div class="card">
        <h3>Cohort-wide totals</h3>
        <div class="h3sub">summed reads across all samples</div>
        <div class="chart" id="chart-region-totals"></div>
      </div>
    </div>
    <div class="card" style="margin-top:1rem;">
      <h3>Region read counts by sample</h3>
      <div class="table-wrap" style="margin-top:0.6rem;">
        <table class="data" id="region-table"><thead></thead><tbody></tbody></table>
      </div>
    </div>
  </section>

  <!-- NORMALIZATION -->
  <section class="section" id="normalization">
    <div class="head">
      <h2>Normalization factors</h2>
      <div class="sub">CPM and spike-in (siCPM) scaling for cross-sample comparison</div>
    </div>
    <details class="note-box">
      <summary>CPM vs siCPM</summary>
      <div class="body">
        <strong>CPM</strong> rescales by sequencing depth (assumes similar global transcription).
        <strong>siCPM</strong> uses an exogenous spike-in to capture global changes in transcription — essential when conditions
        are expected to shift overall output. Within a condition, factors should be consistent across replicates.
      </div>
    </details>
    <div class="chips" id="norm-chips"></div>
    <div class="grid cols-3">
      <div class="card"><h3>CPM factor</h3><div class="chart" id="chart-cpm"></div></div>
      <div class="card"><h3>siCPM factor</h3><div class="chart" id="chart-sicpm"></div></div>
      <div class="card"><h3>CPM vs siCPM</h3><div class="chart" id="chart-cpm-sicpm"></div></div>
    </div>
  </section>

  <!-- SAMPLES -->
  <section class="section" id="samples">
    <div class="head">
      <h2>Sample-level metrics</h2>
      <div class="sub">Sortable, searchable table — click a row to highlight that sample in every chart</div>
    </div>
    <div class="table-tools">
      <input type="text" id="sample-search" placeholder="Search sample, condition, timepoint…">
      <select id="condition-filter"><option value="">All conditions</option></select>
      <div class="spacer"></div>
      <button class="btn" id="export-csv">Export CSV</button>
    </div>
    <div class="table-wrap">
      <table class="data" id="sample-table">
        <thead><tr id="sample-thead"></tr></thead>
        <tbody id="sample-tbody"></tbody>
      </table>
    </div>
  </section>

  <!-- DETAILS -->
  <section class="section" id="details">
    <div class="head"><h2>Methodology &amp; files</h2></div>
    <details class="note-box">
      <summary>How metrics are computed</summary>
      <div class="body" id="methodology"></div>
    </details>
    <details class="note-box">
      <summary>Output files &amp; run command</summary>
      <div class="body">
        <h4 style="margin:0 0 0.4rem;">Run command</h4>
        <pre class="cmd"><code>%%RUN_CMD%%</code></pre>
        <h4 style="margin:0.8rem 0 0.4rem;">File locations</h4>
        <ul>
          <li>Cohort HTML: <code>%%OUT_HTML%%</code></li>
          <li>Cohort TSV: <code>%%OUT_TSV%%</code></li>
          <li>Cohort JSON: <code>%%OUT_JSON%%</code></li>
          <li>Sample reports: <code>11_reports/samples/&lt;sample&gt;/&lt;sample&gt;.report.html</code></li>
          <li>Normalized tracks: <code>05_normalized_tracks/&lt;sample&gt;/*.bw</code></li>
          <li>Divergent TX: <code>06_divergent_tx/&lt;sample&gt;/divergent_transcription.bed</code></li>
          <li>Functional regions: <code>07_functional_regions/&lt;sample&gt;/functional_regions.bed</code></li>
        </ul>
      </div>
    </details>
  </section>

</div>

<script type="application/json" id="payload">%%DATA%%</script>
%%JS%%
</body>
</html>'''

    html = (template
        .replace("%%CSS%%", css)
        .replace("%%JS%%", js)
        .replace("%%DATA%%", data_json)
        .replace("%%PROFILE%%", _html.escape(str(args.profile)))
        .replace("%%RUN_NAME%%", _html.escape(str(args.run_name)))
        .replace("%%DURATION%%", _html.escape(str(args.duration)))
        .replace("%%TIMESTAMP%%", _html.escape(str(timestamp)))
        .replace("%%RUN_CMD%%", _html.escape(str(run_command)))
        .replace("%%OUT_HTML%%", _html.escape(str(args.out_html)))
        .replace("%%OUT_TSV%%", _html.escape(str(args.out_tsv)))
        .replace("%%OUT_JSON%%", _html.escape(str(args.out_json)))
    )

    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    file_size = os.path.getsize(output_path)
    log("HTML", f"Written: {file_size:,} bytes ({file_size / 1024 / 1024:.2f} MB)")
# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main execution function"""
    
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Combine TrackTx per-sample summaries into cohort report",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
        help="Per-sample JSONs (files or directories). Supports *.json, *.json.gz"
    )
    parser.add_argument(
        "--out-tsv",
        required=True,
        help="Output TSV path"
    )
    parser.add_argument(
        "--out-json",
        required=True,
        help="Output JSON path"
    )
    parser.add_argument(
        "--out-html",
        required=True,
        help="Output HTML path"
    )
    parser.add_argument(
        "--out-regions",
        default=None,
        help="Optional TSV of summed functional-region totals"
    )
    parser.add_argument(
        "--pipeline-version",
        default="unknown",
        help="Pipeline version for metadata"
    )
    parser.add_argument(
        "--run-name",
        default="unnamed",
        help="Run name for metadata"
    )
    parser.add_argument(
        "--duration",
        default="unknown",
        help="Run duration for metadata"
    )
    parser.add_argument(
        "--profile",
        default="unknown",
        help="Execution profile for metadata"
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {VERSION}"
    )
    
    args = parser.parse_args()
    
    # Capture run command
    run_command = " ".join(shlex.quote(arg) for arg in sys.argv)
    
    # Start
    log("START", f"combine_reports.py v{VERSION}")
    log("START", f"Input paths: {len(args.inputs)}")
    log("START", f"Pipeline version: {args.pipeline_version}")
    
    # Discover JSON files
    log("═" * 70, "")
    json_files = discover_json_files(args.inputs)
    
    if not json_files:
        log_error("No usable input files found")
        log_error("Expected: *.summary.json, *.report.json, *.json, or *.json.gz")
        return 2
    
    # Load and normalize samples
    log("═" * 70, "")
    log("LOAD", f"Loading {len(json_files)} JSON files...")
    
    samples = []
    skipped_files = []
    
    for i, filepath in enumerate(json_files, 1):
        if i % 10 == 0:
            log_info(f"Loaded {i}/{len(json_files)} files...")
        
        json_data = read_json_file(filepath)
        if not json_data:
            skipped_files.append(filepath)
            continue

        # Schema version check
        schema = json_data.get("schema_version", "unknown")
        if schema not in ("1.0", "1.0.0", "unknown"):
            log_warning(f"Unexpected schema version '{schema}' in {filepath}")

        fallback_name = Path(filepath).stem.split(".")[0]
        sample = normalize_sample_data(json_data, fallback_name)
        if sample and sample.get("sample_id") == fallback_name and fallback_name.startswith("report_"):
            log_warning(f"Sample ID from fallback: {fallback_name} (JSON missing sample field)")
        
        if sample:
            samples.append(sample)
        else:
            skipped_files.append(filepath)
    
    log("LOAD", f"Successfully loaded: {len(samples)} samples")
    if skipped_files:
        log_warning(f"Skipped: {len(skipped_files)} files")
        for path in skipped_files[:5]:  # Show first 5
            log_warning(f"  {path}")
        if len(skipped_files) > 5:
            log_warning(f"  ... and {len(skipped_files) - 5} more")
    
    if not samples:
        log_error("No valid samples loaded")
        return 3

    # Sort by experimental design (condition, timepoint, replicate, sample_id)
    samples.sort(key=lambda s: (
        str(s.get("condition", "")),
        str(s.get("timepoint", "")),
        str(s.get("replicate", "")),
        str(s.get("sample_id", "")),
    ))
    
    # Aggregate region data
    log("═" * 70, "")
    log("AGGREGATE", "Collecting region information...")
    region_keys = collect_region_keys(samples)
    log("AGGREGATE", f"Found {len(region_keys)} unique regions")
    
    region_totals = aggregate_region_totals(samples)
    log("AGGREGATE", f"Computed totals for {len(region_totals)} regions")
    
    # Build DataFrame
    log("═" * 70, "")
    df = build_cohort_dataframe(samples, region_keys)
    
    # Extract rows for JSON
    rows = df.to_dict("records")

    # Aggregate cohort-level summary statistics
    aggregate = compute_aggregate_summary(rows)
    
    # Write outputs
    log("═" * 70, "")
    write_tsv_output(df, args.out_tsv)
    
    write_json_output(
        samples,
        rows,
        region_totals,
        region_keys,
        skipped_files,
        aggregate,
        args.out_json
    )
    
    if args.out_regions:
        write_region_totals_tsv(region_totals, region_keys, args.out_regions)
    
    # Generate HTML with embedded CSS/JS
    # (Keeping existing CSS and JS strings from original)
    log("═" * 70, "")
    data_json = json.dumps(sanitize_data({
        "n_samples": len(samples),
        "rows": rows,
        "samples": samples,
        "region_totals": region_totals,
        "region_keys": region_keys,
        "aggregate": aggregate
    }), ensure_ascii=False)
    
    # Embedded Assets
    CSS = r'''
<style>
:root {
  --bg:#ffffff; --fg:#1a1d23; --muted:#6b7280; --faint:#9aa1ab;
  --card:#ffffff; --panel:#f7f8fa; --line:#e6e8ec; --line-strong:#d4d8de;
  --accent:#2f6df6; --accent-soft:#eaf0fe;
  --ok:#1f9d57; --warn:#d98a00; --fail:#dc3a3a;
  --grid:#eceef1;
  --font:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
  --mono:ui-monospace,SFMono-Regular,"SF Mono",Menlo,Consolas,monospace;
}
@media (prefers-color-scheme: dark) {
  :root {
    --bg:#0e1116; --fg:#e6e8ec; --muted:#9aa1ab; --faint:#6b7280;
    --card:#161a21; --panel:#12161c; --line:#262c35; --line-strong:#333b46;
    --accent:#5b8cff; --accent-soft:#1a2230;
    --ok:#37c172; --warn:#e0a032; --fail:#ef5e5e;
    --grid:#1d222a;
  }
}

* { box-sizing:border-box; }
html { scroll-behavior:smooth; }
body { background:var(--bg); color:var(--fg); font-family:var(--font); margin:0; line-height:1.55; font-size:15px; -webkit-font-smoothing:antialiased; }
.container { max-width:1240px; margin:0 auto; padding:0 1.5rem 5rem; }
a { color:var(--accent); text-decoration:none; }
a:hover { text-decoration:underline; }
code, .mono { font-family:var(--mono); font-size:0.85em; }

/* Header */
.page-header { padding:2.25rem 0 1.25rem; border-bottom:1px solid var(--line); }
.eyebrow { text-transform:uppercase; font-size:0.7rem; font-weight:700; letter-spacing:0.12em; color:var(--accent); }
.page-header h1 { margin:0.35rem 0 0.4rem; font-size:1.7rem; font-weight:700; letter-spacing:-0.01em; }
.page-header .meta { color:var(--muted); font-size:0.875rem; }
.page-header .meta span { white-space:nowrap; }
.page-header .meta .dot { margin:0 0.5rem; color:var(--faint); }

/* KPI strip */
.kpi-strip { display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr)); gap:1px; background:var(--line); border:1px solid var(--line); border-radius:10px; overflow:hidden; margin:1.5rem 0 0.5rem; }
.kpi { background:var(--card); padding:0.9rem 1.1rem; }
.kpi .k-label { font-size:0.72rem; color:var(--muted); font-weight:600; text-transform:uppercase; letter-spacing:0.04em; }
.kpi .k-value { font-size:1.45rem; font-weight:700; margin-top:0.2rem; letter-spacing:-0.01em; }
.kpi .k-sub { font-size:0.72rem; color:var(--faint); margin-top:0.1rem; }

/* Sticky control bar (nav + controls) */
.toolbar { position:sticky; top:0; z-index:50; background:var(--bg); border-bottom:1px solid var(--line); padding:0.6rem 0; margin-bottom:2rem; }
.toolbar-inner { display:flex; align-items:center; gap:1rem; flex-wrap:wrap; }
.nav { display:flex; gap:0.15rem; flex-wrap:wrap; flex:1; }
.nav a { padding:0.35rem 0.7rem; border-radius:7px; color:var(--muted); font-weight:600; font-size:0.82rem; }
.nav a:hover { background:var(--panel); color:var(--fg); text-decoration:none; }
.nav a.active { background:var(--accent-soft); color:var(--accent); }
.controls { display:flex; align-items:center; gap:0.9rem; flex-wrap:wrap; }
.control { display:flex; align-items:center; gap:0.4rem; font-size:0.8rem; color:var(--muted); }
.control label { font-weight:600; }
.seg { display:inline-flex; border:1px solid var(--line-strong); border-radius:7px; overflow:hidden; }
.seg button { border:none; background:var(--card); color:var(--muted); padding:0.3rem 0.6rem; font-size:0.78rem; font-weight:600; cursor:pointer; }
.seg button + button { border-left:1px solid var(--line); }
.seg button.on { background:var(--accent); color:#fff; }
.toggle { display:inline-flex; align-items:center; gap:0.4rem; cursor:pointer; user-select:none; }
.toggle input { accent-color:var(--accent); width:auto; }

/* Sections */
.section { padding-top:2.4rem; margin-top:1rem; scroll-margin-top:64px; border-top:1px solid var(--line); }
.section:first-of-type { border-top:none; }
.section > .head { margin-bottom:1.1rem; }
.section > .head h2 { margin:0; font-size:1.25rem; font-weight:700; letter-spacing:-0.01em; }
.section > .head .sub { color:var(--muted); font-size:0.875rem; margin-top:0.15rem; }
.section > .head .sub a.note { font-size:0.8rem; margin-left:0.5rem; }

/* Collapsible note */
details.note-box { margin:0.4rem 0 1rem; }
details.note-box > summary { cursor:pointer; font-size:0.82rem; color:var(--accent); font-weight:600; list-style:none; }
details.note-box > summary::-webkit-details-marker { display:none; }
details.note-box > summary::before { content:"ⓘ "; }
details.note-box[open] > summary::before { content:"▾ "; }
details.note-box .body { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:0.9rem 1.1rem; margin-top:0.5rem; font-size:0.86rem; color:var(--fg); }
details.note-box .body ul { margin:0.4rem 0 0 1.1rem; padding:0; }
details.note-box .body li { margin:0.25rem 0; }
details.note-box .body strong { color:var(--fg); }

/* Cards & grids */
.grid { display:grid; gap:1rem; }
.cols-2 { grid-template-columns:repeat(2,1fr); }
.cols-3 { grid-template-columns:repeat(3,1fr); }
.cols-4 { grid-template-columns:repeat(4,1fr); }
@media (max-width:980px){ .cols-2,.cols-3,.cols-4 { grid-template-columns:1fr; } }
.card { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:1rem 1.1rem; }
.card > h3 { margin:0 0 0.15rem; font-size:0.95rem; font-weight:700; }
.card > .h3sub { font-size:0.78rem; color:var(--faint); margin:0 0 0.5rem; }
.chart { width:100%; }
.chart svg { display:block; width:100%; height:auto; overflow:visible; }
.empty { padding:2.2rem 1rem; text-align:center; color:var(--faint); font-size:0.85rem; }

/* Stat chips */
.chips { display:flex; gap:0.5rem; flex-wrap:wrap; margin-bottom:1rem; }
.chip { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:0.45rem 0.75rem; }
.chip .c-label { font-size:0.68rem; color:var(--muted); text-transform:uppercase; letter-spacing:0.04em; font-weight:600; }
.chip .c-value { font-size:1.05rem; font-weight:700; }

/* Condition / design cards */
.cond-card { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:0.85rem 1rem; border-left:3px solid var(--cond,var(--accent)); }
.cond-card .cc-name { font-weight:700; font-size:0.95rem; }
.cond-card .cc-row { font-size:0.8rem; color:var(--muted); margin-top:0.25rem; }
.flag { display:inline-block; font-size:0.68rem; font-weight:700; padding:0.05rem 0.4rem; border-radius:5px; vertical-align:middle; }
.flag.ok { background:rgba(31,157,87,0.14); color:var(--ok); }
.flag.warn { background:rgba(217,138,0,0.16); color:var(--warn); }
.flag.fail { background:rgba(220,58,58,0.14); color:var(--fail); }

/* Design grid (condition x timepoint) */
.design-table { border-collapse:collapse; font-size:0.82rem; }
.design-table th, .design-table td { border:1px solid var(--line); padding:0.4rem 0.6rem; text-align:center; }
.design-table th { background:var(--panel); font-weight:700; color:var(--muted); }
.design-table td.has { background:var(--accent-soft); color:var(--accent); font-weight:700; }
.design-table td.empty-cell { color:var(--faint); }

/* Tables */
.table-tools { display:flex; gap:0.75rem; flex-wrap:wrap; align-items:center; margin-bottom:0.9rem; }
.table-tools input[type=text], .table-tools select { padding:0.45rem 0.7rem; border:1px solid var(--line-strong); border-radius:7px; background:var(--card); color:var(--fg); font-size:0.85rem; }
.table-tools input[type=text] { min-width:240px; }
.table-tools .spacer { flex:1; }
.btn { padding:0.45rem 0.8rem; border:1px solid var(--line-strong); border-radius:7px; background:var(--card); color:var(--fg); font-size:0.82rem; font-weight:600; cursor:pointer; }
.btn:hover { background:var(--panel); }
.btn.primary { background:var(--accent); color:#fff; border-color:var(--accent); }
.table-wrap { overflow-x:auto; border:1px solid var(--line); border-radius:10px; }
table.data { width:100%; border-collapse:collapse; font-size:0.84rem; }
table.data th, table.data td { padding:0.55rem 0.75rem; text-align:right; border-bottom:1px solid var(--line); white-space:nowrap; }
table.data th:first-child, table.data td:first-child,
table.data th.lft, table.data td.lft { text-align:left; }
table.data thead th { background:var(--panel); position:sticky; top:0; font-weight:700; color:var(--muted); font-size:0.72rem; text-transform:uppercase; letter-spacing:0.03em; cursor:pointer; user-select:none; z-index:2; }
table.data thead th.sortable:hover { color:var(--fg); }
table.data thead th .arrow { color:var(--faint); font-size:0.85em; }
table.data tbody tr:hover { background:var(--panel); }
table.data tbody tr.hl { background:var(--accent-soft) !important; }
table.data tbody tr:last-child td { border-bottom:none; }
table.data td.num { font-variant-numeric:tabular-nums; }

/* Floating tooltip */
#tt { position:fixed; z-index:9999; pointer-events:none; background:var(--fg); color:var(--bg); padding:0.45rem 0.6rem; border-radius:7px; font-size:0.78rem; line-height:1.35; box-shadow:0 6px 22px rgba(0,0,0,0.25); max-width:280px; opacity:0; transition:opacity 0.08s; }
#tt .tt-title { font-weight:700; margin-bottom:0.15rem; }
#tt .tt-row { display:flex; justify-content:space-between; gap:0.8rem; }
#tt .tt-row .v { font-variant-numeric:tabular-nums; }
#tt .sw { display:inline-block; width:9px; height:9px; border-radius:2px; margin-right:0.35rem; vertical-align:middle; }

/* Legend */
.legend { display:flex; flex-wrap:wrap; gap:0.5rem 1rem; margin-top:0.6rem; font-size:0.78rem; color:var(--muted); }
.legend .item { display:inline-flex; align-items:center; gap:0.35rem; cursor:default; }
.legend .sw { width:11px; height:11px; border-radius:3px; }

/* Methodology */
.method-row { padding:0.7rem 0; border-bottom:1px solid var(--line); }
.method-row:last-child { border-bottom:none; }
.method-row .m-name { font-weight:700; font-size:0.9rem; }
.method-row .m-formula { font-family:var(--mono); font-size:0.8rem; color:var(--accent); margin:0.2rem 0; }
.method-row .m-desc { font-size:0.82rem; color:var(--muted); }
pre.cmd { background:var(--panel); border:1px solid var(--line); border-radius:8px; padding:0.8rem; overflow-x:auto; font-size:0.78rem; }

.hl-banner { display:none; align-items:center; gap:0.6rem; background:var(--accent-soft); border:1px solid var(--accent); color:var(--accent); border-radius:8px; padding:0.45rem 0.8rem; font-size:0.82rem; font-weight:600; margin-bottom:1rem; }
.hl-banner button { margin-left:auto; }
</style>
'''

    JS = r'''
<script>
(function () {
  "use strict";
  var payload = JSON.parse(document.getElementById('payload').textContent);
  var rows = payload.rows || [];
  var region_totals = payload.region_totals || {};
  var region_keys = payload.region_keys || [];
  var aggregate = payload.aggregate || {};

  // ---------- state ----------
  var state = { order: 'timepoint', colorByCondition: false, highlight: null };
  var redraws = [];
  function register(fn) { redraws.push(fn); }
  function redrawAll() { redraws.forEach(function (f) { try { f(); } catch (e) { console.error(e); } }); }

  // ---------- math / format ----------
  function mean(a) { return a.length ? a.reduce(function (x, y) { return x + y; }, 0) / a.length : 0; }
  function median(a) { if (!a.length) return 0; var s = a.slice().sort(function (x, y) { return x - y; }); var m = Math.floor(s.length / 2); return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }
  function std(a) { if (a.length < 2) return 0; var m = mean(a); return Math.sqrt(mean(a.map(function (v) { return (v - m) * (v - m); }))); }
  function cv(a) { if (a.length < 2) return null; var m = mean(a); return m ? std(a) / m * 100 : null; }
  function isNum(v) { return v != null && isFinite(v); }
  function fmt(n) {
    if (n == null || !isFinite(n)) return '–';
    var a = Math.abs(n);
    if (a >= 1e6) return (n / 1e6).toFixed(1) + 'M';
    if (a >= 1e3) return (n / 1e3).toFixed(1) + 'K';
    if (a === 0) return '0';
    if (a < 1) return n.toFixed(3);
    if (a < 10) return n.toFixed(2);
    if (a < 100) return n.toFixed(1);
    return Math.round(n).toLocaleString();
  }
  function fmtInt(n) { return (n == null) ? '–' : Math.round(n).toLocaleString(); }
  function shortName(id) { return String(id || '').replace(/_merged$/, ''); }
  function parseTP(v) { var x = parseFloat(v); return isFinite(x) ? x : null; }

  // ---------- colors ----------
  var ACCENT = cssVar('--accent') || '#2f6df6';
  var COND_PALETTE = ['#2f6df6', '#e8843c', '#1f9d57', '#b23bd4', '#d83a7a', '#0f9bb0', '#8a6d3b', '#6366f1'];
  function cssVar(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
  var conditions = uniq(rows.map(function (r) { return r.condition || 'Unknown'; }));
  var condColor = {};
  conditions.forEach(function (c, i) { condColor[c] = COND_PALETTE[i % COND_PALETTE.length]; });
  function uniq(a) { return a.filter(function (v, i) { return a.indexOf(v) === i; }); }

  var REGION_COLORS = {
    'promoter': '#f38400', 'activepromoter': '#f38400', 'pppol': '#f38400',
    'divergenttx': '#b23bd4', 'divtx': '#b23bd4', 'divergent': '#b23bd4',
    'enhancers': '#5bbf6a', 'enhancer': '#5bbf6a', 'enh': '#5bbf6a',
    'genebody': '#3a3f4a', 'cps': '#3aa7e0', 'cleavagepolyadenylation': '#3aa7e0',
    'terminationwindow': '#ef5470', 'termination': '#ef5470', 'tw': '#ef5470',
    'nonlocalizedpolymerase': '#aeb4bd', 'nonlocalized': '#aeb4bd', 'unlocalized': '#aeb4bd'
  };
  function regionColor(name) {
    var k = String(name || '').toLowerCase().replace(/[\s_-]+/g, '');
    return REGION_COLORS[k] || '#8a93a0';
  }
  function barColor(row) { return state.colorByCondition ? (condColor[row.condition || 'Unknown'] || ACCENT) : ACCENT; }

  // ---------- ordering ----------
  function orderedRows() {
    var r = rows.slice();
    if (state.order === 'name') {
      r.sort(function (a, b) { return String(a.sample_id).localeCompare(String(b.sample_id)); });
    } else {
      r.sort(function (a, b) {
        var ta = parseTP(a.timepoint), tb = parseTP(b.timepoint);
        if (ta != null && tb != null && ta !== tb) return ta - tb;
        var ca = String(a.condition || ''), cb = String(b.condition || '');
        if (ca !== cb) return ca.localeCompare(cb);
        return String(a.sample_id).localeCompare(String(b.sample_id));
      });
    }
    return r;
  }

  // ---------- tooltip ----------
  var tt = document.getElementById('tt');
  function ttShow(html) { tt.innerHTML = html; tt.style.opacity = '1'; }
  function ttMove(e) {
    var pad = 14, w = tt.offsetWidth, h = tt.offsetHeight;
    var x = e.clientX + pad, y = e.clientY + pad;
    if (x + w > window.innerWidth - 8) x = e.clientX - w - pad;
    if (y + h > window.innerHeight - 8) y = e.clientY - h - pad;
    tt.style.left = x + 'px'; tt.style.top = y + 'px';
  }
  function ttHide() { tt.style.opacity = '0'; }
  function rowTip(row, extra) {
    var s = '<div class="tt-title">' + shortName(row.sample_id) + '</div>';
    s += '<div class="tt-row"><span>' + (row.condition || '–') + '</span><span class="v">tp ' + (row.timepoint != null ? row.timepoint : '–') + '</span></div>';
    if (extra) s += extra;
    return s;
  }
  function tipRow(label, val, color) {
    return '<div class="tt-row"><span>' + (color ? '<span class="sw" style="background:' + color + '"></span>' : '') + label + '</span><span class="v">' + val + '</span></div>';
  }

  // ---------- highlight ----------
  function setHighlight(id) {
    state.highlight = (state.highlight === id) ? null : id;
    var b = document.getElementById('hl-banner');
    if (state.highlight) { document.getElementById('hl-name').textContent = shortName(state.highlight); b.style.display = 'flex'; }
    else b.style.display = 'none';
    redrawAll();
    syncTableHl();
  }
  function dim(id) { return state.highlight && state.highlight !== id ? 0.22 : 1; }
  function isHl(id) { return state.highlight === id; }

  // ---------- svg helpers ----------
  function svgEl(w, h) { return '<svg viewBox="0 0 ' + w + ' ' + h + '" preserveAspectRatio="xMidYMid meet" font-family="var(--font)">'; }
  function widthOf(el) { var w = el.getBoundingClientRect().width; if (!w && el.parentElement) w = el.parentElement.getBoundingClientRect().width; return Math.max(260, Math.round(w || 360)); }
  function txt(x, y, s, opts) { opts = opts || {}; return '<text x="' + x + '" y="' + y + '" text-anchor="' + (opts.anchor || 'middle') + '" font-size="' + (opts.size || 10) + '" fill="' + (opts.fill || 'var(--muted)') + '"' + (opts.weight ? ' font-weight="' + opts.weight + '"' : '') + (opts.rotate ? ' transform="rotate(' + opts.rotate + ' ' + x + ' ' + y + ')"' : '') + '>' + s + '</text>'; }
  function niceMax(v) { if (v <= 0) return 1; var mag = Math.pow(10, Math.floor(Math.log10(v))); return Math.ceil(v / mag) * mag; }
  function empty(el, msg) { el.innerHTML = '<div class="empty">' + (msg || 'No data') + '</div>'; }

  // ---------- generic vertical bar chart ----------
  // items: [{id, label, value, row}]
  function drawBars(elId, getItems, opts) {
    opts = opts || {};
    var el = document.getElementById(elId);
    if (!el) return;
    function render() {
      var items = getItems();
      if (!items || !items.length || items.every(function (d) { return !isNum(d.value); })) { empty(el, opts.emptyMsg); return; }
      var W = widthOf(el), H = opts.height || 230;
      var m = { t: 18, r: 10, b: 52, l: 44 };
      var iw = W - m.l - m.r, ih = H - m.t - m.b;
      var vals = items.map(function (d) { return isNum(d.value) ? d.value : 0; });
      var dmax = Math.max.apply(null, vals.concat([0]));
      var yMax = opts.yMax != null ? opts.yMax : (dmax <= 1 ? Math.min(1, dmax * 1.15 || 1) : niceMax(dmax * 1.05));
      if (yMax <= 0) yMax = 1;
      var n = items.length;
      var step = iw / n, bw = Math.min(opts.maxBar || 64, step * 0.66);
      var lblEvery = n <= 18 ? 1 : Math.ceil(n / 18);
      var s = svgEl(W, H);
      // gridlines + y labels
      for (var g = 0; g <= 4; g++) {
        var gy = m.t + ih * g / 4;
        s += '<line x1="' + m.l + '" y1="' + gy + '" x2="' + (W - m.r) + '" y2="' + gy + '" stroke="var(--grid)" stroke-width="1"/>';
        s += txt(m.l - 6, gy + 3, fmt(yMax - yMax * g / 4), { anchor: 'end', size: 9 });
      }
      // threshold band
      if (opts.threshold && opts.threshold.warn != null) {
        var wy = m.t + ih * (1 - Math.min(1, opts.threshold.warn / yMax));
        s += '<line x1="' + m.l + '" y1="' + wy + '" x2="' + (W - m.r) + '" y2="' + wy + '" stroke="var(--warn)" stroke-dasharray="3 3" stroke-width="1" opacity="0.6"/>';
      }
      items.forEach(function (d, i) {
        var v = isNum(d.value) ? d.value : 0;
        var bh = Math.max(v > 0 ? 2 : 0, (v / yMax) * ih);
        var x = m.l + step * i + (step - bw) / 2;
        var y = m.t + ih - bh;
        var col = opts.colorFn ? opts.colorFn(d) : barColor(d.row);
        var op = dim(d.id);
        var stroke = isHl(d.id) ? ' stroke="var(--fg)" stroke-width="1.5"' : '';
        s += '<rect class="bar" data-id="' + esc(d.id) + '" data-i="' + i + '" x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' + bw.toFixed(1) + '" height="' + bh.toFixed(1) + '" rx="2.5" fill="' + col + '" opacity="' + op + '"' + stroke + ' style="cursor:pointer"/>';
        if (bh > 16 && bw > 22) s += txt(x + bw / 2, y - 4, opts.valFmt ? opts.valFmt(v) : fmt(v), { size: 9, fill: 'var(--fg)', weight: '600' });
        if (i % lblEvery === 0) {
          var lab = state.order === 'timepoint' && d.row.timepoint != null && d.row.timepoint !== '' ? String(d.row.timepoint) : shortName(d.id);
          if (lab.length > 9) lab = lab.slice(0, 8) + '…';
          s += txt(x + bw / 2, H - m.b + 12, esc(lab), { size: 9, anchor: 'end', rotate: -42 });
        }
      });
      if (opts.xLabel) s += txt(m.l + iw / 2, H - 4, opts.xLabel, { size: 9, weight: '600' });
      s += '</svg>';
      el.innerHTML = s;
      el.querySelectorAll('rect.bar').forEach(function (r) {
        var d = items[+r.getAttribute('data-i')];
        r.addEventListener('mousemove', function (e) { ttShow(rowTip(d.row, tipRow(opts.tipLabel || 'Value', opts.tipFmt ? opts.tipFmt(d.value) : fmt(d.value)))); ttMove(e); });
        r.addEventListener('mouseleave', ttHide);
        r.addEventListener('click', function () { ttHide(); setHighlight(d.id); });
      });
    }
    register(render); render();
  }

  // ---------- horizontal ranked bars (region totals) ----------
  function drawHBars(elId, items, opts) {
    opts = opts || {};
    var el = document.getElementById(elId);
    if (!el) return;
    function render() {
      if (!items.length) { empty(el); return; }
      var W = widthOf(el), rowH = 26, m = { t: 8, r: 54, b: 8, l: 110 };
      var H = m.t + m.b + items.length * rowH;
      var iw = W - m.l - m.r;
      var dmax = Math.max.apply(null, items.map(function (d) { return d.value; }).concat([1]));
      var s = svgEl(W, H);
      items.forEach(function (d, i) {
        var y = m.t + i * rowH, bw = Math.max(1, d.value / dmax * iw), bh = rowH - 9;
        s += txt(m.l - 8, y + bh / 2 + 3, esc(d.label), { anchor: 'end', size: 10, fill: 'var(--fg)' });
        s += '<rect x="' + m.l + '" y="' + y + '" width="' + bw.toFixed(1) + '" height="' + bh + '" rx="2.5" fill="' + (opts.colorFn ? opts.colorFn(d) : ACCENT) + '"/>';
        s += txt(m.l + bw + 5, y + bh / 2 + 3, fmt(d.value), { anchor: 'start', size: 9, weight: '600', fill: 'var(--muted)' });
      });
      s += '</svg>'; el.innerHTML = s;
    }
    register(render); render();
  }

  // ---------- 100% stacked composition ----------
  function drawStacked(elId, legendId, keys) {
    var el = document.getElementById(elId);
    if (!el) return;
    function render() {
      var items = orderedRows();
      if (!items.length || !keys.length) { empty(el); return; }
      var W = widthOf(el), H = 280, m = { t: 14, r: 10, b: 54, l: 38 };
      var iw = W - m.l - m.r, ih = H - m.t - m.b;
      var n = items.length, step = iw / n, bw = Math.min(58, step * 0.7);
      var lblEvery = n <= 18 ? 1 : Math.ceil(n / 18);
      var s = svgEl(W, H);
      for (var g = 0; g <= 4; g++) { var gy = m.t + ih * g / 4; s += '<line x1="' + m.l + '" y1="' + gy + '" x2="' + (W - m.r) + '" y2="' + gy + '" stroke="var(--grid)"/>'; s += txt(m.l - 6, gy + 3, (100 - 25 * g) + '%', { anchor: 'end', size: 9 }); }
      items.forEach(function (r, i) {
        var tot = keys.reduce(function (a, k) { return a + (r['func_' + k] || 0); }, 0) || 1;
        var x = m.l + step * i + (step - bw) / 2, yA = m.t + ih, op = dim(r.sample_id);
        var hlStroke = isHl(r.sample_id) ? '<rect x="' + (x - 2) + '" y="' + (m.t - 2) + '" width="' + (bw + 4) + '" height="' + (ih + 4) + '" fill="none" stroke="var(--fg)" stroke-width="1.5" rx="3"/>' : '';
        keys.forEach(function (k) {
          var pct = (r['func_' + k] || 0) / tot, sh = pct * ih;
          if (sh <= 0) return;
          var y = yA - sh; yA = y;
          s += '<rect class="seg" data-id="' + esc(r.sample_id) + '" data-k="' + esc(k) + '" x="' + x.toFixed(1) + '" y="' + y.toFixed(1) + '" width="' + bw.toFixed(1) + '" height="' + sh.toFixed(1) + '" fill="' + regionColor(k) + '" opacity="' + op + '" style="cursor:pointer"/>';
        });
        s += hlStroke;
        if (i % lblEvery === 0) {
          var lab = state.order === 'timepoint' && r.timepoint != null && r.timepoint !== '' ? String(r.timepoint) : shortName(r.sample_id);
          if (lab.length > 9) lab = lab.slice(0, 8) + '…';
          s += txt(x + bw / 2, H - m.b + 12, esc(lab), { size: 9, anchor: 'end', rotate: -42 });
        }
      });
      s += '</svg>'; el.innerHTML = s;
      el.querySelectorAll('rect.seg').forEach(function (rc) {
        var id = rc.getAttribute('data-id'), k = rc.getAttribute('data-k');
        var row = rows.find(function (r) { return r.sample_id === id; });
        var tot = keys.reduce(function (a, kk) { return a + (row['func_' + kk] || 0); }, 0) || 1;
        var reads = row['func_' + k] || 0;
        rc.addEventListener('mousemove', function (e) { ttShow(rowTip(row, tipRow(k, (reads / tot * 100).toFixed(1) + '%', regionColor(k)) + tipRow('reads', fmt(reads)))); ttMove(e); });
        rc.addEventListener('mouseleave', ttHide);
        rc.addEventListener('click', function () { ttHide(); setHighlight(id); });
      });
    }
    register(render); render();
    // legend
    var lg = document.getElementById(legendId);
    if (lg) lg.innerHTML = keys.map(function (k) { return '<span class="item"><span class="sw" style="background:' + regionColor(k) + '"></span>' + esc(k) + '</span>'; }).join('');
  }

  // ---------- scatter ----------
  function drawScatter(elId, getPoints, xLabel, yLabel, opts) {
    opts = opts || {};
    var el = document.getElementById(elId);
    if (!el) return;
    function render() {
      var pts = getPoints();
      if (!pts.length) { empty(el, opts.emptyMsg); return; }
      var W = widthOf(el), H = 230, m = { t: 14, r: 14, b: 42, l: 48 };
      var iw = W - m.l - m.r, ih = H - m.t - m.b;
      var xs = pts.map(function (p) { return p.x; }), ys = pts.map(function (p) { return p.y; });
      var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
      var yMin = Math.min.apply(null, ys), yMax = Math.max.apply(null, ys);
      var xr = (xMax - xMin) || 1, yr = (yMax - yMin) || 1;
      xMin -= xr * 0.08; xMax += xr * 0.08; yMin -= yr * 0.08; yMax += yr * 0.08;
      var sx = function (v) { return m.l + (v - xMin) / (xMax - xMin) * iw; };
      var sy = function (v) { return m.t + ih - (v - yMin) / (yMax - yMin) * ih; };
      var s = svgEl(W, H);
      for (var g = 0; g <= 4; g++) {
        var gy = m.t + ih * g / 4; s += '<line x1="' + m.l + '" y1="' + gy + '" x2="' + (W - m.r) + '" y2="' + gy + '" stroke="var(--grid)"/>';
        s += txt(m.l - 6, gy + 3, fmt(yMax - (yMax - yMin) * g / 4), { anchor: 'end', size: 9 });
        var gx = m.l + iw * g / 4; s += txt(gx, H - m.b + 14, fmt(xMin + (xMax - xMin) * g / 4), { size: 9 });
      }
      pts.forEach(function (p) {
        s += '<circle class="pt" data-id="' + esc(p.id) + '" cx="' + sx(p.x).toFixed(1) + '" cy="' + sy(p.y).toFixed(1) + '" r="' + (isHl(p.id) ? 7 : 5) + '" fill="' + barColor(p.row) + '" opacity="' + dim(p.id) + '" stroke="var(--card)" stroke-width="1.5" style="cursor:pointer"/>';
      });
      s += txt(m.l + iw / 2, H - 3, xLabel, { size: 9, weight: '600' });
      s += '<text x="12" y="' + (m.t + ih / 2) + '" text-anchor="middle" font-size="9" font-weight="600" fill="var(--muted)" transform="rotate(-90 12 ' + (m.t + ih / 2) + ')">' + yLabel + '</text>';
      s += '</svg>'; el.innerHTML = s;
      el.querySelectorAll('circle.pt').forEach(function (c) {
        var p = pts.find(function (q) { return q.id === c.getAttribute('data-id'); });
        c.addEventListener('mousemove', function (e) { ttShow(rowTip(p.row, tipRow(xLabel, fmt(p.x)) + tipRow(yLabel, fmt(p.y)))); ttMove(e); });
        c.addEventListener('mouseleave', ttHide);
        c.addEventListener('click', function () { ttHide(); setHighlight(p.id); });
      });
    }
    register(render); render();
  }

  // ---------- line chart (trends) ----------
  function drawTrend() {
    var el = document.getElementById('chart-trend');
    if (!el) return;
    var sel = document.getElementById('trend-metric');
    var metric = TREND_METRICS.find(function (m) { return m.key === sel.value; }) || TREND_METRICS[0];
    // group by condition; x = numeric timepoint
    var series = {};
    rows.forEach(function (r) {
      var tp = parseTP(r.timepoint); if (tp == null) return;
      var v = metric.get(r); if (!isNum(v)) return;
      var c = r.condition || 'Unknown';
      (series[c] = series[c] || []).push({ x: tp, y: v, row: r });
    });
    var names = Object.keys(series).filter(function (c) { return series[c].length; });
    names.forEach(function (c) { series[c].sort(function (a, b) { return a.x - b.x; }); });
    var allPts = names.reduce(function (a, c) { return a.concat(series[c]); }, []);
    if (allPts.length < 2) { empty(el, 'Numeric timepoints required for trends'); document.getElementById('trend-legend').innerHTML = ''; return; }
    var W = widthOf(el), H = 300, m = { t: 16, r: 16, b: 44, l: 52 };
    var iw = W - m.l - m.r, ih = H - m.t - m.b;
    var xs = allPts.map(function (p) { return p.x; }), ys = allPts.map(function (p) { return p.y; });
    var xMin = Math.min.apply(null, xs), xMax = Math.max.apply(null, xs);
    var yMin = Math.min.apply(null, ys.concat(metric.zero ? [0] : [])), yMax = Math.max.apply(null, ys);
    var yr = (yMax - yMin) || 1; yMin -= yr * 0.08; yMax += yr * 0.08;
    if (xMin === xMax) xMax = xMin + 1;
    var sx = function (v) { return m.l + (v - xMin) / (xMax - xMin) * iw; };
    var sy = function (v) { return m.t + ih - (v - yMin) / (yMax - yMin) * ih; };
    var s = svgEl(W, H);
    for (var g = 0; g <= 4; g++) { var gy = m.t + ih * g / 4; s += '<line x1="' + m.l + '" y1="' + gy + '" x2="' + (W - m.r) + '" y2="' + gy + '" stroke="var(--grid)"/>'; s += txt(m.l - 8, gy + 3, fmt(yMax - (yMax - yMin) * g / 4), { anchor: 'end', size: 9 }); }
    // x ticks at actual timepoints
    uniq(xs).sort(function (a, b) { return a - b; }).forEach(function (xv) { s += txt(sx(xv), H - m.b + 14, String(xv), { size: 9 }); });
    s += txt(m.l + iw / 2, H - 3, metric.x || 'Timepoint', { size: 9, weight: '600' });
    s += '<text x="13" y="' + (m.t + ih / 2) + '" text-anchor="middle" font-size="9" font-weight="600" fill="var(--muted)" transform="rotate(-90 13 ' + (m.t + ih / 2) + ')">' + metric.label + '</text>';
    names.forEach(function (c) {
      var col = condColor[c] || ACCENT, pp = series[c];
      var path = pp.map(function (p, i) { return (i ? 'L' : 'M') + sx(p.x).toFixed(1) + ' ' + sy(p.y).toFixed(1); }).join(' ');
      s += '<path d="' + path + '" fill="none" stroke="' + col + '" stroke-width="2" opacity="0.9"/>';
      pp.forEach(function (p) {
        s += '<circle class="tp" data-id="' + esc(p.row.sample_id) + '" cx="' + sx(p.x).toFixed(1) + '" cy="' + sy(p.y).toFixed(1) + '" r="' + (isHl(p.row.sample_id) ? 6.5 : 4) + '" fill="' + col + '" stroke="var(--card)" stroke-width="1.5" opacity="' + dim(p.row.sample_id) + '" style="cursor:pointer"/>';
      });
    });
    s += '</svg>'; el.innerHTML = s;
    el.querySelectorAll('circle.tp').forEach(function (c) {
      var id = c.getAttribute('data-id'), p = allPts.find(function (q) { return q.row.sample_id === id; });
      c.addEventListener('mousemove', function (e) { ttShow(rowTip(p.row, tipRow(metric.label, fmt(p.y), condColor[p.row.condition || 'Unknown']))); ttMove(e); });
      c.addEventListener('mouseleave', ttHide);
      c.addEventListener('click', function () { ttHide(); setHighlight(id); });
    });
    document.getElementById('trend-legend').innerHTML = names.map(function (c) { return '<span class="item"><span class="sw" style="background:' + (condColor[c] || ACCENT) + '"></span>' + esc(c) + '</span>'; }).join('');
    document.getElementById('trend-note').textContent = metric.note || '';
  }

  // ---------- localized region keys + helpers ----------
  var locKeys = region_keys.filter(function (k) { return !/localized/i.test(k); });
  function locSum(r) { return locKeys.reduce(function (a, k) { return a + (r['func_' + k] || 0); }, 0); }
  function dupVal(r) { return r.umi_deduplication_enabled === true ? (isNum(r.umi_deduplication_percent) ? r.umi_deduplication_percent : null) : (isNum(r.duplicate_percent) ? r.duplicate_percent : null); }
  var hasUMI = rows.some(function (r) { return r.umi_deduplication_enabled === true; });

  var TREND_METRICS = [
    { key: 'pi', label: 'Median pausing index', get: function (r) { return r.median_pausing_index; }, note: 'Higher = more promoter-proximal pausing.' },
    { key: 'div', label: 'Divergent loci', get: function (r) { return r.divergent_regions; }, zero: true },
    { key: 'depth', label: 'Read depth (M)', get: function (r) { return isNum(r.input_reads) ? r.input_reads / 1e6 : null; }, zero: true },
    { key: 'prom', label: 'Promoter share (%)', get: function (r) { var t = locSum(r); return t ? (r['func_Promoter'] || 0) / t * 100 : null; } },
    { key: 'body', label: 'Gene-body share (%)', get: function (r) { var t = locSum(r); return t ? (r['func_Gene body'] || 0) / t * 100 : null; } },
    { key: 'unloc', label: 'Unlocalized (%)', get: function (r) { return isNum(r.unlocalized_fraction) ? r.unlocalized_fraction * 100 : null; }, zero: true }
  ];

  function esc(s) { return String(s).replace(/[&<>"]/g, function (c) { return { '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]; }); }

  // ===================== BUILD =====================
  document.getElementById('hdr-samples').textContent = rows.length + ' samples · ' + conditions.length + ' condition' + (conditions.length !== 1 ? 's' : '');

  // KPI strip
  (function () {
    var depths = rows.map(function (r) { return r.input_reads; }).filter(function (v) { return isNum(v) && v > 0; }).map(function (v) { return v / 1e6; });
    var pis = rows.map(function (r) { return r.median_pausing_index; }).filter(function (v) { return isNum(v) && v > 0; });
    var unlocs = rows.map(function (r) { return r.unlocalized_fraction; }).filter(isNum);
    var totalDiv = rows.reduce(function (a, r) { return a + (r.divergent_regions || 0); }, 0);
    var kpis = [
      { l: 'Samples', v: rows.length, s: conditions.length + ' conditions' },
      { l: 'Median depth', v: depths.length ? median(depths).toFixed(1) + 'M' : '–', s: 'usable reads' },
      { l: 'Median pausing', v: pis.length ? median(pis).toFixed(2) : '–', s: 'PI (len-norm)' },
      { l: 'Total divergent', v: fmtInt(totalDiv), s: 'loci (Σ samples)' },
      { l: 'Median unlocalized', v: unlocs.length ? (median(unlocs) * 100).toFixed(1) + '%' : '–', s: 'off-feature' }
    ];
    document.getElementById('kpi-strip').innerHTML = kpis.map(function (k) { return '<div class="kpi"><div class="k-label">' + k.l + '</div><div class="k-value">' + k.v + '</div><div class="k-sub">' + k.s + '</div></div>'; }).join('');
  })();

  // Design matrix
  (function () {
    var tps = uniq(rows.map(function (r) { return r.timepoint; })).sort(function (a, b) { var na = parseTP(a), nb = parseTP(b); if (na != null && nb != null) return na - nb; return String(a).localeCompare(String(b)); });
    var conds = conditions.slice();
    var head = '<tr><th>Condition \\ tp</th>' + tps.map(function (t) { return '<th>' + esc(t == null || t === '' ? '–' : t) + '</th>'; }).join('') + '</tr>';
    var body = conds.map(function (c) {
      return '<tr><th style="text-align:left">' + esc(c) + '</th>' + tps.map(function (t) {
        var n = rows.filter(function (r) { return (r.condition || 'Unknown') === c && r.timepoint === t; }).length;
        return n ? '<td class="has">' + n + '</td>' : '<td class="empty-cell">·</td>';
      }).join('') + '</tr>';
    }).join('');
    document.getElementById('design-table').innerHTML = head + body;
  })();

  // Condition cards
  (function () {
    var html = conditions.map(function (c) {
      var rs = rows.filter(function (r) { return (r.condition || 'Unknown') === c; });
      var d = rs.map(function (r) { return r.input_reads; }).filter(function (v) { return isNum(v) && v > 0; }).map(function (v) { return v / 1e6; });
      var pi = rs.map(function (r) { return r.median_pausing_index; }).filter(function (v) { return isNum(v) && v > 0; });
      var dvar = cv(d);
      var flag = dvar == null ? '' : dvar < 20 ? '<span class="flag ok">consistent</span>' : dvar < 40 ? '<span class="flag warn">moderate</span>' : '<span class="flag fail">variable</span>';
      return '<div class="cond-card" style="--cond:' + (condColor[c] || ACCENT) + '">' +
        '<div class="cc-name">' + esc(c) + ' <span style="color:var(--faint);font-weight:500">· ' + rs.length + '</span></div>' +
        '<div class="cc-row">Depth ' + (d.length ? median(d).toFixed(1) + 'M' : '–') + ' · PI ' + (pi.length ? median(pi).toFixed(2) : '–') + '</div>' +
        '<div class="cc-row">Replicates ' + (dvar == null ? 'n/a' : 'CV ' + dvar.toFixed(0) + '%') + ' ' + flag + '</div>' +
        '</div>';
    }).join('');
    document.getElementById('condition-cards').innerHTML = html;
  })();

  // Trends
  (function () {
    var numericTPs = uniq(rows.map(function (r) { return parseTP(r.timepoint); }).filter(function (v) { return v != null; }));
    if (numericTPs.length < 2) {
      document.getElementById('trends').style.display = 'none';
      var nl = document.querySelector('#nav a[href="#trends"]'); if (nl) nl.style.display = 'none';
      return;
    }
    var sel = document.getElementById('trend-metric');
    sel.innerHTML = TREND_METRICS.map(function (m) { return '<option value="' + m.key + '">' + m.label + '</option>'; }).join('');
    sel.addEventListener('change', drawTrend);
    register(drawTrend); drawTrend();
  })();

  // QC charts
  drawBars('chart-depth', function () { return orderedRows().map(function (r) { return { id: r.sample_id, value: isNum(r.input_reads) ? r.input_reads / 1e6 : null, row: r }; }); },
    { tipLabel: 'Depth', tipFmt: function (v) { return v.toFixed(1) + 'M'; }, valFmt: function (v) { return v.toFixed(0); }, threshold: { warn: 5 },
      colorFn: function (d) { if (state.colorByCondition) return barColor(d.row); var v = d.value; return v == null ? '#ccc' : v >= 10 ? 'var(--ok)' : v >= 5 ? 'var(--warn)' : 'var(--fail)'; } });
  document.getElementById('dup-title').textContent = hasUMI ? 'UMI-dedup %' : 'Duplication %';
  drawBars('chart-dup', function () { return orderedRows().map(function (r) { return { id: r.sample_id, value: dupVal(r), row: r }; }); },
    { tipLabel: hasUMI ? 'UMI-dedup' : 'Dup', tipFmt: function (v) { return v.toFixed(1) + '%'; }, emptyMsg: 'No deduplication performed' });
  drawBars('chart-mm', function () { return orderedRows().map(function (r) { return { id: r.sample_id, value: isNum(r.multimapper_percent) ? r.multimapper_percent : null, row: r }; }); },
    { tipLabel: 'Multimap', tipFmt: function (v) { return v.toFixed(1) + '%'; } });
  drawBars('chart-unloc', function () { return orderedRows().map(function (r) { return { id: r.sample_id, value: isNum(r.unlocalized_fraction) ? r.unlocalized_fraction * 100 : null, row: r }; }); },
    { tipLabel: 'Unlocalized', tipFmt: function (v) { return v.toFixed(1) + '%'; }, threshold: { warn: 20 } });

  // Divergent
  (function () {
    var dv = rows.map(function (r) { return r.divergent_regions || 0; });
    document.getElementById('div-chips').innerHTML = chip('Total', fmtInt(dv.reduce(function (a, b) { return a + b; }, 0))) + chip('Median/sample', fmtInt(median(dv))) + chip('Range', fmtInt(Math.min.apply(null, dv)) + '–' + fmtInt(Math.max.apply(null, dv)));
    drawBars('chart-div', function () { return orderedRows().map(function (r) { return { id: r.sample_id, value: r.divergent_regions || 0, row: r }; }); }, { tipLabel: 'Loci', tipFmt: fmtInt });
    var byc = {}; rows.forEach(function (r) { var c = r.condition || 'Unknown'; (byc[c] = byc[c] || []).push(r.divergent_regions || 0); });
    var citems = Object.keys(byc).map(function (c) { return { id: 'cond:' + c, label: c, value: mean(byc[c]), row: { sample_id: c, condition: c, timepoint: null } }; });
    drawBars('chart-div-cond', function () { return citems; }, { tipLabel: 'Mean', tipFmt: fmtInt, colorFn: function (d) { return condColor[d.label] || ACCENT; } });
  })();

  // Pausing
  (function () {
    var pis = rows.map(function (r) { return r.median_pausing_index; }).filter(function (v) { return isNum(v) && v > 0; });
    if (pis.length) document.getElementById('pi-chips').innerHTML = chip('Median', median(pis).toFixed(2)) + chip('Mean', mean(pis).toFixed(2)) + chip('Range', Math.min.apply(null, pis).toFixed(2) + '–' + Math.max.apply(null, pis).toFixed(2));
    drawBars('chart-pi', function () { return orderedRows().map(function (r) { return { id: r.sample_id, value: isNum(r.median_pausing_index) && r.median_pausing_index > 0 ? r.median_pausing_index : null, row: r }; }); },
      { tipLabel: 'Median PI', tipFmt: function (v) { return v.toFixed(2); }, valFmt: function (v) { return v.toFixed(1); }, emptyMsg: 'No pausing index data' });
    drawScatter('chart-pi-depth', function () { return rows.filter(function (r) { return isNum(r.median_pausing_index) && isNum(r.input_reads) && r.input_reads > 0; }).map(function (r) { return { id: r.sample_id, x: r.input_reads / 1e6, y: r.median_pausing_index, row: r }; }); },
      'Read depth (M)', 'Median PI', { emptyMsg: 'No PI / depth pairs' });
  })();

  // Functional regions
  (function () {
    if (!locKeys.length) { empty(document.getElementById('chart-composition')); empty(document.getElementById('chart-region-totals')); return; }
    drawStacked('chart-composition', 'composition-legend', locKeys);
    var totItems = locKeys.map(function (k) { return { label: k, value: region_totals[k] || 0 }; }).sort(function (a, b) { return b.value - a.value; });
    drawHBars('chart-region-totals', totItems, { colorFn: function (d) { return regionColor(d.label); } });
    // region table
    var thead = document.querySelector('#region-table thead');
    var tbody = document.querySelector('#region-table tbody');
    thead.innerHTML = '<tr><th class="lft">Sample</th>' + locKeys.map(function (k) { return '<th>' + esc(k) + '</th>'; }).join('') + '</tr>';
    function fillRegion() {
      tbody.innerHTML = orderedRows().map(function (r) {
        return '<tr data-id="' + esc(r.sample_id) + '"' + (isHl(r.sample_id) ? ' class="hl"' : '') + '><td class="lft">' + esc(shortName(r.sample_id)) + '</td>' +
          locKeys.map(function (k) { return '<td class="num">' + fmtInt(r['count_' + k] || 0) + '</td>'; }).join('') + '</tr>';
      }).join('');
      bindRowClicks(tbody);
    }
    register(fillRegion); fillRegion();
  })();

  // Normalization
  (function () {
    var cpm = rows.map(function (r) { return r.cpm_factor; }).filter(isNum);
    var si = rows.map(function (r) { return r.crpmsi_factor; }).filter(isNum);
    document.getElementById('norm-chips').innerHTML =
      chip('CPM samples', cpm.length) + chip('siCPM samples', si.length) +
      chip('CPM range', cpm.length ? cpm.reduce(mn).toFixed(3) + '–' + cpm.reduce(mx).toFixed(3) : '–') +
      chip('siCPM range', si.length ? si.reduce(mn).toFixed(3) + '–' + si.reduce(mx).toFixed(3) : '–');
    drawBars('chart-cpm', function () { return orderedRows().map(function (r) { return { id: r.sample_id, value: isNum(r.cpm_factor) ? r.cpm_factor : null, row: r }; }); }, { tipLabel: 'CPM', tipFmt: function (v) { return v.toFixed(4); }, emptyMsg: 'No CPM factors' });
    drawBars('chart-sicpm', function () { return orderedRows().map(function (r) { return { id: r.sample_id, value: isNum(r.crpmsi_factor) ? r.crpmsi_factor : null, row: r }; }); }, { tipLabel: 'siCPM', tipFmt: function (v) { return v.toFixed(4); }, emptyMsg: 'No spike-in factors' });
    drawScatter('chart-cpm-sicpm', function () { return rows.filter(function (r) { return isNum(r.cpm_factor) && isNum(r.crpmsi_factor); }).map(function (r) { return { id: r.sample_id, x: r.cpm_factor, y: r.crpmsi_factor, row: r }; }); }, 'CPM factor', 'siCPM factor', { emptyMsg: 'No CPM/siCPM pairs' });
  })();
  function mn(a, b) { return Math.min(a, b); } function mx(a, b) { return Math.max(a, b); }
  function chip(l, v) { return '<div class="chip"><div class="c-label">' + l + '</div><div class="c-value">' + v + '</div></div>'; }

  // Sample table
  (function () {
    var COLS = [
      { k: 'sample_id', t: 'Sample', lft: true, fmt: function (r) { var u = '../samples/' + encodeURIComponent(r.sample_id) + '/' + encodeURIComponent(r.sample_id) + '.report.html'; return '<a href="' + u + '" target="_blank">' + esc(shortName(r.sample_id)) + '</a>'; }, val: function (r) { return r.sample_id; } },
      { k: 'condition', t: 'Condition', lft: true, val: function (r) { return r.condition || ''; }, fmt: function (r) { return esc(r.condition || '–'); } },
      { k: 'timepoint', t: 'Time', val: function (r) { var n = parseTP(r.timepoint); return n == null ? r.timepoint : n; }, fmt: function (r) { return esc(r.timepoint != null ? r.timepoint : '–'); } },
      { k: 'input_reads', t: 'Depth', val: function (r) { return r.input_reads || 0; }, fmt: function (r) { return ((r.input_reads || 0) / 1e6).toFixed(1) + 'M'; } },
      { k: 'reads_total_functional', t: 'Func reads', val: function (r) { return r.reads_total_functional || 0; }, fmt: function (r) { return fmt(r.reads_total_functional || 0); } },
      { k: 'dup', t: hasUMI ? 'UMI-dd%' : 'Dup%', val: function (r) { return dupVal(r) == null ? -1 : dupVal(r); }, fmt: function (r) { var v = dupVal(r); return v == null ? '–' : v.toFixed(1) + '%'; } },
      { k: 'unlocalized_fraction', t: 'Unloc%', val: function (r) { return r.unlocalized_fraction || 0; }, fmt: function (r) { return ((r.unlocalized_fraction || 0) * 100).toFixed(1) + '%'; } },
      { k: 'multimapper_percent', t: 'MM%', val: function (r) { return isNum(r.multimapper_percent) ? r.multimapper_percent : -1; }, fmt: function (r) { return isNum(r.multimapper_percent) ? r.multimapper_percent.toFixed(1) + '%' : '–'; } },
      { k: 'divergent_regions', t: 'Div loci', val: function (r) { return r.divergent_regions || 0; }, fmt: function (r) { return fmtInt(r.divergent_regions || 0); } },
      { k: 'median_pausing_index', t: 'Med PI', val: function (r) { return isNum(r.median_pausing_index) ? r.median_pausing_index : -1; }, fmt: function (r) { return isNum(r.median_pausing_index) ? r.median_pausing_index.toFixed(2) : '–'; } },
      { k: 'cpm_factor', t: 'CPM', val: function (r) { return isNum(r.cpm_factor) ? r.cpm_factor : -1; }, fmt: function (r) { return isNum(r.cpm_factor) ? r.cpm_factor.toFixed(4) : '–'; } },
      { k: 'crpmsi_factor', t: 'siCPM', val: function (r) { return isNum(r.crpmsi_factor) ? r.crpmsi_factor : -1; }, fmt: function (r) { return isNum(r.crpmsi_factor) ? r.crpmsi_factor.toFixed(4) : '–'; } }
    ];
    var sortK = null, sortDir = 1;
    var thead = document.getElementById('sample-thead');
    thead.innerHTML = COLS.map(function (c) { return '<th class="sortable' + (c.lft ? ' lft' : '') + '" data-k="' + c.k + '">' + c.t + ' <span class="arrow"></span></th>'; }).join('');
    var tbody = document.getElementById('sample-tbody');
    var search = document.getElementById('sample-search');
    var condFilter = document.getElementById('condition-filter');
    conditions.forEach(function (c) { var o = document.createElement('option'); o.value = c; o.textContent = c; condFilter.appendChild(o); });

    function current() {
      var term = search.value.toLowerCase(), cond = condFilter.value;
      var data = rows.filter(function (r) {
        var hay = (r.sample_id + ' ' + (r.condition || '') + ' ' + (r.timepoint || '')).toLowerCase();
        return hay.indexOf(term) >= 0 && (!cond || r.condition === cond);
      });
      if (sortK) { var col = COLS.find(function (c) { return c.k === sortK; }); data.sort(function (a, b) { var va = col.val(a), vb = col.val(b); if (typeof va === 'number' && typeof vb === 'number') return (va - vb) * sortDir; return String(va).localeCompare(String(vb)) * sortDir; }); }
      else { data = orderedRows().filter(function (r) { return data.indexOf(r) >= 0; }); }
      return data;
    }
    function render() {
      tbody.innerHTML = current().map(function (r) {
        return '<tr data-id="' + esc(r.sample_id) + '"' + (isHl(r.sample_id) ? ' class="hl"' : '') + '>' + COLS.map(function (c) { return '<td class="' + (c.lft ? 'lft' : 'num') + '">' + c.fmt(r) + '</td>'; }).join('') + '</tr>';
      }).join('');
      bindRowClicks(tbody);
      thead.querySelectorAll('th').forEach(function (th) { var a = th.querySelector('.arrow'); a.textContent = th.getAttribute('data-k') === sortK ? (sortDir > 0 ? '▲' : '▼') : ''; });
    }
    thead.querySelectorAll('th').forEach(function (th) { th.addEventListener('click', function () { var k = th.getAttribute('data-k'); if (sortK === k) sortDir *= -1; else { sortK = k; sortDir = 1; } render(); }); });
    search.addEventListener('input', render);
    condFilter.addEventListener('change', render);
    register(render); render();

    document.getElementById('export-csv').addEventListener('click', function () {
      var head = COLS.map(function (c) { return c.t; });
      var lines = [head.join(',')].concat(current().map(function (r) { return COLS.map(function (c) { var v = c.val(r); return typeof v === 'string' && v.indexOf(',') >= 0 ? '"' + v + '"' : v; }).join(','); }));
      var blob = new Blob([lines.join('\n')], { type: 'text/csv' });
      var a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = 'tracktx_cohort_metrics.csv'; a.click();
    });
  })();

  function bindRowClicks(tbody) {
    tbody.querySelectorAll('tr[data-id]').forEach(function (tr) { tr.addEventListener('click', function (e) { if (e.target.tagName === 'A') return; setHighlight(tr.getAttribute('data-id')); }); });
  }
  function syncTableHl() {
    document.querySelectorAll('tr[data-id]').forEach(function (tr) { tr.classList.toggle('hl', tr.getAttribute('data-id') === state.highlight); });
  }

  // Methodology
  document.getElementById('methodology').innerHTML = [
    ['Read depth', 'total_reads_raw (qc_pol.json)', 'samtools view -c -F 0x900 on the aligned BAM; falls back to dedup or functional reads if QC parsing failed.'],
    ['Duplication', '100 × duplicate / total reads', 'From samtools stats (flag 0x400). UMI deduplication uses reads_before/after when enabled.'],
    ['Multimapper %', '100 × (1 − unique / mapped)', 'Uniqueness = NH==1 when bowtie2 -k multimapping is active, else MAPQ≥threshold. High values flag repetitive genomes.'],
    ['Divergent loci', 'count from divergent BED', 'Gaussian Mixture Model + FDR control; one value per sample.'],
    ['Pausing index', '(TSS density) / (gene-body density)', 'pi_len_norm = (tss_count/tss_width) / (body_count/body_len); median over genes. >1.5 strong, ≈1 balanced, <1 elongation.'],
    ['Functional signal', 'Σ |pos|+|neg| reads per region', 'Promoter, Gene body, CPS, Enhancers, Termination window, DivergentTx from normalized bedGraphs.'],
    ['Unlocalized fraction', 'unlocalized / (functional + unlocalized)', 'Reads outside annotated features; <20% typical.'],
    ['CPM / siCPM', 'CPM = 1e6/total_mapped; siCPM = CPM × spike_ctrl/spike_sample', 'Convert raw counts to normalized signal in tracks and density tables.']
  ].map(function (m) { return '<div class="method-row"><div class="m-name">' + m[0] + '</div><div class="m-formula">' + esc(m[1]) + '</div><div class="m-desc">' + esc(m[2]) + '</div></div>'; }).join('');

  // controls
  document.getElementById('order-seg').addEventListener('click', function (e) {
    var b = e.target.closest('button'); if (!b) return;
    state.order = b.getAttribute('data-order');
    this.querySelectorAll('button').forEach(function (x) { x.classList.toggle('on', x === b); });
    redrawAll();
  });
  document.getElementById('color-by-cond').addEventListener('change', function () { state.colorByCondition = this.checked; redrawAll(); });
  document.getElementById('hl-clear').addEventListener('click', function () { state.highlight = null; document.getElementById('hl-banner').style.display = 'none'; redrawAll(); syncTableHl(); });

  // nav scroll-spy
  (function () {
    var links = Array.prototype.slice.call(document.querySelectorAll('#nav a'));
    var secs = links.map(function (a) { return document.getElementById(a.getAttribute('href').slice(1)); });
    function spy() {
      var y = window.scrollY + 90, idx = 0;
      secs.forEach(function (s, i) { if (s && s.offsetTop <= y) idx = i; });
      links.forEach(function (a, i) { a.classList.toggle('active', i === idx); });
    }
    window.addEventListener('scroll', spy, { passive: true }); spy();
  })();

  // responsive redraw
  var rt; window.addEventListener('resize', function () { clearTimeout(rt); rt = setTimeout(redrawAll, 180); });

  console.log('TrackTx cohort report ready:', rows.length, 'samples');
})();
</script>
'''
    
    generate_html_report(
        data_json,
        CSS,
        JS,
        args,
        run_command,
        args.out_html
    )
    
    # Success
    log("═" * 70, "")
    log("COMPLETE", f"Cohort report generated successfully")
    log("COMPLETE", f"Samples: {len(samples)}")
    log("COMPLETE", f"Regions: {len(region_keys)}")
    log("COMPLETE", f"Outputs: 3-4 files")
    
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
