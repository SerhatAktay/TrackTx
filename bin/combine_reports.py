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
        except Exception:
            continue
    
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
        "input_reads": qc.get("total_reads_raw") or json_data.get("input_reads"),
        "dedup_reads": qc.get("dedup_reads_mapq_ge") or json_data.get("dedup_reads"),
        "duplicate_percent": dup_percent,
        "umi_deduplication_enabled": qc.get("umi_deduplication_enabled", False),
        "umi_deduplication_percent": qc.get("umi_deduplication_percent"),
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
        "input_reads", "dedup_reads", "duplicate_percent",
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
        "input_reads", "dedup_reads", "duplicate_percent",
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
        "skipped_inputs": skipped_files
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
    """
    Generate single-file HTML report with embedded CSS/JS
    
    Args:
        data_json: JSON data string
        css: CSS string
        js: JavaScript string
        args: Command-line arguments
        run_command: Full command line
        output_path: Output HTML path
    """
    log("HTML", "Generating interactive dashboard...")
    
    # Generate timestamp
    timestamp = datetime.datetime.now().strftime("%B %d, %Y at %H:%M")
    
    # Build comprehensive HTML report
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width,initial-scale=1">
  <title>TrackTx Cohort Report</title>
  {css}
</head>
<body>
<div class="container">
  <!-- Header -->
  <header class="page-header">
    <div class="eyebrow">TrackTx PRO-seq Analysis • Cohort Report</div>
    <h1>Global Summary</h1>
    <p class="muted">Profile: {args.profile} • Run: {args.run_name} • Duration: {args.duration}</p>
    <p class="muted">Generated {timestamp} • <span id="sample-count-header"></span></p>
  </header>

  <!-- Navigation -->
  <nav class="nav-pills">
    <a href="#overview">📊 Overview</a>
    <a href="#qc">✓ Quality Control</a>
    <a href="#divergent">🔀 Divergent TX</a>
    <a href="#pausing">⏸️ Pausing</a>
    <a href="#regions">🎯 Functional Regions</a>
    <a href="#normalization">📏 Normalization</a>
    <a href="#samples">📋 Sample Table</a>
    <a href="#files">📁 Files</a>
  </nav>

  <!-- SECTION: Overview -->
  <section class="layer" id="overview">
    <div class="layer-heading">
      <h2>Cohort Overview</h2>
      <div class="subtitle">High-level summary and experimental design</div>
    </div>
    
    <div class="help-box">
      <strong>What is this report?</strong> This cohort-level dashboard aggregates metrics from all individual sample reports in your TrackTx PRO-seq analysis.
      Use it to assess overall experiment quality, identify outliers, compare conditions, and understand transcriptional dynamics across your cohort.
      <ul>
        <li><strong>QC metrics:</strong> Read depth, duplication rates, mapping stats</li>
        <li><strong>Biological metrics:</strong> Divergent transcription, Pol II pausing, functional region distributions</li>
        <li><strong>Normalization:</strong> CPM and spike-in factors for cross-sample comparisons</li>
      </ul>
    </div>

    <div class="kpi-grid">
      <div class="kpi-card">
        <div class="kpi-label">Total Samples <span class="info-tip" title="Total number of samples successfully processed and included in this cohort report">?</span></div>
        <div class="kpi-value" id="kpi-total">0</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Conditions <span class="info-tip" title="Number of unique experimental conditions in the cohort">?</span></div>
        <div class="kpi-value" id="kpi-conditions">0</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Avg Read Depth <span class="info-tip" title="Average total input reads per sample (millions)">?</span></div>
        <div class="kpi-value" id="kpi-depth">0M</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Total Divergent Loci <span class="info-tip" title="Sum of divergent transcription regions detected across all samples">?</span></div>
        <div class="kpi-value" id="kpi-div-total">0</div>
      </div>
      <div class="kpi-card">
        <div class="kpi-label">Avg Functional Regions <span class="info-tip" title="Average number of functional genomic regions (genes, promoters, etc.) with signal per sample">?</span></div>
        <div class="kpi-value" id="kpi-regions-avg">0</div>
      </div>
    </div>

    <div id="experimental-design-summary" style="margin-top:2rem;"></div>
  </section>

  <!-- SECTION: Quality Control -->
  <section class="layer" id="qc">
    <div class="layer-heading">
      <h2>Quality Control Analysis</h2>
      <div class="subtitle">Sequencing depth, duplication rates, and sample consistency</div>
    </div>
    
    <div class="help-box">
      <strong>What to look for:</strong>
      <ul>
        <li><strong>Read Depth:</strong> PRO-seq typically requires 5-20M reads per sample for good gene coverage. Lower depth may miss lowly-expressed genes.</li>
        <li><strong>Duplication Rate:</strong> <15% is excellent, 15-30% is acceptable, >30% suggests PCR over-amplification or low library complexity.</li>
        <li><strong>Unlocalized Reads:</strong> <20% is typical. Higher values may indicate contamination, rRNA, or incomplete genome annotation.</li>
        <li><strong>Sample Consistency:</strong> Replicates within a condition should cluster together. Large variability suggests technical issues or biological heterogeneity.</li>
      </ul>
    </div>

    <div class="stats-grid">
      <div class="stat-item">
        <div class="stat-label">Median Depth</div>
        <div class="stat-value" id="qc-median-depth">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">Median Duplication</div>
        <div class="stat-value" id="qc-median-dup">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">Median Unlocalized</div>
        <div class="stat-value" id="qc-median-unloc">-</div>
      </div>
    </div>

    <div class="viz-grid">
      <div class="viz-card">
        <h3>Read Depth per Sample <span class="info-tip" title="Total input reads for each sample. HOVER over bars to see sample names and exact values.">?</span></h3>
        <div id="chart-depth-per-sample" class="viz-body"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">💡 Hover over bars to see sample IDs</p>
      </div>
      <div class="viz-card">
        <h3>Duplication Rate per Sample <span class="info-tip" title="PCR duplication rate for each sample. HOVER over bars to see sample names. Lower is better.">?</span></h3>
        <div id="chart-dup-per-sample" class="viz-body"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">💡 Hover over bars to see sample IDs</p>
      </div>
      <div class="viz-card">
        <h3>Read Depth Distribution <span class="info-tip" title="Histogram showing how read depths are distributed. HOVER over bars to see which samples fall in each bin.">?</span></h3>
        <div id="chart-depth-dist" class="viz-body"></div>
        <div class="distrib-summary" id="depth-stats"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">📊 Aggregate distribution - hover to see samples in each bin</p>
      </div>
      <div class="viz-card">
        <h3>Duplication Rate Distribution <span class="info-tip" title="Histogram showing how duplication rates are distributed. HOVER over bars to see which samples fall in each bin.">?</span></h3>
        <div id="chart-dup-dist" class="viz-body"></div>
        <div class="distrib-summary" id="dup-stats"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">📊 Aggregate distribution - hover to see samples in each bin</p>
      </div>
    </div>

  </section>

  <!-- SECTION: Divergent Transcription -->
  <section class="layer" id="divergent">
    <div class="layer-heading">
      <h2>Divergent Transcription Analysis</h2>
      <div class="subtitle">Bidirectional transcription from promoters and enhancers</div>
    </div>
    
    <div class="help-box">
      <strong>What is divergent transcription?</strong> Divergent (bidirectional) transcription occurs when RNA Polymerase II initiates in both directions from a promoter or enhancer.
      It's a hallmark of active regulatory elements. TrackTx uses statistical detection (Gaussian Mixture Models + FDR control) to identify high-confidence divergent regions.
      <ul>
        <li><strong>More regions:</strong> Indicates higher transcriptional activity or more active enhancers</li>
        <li><strong>Condition differences:</strong> Can reflect stimulus-dependent enhancer activation</li>
        <li><strong>Variability:</strong> High variability within replicates may suggest technical noise or biological heterogeneity</li>
      </ul>
    </div>

    <div class="stats-grid">
      <div class="stat-item">
        <div class="stat-label">Total Regions</div>
        <div class="stat-value" id="div-total-regions">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">Mean per Sample</div>
        <div class="stat-value" id="div-mean">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">Median per Sample</div>
        <div class="stat-value" id="div-median">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">Range</div>
        <div class="stat-value" id="div-range">-</div>
      </div>
    </div>

    <div class="viz-grid">
      <div class="viz-card">
        <h3>Divergent Regions per Sample <span class="info-tip" title="Number of high-confidence divergent transcription sites detected in each sample">?</span></h3>
        <div id="chart-div-per-sample" class="viz-body"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">💡 Hover over bars to see sample IDs</p>
      </div>
      <div class="viz-card">
        <h3>Distribution Across Cohort <span class="info-tip" title="Histogram showing how divergent region counts are distributed. HOVER over bars to see which samples fall in each bin.">?</span></h3>
        <div id="chart-div-hist" class="viz-body"></div>
        <div class="distrib-summary" id="div-stats"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">📊 Aggregate distribution - hover to see samples in each bin</p>
      </div>
      <div class="viz-card">
        <h3>By Condition <span class="info-tip" title="Compare divergent transcription levels across experimental conditions">?</span></h3>
        <div id="chart-div-by-condition" class="viz-body"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">💡 Hover over bars to see condition names</p>
      </div>
      <div class="viz-card">
        <h3>Replicate Consistency <span class="info-tip" title="Coefficient of variation within replicate groups. Lower is better (more consistent).">?</span></h3>
        <div id="chart-div-cv" class="viz-body"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">💡 Hover over bars to see condition names</p>
      </div>
    </div>
  </section>

  <!-- SECTION: Pausing Index -->
  <section class="layer" id="pausing">
    <div class="layer-heading">
      <h2>Pol II Pausing Analysis</h2>
      <div class="subtitle">Promoter-proximal pausing and elongation dynamics</div>
    </div>
    
    <div class="help-box">
      <strong>What is Pol II pausing?</strong> After transcription initiation, RNA Polymerase II often pauses 20-60 bp downstream of the TSS before entering productive elongation.
      The <strong>Pausing Index (PI)</strong> quantifies this: PI = (Promoter Signal) / (Gene Body Signal). 
      <ul>
        <li><strong>PI &gt; 2:</strong> Strong pausing - typical for rapidly-induced genes (e.g., heat shock genes)</li>
        <li><strong>PI = 0.5-2:</strong> Moderate pausing - most constitutive genes</li>
        <li><strong>PI &lt; 0.5:</strong> Weak/no pausing - highly elongating genes</li>
        <li><strong>Condition changes:</strong> PI shifts can indicate transcriptional regulation at elongation step vs. initiation</li>
      </ul>
    </div>

    <div class="stats-grid">
      <div class="stat-item">
        <div class="stat-label">Cohort Median PI</div>
        <div class="stat-value" id="pi-cohort-median">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">Mean PI</div>
        <div class="stat-value" id="pi-cohort-mean">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">Std Dev</div>
        <div class="stat-value" id="pi-std">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">Range</div>
        <div class="stat-value" id="pi-range">-</div>
      </div>
    </div>

    <div class="viz-grid">
      <div class="viz-card">
        <h3>Median PI per Sample <span class="info-tip" title="Median pausing index for each sample. Shows overall pausing landscape.">?</span></h3>
        <div id="chart-pi-per-sample" class="viz-body"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">💡 Hover over bars to see sample IDs</p>
      </div>
      <div class="viz-card">
        <h3>PI Distribution <span class="info-tip" title="Histogram showing how pausing indices are distributed. HOVER over bars to see which samples fall in each bin.">?</span></h3>
        <div id="chart-pi-dist" class="viz-body"></div>
        <div class="distrib-summary" id="pi-stats"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">📊 Aggregate distribution - hover to see samples in each bin</p>
      </div>
      <div class="viz-card">
        <h3>PI by Condition <span class="info-tip" title="Compare pausing indices across experimental conditions. Differences may reflect regulation at elongation stage.">?</span></h3>
        <div id="chart-pi-by-condition" class="viz-body"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">💡 Hover over bars to see condition names</p>
      </div>
      <div class="viz-card">
        <h3>PI vs Read Depth <span class="info-tip" title="Pausing index should be independent of sequencing depth. Strong correlation suggests technical bias.">?</span></h3>
        <div id="chart-pi-vs-depth" class="viz-body"></div>
        <p style="text-align:center;font-size:0.85rem;color:var(--muted);margin-top:0.5rem;">💡 Hover over points to see sample IDs</p>
      </div>
    </div>
  </section>

  <!-- SECTION: Functional Regions -->
  <section class="layer" id="regions">
    <div class="layer-heading">
      <h2>Functional Region Composition</h2>
      <div class="subtitle">Read distribution across genomic features</div>
    </div>
    
    <div class="help-box">
      <strong>Functional regions</strong> are genomic features where Pol II signal is measured:
      <ul>
        <li><strong>Promoters:</strong> TSS ± window (default: -500 to +250 bp). High signal indicates active transcription initiation.</li>
        <li><strong>Gene Bodies:</strong> From TSS + offset to TES (default: start +2000 bp). Measures elongating Pol II.</li>
        <li><strong>CPS (Cleavage/PolyA Sites):</strong> 3' end regions. Signal here indicates termination events.</li>
        <li><strong>Enhancers:</strong> Distal regulatory elements (if annotated). Active enhancers produce eRNAs.</li>
        <li><strong>Divergent TX:</strong> Bidirectional transcription sites detected by pipeline.</li>
        <li><strong>Non-localized:</strong> Reads not mapping to any defined feature. Should be <20%.</li>
      </ul>
      <strong>What to look for:</strong> Samples should have similar functional region distributions. Large differences may indicate varying library quality or biological states.
    </div>

    <div class="viz-grid">
      <div class="viz-card">
        <h3>Cohort-wide Region Totals <span class="info-tip" title="Aggregate read counts across all samples for each functional region type">?</span></h3>
        <div id="chart-region-totals" class="viz-body"></div>
      </div>
      <div class="viz-card">
        <h3>Region Composition by Sample <span class="info-tip" title="Stacked bar showing percentage distribution of reads across regions for each sample">?</span></h3>
        <div id="chart-region-composition" class="viz-body"></div>
      </div>
      <div class="viz-card">
        <h3>Promoter Signal per Sample <span class="info-tip" title="Promoter signal for each sample. Hover to see sample names and exact counts.">?</span></h3>
        <div id="chart-promoter-per-sample" class="viz-body"></div>
      </div>
      <div class="viz-card">
        <h3>Gene Body Signal per Sample <span class="info-tip" title="Gene body signal for each sample. Hover to see sample names and exact counts.">?</span></h3>
        <div id="chart-genebody-per-sample" class="viz-body"></div>
      </div>
    </div>

    <div style="margin-top:2rem;">
      <h3 style="margin-bottom:1rem;font-size:1.25rem;">Detailed Region Counts by Sample</h3>
      <div class="table-wrap">
        <table id="region-counts-table">
          <thead id="region-counts-thead"></thead>
          <tbody id="region-counts-tbody"></tbody>
        </table>
      </div>
    </div>
  </section>

  <!-- SECTION: Normalization -->
  <section class="layer" id="normalization">
    <div class="layer-heading">
      <h2>Normalization Factors</h2>
      <div class="subtitle">CPM and spike-in normalization for cross-sample comparisons</div>
    </div>
    
    <div class="help-box">
      <strong>Why normalize?</strong> Raw read counts vary with sequencing depth and library preparation efficiency. Normalization enables quantitative comparison across samples.
      <ul>
        <li><strong>CPM (Counts Per Million):</strong> Simple depth normalization. Assumes total transcription is similar across samples.
          <br>Formula: CPM = (raw_counts / total_reads) × 1,000,000</li>
        <li><strong>siCPM (Spike-in CPM):</strong> Uses exogenous spike-in control (e.g., Drosophila) for absolute quantification. Accounts for global changes in transcription.
          <br>Formula: siCPM = CPM × (spike-in_reads_control / spike-in_reads_sample)</li>
        <li><strong>When to use spike-in:</strong> Essential when comparing samples with expected global transcriptional changes (e.g., stress conditions, differentiation)</li>
        <li><strong>Quality check:</strong> Spike-in factors should be consistent across replicates within a condition. Large variation suggests pipetting errors or contamination.</li>
      </ul>
    </div>

    <div class="stats-grid">
      <div class="stat-item">
        <div class="stat-label">Samples with CPM</div>
        <div class="stat-value" id="norm-cpm-count">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">Samples with siCPM</div>
        <div class="stat-value" id="norm-sicpm-count">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">CPM Factor Range</div>
        <div class="stat-value" id="norm-cpm-range">-</div>
      </div>
      <div class="stat-item">
        <div class="stat-label">siCPM Factor Range</div>
        <div class="stat-value" id="norm-sicpm-range">-</div>
      </div>
    </div>

    <div class="viz-grid">
      <div class="viz-card">
        <h3>CPM Factors <span class="info-tip" title="CPM normalization factors for each sample. Should be inversely proportional to sequencing depth.">?</span></h3>
        <div id="chart-cpm-factors" class="viz-body"></div>
      </div>
      <div class="viz-card">
        <h3>siCPM Factors <span class="info-tip" title="Spike-in normalization factors. Consistent within condition indicates good technical reproducibility.">?</span></h3>
        <div id="chart-sicpm-factors" class="viz-body"></div>
      </div>
      <div class="viz-card">
        <h3>CPM vs siCPM Comparison <span class="info-tip" title="Compare CPM and siCPM factors. Large differences suggest global transcriptional changes between conditions.">?</span></h3>
        <div id="chart-cpm-vs-sicpm" class="viz-body"></div>
      </div>
      <div class="viz-card">
        <h3>Normalization Factor CV by Condition <span class="info-tip" title="Coefficient of variation for normalization factors within each condition. Lower = more consistent.">?</span></h3>
        <div id="chart-norm-cv" class="viz-body"></div>
      </div>
    </div>
  </section>

  <!-- SECTION: Sample Table -->
  <section class="layer" id="samples">
    <div class="layer-heading">
      <h2>Sample-Level Metrics</h2>
      <div class="subtitle">Interactive table with all metrics, filters, and export</div>
    </div>

    <div class="table-tools">
      <div>
        <label>Search</label>
        <input id="sample-search" type="text" placeholder="Search samples, conditions, timepoints...">
      </div>
      <div>
        <label>Condition</label>
        <select id="condition-filter"><option value="">All</option></select>
      </div>
      <div>
        <label>Actions</label>
        <button id="export-csv">Export CSV</button>
      </div>
    </div>
    
    <div class="table-wrap">
      <table id="sample-table">
        <thead>
          <tr>
            <th>Sample</th>
            <th>Condition</th>
            <th>Time</th>
            <th>Rep</th>
            <th>Input Reads</th>
            <th>Functional Reads</th>
            <th>Dup %</th>
            <th>Unloc %</th>
            <th>Div Regions</th>
            <th>Total Regions</th>
            <th>Median PI</th>
            <th>Median Density</th>
            <th>CPM Factor</th>
            <th>siCPM Factor</th>
          </tr>
        </thead>
        <tbody id="sample-table-body"></tbody>
      </table>
    </div>
  </section>

  <!-- SECTION: Files -->
  <section class="layer" id="files">
    <div class="layer-heading">
      <h2>Output Files & Documentation</h2>
      <div class="subtitle">Locate results and reproduce analysis</div>
    </div>

    <div class="help-box">
      <strong>Output Organization:</strong> TrackTx organizes results hierarchically by analysis step and sample:
      <ul>
        <li><strong>11_reports/:</strong> This cohort report and individual sample HTML reports</li>
        <li><strong>05_normalized_tracks/:</strong> BigWig files for genome browser visualization (CPM and siCPM normalized)</li>
        <li><strong>06_divergent_tx/:</strong> Divergent transcription BED files and QC reports</li>
        <li><strong>07_functional_regions/:</strong> Region assignments and read count summaries</li>
        <li><strong>08_pol2_metrics/:</strong> Pausing indices and density calculations per gene/region</li>
        <li><strong>10_qc/:</strong> Quality control JSON files with alignment stats</li>
      </ul>
    </div>

    <div style="background:var(--card);padding:1.5rem;border-radius:0.75rem;border:1px solid var(--line);">
      <h3 style="margin:0 0 1rem 0;font-size:1.125rem;">Run Command</h3>
      <pre style="background:var(--bg);padding:1rem;border-radius:0.5rem;overflow-x:auto;"><code>{run_command}</code></pre>
    </div>

    <div style="background:var(--card);padding:1.5rem;border-radius:0.75rem;border:1px solid var(--line);margin-top:1.5rem;">
      <h3 style="margin:0 0 1rem 0;font-size:1.125rem;">File Locations</h3>
      <ul style="line-height:1.8;">
        <li>Cohort HTML: <code>{args.out_html}</code></li>
        <li>Cohort TSV: <code>{args.out_tsv}</code></li>
        <li>Cohort JSON: <code>{args.out_json}</code></li>
        <li>Sample reports: <code>{{output_dir}}/11_reports/samples/&lt;sample&gt;/&lt;sample&gt;.report.html</code></li>
        <li>Normalized tracks: <code>{{output_dir}}/05_normalized_tracks/&lt;sample&gt;/*.bw</code></li>
        <li>Divergent TX: <code>{{output_dir}}/06_divergent_tx/&lt;sample&gt;/divergent_transcription.bed</code></li>
        <li>Functional regions: <code>{{output_dir}}/07_functional_regions/&lt;sample&gt;/functional_regions.bed</code></li>
        <li>Pol II metrics: <code>{{output_dir}}/08_pol2_metrics/&lt;sample&gt;/*.tsv</code></li>
        <li>QC stats: <code>{{output_dir}}/10_qc/&lt;sample&gt;/qc_summary.json</code></li>
      </ul>
    </div>
  </section>

</div>

<script type="application/json" id="payload">{data_json}</script>
{js}
</body>
</html>
"""
    
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
        
        fallback_name = Path(filepath).stem.split(".")[0]
        sample = normalize_sample_data(json_data, fallback_name)
        
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
    
    # Write outputs
    log("═" * 70, "")
    write_tsv_output(df, args.out_tsv)
    
    write_json_output(
        samples,
        rows,
        region_totals,
        region_keys,
        skipped_files,
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
        "region_keys": region_keys
    }), ensure_ascii=False)
    
    # Embedded Assets
    CSS = """
    <style>
    :root {
      --bg: #ffffff; --fg: #111827; --muted: #6b7280;
      --card: #f9fafb; --line: #e5e7eb; --primary: #3b82f6;
      --ok: #10b981; --warn: #f59e0b; --fail: #ef4444;
      --font: system-ui, -apple-system, sans-serif;
      --accent: #8b5cf6; --accent-light: #c4b5fd;
      --info: #06b6d4; --danger: #ef4444;
    }
    @media (prefers-color-scheme: dark) {
      :root {
        --bg: #0f172a; --fg: #f1f5f9; --muted: #94a3b8;
        --card: #1e293b; --line: #334155; --primary: #60a5fa;
        --accent: #a78bfa; --accent-light: #6d28d9;
        --info: #22d3ee; --danger: #f87171;
      }
    }
    
    * { box-sizing: border-box; }
    body { background: var(--bg); color: var(--fg); font-family: var(--font); margin: 0; line-height: 1.6; }
    .container { max-width: 1400px; margin: 0 auto; padding: 2.5rem; }
    
    /* Header */
    .page-header { 
      background: linear-gradient(135deg, var(--primary) 0%, var(--accent) 100%);
      color: white; padding: 3rem 2rem; border-radius: 1rem; margin-bottom: 3rem;
      box-shadow: 0 10px 40px rgba(0,0,0,0.1);
    }
    .eyebrow { text-transform: uppercase; font-size: 0.75rem; font-weight: 700; opacity: 0.9; letter-spacing: 0.1em; }
    h1 { margin: 0.5rem 0; font-size: 2.5rem; font-weight: 800; }
    .page-header .muted { color: rgba(255,255,255,0.9); font-size: 1rem; }
    
    /* Navigation */
    .nav-pills { 
      display: flex; gap: 0.5rem; flex-wrap: wrap; margin-bottom: 2rem;
      position: sticky; top: 0; background: var(--bg); padding: 1rem 0; z-index: 100;
      border-bottom: 2px solid var(--line);
    }
    .nav-pills a {
      padding: 0.625rem 1.25rem; border-radius: 0.5rem; text-decoration: none;
      color: var(--fg); background: var(--card); border: 1px solid var(--line);
      font-weight: 600; font-size: 0.875rem; transition: all 0.2s;
    }
    .nav-pills a:hover { background: var(--primary); color: white; border-color: var(--primary); }
    
    /* Sections */
    .layer { margin-bottom: 4rem; scroll-margin-top: 80px; }
    .layer-heading { margin-bottom: 1.5rem; }
    .layer-heading h2 { 
      margin: 0; font-size: 1.75rem; font-weight: 700;
      background: linear-gradient(135deg, var(--primary), var(--accent));
      -webkit-background-clip: text; -webkit-text-fill-color: transparent;
      background-clip: text;
    }
    .layer-heading .subtitle { color: var(--muted); font-size: 1rem; margin-top: 0.5rem; }
    
    /* Help boxes */
    .help-box {
      background: linear-gradient(135deg, rgba(59,130,246,0.1), rgba(139,92,246,0.1));
      border-left: 4px solid var(--primary); padding: 1.25rem; border-radius: 0.5rem;
      margin: 1.5rem 0; font-size: 0.95rem; line-height: 1.7;
    }
    .help-box strong { color: var(--primary); }
    .help-box ul { margin: 0.75rem 0 0 1.5rem; }
    .help-box li { margin: 0.5rem 0; }
    
    /* KPI Cards */
    .kpi-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1.25rem; }
    .kpi-card { 
      background: var(--card); padding: 1.5rem; border-radius: 0.75rem;
      border: 1px solid var(--line); transition: all 0.3s;
      box-shadow: 0 2px 8px rgba(0,0,0,0.05);
    }
    .kpi-card:hover { box-shadow: 0 8px 24px rgba(0,0,0,0.12); transform: translateY(-2px); }
    .kpi-label { font-size: 0.875rem; color: var(--muted); margin-bottom: 0.75rem; font-weight: 600; }
    .kpi-value { font-size: 2rem; font-weight: 800; color: var(--fg); }
    .kpi-trend { font-size: 0.8rem; color: var(--muted); margin-top: 0.5rem; }
    .kpi-sparkline { height: 40px; margin-top: 0.75rem; }
    
    /* Visualization Cards */
    .viz-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(450px, 1fr)); gap: 1.5rem; }
    .viz-card { 
      background: var(--card); padding: 1.5rem; border-radius: 0.75rem;
      border: 1px solid var(--line); box-shadow: 0 2px 8px rgba(0,0,0,0.05);
    }
    .viz-card h3 { margin: 0 0 1rem 0; font-size: 1.125rem; font-weight: 700; }
    .viz-body { min-height: 380px; height: 380px; position: relative; width: 100%; }
    
    /* Stats Grid */
    .stats-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin: 1.5rem 0; }
    .stat-item { 
      background: var(--card); padding: 1.25rem; border-radius: 0.5rem;
      border-left: 4px solid var(--primary);
    }
    .stat-label { font-size: 0.8rem; color: var(--muted); font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; }
    .stat-value { font-size: 1.5rem; font-weight: 700; margin-top: 0.25rem; }
    
    /* Tables */
    .table-tools { 
      display: flex; gap: 1rem; margin-bottom: 1.5rem; flex-wrap: wrap; align-items: flex-end;
      background: var(--card); padding: 1.25rem; border-radius: 0.75rem; border: 1px solid var(--line);
    }
    .table-tools > div { flex: 1; min-width: 200px; }
    .table-tools label { display: block; font-size: 0.75rem; font-weight: 700; margin-bottom: 0.5rem; text-transform: uppercase; letter-spacing: 0.05em; }
    input, select, button { 
      width: 100%; padding: 0.625rem 0.875rem; border: 1px solid var(--line); 
      border-radius: 0.5rem; background: var(--bg); color: var(--fg); font-size: 0.875rem;
      transition: all 0.2s;
    }
    input:focus, select:focus { outline: none; border-color: var(--primary); box-shadow: 0 0 0 3px rgba(59,130,246,0.1); }
    button { 
      cursor: pointer; background: var(--primary); color: white; border: none; font-weight: 600;
      box-shadow: 0 2px 8px rgba(59,130,246,0.3);
    }
    button:hover { background: var(--accent); box-shadow: 0 4px 12px rgba(139,92,246,0.4); }
    
    .table-wrap { 
      overflow-x: auto; border: 1px solid var(--line); border-radius: 0.75rem;
      box-shadow: 0 2px 8px rgba(0,0,0,0.05);
    }
    table { width: 100%; border-collapse: collapse; font-size: 0.875rem; }
    th, td { padding: 0.875rem 1.125rem; text-align: left; border-bottom: 1px solid var(--line); }
    th { 
      background: var(--card); font-weight: 700; white-space: nowrap; position: sticky; top: 0;
      text-transform: uppercase; font-size: 0.75rem; letter-spacing: 0.05em; color: var(--muted);
    }
    tbody tr { transition: background 0.2s; }
    tbody tr:hover { background: rgba(59,130,246,0.05); }
    tr:last-child td { border-bottom: none; }
    
    /* Status badges */
    .status { 
      display: inline-flex; align-items: center; gap: 0.5rem; font-weight: 600;
      font-size: 0.75rem; padding: 0.25rem 0.75rem; border-radius: 999px;
      background: var(--card); border: 1px solid var(--line);
    }
    .status::before { content: ''; width: 0.625rem; height: 0.625rem; border-radius: 50%; background: var(--status); }
    
    /* Info tooltips */
    .info-tip { 
      cursor: help; color: var(--primary); border-bottom: 1px dotted var(--primary);
      margin-left: 0.25rem; font-weight: 600; position: relative;
      display: inline-block;
    }
    .info-tip:hover { color: var(--accent); }
    .info-tip:hover::after {
      content: attr(title);
      position: absolute;
      bottom: 100%;
      left: 50%;
      transform: translateX(-50%);
      background: rgba(17, 24, 39, 0.95);
      color: white;
      padding: 0.5rem 0.75rem;
      border-radius: 0.375rem;
      font-size: 0.8125rem;
      white-space: normal;
      width: max-content;
      max-width: 300px;
      z-index: 1000;
      box-shadow: 0 4px 12px rgba(0,0,0,0.3);
      margin-bottom: 0.5rem;
      pointer-events: none;
      line-height: 1.4;
    }
    @media (prefers-color-scheme: dark) {
      .info-tip:hover::after {
        background: rgba(241, 245, 249, 0.95);
        color: #0f172a;
      }
    }
    
    /* Charts */
    .bar-chart { display: flex; align-items: flex-end; height: 100%; gap: 4px; padding: 1rem 0; }
    .bar-col { 
      flex: 1; display: flex; flex-direction: column; justify-content: flex-end;
      height: 100%; position: relative; cursor: pointer;
    }
    .bar { 
      background: linear-gradient(180deg, var(--primary), var(--accent)); 
      width: 100%; border-radius: 4px 4px 0 0; transition: all 0.3s;
      box-shadow: 0 2px 8px rgba(59,130,246,0.3);
    }
    .bar:hover { opacity: 0.8; transform: translateY(-4px); box-shadow: 0 4px 16px rgba(59,130,246,0.5); }
    .bar-label { 
      position: absolute; bottom: -25px; left: 50%; transform: translateX(-50%);
      font-size: 0.7rem; color: var(--muted); white-space: nowrap;
    }
    
    /* Histogram */
    .histogram { display: flex; align-items: flex-end; height: 200px; gap: 2px; }
    .hist-bar { background: var(--primary); flex-grow: 1; border-radius: 2px 2px 0 0; transition: opacity 0.2s; }
    .hist-bar:hover { opacity: 0.7; }
    
    /* Distribution summary */
    .distrib-summary {
      display: grid; grid-template-columns: repeat(5, 1fr); gap: 0.75rem;
      margin-top: 1rem; padding: 1rem; background: rgba(59,130,246,0.05);
      border-radius: 0.5rem;
    }
    .distrib-stat { text-align: center; }
    .distrib-stat .label { font-size: 0.7rem; color: var(--muted); text-transform: uppercase; font-weight: 600; }
    .distrib-stat .value { font-size: 1.125rem; font-weight: 700; margin-top: 0.25rem; }
    
    /* Alert boxes */
    .alert { 
      padding: 1rem 1.25rem; border-radius: 0.5rem; margin: 1rem 0;
      border-left: 4px solid;
    }
    .alert-info { background: rgba(6,182,212,0.1); border-color: var(--info); }
    .alert-warning { background: rgba(245,158,11,0.1); border-color: var(--warn); }
    .alert-danger { background: rgba(239,68,68,0.1); border-color: var(--danger); }
    .alert-success { background: rgba(16,185,129,0.1); border-color: var(--ok); }
    
    /* Responsive */
    @media (max-width: 768px) {
      .container { padding: 1rem; }
      .page-header { padding: 2rem 1.5rem; }
      h1 { font-size: 1.75rem; }
      .kpi-grid, .viz-grid { grid-template-columns: 1fr; }
    }
    </style>
    """

    JS = """
    <script>
    document.addEventListener('DOMContentLoaded', () => {
        const payload = JSON.parse(document.getElementById('payload').textContent);
        const { samples, rows, region_totals, region_keys } = payload;
        
        console.log('Loaded', rows.length, 'samples');
        
        // ===== UTILITY FUNCTIONS =====
        function median(arr) {
            if (!arr.length) return 0;
            const s = [...arr].sort((a, b) => a - b);
            const mid = Math.floor(s.length / 2);
            return s.length % 2 !== 0 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
        }
        
        function mean(arr) {
            if (!arr.length) return 0;
            return arr.reduce((a,b) => a+b, 0) / arr.length;
        }
        
        function stdDev(arr) {
            if (!arr.length) return 0;
            const m = mean(arr);
            const variance = arr.reduce((sum, val) => sum + Math.pow(val - m, 2), 0) / arr.length;
            return Math.sqrt(variance);
        }
        
        
        function formatNumber(n, forceDecimals = false) {
            if (n >= 1e6) return (n/1e6).toFixed(1) + 'M';
            if (n >= 1e3) return (n/1e3).toFixed(1) + 'K';
            // For values < 1, show up to 3 decimal places
            if (n < 1 && n > 0) return n.toFixed(3);
            // For values 1-10, show 2 decimals
            if (n < 10) return n.toFixed(2);
            // For larger values < 1000, show 1 decimal or integer
            if (n < 100) return n.toFixed(1);
            return n.toFixed(0);
        }
        
        function renderDistribSummary(containerId, data) {
            const container = document.getElementById(containerId);
            if (!container) return;
            if (!data || data.length === 0) {
                container.innerHTML = '<div style="padding:0.5rem;text-align:center;color:var(--muted);font-size:0.875rem;">No data</div>';
                return;
            }
            const min = Math.min(...data);
            const max = Math.max(...data);
            const med = median(data);
            const avg = mean(data);
            const std = stdDev(data);
            container.innerHTML = `
                <div class="distrib-stat"><div class="label">Min</div><div class="value">${min.toFixed(1)}</div></div>
                <div class="distrib-stat"><div class="label">Max</div><div class="value">${max.toFixed(1)}</div></div>
                <div class="distrib-stat"><div class="label">Median</div><div class="value">${med.toFixed(1)}</div></div>
                <div class="distrib-stat"><div class="label">Mean</div><div class="value">${avg.toFixed(1)}</div></div>
                <div class="distrib-stat"><div class="label">Std Dev</div><div class="value">${std.toFixed(1)}</div></div>
            `;
        }
        
        function renderBarChart(containerId, data, labels = null) {
            const container = document.getElementById(containerId);
            if (!container) return;
            if (!data || data.length === 0) {
                container.innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No data available</div>';
                return;
            }
            
            // Get container dimensions - use parent or reasonable defaults
            const containerRect = container.getBoundingClientRect();
            const parentWidth = container.parentElement ? container.parentElement.getBoundingClientRect().width : 0;
            const containerWidth = Math.max(containerRect.width || parentWidth || 600, 500);
            const containerHeight = Math.max(containerRect.height || 320, 300);
            
            // Smart Y-axis max: add small padding (5%), but ensure it's tight to data
            const dataMax = Math.max(...data);
            const dataMin = Math.min(...data);
            const dataRange = dataMax - dataMin;
            let yMax = dataMax;
            if (dataRange > 0) {
                // Add 5% padding above max
                yMax = dataMax + (dataRange * 0.05);
                // For small ranges (like PI 0-1), use tighter scaling
                if (dataMax <= 1 && dataMin >= 0) {
                    yMax = Math.min(1, dataMax + 0.05);
                } else {
                    // Round to nice number for larger ranges
                    const magnitude = Math.pow(10, Math.floor(Math.log10(yMax)));
                    yMax = Math.ceil(yMax / magnitude) * magnitude;
                }
            } else {
                // Single value or no range - add small padding
                yMax = dataMax > 0 ? dataMax * 1.05 : 1;
            }
            // Ensure minimum height for visibility
            if (yMax <= dataMax) yMax = dataMax * 1.05;
            
            // Adaptive bar sizing to fill container
            const barCount = data.length;
            const margin = {top: 40, right: 20, bottom: 80, left: 70};
            const chartWidth = containerWidth - margin.left - margin.right;
            const chartHeight = containerHeight - margin.top - margin.bottom;
            const minBarWidth = 30;
            const maxBarWidth = 80;
            const barSpacing = 8;
            const totalBarSpace = chartWidth - (barSpacing * (barCount - 1));
            const barWidth = Math.max(minBarWidth, Math.min(maxBarWidth, totalBarSpace / barCount));
            
            // Calculate actual width needed (only for bars that exist)
            const actualChartWidth = (barWidth * barCount) + (barSpacing * (barCount - 1));
            const width = margin.left + actualChartWidth + margin.right;
            const height = containerHeight;
            
            // Color palette
            const colorPalette = ['#3b82f6', '#8b5cf6', '#ec4899', '#14b8a6', '#f59e0b', '#10b981', '#6366f1', '#06b6d4'];
            
            // Create SVG
            let svg = `<svg width="100%" height="${height}" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" style="max-width:${width}px;font-family:var(--font);">`;
            
            // Gradients
            svg += '<defs>';
            colorPalette.forEach((color, i) => {
                svg += `<linearGradient id="grad${i}" x1="0%" y1="0%" x2="0%" y2="100%">
                    <stop offset="0%" style="stop-color:${color};stop-opacity:0.9" />
                    <stop offset="100%" style="stop-color:${color};stop-opacity:0.7" />
                </linearGradient>
                <filter id="shadow${i}">
                    <feDropShadow dx="0" dy="2" stdDeviation="3" flood-opacity="0.3"/>
                </filter>`;
            });
            svg += '</defs>';
            
            // Grid lines
            const ySteps = 5;
            for (let i = 0; i <= ySteps; i++) {
                const y = margin.top + (chartHeight * i / ySteps);
                svg += `<line x1="${margin.left}" y1="${y}" x2="${width - margin.right}" y2="${y}" stroke="var(--line)" stroke-width="0.5" opacity="0.4"/>`;
            }
            
            // Y-axis
            svg += `<line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${height - margin.bottom}" stroke="var(--muted)" stroke-width="1.5" opacity="0.6"/>`;
            
            // Y-axis labels
            for (let i = 0; i <= ySteps; i++) {
                const val = yMax - (yMax * i / ySteps);
                const y = margin.top + (chartHeight * i / ySteps);
                svg += `<text x="${margin.left - 12}" y="${y + 4}" text-anchor="end" font-size="11" fill="var(--muted)" font-weight="500">${formatNumber(val)}</text>`;
            }
            
            // Bars
            data.forEach((v, i) => {
                const x = margin.left + (i * (barWidth + barSpacing));
                const barHeight = Math.max(2, (v / yMax) * chartHeight);
                const y = height - margin.bottom - barHeight;
                const colorIdx = i % colorPalette.length;
                const label = labels ? labels[i] : `Item ${i+1}`;
                
                // Bar shadow
                svg += `<rect x="${x}" y="${y+2}" width="${barWidth}" height="${barHeight}" fill="black" opacity="0.1" rx="3"/>`;
                
                // Bar
                const formattedValue = v < 1 ? v.toFixed(3) : (v < 10 ? v.toFixed(2) : v.toLocaleString());
                svg += `<rect x="${x}" y="${y}" width="${barWidth}" height="${barHeight}" fill="url(#grad${colorIdx})" rx="3" filter="url(#shadow${colorIdx})" style="cursor:pointer;transition:opacity 0.15s;" onmouseover="evt.target.style.opacity='0.7'" onmouseout="evt.target.style.opacity='1'">
                    <title>${label}: ${formattedValue}</title>
                </rect>`;
                
                // Value labels on top
                if (barHeight > 25) {
                    svg += `<text x="${x + barWidth/2}" y="${y - 8}" text-anchor="middle" font-size="10" fill="var(--fg)" font-weight="700" opacity="0.8">${formatNumber(v)}</text>`;
                }
            });
            
            // X-axis (only as wide as needed)
            const xAxisEnd = margin.left + actualChartWidth;
            svg += `<line x1="${margin.left}" y1="${height - margin.bottom}" x2="${xAxisEnd}" y2="${height - margin.bottom}" stroke="var(--muted)" stroke-width="1.5" opacity="0.6"/>`;
            
            // X-axis labels
            data.forEach((v, i) => {
                const x = margin.left + (i * (barWidth + barSpacing)) + barWidth/2;
                const label = labels ? labels[i] : `Item ${i+1}`;
                const displayLabel = label.length > 10 ? label.substring(0,9) + '..' : label;
                
                svg += `<text x="${x}" y="${height - margin.bottom + 18}" text-anchor="end" font-size="10" fill="var(--muted)" font-weight="500" transform="rotate(-45 ${x} ${height - margin.bottom + 18})"><title>${label}</title>${displayLabel}</text>`;
            });
            
            // X-axis label
            svg += `<text x="${width/2}" y="${height - 8}" text-anchor="middle" font-size="11" fill="var(--muted)" font-weight="600">Sample</text>`;
            
            svg += '</svg>';
            container.innerHTML = svg;
        }
        
        function renderHistogram(containerId, data, bins = 20, xAxisLabel = 'Value', labels = null) {
            const container = document.getElementById(containerId);
            if (!container) return;
            if (!data || data.length === 0) {
                container.innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No data available</div>';
                return;
            }
            if (data.length === 1) {
                const sampleInfo = labels && labels[0] ? ` (${labels[0]})` : '';
                container.innerHTML = `<div style="padding:2rem;text-align:center;color:var(--muted)">Single sample${sampleInfo}: <strong>${data[0].toFixed(2)}</strong></div>`;
                return;
            }
            
            const min = Math.min(...data);
            const max = Math.max(...data);
            const range = max - min;
            
            if (range === 0) {
                container.innerHTML = `<div style="padding:2rem;text-align:center;color:var(--muted)">All samples: <strong>${min.toFixed(2)}</strong></div>`;
                return;
            }
            
            // Get container dimensions - use parent or reasonable defaults
            const containerRect = container.getBoundingClientRect();
            const parentWidth = container.parentElement ? container.parentElement.getBoundingClientRect().width : 0;
            const containerWidth = Math.max(containerRect.width || parentWidth || 600, 500);
            const containerHeight = Math.max(containerRect.height || 360, 340);
            
            // Create histogram bins with sample tracking
            const binSize = range / bins;
            const histogram = new Array(bins).fill(0).map(() => ({ count: 0, samples: [] }));
            data.forEach((v, idx) => {
                const binIndex = Math.min(Math.floor((v - min) / binSize), bins - 1);
                histogram[binIndex].count++;
                if (labels && labels[idx]) {
                    histogram[binIndex].samples.push({ name: labels[idx], value: v });
                }
            });
            const maxCount = Math.max(...histogram.map(h => h.count), 1);
            
            // Smart Y-axis max: add 10% padding
            const yMax = Math.ceil(maxCount * 1.1);
            
            // Chart dimensions - fill container
            const margin = {top: 30, right: 20, bottom: 50, left: 60};
            const width = containerWidth;
            const height = containerHeight;
            const chartWidth = width - margin.left - margin.right;
            const chartHeight = height - margin.top - margin.bottom;
            const barWidth = chartWidth / bins;
            
            // Create SVG
            let svg = `<svg width="100%" height="${height}" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" style="max-width:${width}px;font-family:var(--font);">`;
            
            // Gradients
            svg += `<defs>
                <linearGradient id="histGrad" x1="0%" y1="0%" x2="0%" y2="100%">
                    <stop offset="0%" style="stop-color:#3b82f6;stop-opacity:0.9" />
                    <stop offset="100%" style="stop-color:#8b5cf6;stop-opacity:0.7" />
                </linearGradient>
                <filter id="histShadow">
                    <feGaussianBlur in="SourceAlpha" stdDeviation="2"/>
                    <feOffset dx="0" dy="1" result="offsetblur"/>
                    <feComponentTransfer>
                        <feFuncA type="linear" slope="0.3"/>
                    </feComponentTransfer>
                    <feMerge>
                        <feMergeNode/>
                        <feMergeNode in="SourceGraphic"/>
                    </feMerge>
                </filter>
            </defs>`;
            
            // Grid lines
            const ySteps = 4;
            for (let i = 0; i <= ySteps; i++) {
                const y = margin.top + (chartHeight * i / ySteps);
                svg += `<line x1="${margin.left}" y1="${y}" x2="${width - margin.right}" y2="${y}" stroke="var(--line)" stroke-width="0.5" opacity="0.3"/>`;
            }
            
            // Y-axis
            svg += `<line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${height - margin.bottom}" stroke="var(--muted)" stroke-width="1.5" opacity="0.6"/>`;
            
            // Y-axis labels (count)
            for (let i = 0; i <= ySteps; i++) {
                const val = Math.round(yMax - (yMax * i / ySteps));
                const y = margin.top + (chartHeight * i / ySteps);
                svg += `<text x="${margin.left - 10}" y="${y + 4}" text-anchor="end" font-size="11" fill="var(--muted)" font-weight="500">${val}</text>`;
            }
            
            // Y-axis label
            svg += `<text x="${15}" y="${height/2}" text-anchor="middle" font-size="11" fill="var(--muted)" font-weight="600" transform="rotate(-90 15 ${height/2})">Count</text>`;
            
            // Histogram bars
            histogram.forEach((bin, i) => {
                const x = margin.left + (i * barWidth);
                const barHeight = (bin.count / yMax) * chartHeight;
                const y = height - margin.bottom - barHeight;
                const binStart = min + (i * binSize);
                const binEnd = binStart + binSize;
                
                // Format bin range based on magnitude
                const formatBinValue = (v) => v < 1 ? v.toFixed(3) : (v < 10 ? v.toFixed(2) : v.toFixed(1));
                
                // Build tooltip with sample names if available
                let tooltip = `${bin.count} sample${bin.count !== 1 ? 's' : ''} in range ${formatBinValue(binStart)} - ${formatBinValue(binEnd)}`;
                if (bin.samples.length > 0 && bin.samples.length <= 8) {
                    tooltip += ':\\n' + bin.samples.map(s => `  ${s.name}: ${formatBinValue(s.value)}`).join('\\n');
                } else if (bin.samples.length > 8) {
                    tooltip += ':\\n' + bin.samples.slice(0, 8).map(s => `  ${s.name}: ${formatBinValue(s.value)}`).join('\\n') + `\\n  ... and ${bin.samples.length - 8} more`;
                }
                
                // Bar
                svg += `<rect x="${x}" y="${y}" width="${barWidth - 1}" height="${barHeight}" fill="url(#histGrad)" filter="url(#histShadow)" rx="2" style="cursor:pointer;transition:opacity 0.15s;" onmouseover="evt.target.style.opacity='0.7'" onmouseout="evt.target.style.opacity='1'">
                    <title>${tooltip}</title>
                </rect>`;
            });
            
            // X-axis
            svg += `<line x1="${margin.left}" y1="${height - margin.bottom}" x2="${width - margin.right}" y2="${height - margin.bottom}" stroke="var(--muted)" stroke-width="1.5" opacity="0.6"/>`;
            
            // X-axis labels (value range)
            const xLabelCount = Math.min(5, bins);
            for (let i = 0; i <= xLabelCount; i++) {
                const val = min + (range * i / xLabelCount);
                const x = margin.left + (chartWidth * i / xLabelCount);
                const formattedVal = val < 1 ? val.toFixed(3) : (val < 10 ? val.toFixed(2) : val.toFixed(1));
                svg += `<text x="${x}" y="${height - margin.bottom + 20}" text-anchor="middle" font-size="10" fill="var(--muted)" font-weight="500">${formattedVal}</text>`;
            }
            
            // X-axis label
            svg += `<text x="${width/2}" y="${height - 8}" text-anchor="middle" font-size="11" fill="var(--muted)" font-weight="600">${xAxisLabel}</text>`;
            
            svg += '</svg>';
            container.innerHTML = svg;
        }
        
        function renderScatterPlot(containerId, xData, yData, xLabel, yLabel, labels = null) {
            const container = document.getElementById(containerId);
            if (!container) return;
            if (!xData || !yData || xData.length === 0 || yData.length === 0) {
                container.innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No data available</div>';
                return;
            }
            
            // Filter valid pairs
            const points = [];
            for (let i = 0; i < Math.min(xData.length, yData.length); i++) {
                if (xData[i] != null && yData[i] != null && isFinite(xData[i]) && isFinite(yData[i])) {
                    points.push({x: xData[i], y: yData[i], label: labels ? labels[i] : `Sample ${i+1}`});
                }
            }
            
            if (points.length === 0) {
                container.innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No valid data points</div>';
                return;
            }
            
            const containerRect = container.getBoundingClientRect();
            const parentWidth = container.parentElement ? container.parentElement.getBoundingClientRect().width : 0;
            const containerWidth = Math.max(containerRect.width || parentWidth || 600, 500);
            const containerHeight = Math.max(containerRect.height || 320, 300);
            
            const xMin = Math.min(...points.map(p => p.x));
            const xMax = Math.max(...points.map(p => p.x));
            const yMin = Math.min(...points.map(p => p.y));
            const yMax = Math.max(...points.map(p => p.y));
            const xRange = xMax - xMin || 1;
            const yRange = yMax - yMin || 1;
            
            // Add padding
            const xPadding = xRange * 0.1;
            const yPadding = yRange * 0.1;
            const xScaleMin = xMin - xPadding;
            const xScaleMax = xMax + xPadding;
            const yScaleMin = yMin - yPadding;
            const yScaleMax = yMax + yPadding;
            
            const margin = {top: 30, right: 20, bottom: 60, left: 70};
            const width = containerWidth;
            const height = containerHeight;
            const chartWidth = width - margin.left - margin.right;
            const chartHeight = height - margin.top - margin.bottom;
            
            const scaleX = (val) => margin.left + ((val - xScaleMin) / (xScaleMax - xScaleMin)) * chartWidth;
            const scaleY = (val) => height - margin.bottom - ((val - yScaleMin) / (yScaleMax - yScaleMin)) * chartHeight;
            
            let svg = `<svg width="100%" height="${height}" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" style="font-family:var(--font);">`;
            
            // Grid lines
            const xSteps = 5;
            const ySteps = 5;
            for (let i = 0; i <= xSteps; i++) {
                const x = margin.left + (chartWidth * i / xSteps);
                svg += `<line x1="${x}" y1="${margin.top}" x2="${x}" y2="${height - margin.bottom}" stroke="var(--line)" stroke-width="0.5" opacity="0.3"/>`;
            }
            for (let i = 0; i <= ySteps; i++) {
                const y = margin.top + (chartHeight * i / ySteps);
                svg += `<line x1="${margin.left}" y1="${y}" x2="${width - margin.right}" y2="${y}" stroke="var(--line)" stroke-width="0.5" opacity="0.3"/>`;
            }
            
            // Axes
            svg += `<line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${height - margin.bottom}" stroke="var(--muted)" stroke-width="1.5" opacity="0.6"/>`;
            svg += `<line x1="${margin.left}" y1="${height - margin.bottom}" x2="${width - margin.right}" y2="${height - margin.bottom}" stroke="var(--muted)" stroke-width="1.5" opacity="0.6"/>`;
            
            // Axis labels
            for (let i = 0; i <= xSteps; i++) {
                const val = xScaleMin + (xScaleMax - xScaleMin) * i / xSteps;
                const x = margin.left + (chartWidth * i / xSteps);
                const formattedVal = val < 1 ? val.toFixed(3) : (val < 10 ? val.toFixed(2) : val.toFixed(1));
                svg += `<text x="${x}" y="${height - margin.bottom + 20}" text-anchor="middle" font-size="10" fill="var(--muted)" font-weight="500">${formattedVal}</text>`;
            }
            for (let i = 0; i <= ySteps; i++) {
                const val = yScaleMax - (yScaleMax - yScaleMin) * i / ySteps;
                const y = margin.top + (chartHeight * i / ySteps);
                const formattedVal = val < 1 ? val.toFixed(3) : (val < 10 ? val.toFixed(2) : val.toFixed(1));
                svg += `<text x="${margin.left - 10}" y="${y + 4}" text-anchor="end" font-size="10" fill="var(--muted)" font-weight="500">${formattedVal}</text>`;
            }
            
            // Axis titles
            svg += `<text x="${width/2}" y="${height - 10}" text-anchor="middle" font-size="11" fill="var(--muted)" font-weight="600">${xLabel}</text>`;
            svg += `<text x="${15}" y="${height/2}" text-anchor="middle" font-size="11" fill="var(--muted)" font-weight="600" transform="rotate(-90 15 ${height/2})">${yLabel}</text>`;
            
            // Color palette for scatter points
            const colorPalette = ['#3b82f6', '#8b5cf6', '#ec4899', '#14b8a6', '#f59e0b', '#10b981', '#6366f1', '#06b6d4'];
            
            // Points with different colors per sample
            points.forEach((p, i) => {
                const x = scaleX(p.x);
                const y = scaleY(p.y);
                const colorIdx = i % colorPalette.length;
                const color = colorPalette[colorIdx];
                const formatScatterVal = (v) => v < 1 ? v.toFixed(4) : (v < 10 ? v.toFixed(3) : v.toFixed(2));
                svg += `<circle cx="${x}" cy="${y}" r="5" fill="${color}" opacity="0.7" stroke="white" stroke-width="1.5" style="cursor:pointer;" onmouseover="evt.target.setAttribute('r', '7');evt.target.setAttribute('opacity', '1');" onmouseout="evt.target.setAttribute('r', '5');evt.target.setAttribute('opacity', '0.7');">
                    <title>${p.label}\n${xLabel}: ${formatScatterVal(p.x)}\n${yLabel}: ${formatScatterVal(p.y)}</title>
                </circle>`;
            });
            
            svg += '</svg>';
            container.innerHTML = svg;
        }
        
        function renderStackedBarChart(containerId, dataBySample, regionKeys, sampleLabels) {
            const container = document.getElementById(containerId);
            if (!container) return;
            if (!dataBySample || dataBySample.length === 0 || !regionKeys || regionKeys.length === 0) {
                container.innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No data available</div>';
                return;
            }
            
            const containerRect = container.getBoundingClientRect();
            const parentWidth = container.parentElement ? container.parentElement.getBoundingClientRect().width : 0;
            const containerWidth = Math.max(containerRect.width || parentWidth || 600, 500);
            const containerHeight = Math.max(containerRect.height || 320, 300);
            
            // Calculate percentages for each sample
            const percentages = dataBySample.map(sample => {
                const total = Object.values(sample).reduce((sum, v) => sum + (v || 0), 0);
                const pct = {};
                regionKeys.forEach(key => {
                    pct[key] = total > 0 ? ((sample[key] || 0) / total) * 100 : 0;
                });
                return pct;
            });
            
            const margin = {top: 40, right: 20, bottom: 80, left: 70};
            const width = containerWidth;
            const height = containerHeight;
            const chartWidth = width - margin.left - margin.right;
            const chartHeight = height - margin.top - margin.bottom;
            
            const barCount = dataBySample.length;
            const barSpacing = 8;
            const barWidth = Math.max(30, (chartWidth - (barSpacing * (barCount - 1))) / barCount);
            
            const colorPalette = ['#3b82f6', '#8b5cf6', '#ec4899', '#14b8a6', '#f59e0b', '#10b981', '#6366f1', '#06b6d4'];
            
            let svg = `<svg width="100%" height="${height}" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" style="font-family:var(--font);">`;
            
            // Gradients
            svg += '<defs>';
            colorPalette.forEach((color, i) => {
                svg += `<linearGradient id="stackGrad${i}" x1="0%" y1="0%" x2="0%" y2="100%">
                    <stop offset="0%" style="stop-color:${color};stop-opacity:0.9" />
                    <stop offset="100%" style="stop-color:${color};stop-opacity:0.7" />
                </linearGradient>`;
            });
            svg += '</defs>';
            
            // Grid and axes
            const ySteps = 5;
            for (let i = 0; i <= ySteps; i++) {
                const y = margin.top + (chartHeight * i / ySteps);
                svg += `<line x1="${margin.left}" y1="${y}" x2="${width - margin.right}" y2="${y}" stroke="var(--line)" stroke-width="0.5" opacity="0.4"/>`;
            }
            svg += `<line x1="${margin.left}" y1="${margin.top}" x2="${margin.left}" y2="${height - margin.bottom}" stroke="var(--muted)" stroke-width="1.5" opacity="0.6"/>`;
            svg += `<line x1="${margin.left}" y1="${height - margin.bottom}" x2="${width - margin.right}" y2="${height - margin.bottom}" stroke="var(--muted)" stroke-width="1.5" opacity="0.6"/>`;
            
            // Y-axis labels (0-100%)
            for (let i = 0; i <= ySteps; i++) {
                const val = 100 - (100 * i / ySteps);
                const y = margin.top + (chartHeight * i / ySteps);
                svg += `<text x="${margin.left - 12}" y="${y + 4}" text-anchor="end" font-size="11" fill="var(--muted)" font-weight="500">${val}%</text>`;
            }
            
            // Stacked bars
            percentages.forEach((pct, sampleIdx) => {
                const x = margin.left + (sampleIdx * (barWidth + barSpacing));
                let yAccum = height - margin.bottom;
                
                regionKeys.forEach((key, regionIdx) => {
                    const segmentHeight = (pct[key] / 100) * chartHeight;
                    const colorIdx = regionIdx % colorPalette.length;
                    const y = yAccum - segmentHeight;
                    
                    svg += `<rect x="${x}" y="${y}" width="${barWidth}" height="${segmentHeight}" fill="url(#stackGrad${colorIdx})" rx="0" style="cursor:pointer;" onmouseover="evt.target.style.opacity='0.8'" onmouseout="evt.target.style.opacity='1'">
                        <title>${sampleLabels[sampleIdx]}: ${key} = ${pct[key].toFixed(1)}%</title>
                    </rect>`;
                    
                    yAccum = y;
                });
            });
            
            // X-axis labels
            percentages.forEach((pct, sampleIdx) => {
                const x = margin.left + (sampleIdx * (barWidth + barSpacing)) + barWidth/2;
                const label = sampleLabels[sampleIdx] || `Sample ${sampleIdx+1}`;
                const displayLabel = label.length > 10 ? label.substring(0,9) + '..' : label;
                svg += `<text x="${x}" y="${height - margin.bottom + 18}" text-anchor="end" font-size="10" fill="var(--muted)" font-weight="500" transform="rotate(-45 ${x} ${height - margin.bottom + 18})"><title>${label}</title>${displayLabel}</text>`;
            });
            
            // Legend
            const legendX = width - margin.right - 150;
            const legendY = margin.top + 20;
            regionKeys.forEach((key, i) => {
                const colorIdx = i % colorPalette.length;
                const y = legendY + (i * 18);
                svg += `<rect x="${legendX}" y="${y - 8}" width="12" height="12" fill="url(#stackGrad${colorIdx})" rx="2"/>`;
                svg += `<text x="${legendX + 18}" y="${y}" font-size="10" fill="var(--fg)" font-weight="500">${key}</text>`;
            });
            
            svg += '</svg>';
            container.innerHTML = svg;
        }
        
        // ===== OVERVIEW SECTION =====
        document.getElementById('sample-count-header').textContent = rows.length + ' samples analyzed';
        document.getElementById('kpi-total').textContent = rows.length;
        
        const conditions = [...new Set(rows.map(r => r.condition).filter(Boolean))];
        document.getElementById('kpi-conditions').textContent = conditions.length;
        
        const depthValues = rows.map(r => r.input_reads || 0).filter(v => v > 0);
        const avgDepth = depthValues.length > 0 ? mean(depthValues) / 1e6 : 0;
        document.getElementById('kpi-depth').textContent = avgDepth > 0 ? avgDepth.toFixed(1) + 'M' : '-';
        
        const totalDiv = rows.reduce((sum, r) => sum + (r.divergent_regions || 0), 0);
        document.getElementById('kpi-div-total').textContent = totalDiv.toLocaleString();
        
        const avgRegions = mean(rows.map(r => r.total_regions || 0));
        document.getElementById('kpi-regions-avg').textContent = avgRegions.toFixed(0);
        
        // ===== QC SECTION =====
        // Extract read depths (filter out null/undefined values)
        const depths = rows.map(r => r.input_reads).filter(v => v != null && v > 0).map(v => v / 1e6);
        
        // Extract duplication rates - check if UMI deduplication was used
        const hasUMIDedup = rows.some(r => r.umi_deduplication_enabled === true);
        const dups = rows.map(r => {
            if (r.umi_deduplication_enabled === true) {
                return r.umi_deduplication_percent != null ? r.umi_deduplication_percent : null;
            } else {
                return r.duplicate_percent != null ? r.duplicate_percent : null;
            }
        }).filter(v => v != null);
        
        // Display median depth
        if (depths.length > 0) {
            document.getElementById('qc-median-depth').textContent = median(depths).toFixed(1) + 'M';
        } else {
            document.getElementById('qc-median-depth').textContent = '-';
        }
        
        // Display median duplication with special handling for no UMI
        if (dups.length > 0) {
            document.getElementById('qc-median-dup').textContent = median(dups).toFixed(1) + '%';
        } else if (!hasUMIDedup) {
            document.getElementById('qc-median-dup').textContent = 'No deduplication done';
        } else {
            document.getElementById('qc-median-dup').textContent = '-';
        }
        
        // Render per-sample depth bar chart
        if (depths.length > 0) {
            const depthSampleIds = rows.filter(r => r.input_reads != null && r.input_reads > 0).map(r => r.sample_id);
            renderBarChart('chart-depth-per-sample', depths, depthSampleIds);
        } else {
            document.getElementById('chart-depth-per-sample').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No read depth data available</div>';
        }
        
        // Render per-sample duplication bar chart
        if (dups.length > 0) {
            const dupSampleIds = rows.filter(r => {
                if (r.umi_deduplication_enabled === true) {
                    return r.umi_deduplication_percent != null;
                } else {
                    return r.duplicate_percent != null;
                }
            }).map(r => r.sample_id);
            renderBarChart('chart-dup-per-sample', dups, dupSampleIds);
        } else if (!hasUMIDedup) {
            document.getElementById('chart-dup-per-sample').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No deduplication done</div>';
        } else {
            document.getElementById('chart-dup-per-sample').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No duplication data available</div>';
        }
        
        // Render depth distribution histogram (only for larger cohorts)
        if (depths.length > 0) {
            if (depths.length >= 10) {
                const depthSampleIds = rows.filter(r => r.input_reads != null && r.input_reads > 0).map(r => r.sample_id);
                renderHistogram('chart-depth-dist', depths, 15, 'Read Depth (M)', depthSampleIds);
                renderDistribSummary('depth-stats', depths);
            } else {
                // For small cohorts, show summary stats instead of histogram
                document.getElementById('chart-depth-dist').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">Distribution histograms require ≥10 samples<br>(current: ' + depths.length + ' samples)</div>';
                renderDistribSummary('depth-stats', depths);
            }
        } else {
            document.getElementById('chart-depth-dist').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No read depth data available</div>';
        }
        
        // Render duplication distribution histogram (only for larger cohorts)
        if (dups.length > 0) {
            if (dups.length >= 10) {
                const dupSampleIds = rows.filter(r => {
                    if (r.umi_deduplication_enabled === true) {
                        return r.umi_deduplication_percent != null;
                    } else {
                        return r.duplicate_percent != null;
                    }
                }).map(r => r.sample_id);
                renderHistogram('chart-dup-dist', dups, 15, 'Duplication Rate (%)', dupSampleIds);
                renderDistribSummary('dup-stats', dups);
            } else {
                document.getElementById('chart-dup-dist').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">Distribution histograms require ≥10 samples<br>(current: ' + dups.length + ' samples)</div>';
                renderDistribSummary('dup-stats', dups);
            }
        } else if (!hasUMIDedup) {
            document.getElementById('chart-dup-dist').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No deduplication done</div>';
        } else {
            document.getElementById('chart-dup-dist').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No duplication data available</div>';
        }
        
        // ===== DIVERGENT SECTION =====
        const divRegions = rows.map(r => r.divergent_regions || 0);
        const divMean = mean(divRegions);
        const divMedian = median(divRegions);
        const divMin = Math.min(...divRegions);
        const divMax = Math.max(...divRegions);
        
        document.getElementById('div-total-regions').textContent = totalDiv.toLocaleString();
        document.getElementById('div-mean').textContent = divMean.toFixed(0);
        document.getElementById('div-median').textContent = divMedian.toFixed(0);
        document.getElementById('div-range').textContent = `${divMin}-${divMax}`;
        
        renderBarChart('chart-div-per-sample', divRegions, rows.map(r => r.sample_id));
        
        // Only show histogram for larger cohorts
        if (divRegions.length >= 10) {
            renderHistogram('chart-div-hist', divRegions, 15, 'Divergent Regions per Sample', rows.map(r => r.sample_id));
            renderDistribSummary('div-stats', divRegions);
        } else {
            document.getElementById('chart-div-hist').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">Distribution histograms require ≥10 samples<br>(current: ' + divRegions.length + ' samples)</div>';
            renderDistribSummary('div-stats', divRegions);
        }
        
        // Divergent by condition
        const divByCondition = {};
        rows.forEach(r => {
            const cond = r.condition || 'Unknown';
            if (!divByCondition[cond]) divByCondition[cond] = [];
            divByCondition[cond].push(r.divergent_regions || 0);
        });
        const divCondLabels = Object.keys(divByCondition);
        const divCondMeans = divCondLabels.map(c => mean(divByCondition[c]));
        renderBarChart('chart-div-by-condition', divCondMeans, divCondLabels);
        
        // Replicate CV
        const divCVs = divCondLabels.map(c => {
            const data = divByCondition[c];
            return data.length > 1 ? (stdDev(data) / mean(data)) * 100 : 0;
        });
        renderBarChart('chart-div-cv', divCVs, divCondLabels);
        
        // ===== PAUSING SECTION =====
        const pausingIndices = rows.map(r => r.median_pausing_index).filter(p => p != null && p > 0 && isFinite(p));
        if (pausingIndices.length > 0) {
            const piMedian = median(pausingIndices);
            const piMean = mean(pausingIndices);
            const piStd = stdDev(pausingIndices);
            const piMin = Math.min(...pausingIndices);
            const piMax = Math.max(...pausingIndices);
            
            document.getElementById('pi-cohort-median').textContent = piMedian.toFixed(2);
            document.getElementById('pi-cohort-mean').textContent = piMean.toFixed(2);
            document.getElementById('pi-std').textContent = piStd.toFixed(2);
            document.getElementById('pi-range').textContent = `${piMin.toFixed(2)}-${piMax.toFixed(2)}`;
            
            const piSampleIds = rows.filter(r => r.median_pausing_index != null && r.median_pausing_index > 0).map(r => r.sample_id);
            renderBarChart('chart-pi-per-sample', pausingIndices, piSampleIds);
            
            // Only show histogram for larger cohorts
            if (pausingIndices.length >= 10) {
                renderHistogram('chart-pi-dist', pausingIndices, 15, 'Median Pausing Index', piSampleIds);
                renderDistribSummary('pi-stats', pausingIndices);
            } else {
                document.getElementById('chart-pi-dist').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">Distribution histograms require ≥10 samples<br>(current: ' + pausingIndices.length + ' samples)</div>';
                renderDistribSummary('pi-stats', pausingIndices);
            }
            
            // PI by condition
            const piByCondition = {};
            rows.forEach(r => {
                if (r.median_pausing_index != null && r.median_pausing_index > 0 && isFinite(r.median_pausing_index)) {
                    const cond = r.condition || 'Unknown';
                    if (!piByCondition[cond]) piByCondition[cond] = [];
                    piByCondition[cond].push(r.median_pausing_index);
                }
            });
            const piCondLabels = Object.keys(piByCondition);
            const piCondMeans = piCondLabels.map(c => mean(piByCondition[c]));
            renderBarChart('chart-pi-by-condition', piCondMeans, piCondLabels);
        } else {
            document.getElementById('pi-cohort-median').textContent = '-';
            document.getElementById('pi-cohort-mean').textContent = '-';
            document.getElementById('pi-std').textContent = '-';
            document.getElementById('pi-range').textContent = '-';
            document.getElementById('chart-pi-per-sample').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No pausing index data available</div>';
            document.getElementById('chart-pi-dist').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No data</div>';
            document.getElementById('chart-pi-by-condition').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No data</div>';
        }
        
        // ===== FUNCTIONAL REGIONS SECTION =====
        if (region_totals && region_keys && region_keys.length > 0) {
            const regionData = region_keys.filter(k => !k.toLowerCase().includes('localized')).map(k => region_totals[k] || 0);
            const regionLabels = region_keys.filter(k => !k.toLowerCase().includes('localized'));
            
            if (regionLabels.length > 0 && regionData.some(d => d > 0)) {
                renderBarChart('chart-region-totals', regionData, regionLabels);
                
                // Region counts table
                const regionTableHead = document.getElementById('region-counts-thead');
                const regionTableBody = document.getElementById('region-counts-tbody');
                regionTableHead.innerHTML = '<tr><th>Sample</th>' + regionLabels.map(k => `<th>${k}</th>`).join('') + '</tr>';
                regionTableBody.innerHTML = rows.map(r => {
                    return '<tr><td>' + r.sample_id + '</td>' +
                        regionLabels.map(k => {
                            const count = r['count_' + k] || 0;
                            return `<td>${count.toLocaleString()}</td>`;
                        }).join('') + '</tr>';
                }).join('');
            } else {
                document.getElementById('chart-region-totals').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No functional region data available</div>';
            }
        } else {
            document.getElementById('chart-region-totals').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No region data available</div>';
        }
        
        // Region composition by sample - Stacked bar chart
        if (region_keys && region_keys.length > 0) {
            const filteredRegionKeys = region_keys.filter(k => !k.toLowerCase().includes('localized'));
            if (filteredRegionKeys.length > 0) {
                const regionDataBySample = rows.map(r => {
                    const sampleData = {};
                    filteredRegionKeys.forEach(key => {
                        sampleData[key] = r['func_' + key] || 0;
                    });
                    return sampleData;
                });
                const sampleLabels = rows.map(r => r.sample_id);
                renderStackedBarChart('chart-region-composition', regionDataBySample, filteredRegionKeys, sampleLabels);
            } else {
                document.getElementById('chart-region-composition').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No region data available</div>';
            }
        } else {
            document.getElementById('chart-region-composition').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No region data available</div>';
        }
        
        // Promoter and gene body per-sample signals
        const promoterSignals = rows.map(r => r.func_Promoter || 0).filter(v => v > 0);
        const geneBodySignals = rows.map(r => r['func_Gene body'] || 0).filter(v => v > 0);
        
        // Render per-sample promoter bar chart
        if (promoterSignals.length > 0) {
            const promoterSampleIds = rows.filter(r => (r.func_Promoter || 0) > 0).map(r => r.sample_id);
            renderBarChart('chart-promoter-per-sample', promoterSignals, promoterSampleIds);
        } else {
            document.getElementById('chart-promoter-per-sample').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No promoter data</div>';
        }
        
        // Render per-sample gene body bar chart
        if (geneBodySignals.length > 0) {
            const geneBodySampleIds = rows.filter(r => (r['func_Gene body'] || 0) > 0).map(r => r.sample_id);
            renderBarChart('chart-genebody-per-sample', geneBodySignals, geneBodySampleIds);
        } else {
            document.getElementById('chart-genebody-per-sample').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No gene body data</div>';
        }
        
        // ===== NORMALIZATION SECTION =====
        const cpmFactors = rows.map(r => r.cpm_factor).filter(f => f != null && isFinite(f));
        const sicpmFactors = rows.map(r => r.crpmsi_factor).filter(f => f != null && isFinite(f));
        
        document.getElementById('norm-cpm-count').textContent = cpmFactors.length;
        document.getElementById('norm-sicpm-count').textContent = sicpmFactors.length;
        
        if (cpmFactors.length > 0) {
            const cpmMin = Math.min(...cpmFactors);
            const cpmMax = Math.max(...cpmFactors);
            document.getElementById('norm-cpm-range').textContent = `${cpmMin.toFixed(2)}-${cpmMax.toFixed(2)}`;
            renderBarChart('chart-cpm-factors', cpmFactors, rows.filter(r => r.cpm_factor != null).map(r => r.sample_id));
        } else {
            document.getElementById('norm-cpm-range').textContent = 'No data';
            document.getElementById('chart-cpm-factors').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No CPM normalization factors available</div>';
        }
        
        if (sicpmFactors.length > 0) {
            const sicpmMin = Math.min(...sicpmFactors);
            const sicpmMax = Math.max(...sicpmFactors);
            document.getElementById('norm-sicpm-range').textContent = `${sicpmMin.toFixed(2)}-${sicpmMax.toFixed(2)}`;
            renderBarChart('chart-sicpm-factors', sicpmFactors, rows.filter(r => r.crpmsi_factor != null).map(r => r.sample_id));
        } else {
            document.getElementById('norm-sicpm-range').textContent = 'No data';
            document.getElementById('chart-sicpm-factors').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No spike-in normalization available</div>';
        }
        
        // CPM vs siCPM comparison scatter plot
        const cpmVsSicpmPairs = rows.filter(r => r.cpm_factor != null && r.crpmsi_factor != null && isFinite(r.cpm_factor) && isFinite(r.crpmsi_factor));
        if (cpmVsSicpmPairs.length > 0) {
            const cpmValues = cpmVsSicpmPairs.map(r => r.cpm_factor);
            const sicpmValues = cpmVsSicpmPairs.map(r => r.crpmsi_factor);
            const sampleLabels = cpmVsSicpmPairs.map(r => r.sample_id);
            renderScatterPlot('chart-cpm-vs-sicpm', cpmValues, sicpmValues, 'CPM Factor', 'siCPM Factor', sampleLabels);
        } else {
            document.getElementById('chart-cpm-vs-sicpm').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No CPM/siCPM data available</div>';
        }
        
        // Normalization CV by condition
        const normCVByCondition = {};
        conditions.forEach(cond => {
            const condRows = rows.filter(r => r.condition === cond);
            const cpmVals = condRows.map(r => r.cpm_factor).filter(f => f != null && isFinite(f));
            const sicpmVals = condRows.map(r => r.crpmsi_factor).filter(f => f != null && isFinite(f));
            
            if (cpmVals.length > 1) {
                const cpmCV = (stdDev(cpmVals) / mean(cpmVals)) * 100;
                normCVByCondition[cond] = normCVByCondition[cond] || {};
                normCVByCondition[cond].cpm = cpmCV;
            }
            if (sicpmVals.length > 1) {
                const sicpmCV = (stdDev(sicpmVals) / mean(sicpmVals)) * 100;
                normCVByCondition[cond] = normCVByCondition[cond] || {};
                normCVByCondition[cond].sicpm = sicpmCV;
            }
        });
        
        if (Object.keys(normCVByCondition).length > 0) {
            const condLabels = Object.keys(normCVByCondition);
            const cpmCVs = condLabels.map(c => normCVByCondition[c].cpm || 0);
            const sicpmCVs = condLabels.map(c => normCVByCondition[c].sicpm || 0);
            
            // Show both CPM and siCPM CVs as grouped bars
            if (cpmCVs.some(v => v > 0) || sicpmCVs.some(v => v > 0)) {
                // For simplicity, show CPM CV (can be enhanced to show both)
                const cvData = condLabels.map(c => normCVByCondition[c].cpm || normCVByCondition[c].sicpm || 0);
                renderBarChart('chart-norm-cv', cvData, condLabels);
            } else {
                document.getElementById('chart-norm-cv').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">Insufficient data for CV calculation</div>';
            }
        } else {
            document.getElementById('chart-norm-cv').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No normalization data available</div>';
        }
        
        // PI vs depth scatter plot
        const piVsDepthPairs = rows.filter(r => r.median_pausing_index != null && r.input_reads != null && isFinite(r.median_pausing_index) && r.input_reads > 0);
        if (piVsDepthPairs.length > 0) {
            const depths = piVsDepthPairs.map(r => (r.input_reads || 0) / 1e6); // Convert to millions
            const pis = piVsDepthPairs.map(r => r.median_pausing_index);
            const sampleLabels = piVsDepthPairs.map(r => r.sample_id);
            renderScatterPlot('chart-pi-vs-depth', depths, pis, 'Read Depth (M)', 'Median Pausing Index', sampleLabels);
        } else {
            document.getElementById('chart-pi-vs-depth').innerHTML = '<div style="padding:2rem;text-align:center;color:var(--muted)">No pausing index or depth data available</div>';
        }
        
        // ===== SAMPLE TABLE =====
        const tbody = document.getElementById('sample-table-body');
        const searchInput = document.getElementById('sample-search');
        const condFilter = document.getElementById('condition-filter');
        
        // Populate condition filter
        conditions.forEach(c => {
            const opt = document.createElement('option');
            opt.value = c;
            opt.textContent = c;
            condFilter.appendChild(opt);
        });
        
        function renderTable(data) {
            tbody.innerHTML = data.map(row => {
                // Format duplication percentage
                let dupDisplay = '-';
                if (row.umi_deduplication_enabled === true) {
                    dupDisplay = row.umi_deduplication_percent != null ? row.umi_deduplication_percent.toFixed(1) + '%' : '-';
                } else if (row.duplicate_percent != null) {
                    dupDisplay = row.duplicate_percent.toFixed(1) + '%';
                } else {
                    dupDisplay = 'No deduplication done';
                }
                
                return `<tr>
                    <td><strong>${row.sample_id}</strong></td>
                    <td>${row.condition || '-'}</td>
                    <td>${row.timepoint || '-'}</td>
                    <td>${row.replicate || '-'}</td>
                    <td>${((row.input_reads || 0) / 1e6).toFixed(1)}M</td>
                    <td>${formatNumber(row.reads_total_functional || 0)}</td>
                    <td>${dupDisplay}</td>
                    <td>${((row.unlocalized_fraction || 0) * 100).toFixed(1)}%</td>
                    <td>${(row.divergent_regions || 0).toLocaleString()}</td>
                    <td>${(row.total_regions || 0).toLocaleString()}</td>
                    <td>${row.median_pausing_index ? row.median_pausing_index.toFixed(2) : '-'}</td>
                    <td>${row.median_density ? row.median_density.toFixed(2) : '-'}</td>
                    <td>${row.cpm_factor ? row.cpm_factor.toFixed(2) : '-'}</td>
                    <td>${row.crpmsi_factor ? row.crpmsi_factor.toFixed(2) : '-'}</td>
                </tr>`;
            }).join('');
        }
        
        function filterTable() {
            const term = searchInput.value.toLowerCase();
            const cond = condFilter.value;
            
            const filtered = rows.filter(r => {
                const matchesSearch = (r.sample_id + (r.condition || '') + (r.timepoint || '')).toLowerCase().includes(term);
                const matchesCond = !cond || r.condition === cond;
                return matchesSearch && matchesCond;
            });
            renderTable(filtered);
        }
        
        searchInput.addEventListener('input', filterTable);
        condFilter.addEventListener('change', filterTable);
        
        // Export CSV
        document.getElementById('export-csv').addEventListener('click', () => {
            const headers = ['Sample','Condition','Timepoint','Replicate','Input_Reads','Functional_Reads','Dup_Pct','Unloc_Pct','Div_Regions','Total_Regions','Median_PI','Median_Density','CPM_Factor','siCPM_Factor'];
            const csv = [headers.join(',')].concat(rows.map(r => [
                r.sample_id, r.condition || '', r.timepoint || '', r.replicate || '',
                r.input_reads || 0, r.reads_total_functional || 0,
                (r.duplicate_percent || 0).toFixed(2),
                ((r.unlocalized_fraction || 0) * 100).toFixed(2),
                r.divergent_regions || 0, r.total_regions || 0,
                r.median_pausing_index || '', r.median_density || '',
                r.cpm_factor || '', r.crpmsi_factor || ''
            ].join(','))).join('\\n');
            const blob = new Blob([csv], { type: 'text/csv' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = 'tracktx_cohort_metrics.csv';
            a.click();
        });
        
        renderTable(rows);
        
        // Smooth scroll for navigation
        document.querySelectorAll('.nav-pills a').forEach(link => {
            link.addEventListener('click', (e) => {
                e.preventDefault();
                const targetId = link.getAttribute('href').substring(1);
                const target = document.getElementById(targetId);
                if (target) {
                    target.scrollIntoView({ behavior: 'smooth', block: 'start' });
                }
            });
        });
        
        console.log('Report fully loaded');
    });
    </script>
    """
    
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