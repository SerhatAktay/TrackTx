#!/usr/bin/env python3
# =============================================================================
# calculate_pol2_metrics.py — Per-Sample Pol-II Pausing and Gene Metrics
# =============================================================================
#
# Purpose:
#   Calculates RNA Polymerase II pausing index and gene-level metrics from
#   aligned BAM files and gene annotations.
#
# Features:
#   • TSS window counting (promoter occupancy)
#   • Gene body counting (elongation)
#   • Pausing index calculation (TSS/body ratio)
#   • CPM normalization
#   • Per-base-pair density calculation
#   • Length-normalized pausing index
#   • Memory-efficient processing with pysam
#
# Method:
#   1. Parse GTF to extract gene coordinates
#   2. Build TSS windows (TSS ± window_size)
#   3. Build gene body regions (TSS + offset to gene end)
#   4. Count reads in each region using pysam or bedtools
#   5. Calculate pausing indices and densities
#   6. Generate QC summary
#
# Outputs:
#   • pausing_index.tsv: Lean pausing index table
#   • pol2_gene_metrics.tsv: Comprehensive gene metrics
#   • pol2_qc.json: Quality control summary
#
# =============================================================================

from __future__ import annotations
import argparse
import csv
import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# =============================================================================
# CONSTANTS
# =============================================================================

VERSION = "2.0.0"
LOG_PREFIX = "[POL2_CALC]"

# =============================================================================
# LOGGING UTILITIES
# =============================================================================

def log(section: str, message: str, flush: bool = True):
    """Consistent logging format: [POL2_CALC] SECTION | message"""
    timestamp = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    print(f"{LOG_PREFIX} {section} | {message} | ts={timestamp}", flush=flush)

def log_info(message: str):
    """Log informational message"""
    print(f"{LOG_PREFIX} INFO | {message}", flush=True)

def log_error(message: str):
    """Log error message"""
    print(f"{LOG_PREFIX} ERROR | {message}", file=sys.stderr, flush=True)

def log_warning(message: str):
    """Log warning message"""
    print(f"{LOG_PREFIX} WARNING | {message}", flush=True)

def log_progress(section: str, current: int, total: int):
    """Log progress indicator"""
    if total > 0:
        percent = (current / total) * 100
        print(f"{LOG_PREFIX} {section} | Progress: {current}/{total} ({percent:.1f}%)", flush=True)

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def run_command(cmd: List[str], description: str = "") -> subprocess.CompletedProcess:
    """
    Run shell command with error handling
    
    Args:
        cmd: Command and arguments as list
        description: Description for logging
        
    Returns:
        CompletedProcess with stdout/stderr
        
    Raises:
        SystemExit on command failure
    """
    try:
        if description:
            log_info(f"Running: {description}")
        
        result = subprocess.run(
            cmd,
            text=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        return result
    except subprocess.CalledProcessError as e:
        log_error(f"Command failed: {' '.join(cmd)}")
        log_error(f"Exit code: {e.returncode}")
        log_error(f"Stderr: {e.stderr}")
        raise SystemExit(1)
    except FileNotFoundError:
        log_error(f"Command not found: {cmd[0]}")
        log_error("Please ensure the tool is installed and in PATH")
        raise SystemExit(1)

def is_gzipped(path: str) -> bool:
    """
    Check if file is gzipped
    
    Args:
        path: File path to check
        
    Returns:
        True if gzipped, False otherwise
    """
    if str(path).endswith(".gz"):
        return True
    
    try:
        import gzip
        with gzip.open(path, "rb") as f:
            f.read(1)
        return True
    except Exception:
        return False

def open_text_file(path: str):
    """
    Open text file, handling gzip compression
    
    Args:
        path: File path to open
        
    Returns:
        File handle (text mode)
    """
    if is_gzipped(path):
        import gzip
        import io
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    else:
        return open(path, "r", encoding="utf-8", errors="replace")

def parse_gtf_attributes(attr_string: str) -> Dict[str, str]:
    """
    Parse GTF attribute column
    
    Handles both GTF format (key "value") and GFF3 format (key=value)
    
    Args:
        attr_string: Attribute column string
        
    Returns:
        Dictionary of parsed attributes
    """
    attrs = {}
    
    # Parse common attributes
    for key in ["gene_id", "gene_name", "ID", "Name"]:
        # Try GTF format: key "value"
        match = re.search(rf'{key}\s+"([^"]+)"', attr_string)
        if not match:
            # Try GFF3 format: key=value
            match = re.search(rf'{key}=([^;]+)', attr_string)
        
        if match:
            attrs[key] = match.group(1)
    
    # Set defaults
    attrs.setdefault("gene_id", attrs.get("ID", "NA"))
    attrs.setdefault("gene_name", attrs.get("Name", attrs.get("gene_id", "NA")))
    
    return attrs

# =============================================================================
# READ COUNTING
# =============================================================================

def get_mapped_read_count(bam_path: str) -> int:
    """
    Get total mapped read count from BAM
    
    Uses samtools idxstats for efficiency, falls back to samtools view -c
    
    Args:
        bam_path: Path to BAM file
        
    Returns:
        Number of mapped reads
    """
    log_info("Counting total mapped reads...")
    
    try:
        # Primary method: samtools idxstats (faster)
        result = run_command(
            ["samtools", "idxstats", bam_path],
            "Extracting read counts from BAM index"
        )
        
        total = 0
        for line in result.stdout.splitlines():
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[0] != "*":
                try:
                    total += int(parts[2])
                except ValueError:
                    continue
        
        log_info(f"Total mapped reads: {total:,}")
        return total
        
    except Exception:
        # Fallback method: samtools view -c (slower but more reliable)
        log_warning("idxstats failed, using samtools view -c")
        result = run_command(
            ["samtools", "view", "-c", "-F", "4", bam_path],
            "Counting mapped reads"
        )
        total = int(result.stdout.strip())
        log_info(f"Total mapped reads: {total:,}")
        return total

def count_reads_pysam(bed_path: Path, bam_path: str, region_type: str) -> Dict[str, int]:
    """
    Count reads in BED regions using pysam (memory efficient)
    
    Args:
        bed_path: Path to BED file with regions
        region_type: Description of regions (for logging)
        bam_path: Path to BAM file
        
    Returns:
        Dictionary mapping gene_id to read count
    """
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        log_warning(f"Empty {region_type} BED file")
        return {}
    
    # Try pysam first (preferred)
    try:
        import pysam
    except ImportError:
        log_warning("pysam not available, using bedtools intersect")
        return count_reads_bedtools(bed_path, bam_path, region_type)
    
    log_info(f"Counting {region_type} reads with pysam...")
    
    try:
        bamfile = pysam.AlignmentFile(bam_path, "rb")
        counts = {}
        
        with open(bed_path) as f:
            regions = [line for line in f if line.strip() and not line.startswith(("#", "track", "browser"))]
        
        total_regions = len(regions)
        log_info(f"Processing {total_regions:,} {region_type} regions...")
        
        for i, line in enumerate(regions, 1):
            fields = line.strip().split("\t")
            if len(fields) < 4:
                continue
            
            chrom, start, end, gene_id = fields[0], int(fields[1]), int(fields[2]), fields[3]
            
            # Progress indicator every 5000 regions
            if i % 5000 == 0:
                log_progress(region_type.upper(), i, total_regions)
            
            try:
                count = bamfile.count(contig=chrom, start=start, stop=end)
                counts[gene_id] = counts.get(gene_id, 0) + count
            except Exception:
                # Chromosome not in BAM, skip silently
                continue
        
        bamfile.close()
        log_info(f"Counted reads for {len(counts):,} genes in {region_type}")
        return counts
        
    except Exception as e:
        log_error(f"pysam counting failed: {e}")
        log_info("Falling back to bedtools intersect")
        return count_reads_bedtools(bed_path, bam_path, region_type)

def count_reads_bedtools(bed_path: Path, bam_path: str, region_type: str) -> Dict[str, int]:
    """
    Count reads using bedtools intersect (fallback method)
    
    Args:
        bed_path: Path to BED file
        bam_path: Path to BAM file
        region_type: Description of regions
        
    Returns:
        Dictionary mapping gene_id to read count
    """
    log_info(f"Counting {region_type} reads with bedtools...")
    
    result = run_command(
        ["bedtools", "intersect", "-c", "-a", str(bed_path), "-b", bam_path],
        f"Counting {region_type} overlaps"
    )
    
    counts = {}
    for line in result.stdout.splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        
        fields = line.split("\t")
        if len(fields) < 7:
            continue
        
        gene_id = fields[3]
        count = int(fields[-1])
        counts[gene_id] = counts.get(gene_id, 0) + count
    
    log_info(f"Counted reads for {len(counts):,} genes in {region_type}")
    return counts

# =============================================================================
# GTF PARSING
# =============================================================================

def parse_gtf_file(
    gtf_path: str,
    feature_types: set,
    tss_window: int,
    body_offset_min: int,
    body_offset_frac: float
) -> List[Tuple]:
    """
    Parse GTF file and extract gene coordinates
    
    Args:
        gtf_path: Path to GTF file
        feature_types: Set of acceptable feature types
        tss_window: TSS window size (±bp)
        body_offset_min: Minimum body offset (bp)
        body_offset_frac: Body offset as fraction of gene length
        
    Returns:
        List of tuples: (gene_id, gene_name, chrom, strand, 
                        tss_lo, tss_hi, body_lo, body_hi, body_len)
    """
    log("PARSE", f"Reading GTF: {gtf_path}")
    
    # Store transcript data per gene
    gene_data = {}  # gene_id -> {gname, chrom, strand, transcripts[]}
    
    line_count = 0
    feature_count = 0
    
    with open_text_file(gtf_path) as f:
        for line in f:
            line_count += 1
            
            # Progress indicator
            if line_count % 100000 == 0:
                log_info(f"Parsed {line_count:,} GTF lines...")
            
            # Skip comments and empty lines
            if not line.strip() or line.startswith("#"):
                continue
            
            # Parse line
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            
            # Check feature type
            feature_type = fields[2]
            if feature_type not in feature_types:
                continue
            
            feature_count += 1
            
            # Extract coordinates
            chrom = fields[0]
            strand = fields[6]
            
            try:
                start = int(fields[3])
                end = int(fields[4])
            except ValueError:
                log_warning(f"Invalid coordinates at line {line_count}")
                continue
            
            if end <= start:
                log_warning(f"Invalid region (end <= start) at line {line_count}")
                continue
            
            # Parse attributes
            attrs = parse_gtf_attributes(fields[8])
            gene_id = attrs["gene_id"]
            gene_name = attrs["gene_name"]
            
            # Calculate TSS and body coordinates
            tss = start if strand == "+" else end
            gene_length = max(1, end - start)
            offset = max(body_offset_min, int(gene_length * body_offset_frac))
            
            if strand == "+":
                body_lo = tss + offset
                body_hi = end
            else:
                body_lo = start
                body_hi = tss - offset
            
            body_lo = max(0, body_lo)
            body_len = max(0, body_hi - body_lo)
            
            tss_lo = max(0, tss - tss_window)
            tss_hi = tss + tss_window
            
            # Store transcript
            if gene_id not in gene_data:
                gene_data[gene_id] = {
                    'gname': gene_name,
                    'chrom': chrom,
                    'strand': strand,
                    'transcripts': []
                }
            
            gene_data[gene_id]['transcripts'].append((tss_lo, tss_hi, body_lo, body_hi, body_len))
    
    log("PARSE", f"Processed {line_count:,} GTF lines")
    log("PARSE", f"Found {feature_count:,} matching features")
    
    # Aggregate transcripts per gene
    # Use LONGEST transcript per gene (by body length) instead of union to avoid
    # huge bogus spans from genes with dispersed transcripts (e.g. chrY PAR).
    # Union of distant transcripts produced 90+ Mb TSS windows and wrong PIs.
    log("PARSE", "Aggregating transcript coordinates per gene (longest transcript)...")
    genes = []
    max_tss_span = 1000   # Reject genes with TSS window > 1 kb (indicates bad aggregation)
    max_body_span = 500_000  # Cap body at 500 kb; longer suggests multi-locus gene
    
    for gene_id, data in gene_data.items():
        transcripts = data['transcripts']
        if not transcripts:
            continue

        # Pick transcript with longest body (most representative for pausing)
        best = max(transcripts, key=lambda t: t[4])  # t[4] = body_len
        tss_lo, tss_hi, body_lo, body_hi, body_len = best
        
        # Sanity: reject genes with bogus TSS span (should be ~2*tss_window)
        tss_span = tss_hi - tss_lo
        if tss_span > max_tss_span:
            log_warning(f"Gene {gene_id}: TSS span {tss_span} bp > {max_tss_span}, skipping")
            continue
        
        # Cap body length to avoid multi-locus genes inflating body counts
        if body_len > max_body_span:
            log_warning(f"Gene {gene_id}: body {body_len} bp truncated to {max_body_span}")
            if data['strand'] == "+":
                body_hi = body_lo + max_body_span
            else:
                body_lo = max(0, body_hi - max_body_span)
            body_len = body_hi - body_lo  # Ensure consistency after coord adjustment
        else:
            body_len = body_hi - body_lo  # Recompute in case of float rounding

        genes.append((
            gene_id,
            data['gname'],
            data['chrom'],
            data['strand'],
            tss_lo,
            tss_hi,
            body_lo,
            body_hi,
            body_len
        ))
    
    log("PARSE", f"Extracted {len(genes):,} unique genes")
    return genes

# =============================================================================
# BED FILE OPERATIONS
# =============================================================================

def write_bed_files(
    genes: List[Tuple],
    tss_bed_path: Path,
    body_bed_path: Path
):
    """
    Write TSS and body BED files
    
    Args:
        genes: List of gene tuples
        tss_bed_path: Output path for TSS BED
        body_bed_path: Output path for body BED
    """
    log("BED", "Writing BED files...")
    
    tss_count = 0
    body_count = 0
    
    with open(tss_bed_path, "w") as tss_f, open(body_bed_path, "w") as body_f:
        for (gene_id, gene_name, chrom, strand, tss_lo, tss_hi, 
             body_lo, body_hi, body_len) in genes:
            
            # TSS window (always write)
            tss_f.write(f"{chrom}\t{tss_lo}\t{tss_hi}\t{gene_id}\t0\t{strand}\n")
            tss_count += 1
            
            # Body region (only if length > 0)
            if body_len > 0:
                body_f.write(f"{chrom}\t{body_lo}\t{body_hi}\t{gene_id}\t0\t{strand}\n")
                body_count += 1
    
    log("BED", f"TSS regions: {tss_count:,}")
    log("BED", f"Body regions: {body_count:,}")

def extract_genome_file_from_bam(bam_path: str, output_path: Path):
    """
    Extract chromosome names and sizes from BAM header
    
    Args:
        bam_path: Path to BAM file
        output_path: Output path for genome file
    """
    log("GENOME", "Extracting chromosome info from BAM...")
    
    result = run_command(
        ["samtools", "idxstats", bam_path],
        "Reading BAM chromosome info"
    )
    
    lines = []
    for line in result.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0] != "*":
            lines.append(f"{parts[0]}\t{parts[1]}")
    
    output_path.write_text("\n".join(lines) + "\n")
    log("GENOME", f"Wrote {len(lines)} chromosomes to genome file")

def sort_bed_file(bed_path: Path, genome_path: Path):
    """
    Sort BED file by genomic coordinates
    
    Args:
        bed_path: BED file to sort
        genome_path: Genome file for chromosome order
    """
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        return
    
    log("SORT", f"Sorting {bed_path.name}...")
    
    try:
        # Try bedtools sort with genome file
        result = subprocess.run(
            ["bedtools", "sort", "-g", str(genome_path), "-i", str(bed_path)],
            text=True,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        bed_path.write_text(result.stdout)
        log("SORT", f"Sorted {bed_path.name} with bedtools")
        
    except (subprocess.CalledProcessError, FileNotFoundError):
        # Fallback: POSIX sort
        log_warning("bedtools sort failed, using LC_ALL=C sort")
        subprocess.run(
            ["bash", "-c", f"LC_ALL=C sort -k1,1 -k2,2n -k3,3n '{bed_path}' -o '{bed_path}'"],
            check=True
        )
        log("SORT", f"Sorted {bed_path.name} with POSIX sort")

# =============================================================================
# OUTPUT GENERATION
# =============================================================================

def write_output_files(
    genes: List[Tuple],
    tss_counts: Dict[str, int],
    body_counts: Dict[str, int],
    mapped_reads: int,
    pausing_output: str,
    genes_output: str
):
    """
    Write pausing index and gene metrics output files
    
    Args:
        genes: List of gene tuples
        tss_counts: TSS read counts per gene
        body_counts: Body read counts per gene
        mapped_reads: Total mapped reads for CPM
        pausing_output: Output path for pausing index
        genes_output: Output path for gene metrics
    """
    log("OUTPUT", "Writing output tables...")
    
    # CPM denominator
    cpm_denom = (mapped_reads / 1_000_000.0) if mapped_reads > 0 else 1e-9
    
    with open(pausing_output, "w") as p_out, open(genes_output, "w") as g_out:
        # Writers
        p_writer = csv.writer(p_out, delimiter="\t", lineterminator="\n")
        g_writer = csv.writer(g_out, delimiter="\t", lineterminator="\n")
        
        # Headers (pi_len_norm preferred for reporting - corrects for TSS vs body length)
        p_writer.writerow([
            "gene_id", "chrom", "strand",
            "tss_count", "gene_body_count",
            "pi_raw", "pi_len_norm", "is_truncated"
        ])
        
        g_writer.writerow([
            "gene_id", "gene_name", "chrom", "strand",
            "tss_lo", "tss_hi", "tss_width",
            "body_lo", "body_hi", "body_len",
            "tss_count", "tss_cpm", "tss_density_per_bp",
            "body_count", "body_cpm", "body_density_per_bp",
            "pi_raw", "pi_len_norm", "is_truncated"
        ])
        
        # Write data
        for (gene_id, gene_name, chrom, strand, tss_lo, tss_hi,
             body_lo, body_hi, body_len) in genes:
            
            # Get counts
            tss_count = int(tss_counts.get(gene_id, 0))
            body_count = int(body_counts.get(gene_id, 0))
            
            # Calculate metrics
            tss_width = max(1, tss_hi - tss_lo)
            
            # Pausing indices
            pi_raw = (tss_count / body_count) if body_count > 0 else float("nan")
            pi_len_norm = (
                (tss_count / tss_width) / (body_count / body_len)
                if (body_count > 0 and body_len > 0)
                else float("nan")
            )
            
            # Truncation flag
            is_truncated = int(body_count == 0 or body_len == 0)
            
            # CPM and densities
            tss_cpm = tss_count / cpm_denom
            body_cpm = body_count / cpm_denom
            tss_density = tss_count / max(1, tss_width)
            body_density = (body_count / body_len) if body_len > 0 else 0.0
            
            # Write pausing index (pi_len_norm = (TSS_density)/(body_density) for proper comparison)
            p_writer.writerow([
                gene_id, chrom, strand,
                tss_count, body_count,
                pi_raw, pi_len_norm, is_truncated
            ])
            
            # Write gene metrics
            g_writer.writerow([
                gene_id, gene_name, chrom, strand,
                tss_lo, tss_hi, tss_width,
                body_lo, body_hi, body_len,
                tss_count, tss_cpm, tss_density,
                body_count, body_cpm, body_density,
                pi_raw, pi_len_norm, is_truncated
            ])
    
    log("OUTPUT", f"Wrote metrics for {len(genes):,} genes")

def write_qc_json(
    output_path: str,
    mapped_reads: int,
    gene_count: int,
    tss_window: int,
    body_offset_min: int,
    body_offset_frac: float
):
    """
    Write QC JSON summary
    
    Args:
        output_path: Output JSON path
        mapped_reads: Total mapped reads
        gene_count: Number of genes processed
        tss_window: TSS window size
        body_offset_min: Body offset minimum
        body_offset_frac: Body offset fraction
    """
    qc_data = {
        "total_mapped": int(mapped_reads),
        "genes_seen": int(gene_count),
        "tss_window_bp": int(tss_window),
        "body_offset_min_bp": int(body_offset_min),
        "body_offset_frac": float(body_offset_frac),
    }
    
    Path(output_path).write_text(json.dumps(qc_data, indent=2))
    log("QC", f"Wrote QC JSON: {output_path}")

# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main execution function"""
    
    # Parse arguments
    parser = argparse.ArgumentParser(
        description="Calculate Pol-II pausing index and gene metrics",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    parser.add_argument("--bam", required=True, help="Input BAM file")
    parser.add_argument("--gtf", required=True, help="Gene annotation GTF file")
    parser.add_argument("--tss-win", type=int, default=50, 
                       help="TSS window size (±bp) [default: 50]")
    parser.add_argument("--body-offset-min", type=int, default=2000,
                       help="Minimum body offset (bp) [default: 2000]")
    parser.add_argument("--body-offset-frac", type=float, default=0.10,
                       help="Body offset fraction [default: 0.10]")
    parser.add_argument("--feature-types", default="gene,transcript",
                       help="Comma-separated feature types [default: gene,transcript]")
    parser.add_argument("--out-pausing", required=True, 
                       help="Output pausing index TSV")
    parser.add_argument("--out-genes", required=True,
                       help="Output gene metrics TSV")
    parser.add_argument("--out-qc", required=True,
                       help="Output QC JSON")
    parser.add_argument("--threads", type=int, default=1,
                       help="Number of threads [default: 1]")
    parser.add_argument("--fail-if-empty", default="false",
                       help="Fail if no genes parsed [default: false]")
    parser.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    
    args = parser.parse_args()
    
    # Start
    log("START", f"calculate_pol2_metrics.py v{VERSION}")
    log("START", f"BAM: {args.bam}")
    log("START", f"GTF: {args.gtf}")
    
    # Configuration
    log("CONFIG", f"TSS window: ±{args.tss_win} bp")
    log("CONFIG", f"Body offset min: {args.body_offset_min} bp")
    log("CONFIG", f"Body offset fraction: {args.body_offset_frac}")
    log("CONFIG", f"Feature types: {args.feature_types}")
    
    # Parse feature types
    feature_types = {ft.strip() for ft in args.feature_types.split(",") if ft.strip()}
    log("CONFIG", f"Parsed {len(feature_types)} feature types")
    
    # Create temporary directory
    tmpdir = Path(tempfile.mkdtemp(prefix=".pol2_calc_", dir=".")).resolve()
    log("SETUP", f"Temporary directory: {tmpdir}")
    
    try:
        # Define paths
        genome_file = tmpdir / "genome.tsv"
        tss_bed = tmpdir / "tss.bed"
        body_bed = tmpdir / "body.bed"
        
        # Parse GTF
        log("═" * 70, "")
        genes = parse_gtf_file(
            args.gtf,
            feature_types,
            args.tss_win,
            args.body_offset_min,
            args.body_offset_frac
        )
        
        # Check if genes were found
        if not genes:
            log_error("No genes parsed from GTF")
            
            if args.fail_if_empty.lower() in ("true", "1", "yes"):
                log_error("fail-if-empty is enabled, exiting with error")
                raise SystemExit(1)
            
            # Create empty outputs
            log_warning("Creating empty output files")
            Path(args.out_pausing).write_text(
                "gene_id\tchrom\tstrand\ttss_count\tgene_body_count\tpausing_index\tis_truncated\n"
            )
            Path(args.out_genes).write_text(
                "gene_id\tgene_name\tchrom\tstrand\ttss_lo\ttss_hi\ttss_width\t"
                "body_lo\tbody_hi\tbody_len\ttss_count\ttss_cpm\ttss_density_per_bp\t"
                "body_count\tbody_cpm\tbody_density_per_bp\tpi_raw\tpi_len_norm\tis_truncated\n"
            )
            write_qc_json(args.out_qc, 0, 0, args.tss_win, 
                         args.body_offset_min, args.body_offset_frac)
            log("COMPLETE", "Empty outputs created")
            return 0
        
        # Write BED files
        log("═" * 70, "")
        write_bed_files(genes, tss_bed, body_bed)
        
        # Extract genome file and sort BEDs
        log("═" * 70, "")
        extract_genome_file_from_bam(args.bam, genome_file)
        sort_bed_file(tss_bed, genome_file)
        sort_bed_file(body_bed, genome_file)
        
        # Count reads
        log("═" * 70, "")
        tss_counts = count_reads_pysam(tss_bed, args.bam, "TSS")
        
        log("═" * 70, "")
        body_counts = count_reads_pysam(body_bed, args.bam, "body")
        
        # Get total mapped reads
        log("═" * 70, "")
        mapped = get_mapped_read_count(args.bam)
        
        # Write outputs
        log("═" * 70, "")
        write_output_files(
            genes,
            tss_counts,
            body_counts,
            mapped,
            args.out_pausing,
            args.out_genes
        )
        
        # Write QC JSON
        write_qc_json(
            args.out_qc,
            mapped,
            len(genes),
            args.tss_win,
            args.body_offset_min,
            args.body_offset_frac
        )
        
        # Success
        log("═" * 70, "")
        log("COMPLETE", f"Successfully processed {len(genes):,} genes")
        log("COMPLETE", f"Total mapped reads: {mapped:,}")
        log("COMPLETE", "All outputs written")
        
    finally:
        # Cleanup
        try:
            shutil.rmtree(tmpdir, ignore_errors=True)
            log("CLEANUP", "Temporary files removed")
        except Exception as e:
            log_warning(f"Could not remove temp directory: {e}")
    
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