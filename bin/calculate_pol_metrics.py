#!/usr/bin/env python3
# =============================================================================
# calculate_pol_metrics.py — Per-Sample Pol-II Pausing and Gene Metrics
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
#   • pol_gene_metrics.tsv: Comprehensive gene metrics
#   • pol_qc.json: Quality control summary
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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import make_logger, run_main, AtomicFileWriter

# =============================================================================
# CONSTANTS
# =============================================================================

VERSION = "2.0.0"
LOG_PREFIX = "[POL_CALC]"

# =============================================================================
# LOGGING UTILITIES
# =============================================================================

log_info, log_warning, log_error, log_progress = make_logger("POL_CALC")

def log(section: str, message: str, flush: bool = True):
    """Consistent logging format: [POL_CALC] SECTION | message"""
    timestamp = datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    print(f"{LOG_PREFIX} {section} | {message} | ts={timestamp}", flush=flush)

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
    # newline="\n": only split lines on LF. With the default (universal newline)
    # mode, a stray CR mid-line — e.g. "chr6\r\tBestRefSeq..." from a GTF whose
    # chromosomes were renamed using a CRLF NCBI assembly report — would be
    # treated as a line break, leaving the data line starting with a tab and an
    # empty chromosome field. Callers strip the remaining CR per line.
    if is_gzipped(path):
        import gzip
        import io
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace", newline="\n")
    else:
        return open(path, "r", encoding="utf-8", errors="replace", newline="\n")

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

def count_reads_pysam_multi(bed_paths: Dict[str, Path], bam_path: str) -> Dict[str, Dict[str, int]]:
    """
    Count strand-specific reads for multiple BED files (e.g. TSS and body)
    against one BAM, opening the BAM and resolving contigs only once instead
    of once per region type.

    Each BED file must have a strand column (col 6). Only reads whose mapping
    strand matches the gene strand are counted, eliminating contamination from
    antisense transcription at convergently-oriented loci or nearby enhancers.

    Strand logic for PRO-seq (after RC(R1) alignment):
      Gene +  →  count reads that are NOT reverse (forward-strand reads)
      Gene -  →  count reads that ARE  reverse  (reverse-strand reads)

    Args:
        bed_paths: mapping of region_type -> BED path (region_type used for logging)
        bam_path: Path to BAM file

    Returns:
        Mapping of region_type -> {gene_id: read count}
    """
    # Try pysam first (preferred)
    try:
        import pysam
    except ImportError:
        log_warning("pysam not available, using bedtools intersect")
        return {rt: count_reads_bedtools(p, bam_path, rt) for rt, p in bed_paths.items()}

    # Verify BAM index exists — pysam.fetch() needs an index; without it every
    # region returns 0 with no error raised.
    bam_p = Path(bam_path)
    if not (bam_p.with_suffix(".bai").exists() or Path(bam_path + ".bai").exists()):
        log_warning(f"BAM index (.bai) not found for {bam_path} — falling back to bedtools")
        return {rt: count_reads_bedtools(p, bam_path, rt) for rt, p in bed_paths.items()}

    try:
        bamfile = pysam.AlignmentFile(bam_path, "rb")

        # Resolve BED chrom names against what the BAM header actually has.
        # A gene catalog and its BAM can disagree on contig naming (e.g. the
        # catalog carries an alt-haplotype/patch contig, or one file uses
        # "chr1" and the other "1") even when the underlying assembly is the
        # same. Building this set once (shared across all region types, not
        # rebuilt per BED file) and normalizing per-region avoids spurious
        # fetch failures for every affected region up front, rather than only
        # discovering the mismatch region-by-region below.
        bam_contigs = set(bamfile.references)

        def _resolve_contig(chrom: str):
            if chrom in bam_contigs:
                return chrom
            alt = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            if alt in bam_contigs:
                return alt
            return None

        results: Dict[str, Dict[str, int]] = {}
        for region_type, bed_path in bed_paths.items():
            results[region_type] = _count_reads_from_open_bam(bamfile, _resolve_contig, bed_path, region_type)

        bamfile.close()
        return results

    except Exception as e:
        log_error(f"pysam counting failed: {e}")
        log_info("Falling back to bedtools intersect")
        return {rt: count_reads_bedtools(p, bam_path, rt) for rt, p in bed_paths.items()}


def _count_reads_from_open_bam(bamfile, _resolve_contig, bed_path: Path, region_type: str) -> Dict[str, int]:
    """Count strand-specific, 3'-end-based reads for one BED file against an already-open BAM."""
    if not bed_path.exists() or bed_path.stat().st_size == 0:
        log_warning(f"Empty {region_type} BED file")
        return {}

    log_info(f"Counting {region_type} reads with pysam (strand-specific)...")

    counts = {}
    with open(bed_path) as f:
        regions = [line for line in f if line.strip() and not line.startswith(("#", "track", "browser"))]

    total_regions = len(regions)
    log_info(f"Processing {total_regions:,} {region_type} regions (strand-aware)...")

    # Every region that fails to count is tracked here with WHY, and
    # every such gene still gets an explicit 0 in `counts` (never a
    # silently-missing key) so downstream code can't mistake "no data
    # produced" for "gene genuinely has zero reads" without a trace.
    failed_regions: list = []
    for i, line in enumerate(regions, 1):
        fields = line.replace("\r", "").rstrip("\n").split("\t")
        if len(fields) < 4:
            continue

        # Guard coordinate parsing: a single malformed line must not abort
        # the whole (memory-efficient) pysam path into the bedtools fallback,
        # which loads the entire BAM and can OOM under process concurrency.
        try:
            chrom    = fields[0]
            start    = int(fields[1])
            end      = int(fields[2])
        except (ValueError, IndexError):
            log_warning(f"Skipping malformed {region_type} BED line: {line.rstrip()!r}")
            continue
        gene_id  = fields[3]
        # BED strand column is col 6 (index 5); default to '+' if absent
        strand   = fields[5] if len(fields) >= 6 else "+"

        # Progress indicator every 5000 regions
        if i % 5000 == 0:
            log_progress(region_type.upper(), i, total_regions)

        counts.setdefault(gene_id, 0)

        fetch_chrom = _resolve_contig(chrom)
        if fetch_chrom is None:
            failed_regions.append((gene_id, chrom, "contig_not_in_bam"))
            continue

        try:
            count = 0
            for read in bamfile.fetch(contig=fetch_chrom, start=start, stop=end):
                # Skip unmapped, secondary, and supplementary
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue
                # Strand match: PRO-seq RC(R1) alignment
                #   + strand gene → read maps to forward strand (not reverse)
                #   - strand gene → read maps to reverse strand
                if strand == "+" and not read.is_reverse:
                    pass
                elif strand == "-" and read.is_reverse:
                    pass
                else:
                    continue
                # Count by the Pol II active-site (3'-end) position, not
                # full-read overlap -- matches the -3 convention used for
                # every coverage track/TSN call elsewhere in the pipeline.
                # fetch() returns any read whose alignment span overlaps
                # [start, end), so a read entering/exiting the window
                # without its 3' end inside it must be excluded here,
                # otherwise reads are double-counted or misattributed
                # across the TSS/body boundary.
                three_prime = (read.reference_end - 1) if not read.is_reverse else read.reference_start
                if start <= three_prime < end:
                    count += 1
            counts[gene_id] += count
        except (ValueError, OverflowError) as e:
            # pysam raises these for a region pysam itself considers
            # invalid (bad coordinates, contig present but out of range,
            # etc.) — a per-region data problem, safe to record as a
            # failure and continue. Anything else is NOT caught here:
            # an unexpected exception type propagates up to the outer
            # handler instead of being silently absorbed as a "skip".
            failed_regions.append((gene_id, chrom, str(e)))
            continue

    log_info(f"Counted reads for {len(counts):,} genes in {region_type}")
    if failed_regions:
        by_reason: Dict[str, int] = {}
        for _gid, _chrom, reason in failed_regions:
            by_reason[reason] = by_reason.get(reason, 0) + 1
        log_warning(
            f"{len(failed_regions):,}/{total_regions:,} {region_type} regions "
            f"could not be fetched and were counted as 0 — breakdown: {by_reason}. "
            f"First 10 affected gene_ids: {[g for g, _, _ in failed_regions[:10]]}"
        )
    return counts

def count_reads_bedtools(bed_path: Path, bam_path: str, region_type: str) -> Dict[str, int]:
    """
    Count strand-specific reads using bedtools intersect (fallback method,
    only reached when pysam is unavailable or the BAM lacks a .bai index).

    Counts by each read's 3'-end (Pol II active-site) position, matching the
    pysam path above -- not full-read overlap, which would double-count or
    misattribute reads spanning a region boundary (e.g. the TSS/body split).

    Args:
        bed_path: Path to BED file (must include strand col 6)
        bam_path: Path to BAM file
        region_type: Description of regions

    Returns:
        Dictionary mapping gene_id to read count
    """
    log_info(f"Counting {region_type} reads with bedtools (strand-specific, 3'-end)...")

    # Collapse each primary alignment (-F 0x904 drops unmapped/secondary/
    # supplementary, matching is_unmapped/is_secondary/is_supplementary in
    # the pysam path) to its single-base 3'-end position:
    #   forward read (bamtobed strand '+') -> end-1 (0-based, last covered base)
    #   reverse read (bamtobed strand '-') -> start
    slug = re.sub(r"[^A-Za-z0-9_]", "_", region_type.lower())
    threeprime_bed = bed_path.parent / f"{slug}_3prime_reads.bed"
    awk_3prime = (
        r"""awk -v OFS='\t' '{ if ($6 == "+") print $1, $3-1, $3, $4, $5, $6; """
        r"""else print $1, $2, $2+1, $4, $5, $6 }'"""
    )
    cmd = (
        f"samtools view -b -F 0x904 '{bam_path}' | "
        f"bedtools bamtobed -i stdin | "
        f"{awk_3prime} > '{threeprime_bed}'"
    )
    run_command(["bash", "-c", cmd], f"Converting {region_type} reads to 3'-end positions")

    result = run_command(
        ["bedtools", "intersect", "-c", "-s", "-a", str(bed_path), "-b", str(threeprime_bed)],
        f"Counting {region_type} overlaps (strand-specific, 3'-end)"
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

def clip_windows_to_neighbors(genes: List[Tuple]) -> List[Tuple]:
    """
    Clip each gene's TSS/body window at the overlap midpoint with its nearest
    same-chromosome, same-strand neighbor.

    In compact genomes with short intergenic distances (bacterial operons,
    dense plant/insect genomes), one gene's body/TSS window can otherwise
    extend into a neighboring gene's own TSS/body region, letting antisense/
    adjacent-gene signal bleed into the pausing index. Cross-strand neighbors
    (divergent/convergent pairs) are left alone -- that's handled elsewhere
    in the pipeline (module 09).
    """
    if len(genes) < 2:
        return genes

    groups: Dict[Tuple[str, str], List[int]] = {}
    for idx, g in enumerate(genes):
        groups.setdefault((g[2], g[3]), []).append(idx)  # (chrom, strand)

    genes = list(genes)
    clipped = 0
    for (_chrom, _strand), idxs in groups.items():
        if len(idxs) < 2:
            continue
        # Sort by each gene's leftmost occupied coordinate (min of its own
        # TSS/body window bounds -- robust regardless of strand, since
        # body_lo<=body_hi and tss_lo<=tss_hi are already guaranteed above).
        idxs.sort(key=lambda i: min(genes[i][4], genes[i][6]))
        for pos in range(len(idxs) - 1):
            li, ri = idxs[pos], idxs[pos + 1]
            l = genes[li]
            r = genes[ri]
            l_full_hi = max(l[5], l[6], l[7])   # max(tss_hi, body_lo, body_hi)
            r_full_lo = min(r[4], r[6], r[7])   # min(tss_lo, body_lo, body_hi)
            if l_full_hi <= r_full_lo:
                continue  # no overlap, nothing to clip

            boundary = (l_full_hi + r_full_lo) // 2

            gene_id, gname, gchrom, gstrand, tss_lo, tss_hi, body_lo, body_hi, body_len, glen = l
            new_tss_hi  = min(tss_hi, max(tss_lo, boundary))
            new_body_hi = min(body_hi, max(body_lo, boundary))
            if new_tss_hi != tss_hi or new_body_hi != body_hi:
                genes[li] = (gene_id, gname, gchrom, gstrand, tss_lo, new_tss_hi,
                             body_lo, new_body_hi, max(0, new_body_hi - body_lo), glen)
                clipped += 1

            gene_id, gname, gchrom, gstrand, tss_lo, tss_hi, body_lo, body_hi, body_len, glen = r
            new_tss_lo  = max(tss_lo, min(tss_hi, boundary))
            new_body_lo = max(body_lo, min(body_hi, boundary))
            if new_tss_lo != tss_lo or new_body_lo != body_lo:
                genes[ri] = (gene_id, gname, gchrom, gstrand, new_tss_lo, tss_hi,
                             new_body_lo, body_hi, max(0, body_hi - new_body_lo), glen)
                clipped += 1

    if clipped:
        log_info(f"Clipped {clipped:,} gene TSS/body windows overlapping a same-strand neighbor")
    return genes


def auto_body_offset(
    gene_lengths: List[int],
    user_offset_min: int,
    body_offset_frac: float
) -> int:
    """
    Choose an organism-aware body offset from the gene-length distribution.

    Strategy (median-gated -- only shrink for genuinely compact genomes):
      - Look at the MEDIAN gene length across parsed features.
      - If median >= COMPACT_MEDIAN_BP (mammalian-scale, e.g. human/mouse),
        keep the user's --body-offset-min unchanged.
      - For compact genomes (small median), scale the offset to 20% of P25
        with a 200 bp hard floor, but never above the user value.
      - The caller still applies body_offset_frac on top of this minimum.

    See CHANGELOG.md for why this is median-gated rather than a flat
    P25-based shrink, and how COMPACT_MEDIAN_BP was chosen.

    Args:
        gene_lengths: List of gene lengths (bp) from the GTF
        user_offset_min: Value passed via --body-offset-min (default 2000)
        body_offset_frac: Fraction of gene length also used as offset floor

    Returns:
        Effective body_offset_min (int bp)
    """
    if not gene_lengths:
        return user_offset_min

    # Genomes whose median gene length is at least this are treated as "large"
    # and keep the user offset unchanged. This is the median over ALL genes in
    # the catalog, which is pulled down by the many small ncRNAs (observed
    # ~5.9 kb for human/T2T RefSeq, not the ~24 kb protein-coding median);
    # truly compact genomes (Drosophila / C. elegans) sit at ~2-3 kb. 4 kb
    # keeps human/mouse (~5-6 kb) at the user offset while still shrinking for
    # fly/worm. body_offset_min is overridable per-organism if a particular
    # annotation's all-gene median falls on the wrong side of this heuristic.
    COMPACT_MEDIAN_BP = 4_000

    arr = sorted(gene_lengths)
    p25_idx = max(0, int(len(arr) * 0.25) - 1)
    p25 = arr[p25_idx]
    median_idx = len(arr) // 2
    median = arr[median_idx]

    if median >= COMPACT_MEDIAN_BP:
        log_info(
            f"Gene length distribution: P25={p25:,} bp, median={median:,} bp "
            f"(>= {COMPACT_MEDIAN_BP:,} bp) → large genome; keeping user "
            f"body_offset_min={user_offset_min:,} bp"
        )
        return user_offset_min

    # Compact genome: scale offset to 20% of P25 (>=200 bp), never above user.
    auto_min = max(200, int(p25 * 0.20))
    effective = min(user_offset_min, auto_min)
    log_info(
        f"Gene length distribution: P25={p25:,} bp, median={median:,} bp "
        f"(< {COMPACT_MEDIAN_BP:,} bp) → compact genome; auto "
        f"body_offset_min={auto_min:,} bp "
        f"(user requested {user_offset_min:,}; using {effective:,})"
    )
    return effective


def parse_gtf_file(
    gtf_path: str,
    feature_types: set,
    tss_window: int,
    body_offset_min: int,
    body_offset_frac: float
) -> List[Tuple]:
    """
    Parse GTF file and extract gene coordinates.

    Body offset is auto-calibrated from the gene-length distribution so the
    pipeline works correctly for compact genomes (Drosophila, C. elegans) as
    well as human/mouse without requiring organism-specific parameter tuning.

    Args:
        gtf_path: Path to GTF file
        feature_types: Set of acceptable feature types
        tss_window: TSS window size (±bp)
        body_offset_min: Minimum body offset requested by user (bp); may be
                         reduced automatically for small-genome organisms
        body_offset_frac: Body offset as fraction of gene length

    Returns:
        List of tuples: (gene_id, gene_name, chrom, strand,
                        tss_lo, tss_hi, body_lo, body_hi, body_len, gene_length)
    """
    log("PARSE", f"Reading GTF: {gtf_path}")

    # ── Pass 1: collect raw feature records ──────────────────────────────────
    # We store raw coordinates first so we can auto-calibrate the body offset
    # from the gene-length distribution before committing to final windows.

    # gene_id -> {gname, chrom, strand, raw_transcripts: [(start,end), ...]}
    gene_data: Dict = {}
    line_count = 0
    feature_count = 0

    with open_text_file(gtf_path) as f:
        for line in f:
            line_count += 1

            if line_count % 100000 == 0:
                log_info(f"Parsed {line_count:,} GTF lines...")

            # Strip any stray CR (from CRLF-derived GTFs) so the chromosome
            # field stays clean ("chr6", not "chr6\r") and matches the BAM.
            line = line.replace("\r", "")

            if not line.strip() or line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type not in feature_types:
                continue

            feature_count += 1

            chrom  = fields[0]
            strand = fields[6]

            try:
                start = int(fields[3])
                end   = int(fields[4])
            except ValueError:
                log_warning(f"Invalid coordinates at line {line_count}")
                continue

            if end <= start:
                log_warning(f"Invalid region (end <= start) at line {line_count}")
                continue

            attrs     = parse_gtf_attributes(fields[8])
            gene_id   = attrs["gene_id"]
            gene_name = attrs["gene_name"]

            if gene_id not in gene_data:
                gene_data[gene_id] = {
                    "gname": gene_name, "chrom": chrom, "strand": strand,
                    "raw": []
                }
            gene_data[gene_id]["raw"].append((start, end))

    log("PARSE", f"Processed {line_count:,} GTF lines; {feature_count:,} features for {len(gene_data):,} genes")

    # ── Auto-calibrate body_offset_min from gene-length distribution ─────────
    all_lengths = [max(1, e - s) for d in gene_data.values() for s, e in d["raw"]]
    effective_offset_min = auto_body_offset(all_lengths, body_offset_min, body_offset_frac)

    # ── Pass 2: compute TSS/body windows with calibrated offset ──────────────
    # (Re-uses gene_data populated above; no second file read needed.)
    gene_data_with_windows: Dict = {}
    for gene_id, d in gene_data.items():
        gene_data_with_windows[gene_id] = {
            "gname": d["gname"], "chrom": d["chrom"], "strand": d["strand"],
            "transcripts": []
        }
        for start, end in d["raw"]:
            strand = d["strand"]
            tss         = start if strand == "+" else end
            gene_length = max(1, end - start)
            offset      = max(effective_offset_min, int(gene_length * body_offset_frac))

            if strand == "+":
                body_lo = min(tss + offset, end)
                body_hi = end
            else:
                body_lo = start
                body_hi = max(start, tss - offset)

            body_lo  = max(0, body_lo)
            body_len = max(0, body_hi - body_lo)
            tss_lo   = max(0, tss - tss_window)
            tss_hi   = tss + tss_window

            gene_data_with_windows[gene_id]["transcripts"].append(
                (tss_lo, tss_hi, body_lo, body_hi, body_len, gene_length)
            )

    # Swap gene_data to the windowed version for the aggregation step below
    gene_data = gene_data_with_windows

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
        tss_lo, tss_hi, body_lo, body_hi, body_len, gene_length = best
        
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
            body_len,
            gene_length
        ))

    log("PARSE", f"Extracted {len(genes):,} unique genes")
    genes = clip_windows_to_neighbors(genes)
    return genes, effective_offset_min


def parse_catalog_file(
    catalog_path: str,
    tss_window: int,
    body_offset_min: int,
    body_offset_frac: float
) -> List[Tuple]:
    """
    Parse the gtf_to_catalog genes.tsv (one row per gene) and build the same
    (gene_id, gene_name, chrom, strand, tss_lo, tss_hi, body_lo, body_hi,
    body_len, gene_length) tuples that parse_gtf_file produces.

    Using this catalog (instead of re-parsing the raw GTF and picking the longest
    transcript) makes the TSS/TES used for the pausing index IDENTICAL to the
    TSS/TES used by functional-region calling (module 10), which also consumes
    this catalog. Columns: gene_id, gene_name, chr, strand, start, end, tss, tes,
    biotype (tss/tes already strand-resolved by gtf_to_catalog: for '-' genes
    tss = txEnd, tes = txStart).
    """
    log("PARSE", f"Reading gene catalog: {catalog_path}")

    rows: List[Tuple] = []
    with open_text_file(catalog_path) as f:
        header = f.readline()
        cols = [c.strip().lstrip("﻿").lower() for c in header.rstrip("\n").split("\t")]

        def cidx(*names):
            for n in names:
                if n in cols:
                    return cols.index(n)
            return None

        i_id     = cidx("gene_id", "id")
        i_name   = cidx("gene_name", "name", "symbol")
        i_chr    = cidx("chr", "chrom", "chromosome", "seqname")
        i_strand = cidx("strand", "orientation")
        i_start  = cidx("start", "gene_start")
        i_end    = cidx("end", "gene_end")
        i_tss    = cidx("tss", "tx_start", "txstart")
        i_tes    = cidx("tes", "tx_end", "txend")

        if None in (i_id, i_chr, i_strand, i_start, i_end):
            log_error("Catalog missing required columns (need gene_id, chr, strand, start, end)")
            return [], 0

        for line in f:
            line = line.replace("\r", "")
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            try:
                gid    = p[i_id]
                gname  = p[i_name] if (i_name is not None and i_name < len(p) and p[i_name]) else gid
                chrom  = p[i_chr]
                strand = p[i_strand] if p[i_strand] in ("+", "-") else "+"
                gstart = int(float(p[i_start]))
                gend   = int(float(p[i_end]))
                tss    = int(float(p[i_tss])) if (i_tss is not None and p[i_tss] != "") else (gstart if strand == "+" else gend)
                tes    = int(float(p[i_tes])) if (i_tes is not None and p[i_tes] != "") else (gend if strand == "+" else gstart)
            except (ValueError, IndexError):
                continue
            if gend <= gstart:
                continue
            rows.append((gid, gname, chrom, strand, gstart, gend, tss, tes))

    if not rows:
        log_warning("No usable rows parsed from gene catalog")
        return [], 0

    all_lengths = [max(1, e - s) for (_, _, _, _, s, e, _, _) in rows]
    effective_offset_min = auto_body_offset(all_lengths, body_offset_min, body_offset_frac)

    genes: List[Tuple] = []
    max_body_span = 500_000
    for (gid, gname, chrom, strand, gstart, gend, tss, tes) in rows:
        gene_length = max(1, gend - gstart)
        offset = max(effective_offset_min, int(gene_length * body_offset_frac))

        if strand == "+":
            body_lo = min(tss + offset, tes)
            body_hi = tes
        else:
            body_lo = tes
            body_hi = max(tes, tss - offset)

        if body_hi < body_lo:
            body_lo, body_hi = body_hi, body_lo
        body_lo  = max(0, body_lo)
        body_len = max(0, body_hi - body_lo)
        tss_lo   = max(0, tss - tss_window)
        tss_hi   = tss + tss_window

        if body_len > max_body_span:
            if strand == "+":
                body_hi = body_lo + max_body_span
            else:
                body_lo = max(0, body_hi - max_body_span)
            body_len = body_hi - body_lo

        genes.append((gid, gname, chrom, strand, tss_lo, tss_hi,
                      body_lo, body_hi, body_len, gene_length))

    log("PARSE", f"Extracted {len(genes):,} genes from catalog")
    genes = clip_windows_to_neighbors(genes)
    return genes, effective_offset_min

# =============================================================================
# BED FILE OPERATIONS
# =============================================================================

def required_body_len(gene_length: int, min_body_frac: float, min_body_len: int,
                       effective_offset_min: int = 0) -> int:
    """
    Minimum gene-body window (bp) required for a stable pausing index.

    The threshold scales with gene length: the body must retain at least
    `min_body_frac` of the gene. Because the body offset is max(2000 bp,
    10% of L), the body only collapses for short genes, so a fractional
    floor targets exactly those cases and auto-scales for everything else.
    `min_body_len` is an optional absolute floor (bp) applied on top.

    Also floored by half of `effective_offset_min` -- the same per-organism
    body-offset calibration `auto_body_offset` already derives from the
    gene-length distribution (smaller for compact genomes) -- so this
    truncation threshold shares that calibration instead of being a second,
    uncoordinated heuristic that stays mammalian-scale for every organism.
    """
    return max(int(min_body_len), int(min_body_frac * gene_length), int(effective_offset_min) // 2)


def body_window_ok(body_len: int, gene_length: int, min_body_frac: float, min_body_len: int,
                    effective_offset_min: int = 0) -> bool:
    """True if the body window is non-empty and long enough for a stable density."""
    return body_len > 0 and body_len >= required_body_len(
        gene_length, min_body_frac, min_body_len, effective_offset_min
    )

def write_bed_files(
    genes: List[Tuple],
    tss_bed_path: Path,
    body_bed_path: Path,
    min_body_frac: float = 0.0,
    min_body_len: int = 1,
    effective_offset_min: int = 0
):
    """
    Write TSS and body BED files

    Args:
        genes: List of gene tuples
        tss_bed_path: Output path for TSS BED
        body_bed_path: Output path for body BED
        min_body_frac: Body must be >= this fraction of gene length
        min_body_len: Absolute body-length floor (bp)
        effective_offset_min: Per-organism body-offset calibration (see auto_body_offset)
    """
    log("BED", "Writing BED files...")

    tss_count = 0
    body_count = 0

    with AtomicFileWriter(tss_bed_path) as tss_f, AtomicFileWriter(body_bed_path) as body_f:
        for (gene_id, gene_name, chrom, strand, tss_lo, tss_hi,
             body_lo, body_hi, body_len, gene_length) in genes:

            # TSS window (always write)
            tss_f.write(f"{chrom}\t{tss_lo}\t{tss_hi}\t{gene_id}\t0\t{strand}\n")
            tss_count += 1

            # Body region (only if window is a meaningful fraction of the gene;
            # tiny windows on short genes inflate the pausing index)
            if body_window_ok(body_len, gene_length, min_body_frac, min_body_len, effective_offset_min):
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
    genes_output: str,
    min_body_frac: float = 0.0,
    min_body_len: int = 1,
    effective_offset_min: int = 0
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
        effective_offset_min: Per-organism body-offset calibration (see auto_body_offset)
    """
    log("OUTPUT", "Writing output tables...")
    
    # CPM denominator
    cpm_denom = (mapped_reads / 1_000_000.0) if mapped_reads > 0 else 1e-9
    
    with AtomicFileWriter(pausing_output) as p_out, AtomicFileWriter(genes_output) as g_out:
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
             body_lo, body_hi, body_len, gene_length) in genes:

            # Get counts
            tss_count = int(tss_counts.get(gene_id, 0))
            body_count = int(body_counts.get(gene_id, 0))
            
            # Calculate metrics
            tss_width = max(1, tss_hi - tss_lo)
            
            # A body window that is too short (or empty) relative to the gene
            # gives an unstable elongation density, so treat it as truncated
            # (NaN PI) rather than dividing by a near-zero/zero body.
            # Threshold scales with gene length (see required_body_len).
            body_ok = body_window_ok(body_len, gene_length, min_body_frac, min_body_len, effective_offset_min)

            # Pausing indices
            pi_raw = (tss_count / body_count) if (body_count > 0 and body_ok) else float("nan")
            pi_len_norm = (
                (tss_count / tss_width) / (body_count / body_len)
                if (body_count > 0 and body_ok)
                else float("nan")
            )

            # Truncation flag
            is_truncated = int(body_count == 0 or not body_ok)
            
            # CPM and densities
            tss_cpm = tss_count / cpm_denom
            body_cpm = body_count / cpm_denom
            tss_density = tss_count / max(1, tss_width)
            body_density = (body_count / body_len) if body_len > 0 else float("nan")
            
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
    body_offset_frac: float,
    effective_offset_min: int = 0
):
    """
    Write QC JSON summary

    Args:
        output_path: Output JSON path
        mapped_reads: Total mapped reads
        gene_count: Number of genes processed
        tss_window: TSS window size
        body_offset_min: Body offset minimum (as requested by user)
        body_offset_frac: Body offset fraction
        effective_offset_min: Per-organism calibrated body offset actually used
                               (see auto_body_offset) -- may be smaller than
                               body_offset_min for compact genomes
    """
    qc_data = {
        "total_mapped": int(mapped_reads),
        "genes_seen": int(gene_count),
        "tss_window_bp": int(tss_window),
        "body_offset_min_bp": int(body_offset_min),
        "body_offset_frac": float(body_offset_frac),
        "effective_body_offset_min_bp": int(effective_offset_min),
    }
    
    with AtomicFileWriter(output_path) as f:
        f.write(json.dumps(qc_data, indent=2))
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
    parser.add_argument("--gtf", required=False, default=None,
                       help="Gene annotation GTF file (used only if --genes is not given)")
    parser.add_argument("--genes", required=False, default=None,
                       help="gtf_to_catalog genes.tsv catalog. PREFERRED: makes TSS/TES "
                            "identical to functional-region calling. If given, --gtf is ignored.")
    parser.add_argument("--tss-win", type=int, default=50, 
                       help="TSS window size (±bp) [default: 50]")
    parser.add_argument("--body-offset-min", type=int, default=2000,
                       help="Minimum body offset (bp) [default: 2000]")
    parser.add_argument("--body-offset-frac", type=float, default=0.10,
                       help="Body offset fraction [default: 0.10]")
    parser.add_argument("--min-body-frac", type=float, default=0.10,
                       help="Body window must be >= this fraction of gene "
                            "length for a valid pausing index; shorter bodies "
                            "are flagged truncated and report NaN PI. Scales "
                            "with gene length [default: 0.10]")
    parser.add_argument("--min-body-len", type=int, default=0,
                       help="Optional absolute body-length floor (bp) applied "
                            "on top of --min-body-frac [default: 0 = off]")
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
    log("START", f"calculate_pol_metrics.py v{VERSION}")
    log("START", f"BAM: {args.bam}")
    log("START", f"Gene source: {args.genes or args.gtf or '(none)'}")
    
    # Configuration
    log("CONFIG", f"TSS window: ±{args.tss_win} bp")
    log("CONFIG", f"Body offset min: {args.body_offset_min} bp")
    log("CONFIG", f"Body offset fraction: {args.body_offset_frac}")
    log("CONFIG", f"Min body fraction: {args.min_body_frac} of gene length")
    log("CONFIG", f"Min body length floor: {args.min_body_len} bp")
    log("CONFIG", f"Feature types: {args.feature_types}")
    
    # Parse feature types
    feature_types = {ft.strip() for ft in args.feature_types.split(",") if ft.strip()}
    log("CONFIG", f"Parsed {len(feature_types)} feature types")
    
    # Create temporary directory
    tmpdir = Path(tempfile.mkdtemp(prefix=".pol_calc_", dir=".")).resolve()
    log("SETUP", f"Temporary directory: {tmpdir}")
    
    try:
        # Define paths
        genome_file = tmpdir / "genome.tsv"
        tss_bed = tmpdir / "tss.bed"
        body_bed = tmpdir / "body.bed"
        
        # Parse gene model. Prefer the gtf_to_catalog genes.tsv catalog (so the
        # TSS/TES match functional-region calling exactly); fall back to raw GTF.
        log("═" * 70, "")
        if args.genes:
            log("CONFIG", f"Gene source: catalog ({args.genes})")
            genes, effective_offset_min = parse_catalog_file(
                args.genes,
                args.tss_win,
                args.body_offset_min,
                args.body_offset_frac
            )
        elif args.gtf:
            log("CONFIG", f"Gene source: GTF ({args.gtf})")
            genes, effective_offset_min = parse_gtf_file(
                args.gtf,
                feature_types,
                args.tss_win,
                args.body_offset_min,
                args.body_offset_frac
            )
        else:
            log_error("No gene source provided: pass --genes (catalog) or --gtf")
            raise SystemExit(2)

        # Check if genes were found
        if not genes:
            log_error("No genes parsed from GTF")

            if args.fail_if_empty.lower() in ("true", "1", "yes"):
                log_error("fail-if-empty is enabled, exiting with error")
                raise SystemExit(1)

            # Create empty outputs
            log_warning("Creating empty output files")
            Path(args.out_pausing).write_text(
                "gene_id\tchrom\tstrand\ttss_count\tgene_body_count\tpi_raw\tpi_len_norm\tis_truncated\n"
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
        write_bed_files(genes, tss_bed, body_bed, args.min_body_frac, args.min_body_len,
                         effective_offset_min)

        # Extract genome file and sort BEDs
        log("═" * 70, "")
        extract_genome_file_from_bam(args.bam, genome_file)
        sort_bed_file(tss_bed, genome_file)
        sort_bed_file(body_bed, genome_file)
        
        # Count reads (TSS and body share one BAM open + one contig resolution)
        log("═" * 70, "")
        counts = count_reads_pysam_multi({"TSS": tss_bed, "body": body_bed}, args.bam)
        tss_counts = counts["TSS"]
        body_counts = counts["body"]

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
            args.out_genes,
            args.min_body_frac,
            args.min_body_len,
            effective_offset_min
        )

        # Write QC JSON
        write_qc_json(
            args.out_qc,
            mapped,
            len(genes),
            args.tss_win,
            args.body_offset_min,
            args.body_offset_frac,
            effective_offset_min
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
# SELF-CHECK (pure-function sanity check, runs on every invocation -- no
# test framework wired into this pipeline, so this is the cheapest way to
# catch a regression in the body-length/window-clipping logic above)
# =============================================================================

def _selftest():
    assert required_body_len(1000, 0.10, 0, effective_offset_min=0) == 100
    assert required_body_len(1000, 0.10, 0, effective_offset_min=400) == 200
    assert not body_window_ok(0, 1000, 0.10, 0)     # zero body_len always fails
    assert body_window_ok(150, 1000, 0.10, 0)
    assert not body_window_ok(50, 1000, 0.10, 0)

    # clip_windows_to_neighbors: overlapping same-strand genes get clipped;
    # a different-strand gene is left untouched.
    genes = [
        ("gA", "gA", "chr1", "+", 100, 200, 200, 1000, 800, 900),
        ("gB", "gB", "chr1", "+", 900, 1000, 1000, 2000, 1000, 1100),
        ("gC", "gC", "chr1", "-", 3000, 3100, 2500, 3000, 500, 600),
    ]
    clipped = clip_windows_to_neighbors(genes)
    a = next(g for g in clipped if g[0] == "gA")
    b = next(g for g in clipped if g[0] == "gB")
    c = next(g for g in clipped if g[0] == "gC")
    assert a[7] <= b[4], f"gA body_hi ({a[7]}) should not exceed gB tss_lo ({b[4]}) after clipping"
    assert c == genes[2], "different-strand gene must be left untouched"


# =============================================================================
# ENTRY POINT
# =============================================================================

def _run():
    _selftest()
    return main()

if __name__ == "__main__":
    sys.exit(run_main(_run, log_error))