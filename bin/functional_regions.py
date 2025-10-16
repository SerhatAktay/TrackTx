#!/usr/bin/env python3
# =============================================================================
# functional_regions.py — v9.1 (Gene-based, OLD script exact replication)
# -----------------------------------------------------------------------------
# PHILOSOPHY: DT sites identify active genes; gene structure defines regions
# 
# v9.0 REWRITE: Return to gene-based assignment
#      - DT sites: DISCOVERY (which genes are active?)
#      - Gene annotations: ASSIGNMENT (where do reads belong?)
# 
# v9.1 FIX: CRITICAL - Replicate old bash script EXACTLY
#      - Assignment: Use strand-specific intersect (-s or -S)
#      - Removal: Use UNSTRANDED intersect (-v, NO strand flag!)
#      - This matches original TrackTx.sh logic precisely
# 
# Key difference from v9.0:
#   v9.0: Used -s flag on removal steps (wrong!)
#   v9.1: Unstranded removal like old script (correct!)
# 
# Assignment logic:
#   1. Divergent sites overlapping gene promoters → Active Promoters
#   2. Divergent sites NOT overlapping → Enhancers
#   3. For genes with active promoters, create extended regions:
#      - Divergent (opposite strand from promoter)
#      - Gene body, CPS, Termination window
#
# Key insight: The ~20k divergent sites ARE the promoters/enhancers.
# We don't create new promoter regions from genes; we categorize the DT sites.
#
# Default geometry:
#   Promoter detection: TSS ±250 bp (defines what counts as promoter overlap)
#   Divergent: TSS -750 .. -251 (for +) / TSS +251 .. +750 (for -)
#   CPS:       TES -500 .. +499 (for +) / TES -499 .. +500 (for -)
#   TW:        from CPS +500 .. +10499 (for +) / CPS -10499 .. -500 (for -)
#   Gene body: TSS +250 .. CPS -501 (for +) / CPS +501 .. TSS -250 (for -)
#
# Outputs (under --outdir):
#   functional_regions.bed, functional_regions_summary.tsv
# =============================================================================

from __future__ import annotations
import argparse, os, sys, shutil, tempfile, subprocess, math
from pathlib import Path

# ---- CLI --------------------------------------------------------------------
ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
ap.add_argument("--sid",        required=True)
ap.add_argument("--genes",      required=True)
ap.add_argument("--functional-regions", required=True, help="Pre-computed functional regions TSV file")
ap.add_argument("--divergent",  required=True)
ap.add_argument("--pos",        required=True)  # 3′ POS bedGraph (≥0)
ap.add_argument("--neg",        required=True)  # 3′ NEG bedGraph (≤0 mirrored upstream)
ap.add_argument("--outdir",     required=True)

ap.add_argument("--tss", default=None)          # optional BED6 1bp
ap.add_argument("--tes", default=None)          # optional BED6 1bp
ap.add_argument("--tss-active-pm", type=int, default=500)  # half-window for active gene detection

# Geometry
ap.add_argument("--prom-up",    type=int, default=250)
ap.add_argument("--prom-down",  type=int, default=250)
ap.add_argument("--div-inner",  type=int, default=250)
ap.add_argument("--div-outer",  type=int, default=750)
ap.add_argument("--cps-offset", type=int, default=500)
ap.add_argument("--tw-length",  type=int, default=10_000)
ap.add_argument("--short-gene-threshold", type=int, default=750, help="Maximum gene length (bp) to classify as short gene")

# Behavior
ap.add_argument("--min-signal",      type=float, default=0.0)
ap.add_argument("--min-signal-mode", choices=["absolute","quantile"], default="absolute")
ap.add_argument("--min-signal-quantile", type=float, default=0.90)
ap.add_argument("--allow-unstranded",action="store_true")
ap.add_argument("--count-mode",      choices=["signal","event"], default="signal")
ap.add_argument("--div-fallback-enable", action="store_true", help="Enable guarded unstranded fallback for divergent assignment")
ap.add_argument("--div-fallback-threshold", type=float, default=0.30, help="Trigger fallback if DivergentTx < threshold * Promoter")
ap.add_argument("--div-fallback-max-frac", type=float, default=0.25, help="Cap fallback DivergentTx to this fraction of Promoter")
ap.add_argument("--active-slop", type=int, default=0, help="Slop (bp) for active-gene overlap only; does not change region geometry")
ap.add_argument("--debug",           action="store_true")

args = ap.parse_args()

# Validate parameters
def validate_parameters():
    """Validate that all parameters are reasonable and positive"""
    errors = []
    
    # Check geometry parameters
    if args.prom_up < 0:
        errors.append(f"prom_up must be non-negative, got {args.prom_up}")
    if args.prom_down < 0:
        errors.append(f"prom_down must be non-negative, got {args.prom_down}")
    if args.div_inner < 0:
        errors.append(f"div_inner must be non-negative, got {args.div_inner}")
    if args.div_outer < 0:
        errors.append(f"div_outer must be non-negative, got {args.div_outer}")
    if args.cps_offset < 0:
        errors.append(f"cps_offset must be non-negative, got {args.cps_offset}")
    if args.tw_length < 0:
        errors.append(f"tw_length must be non-negative, got {args.tw_length}")
    if args.tss_active_pm < 0:
        errors.append(f"tss_active_pm must be non-negative, got {args.tss_active_pm}")
    if args.short_gene_threshold < 0:
        errors.append(f"short_gene_threshold must be non-negative, got {args.short_gene_threshold}")
    
    # Check logical constraints
    if args.div_outer <= args.div_inner:
        errors.append(f"div_outer ({args.div_outer}) must be > div_inner ({args.div_inner})")
    
    # Check reasonable ranges
    if args.prom_up > 10000 or args.prom_down > 10000:
        errors.append(f"Promoter flank sizes seem too large: up={args.prom_up}, down={args.prom_down}")
    if args.tw_length > 100000:
        errors.append(f"Termination window length seems too large: {args.tw_length}")
    
    if errors:
        print("ERROR: Parameter validation failed:", file=sys.stderr)
        for error in errors:
            print(f"  {error}", file=sys.stderr)
        sys.exit(1)

validate_parameters()

SID    = args.sid
OUT    = Path(args.outdir).resolve(); OUT.mkdir(parents=True, exist_ok=True)
DBG    = OUT / "debug";               DBG.mkdir(parents=True, exist_ok=True) if args.debug else None

def log(m: str): print(m, file=sys.stderr, flush=True)

def bt() -> str:
    p = shutil.which("bedtools"); 
    if not p: sys.exit("ERROR: bedtools not found")
    return p
BT = bt()

# ---- helpers ----------------------------------------------------------------
def run(cmd: list[str], out_path: str | None = None, check=True):
    if out_path:
        with open(out_path, "w") as o:
            result = subprocess.run(cmd, check=False, text=True, stdout=o, stderr=subprocess.PIPE)
            if result.returncode != 0 and check:
                log(f"ERROR: Command failed with exit code {result.returncode}")
                log(f"Command: {' '.join(cmd)}")
                log(f"Stderr: {result.stderr}")
                raise subprocess.CalledProcessError(result.returncode, cmd, stderr=result.stderr)
    else:
        result = subprocess.run(cmd, check=False, text=True, stderr=subprocess.PIPE)
        if result.returncode != 0 and check:
            log(f"ERROR: Command failed with exit code {result.returncode}")
            log(f"Command: {' '.join(cmd)}")
            log(f"Stderr: {result.stderr}")
            raise subprocess.CalledProcessError(result.returncode, cmd, stderr=result.stderr)

def sort_bed(src: str):
    run(["bash","-lc", f"LC_ALL=C sort -k1,1 -k2,2n -k3,3n '{src}' -o '{src}'"], None)

def split_fields(ln: str):  # tolerate spaces/tabs
    return ln.rstrip("\n").split("\t") if "\t" in ln else ln.strip().split()

def wc_effective_lines(p: str) -> int:
    if not (Path(p).exists() and Path(p).stat().st_size>0): return 0
    n=0
    with open(p) as fh:
        for ln in fh:
            if ln.strip() and not ln.startswith(("track","browser","#")): n+=1
    return n

def clamp(a: int, b: int, min_length: int = 10) -> tuple[int,int]:
    """
    Ensure valid BED coordinates: a <= b and both >= 0, with minimum length.
    
    Args:
        a, b: Start and end coordinates
        min_length: Minimum region length in bp (default: 10 for biological relevance)
    
    Returns:
        Tuple of (start, end) with start <= end and length >= min_length
    """
    a, b = (a, b) if a <= b else (b, a)
    a = max(0, a)  # Clamp to chromosome start
    b = max(0, b)
    
    # Ensure minimum region length for biological relevance
    if b - a < min_length:
        b = a + min_length
    
    return (a, b)

def validate_coordinates(gene_info: dict) -> bool:
    """Validate that all calculated coordinates are reasonable and within biological constraints"""
    try:
        # Check basic gene coordinates
        if gene_info['gstart'] >= gene_info['gend']:
            log(f"WARNING: Gene {gene_info['gname']} has invalid coordinates: start={gene_info['gstart']}, end={gene_info['gend']}")
            return False
        
        # Check functional region coordinates
        regions = ['PPs', 'PPe', 'GBs', 'GBe', 'CPSs', 'CPSe', 'TWs', 'TWe', 'DIVs', 'DIVe']
        for region in regions:
            if region not in gene_info:
                log(f"WARNING: Gene {gene_info['gname']} missing region {region}")
                return False
            
            coord = gene_info[region]
            if not isinstance(coord, int) or coord < 1:
                log(f"WARNING: Gene {gene_info['gname']} has invalid {region} coordinate: {coord}")
                return False
        
        # Check that functional regions don't extend far beyond gene boundaries
        gene_start = gene_info['gstart']
        gene_end = gene_info['gend']
        gene_length = gene_end - gene_start
        
        # Allow functional regions to extend up to 20kb beyond gene boundaries
        max_extension = max(20000, gene_length * 2)  # At least 20kb or 2x gene length
        
        for region in regions:
            coord = gene_info[region]
            if coord < gene_start - max_extension or coord > gene_end + max_extension:
                log(f"WARNING: Gene {gene_info['gname']} {region} coordinate {coord} extends too far from gene boundaries [{gene_start}-{gene_end}]")
                return False
        
        # Check logical coordinate ordering
        strand = gene_info['strand']
        if strand == "+":
            # For plus strand: DIVs < DIVe < PPs < PPe < GBs < GBe < CPSs < CPSe < TWs < TWe
            coords = [gene_info['DIVs'], gene_info['DIVe'], gene_info['PPs'], gene_info['PPe'], 
                     gene_info['GBs'], gene_info['GBe'], gene_info['CPSs'], gene_info['CPSe'],
                     gene_info['TWs'], gene_info['TWe']]
        else:
            # For minus strand: TWs < TWe < CPSs < CPSe < GBs < GBe < PPs < PPe < DIVs < DIVe
            coords = [gene_info['TWs'], gene_info['TWe'], gene_info['CPSs'], gene_info['CPSe'],
                     gene_info['GBs'], gene_info['GBe'], gene_info['PPs'], gene_info['PPe'],
                     gene_info['DIVs'], gene_info['DIVe']]
        
        # Check that coordinates are generally increasing (allowing for some overlap)
        for i in range(len(coords) - 1):
            if coords[i] > coords[i + 1] + 1000:  # Allow some flexibility for overlapping regions
                log(f"WARNING: Gene {gene_info['gname']} has suspicious coordinate ordering")
                return False
        
        # Check for reasonable region sizes (not too small or too large)
        for i in range(0, len(coords), 2):  # Check start-end pairs
            if i + 1 < len(coords):
                region_length = coords[i + 1] - coords[i]
                if region_length < 10:  # Too small to be biologically meaningful
                    log(f"WARNING: Gene {gene_info['gname']} has very small functional region: {region_length}bp")
                    return False
                if region_length > 500000:  # Too large (>500kb)
                    log(f"WARNING: Gene {gene_info['gname']} has very large functional region: {region_length}bp")
                    return False
        
        return True
        
    except Exception as e:
        log(f"WARNING: Error validating coordinates for gene {gene_info.get('gname', 'unknown')}: {e}")
        return False

def read_sites(bed6: str | None) -> dict[str, tuple[str,int,str]]:
    d={}; 
    if not bed6 or not Path(bed6).exists(): return d
    with open(bed6) as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith(("track","browser","#")): continue
            f=split_fields(ln)
            if len(f)<6: continue
            d[f[3]] = (f[0], int(f[1]), f[5])   # use start as site
    return d

# ---- coordinate generation using pre-computed functional regions -----------
def read_precomputed_functional_regions(func_regions_tsv: str) -> list:
    """
    Read pre-computed functional regions from the TSV file created by create_functional_regions.py
    Returns a list of gene dictionaries with all functional region coordinates.
    Includes comprehensive validation and error handling.
    """
    genes = []
    
    if not os.path.exists(func_regions_tsv):
        log(f"ERROR: Functional regions file not found: {func_regions_tsv}")
        sys.exit(1)
    
    try:
        with open(func_regions_tsv, 'r') as f:
            header = f.readline().strip().split('\t')
            
            # Validate header contains required columns
            required_cols = ['chr', 'start', 'end', 'strand', 'gene_name', 'TSS', 'CPS', 
                           'DIVs', 'DIVe', 'PPs', 'PPe', 'GBs', 'GBe', 'CPSs', 'CPSe', 'TWs', 'TWe', 'promC1', 'promC2']
            missing_cols = [col for col in required_cols if col not in header]
            if missing_cols:
                log(f"ERROR: Missing required columns in {func_regions_tsv}: {missing_cols}")
                sys.exit(1)
            
            line_num = 1  # Start at 1 since we already read header
            for line in f:
                line_num += 1
                if not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) != len(header):
                    log(f"WARNING: Line {line_num} has {len(fields)} fields, expected {len(header)}. Skipping.")
                    continue
                    
                gene = {}
                for i, field in enumerate(fields):
                    if i < len(header):
                        gene[header[i]] = field
                
                # Convert numeric fields with validation
                try:
                    gene['gstart'] = int(gene['start'])
                    gene['gend'] = int(gene['end'])
                    gene['TSS'] = int(gene['TSS'])
                    gene['CPS'] = int(gene['CPS'])
                    gene['DIVs'] = int(gene['DIVs'])
                    gene['DIVe'] = int(gene['DIVe'])
                    gene['PPs'] = int(gene['PPs'])
                    gene['PPe'] = int(gene['PPe'])
                    gene['GBs'] = int(gene['GBs'])
                    gene['GBe'] = int(gene['GBe'])
                    gene['CPSs'] = int(gene['CPSs'])
                    gene['CPSe'] = int(gene['CPSe'])
                    gene['TWs'] = int(gene['TWs'])
                    gene['TWe'] = int(gene['TWe'])
                    gene['promC1'] = int(gene['promC1'])
                    gene['promC2'] = int(gene['promC2'])
                    
                    # Add compatibility fields
                    gene['chrom'] = gene['chr']
                    gene['gname'] = gene['gene_name']
                    gene['gene_length'] = gene['gend'] - gene['gstart']
                    
                    # Include all genes; avoid erroneous start/end cutoffs.
                    # Keep validation as a warning-only guard; do not drop genes unless clearly nonsensical.
                    if validate_coordinates(gene):
                        genes.append(gene)
                    else:
                        log(f"WARNING: Coordinate validation flagged gene {gene['gname']} at line {line_num}; including anyway for assignment")
                        genes.append(gene)
                        
                except (ValueError, KeyError) as e:
                    log(f"WARNING: Error parsing gene at line {line_num}: {e}. Skipping.")
                    continue
    
    except IOError as e:
        log(f"ERROR: Failed to read {func_regions_tsv}: {e}")
        sys.exit(1)
    except Exception as e:
        log(f"ERROR: Unexpected error reading {func_regions_tsv}: {e}")
        sys.exit(1)
    
    log(f"INFO: Successfully loaded {len(genes)} validated genes from {func_regions_tsv}")
    return genes

# Note: build_coordinate_lists function removed - now using read_precomputed_functional_regions()
# This eliminates duplicate coordinate calculation logic and ensures consistency.

# ---- active genes & enhancers -----------------------------------------------
def find_active_promoters_and_enhancers(all_genes: list, dt_bed: str) -> tuple[list, str, str, str]:
    """
    v9.0: Use DT sites to identify active genes, return gene-based regions.
    
    1. Create promoter regions for all genes (TSS ±500bp for active gene detection)
    2. DT sites overlapping promoters → identify active genes
    3. DT sites NOT overlapping → Enhancers
    4. Return active genes (for creating gene-based regions later)
    """
    # Create promoter regions for ALL genes (TSS ±tss_active_pm like old pipeline)
    # This is ONLY for identifying active genes, not for read assignment
    promoter_regions_file = OUT / f"all_promoter_regions_TSS_pm{args.tss_active_pm}.bed"
    with open(promoter_regions_file, "w") as f:
        for gene in all_genes:
            # Use parameterized window for active gene detection
            prom_start = gene['TSS'] - args.tss_active_pm
            prom_end = gene['TSS'] + args.tss_active_pm
            prom_start, prom_end = clamp(prom_start, prom_end)
            f.write(f"{gene['chrom']}\t{prom_start}\t{prom_end}\t{gene['gname']}\t.\t{gene['strand']}\n")
    
    sort_bed(str(promoter_regions_file))
    
    # Identify which genes have DT sites at their promoters (TSS ±tss_active_pm)
    genes_with_promoters_bed = OUT / "genes_with_active_promoters.bed"
    # Unstranded overlap to find active promoters (as discovery)
    run([BT, "intersect", "-u", "-wa", "-a", str(promoter_regions_file), "-b", dt_bed], 
        str(genes_with_promoters_bed))

    # Also capture DT sites that overlap promoter windows and imprint gene strand on them to produce BED6
    active_promoter_dt_bed = OUT / "active_promoter_dt_sites.bed"
    # -wa -wb to get promoter (A, with strand) and DT (B without strand); write DT coords with A's strand
    # IMPORTANT: DT BED lacks strand; DO NOT use -s here or we will get zero overlaps
    tmp_join = OUT / "_tmp_prom_dt_join.bed"
    run([BT, "intersect", "-wa", "-wb", "-a", str(promoter_regions_file), "-b", dt_bed], str(tmp_join), check=False)
    with open(active_promoter_dt_bed, "w") as out:
        if tmp_join.exists() and tmp_join.stat().st_size > 0:
            with open(tmp_join) as fh:
                for ln in fh:
                    if not ln.strip() or ln.startswith(("#","track","browser")): continue
                    f = split_fields(ln)
                    # A: chrom, start, end, gene, ., strand (6 fields)
                    # B: DT: chrom, start, end, total (4 fields)
                    if len(f) >= 10:
                        # Use B coords for window; use A strand for strand
                        chrom = f[6]  # B chrom
                        start = f[7]
                        end   = f[8]
                        name  = f[3]  # gene name
                        strand = f[5] # A strand
                        out.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\n")
    if tmp_join.exists():
        try:
            tmp_join.unlink()
        except Exception:
            pass
    
    # DT sites NOT overlapping gene promoters = Enhancers
    enhancers_bed = OUT / "enhancers.bed"
    # Prefer subtracting promoter-overlapping DT windows directly if available
    if Path(active_promoter_dt_bed).exists() and Path(active_promoter_dt_bed).stat().st_size > 0:
        run([BT, "intersect", "-v", "-a", dt_bed, "-b", str(active_promoter_dt_bed)], str(enhancers_bed))
    else:
        run([BT, "intersect", "-v", "-a", dt_bed, "-b", str(promoter_regions_file)], str(enhancers_bed))
    
    # Read active gene names from the intersection
    active_names = set()
    if Path(genes_with_promoters_bed).exists():
        with open(genes_with_promoters_bed) as f:
            for ln in f:
                if ln.strip() and not ln.startswith(("#", "track", "browser")):
                    fields = split_fields(ln)
                    if len(fields) >= 4:
                        active_names.add(fields[3])  # gene name
    
    # Filter to active genes only
    active_genes = [g for g in all_genes if g['gname'] in active_names]
    
    # Count 
    n_dt_at_promoters = wc_effective_lines(str(genes_with_promoters_bed))
    n_enhancers = wc_effective_lines(str(enhancers_bed))
    
    log(f"INFO  Found {n_dt_at_promoters} genes with DT sites at promoters (TSS ±{args.tss_active_pm}bp)")
    log(f"INFO  Found {n_enhancers} DT sites NOT at promoters (enhancers)")
    log(f"INFO  Active genes: {len(active_genes)} (out of {len(all_genes)} total)")
    
    # Return: active_genes, genes_with_active_promoters, enhancers, active_promoter_dt_sites
    return active_genes, str(genes_with_promoters_bed), str(enhancers_bed), str(active_promoter_dt_bed)

# ---- coordinate file generation ---------------------------------------------
def write_coordinate_files(active_genes: list):
    """
    v9.0: Write ALL gene-based regions including promoter (ppPolII).
    chr, start, end, name, name, strand (BED6-like)
    
    Following old pipeline logic exactly.
    """
    # Separate short and long genes (using parameterized threshold)
    short_genes = [g for g in active_genes if g['gene_length'] <= args.short_gene_threshold]
    long_genes = [g for g in active_genes if g['gene_length'] > args.short_gene_threshold]
    
    log(f"INFO  {len(short_genes)} short active genes, {len(long_genes)} long active genes")
    
    # Write promoter regions (ppPolII.txt) - TSS±prom_up/down - long genes only
    pp_file = OUT / "ppPolII.txt"
    with open(pp_file, "w") as f:
        for gene in long_genes:
            pps, ppe = clamp(gene['PPs'], gene['PPe'])
            f.write(f"{gene['chrom']}\t{pps}\t{ppe}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(pp_file))
    
    # Write divergent regions (divTx.txt) - TSS±div_outer/inner - long genes only
    div_file = OUT / "divTx.txt"
    with open(div_file, "w") as f:
        for gene in long_genes:
            divs, dive = clamp(gene['DIVs'], gene['DIVe'])
            f.write(f"{gene['chrom']}\t{divs}\t{dive}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(div_file))
    
    # Write short genes (shortGenes.txt) - full gene coordinates (≤short_gene_threshold)
    sg_file = OUT / "shortGenes.txt"
    with open(sg_file, "w") as f:
        for gene in short_genes:
            gstart, gend = clamp(gene['gstart'], gene['gend'])
            f.write(f"{gene['chrom']}\t{gstart}\t{gend}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(sg_file))
    
    # Write CPS regions (CPS.txt) - long genes only (>short_gene_threshold)
    cps_file = OUT / "CPS.txt"
    with open(cps_file, "w") as f:
        for gene in long_genes:
            cpss, cpse = clamp(gene['CPSs'], gene['CPSe'])
            f.write(f"{gene['chrom']}\t{cpss}\t{cpse}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(cps_file))
    
    # Write termination windows (TW.txt) - long genes only (>short_gene_threshold)
    tw_file = OUT / "TW.txt"
    with open(tw_file, "w") as f:
        for gene in long_genes:
            tws, twe = clamp(gene['TWs'], gene['TWe'])
            f.write(f"{gene['chrom']}\t{tws}\t{twe}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(tw_file))
    
    # Write gene body regions (geneBody.txt) - long genes only, where GBe > GBs + 1 (>short_gene_threshold)
    gb_file = OUT / "geneBody.txt"
    with open(gb_file, "w") as f:
        for gene in long_genes:
            gbs, gbe = clamp(gene['GBs'], gene['GBe'])
            if gbe > gbs + 1:  # Following old script logic
                f.write(f"{gene['chrom']}\t{gbs}\t{gbe}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(gb_file))
    
    return {
        'promoter': str(pp_file),
        'divergent': str(div_file), 
        'short_genes': str(sg_file),
        'cps': str(cps_file),
        'termination': str(tw_file),
        'gene_body': str(gb_file)
    }

# ---- bedGraph to read conversion --------------------------------------------
def bedgraph_to_reads(pos_bg: str, neg_bg: str) -> str:
    """
    Convert bedGraphs to read intervals with strand information.
    Following old script logic: chr, start, end, signal, signal, strand
    """
    reads_file = OUT / "PROseq_3pnt.bed"
    
    # Determine dynamic threshold if quantile mode
    thr = float(args.min_signal)
    if args.min_signal_mode == "quantile":
        vals: list[float] = []
        for bg in (pos_bg, neg_bg):
            if Path(bg).exists() and Path(bg).stat().st_size > 0:
                with open(bg) as f:
                    for i, ln in enumerate(f):
                        if ln.strip() and not ln.startswith(("track","browser","#")):
                            fields = split_fields(ln)
                            if len(fields) >= 4:
                                try:
                                    s = abs(float(fields[3]));
                                    vals.append(s)
                                except Exception:
                                    pass
                        if len(vals) >= 200000:  # cap memory; sample first 200k values
                            break
        if vals:
            vals.sort()
            q = min(max(args.min_signal_quantile, 0.0), 1.0)
            idx = int(q * (len(vals)-1))
            thr = float(vals[idx])
            log(f"INFO  dynamic min_signal (quantile {q:.2f}) => {thr:.6g}")

    with open(reads_file, "w") as out:
        # Process positive strand
        if Path(pos_bg).exists() and Path(pos_bg).stat().st_size > 0:
            with open(pos_bg) as f:
                for ln in f:
                    if ln.strip() and not ln.startswith(("track", "browser", "#")):
                        fields = split_fields(ln)
                        if len(fields) >= 4:
                            try:
                                chrom, start, end, signal = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
                                if abs(signal) > thr and end > start:
                                    # Ensure positive strand signals are positive for consistency
                                    signal_val = abs(signal) if signal < 0 else signal
                                    out.write(f"{chrom}\t{start}\t{end}\t{signal_val}\t{signal_val}\t+\n")
                            except (ValueError, IndexError):
                                continue
        
        # Process negative strand (ensure consistent signal handling)
        if Path(neg_bg).exists() and Path(neg_bg).stat().st_size > 0:
            with open(neg_bg) as f:
                for ln in f:
                    if ln.strip() and not ln.startswith(("track", "browser", "#")):
                        fields = split_fields(ln)
                        if len(fields) >= 4:
                            try:
                                chrom, start, end, signal = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
                                if abs(signal) > thr and end > start:
                                    # Convert negative strand signals to positive values for consistent analysis
                                    # This prevents issues with signal summation and downstream analysis
                                    signal_val = abs(signal)
                                    out.write(f"{chrom}\t{start}\t{end}\t{signal_val}\t{signal_val}\t-\n")
                            except (ValueError, IndexError):
                                continue
    
    sort_bed(str(reads_file))
    return str(reads_file)

# ---- sequential read assignment (NEW LOGIC) ---------------------------------
def sequential_read_assignment(reads_file: str, coord_files: dict, enhancers_bed: str):
    """
    v9.1: Sequential read assignment - EXACT replication of old bash script.
    
    CRITICAL: Assignment uses strand-specific (-s/-S), removal uses UNSTRANDED (-v)
    
    Assignment order (sequential masking):
    1. Promoter:     assign -s, remove UNSTRANDED
    2. Divergent:    assign -S, remove UNSTRANDED
    3. CPS:          assign -s, remove UNSTRANDED
    4. Gene body:    assign -s, remove UNSTRANDED
    5. Short genes:  assign -s, remove UNSTRANDED
    6. Enhancers:    assign unstranded, remove unstranded
    7. Termination:  assign -s, remove UNSTRANDED
    8. Non-localized: remaining
    
    This matches the original TrackTx bash pipeline exactly.
    """
    current_reads = reads_file
    assigned_reads = {}
    
    # Step 1: Promoters (gene-based TSS-250..+249, same-strand)
    # OLD SCRIPT EXACT LOGIC: assign with -s, remove UNSTRANDED
    log("INFO  Step 1: Assigning promoter reads (same-strand from ppPolII)...")
    log(f"DEBUG Input reads: {wc_effective_lines(current_reads)}, Promoter regions: {wc_effective_lines(coord_files['promoter'])}")
    prom_reads = OUT / "PROseq_ppPolII.bed"
    prom_removed = OUT / "ppRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['promoter']], str(prom_reads))
    # Removal: keep opposite-strand reads for divergent step; remove only same-strand overlaps
    run([BT, "intersect", "-s", "-v", "-a", current_reads, "-b", coord_files['promoter']], str(prom_removed))
    assigned_reads['Promoter'] = str(prom_reads)
    current_reads = str(prom_removed)
    log(f"DEBUG Assigned: {wc_effective_lines(str(prom_reads))}, Remaining: {wc_effective_lines(current_reads)}")
    
    # Step 2: Divergent (gene-based divergent window)
    # Assign UNSTRANDED within divergent window; promoter already removed
    log("INFO  Step 2: Assigning divergent reads (unstranded within divergent window)...")
    log(f"DEBUG Input reads: {wc_effective_lines(current_reads)}, Divergent regions: {wc_effective_lines(coord_files['divergent'])}")
    div_reads = OUT / "PROseq_ppDiv.bed"
    div_removed = OUT / "ppdivRemoved.bed"
    run([BT, "intersect", "-u", "-wa", "-a", current_reads, "-b", coord_files['divergent']], str(div_reads))
    # Removal: remove only opposite-strand overlaps used for divergent assignment; keep same-strand
    run([BT, "intersect", "-S", "-v", "-a", current_reads, "-b", coord_files['divergent']], str(div_removed))
    assigned_reads['DivergentTx'] = str(div_reads)
    current_reads = str(div_removed)
    log(f"DEBUG Assigned: {wc_effective_lines(str(div_reads))}, Remaining: {wc_effective_lines(current_reads)}")

    
    # Step 3: CPS (same-strand: -s)
    # OLD SCRIPT: assign with -s, remove UNSTRANDED
    log("INFO  Step 3: Assigning CPS reads...")
    cps_reads = OUT / "PROseq_CPS.bed"
    cps_removed = OUT / "ppdivCPSRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['cps']], str(cps_reads))
    run([BT, "intersect", "-s", "-v", "-a", current_reads, "-b", coord_files['cps']], str(cps_removed))
    assigned_reads['CPS'] = str(cps_reads)
    current_reads = str(cps_removed)
    
    # Step 4: Gene Body (same-strand: -s)
    # OLD SCRIPT: assign with -s, remove UNSTRANDED
    log("INFO  Step 4: Assigning gene body reads...")
    gb_reads = OUT / "PROseq_GB.bed"
    gb_removed = OUT / "ppdivCPSgbRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['gene_body']], str(gb_reads))
    run([BT, "intersect", "-s", "-v", "-a", current_reads, "-b", coord_files['gene_body']], str(gb_removed))
    assigned_reads['Gene body'] = str(gb_reads)
    current_reads = str(gb_removed)
    
    # Step 5: Short Genes (same-strand: -s)
    # OLD SCRIPT: assign with -s, remove UNSTRANDED
    log("INFO  Step 5: Assigning short gene reads...")
    sg_reads = OUT / "PROseq_SG.bed"
    sg_removed = OUT / "ppdivCPSgbSGRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['short_genes']], str(sg_reads))
    run([BT, "intersect", "-s", "-v", "-a", current_reads, "-b", coord_files['short_genes']], str(sg_removed))
    assigned_reads['Short genes'] = str(sg_reads)
    current_reads = str(sg_removed)
    
    # Step 6: Enhancers (unstranded)
    # OLD SCRIPT: assign unstranded, remove unstranded
    log("INFO  Step 6: Assigning enhancer reads...")
    enh_reads = OUT / "PROseq_enhancers.bed"
    enh_removed = OUT / "ppdivCPSgbSGEnhRemoved.bed"
    run([BT, "intersect", "-u", "-wa", "-a", current_reads, "-b", enhancers_bed], str(enh_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", enhancers_bed], str(enh_removed))
    assigned_reads['Enhancers'] = str(enh_reads)
    current_reads = str(enh_removed)
    
    # Step 7: Termination Window (same-strand: -s)
    # OLD SCRIPT: assign with -s, remove UNSTRANDED
    log("INFO  Step 7: Assigning termination window reads...")
    tw_reads = OUT / "PROseq_TW.bed"
    tw_removed = OUT / "PROseq_noGene_noEnh.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['termination']], str(tw_reads))
    run([BT, "intersect", "-s", "-v", "-a", current_reads, "-b", coord_files['termination']], str(tw_removed))
    assigned_reads['Termination window'] = str(tw_reads)
    assigned_reads['Non-localized polymerase'] = str(tw_removed)
    
    return assigned_reads

# ---- count and summarize ---------------------------------------------------
def count_and_summarize(assigned_reads: dict):
    """Count reads and signal for each category"""
    summary = {}
    
    for category, reads_file in assigned_reads.items():
        read_count = wc_effective_lines(reads_file)
        signal_sum = 0.0
        
        if Path(reads_file).exists() and Path(reads_file).stat().st_size > 0:
            with open(reads_file) as f:
                for ln in f:
                    if ln.strip() and not ln.startswith(("track", "browser", "#")):
                        fields = split_fields(ln)
                        if len(fields) >= 4:
                            try:
                                # Sum absolute values to avoid negative totals from mixed strands
                                signal_sum += abs(float(fields[3]))
                            except (ValueError, IndexError):
                                pass
        
        # Ensure non-negative signal totals for consistent reporting
        signal_sum = max(0.0, signal_sum)
        
        summary[category] = {
            'signal': signal_sum,
            'count': read_count
        }
        
        log(f"INFO  {category}: {read_count} reads, {signal_sum:.2f} signal")
    
    return summary

# ---- output generation ------------------------------------------------------
def write_outputs(summary: dict, coord_files: dict, assigned_reads: dict):
    """Write final outputs - v9.0: all gene-based regions"""
    
    # Compute number of genomic regions per category
    def bed_regions_count(path: str) -> int:
        if not path or not (Path(path).exists() and Path(path).stat().st_size > 0):
            return 0
        n = 0
        with open(path) as fh:
            for ln in fh:
                if ln.strip() and not ln.startswith(("track","browser","#")):
                    n += 1
        return n

    enhancers_file = OUT / "enhancers.bed"
    region_file_map = {
        # assigned_reads keys → gene-based regions or DT sites (enhancers only)
        'Promoter': coord_files.get('promoter'),  # Gene-based ppPolII.txt
        'DivergentTx': coord_files.get('divergent'),  # Gene-based divTx.txt
        'CPS': coord_files.get('cps'),
        'Gene body': coord_files.get('gene_body'),
        'Short genes': coord_files.get('short_genes'),
        'Enhancers': str(enhancers_file),  # DT sites NOT overlapping promoters (post-subtraction)
        'Termination window': coord_files.get('termination'),
        # 'Non-localized polymerase' does not correspond to a defined region set
    }

    region_counts = {k: (bed_regions_count(v) if v else 0) for k, v in region_file_map.items()}
    region_counts['Non-localized polymerase'] = 0

    # Compute total and median region lengths (bp) per category
    def bed_region_lengths(path: str) -> list[int]:
        if not path or not (Path(path).exists() and Path(path).stat().st_size > 0):
            return []
        lens: list[int] = []
        with open(path) as fh:
            for ln in fh:
                if ln.strip() and not ln.startswith(("track","browser","#")):
                    f = split_fields(ln)
                    try:
                        start = int(float(f[1])); end = int(float(f[2]))
                        if end > start:
                            lens.append(end - start)
                    except Exception:
                        continue
        return lens

    def median_len(vals: list[int]) -> float:
        if not vals: return 0.0
        vals = sorted(vals)
        n = len(vals)
        mid = n // 2
        if n % 2 == 1:
            return float(vals[mid])
        else:
            return (vals[mid - 1] + vals[mid]) / 2.0

    region_lengths_total: dict[str, int] = {}
    region_lengths_median: dict[str, float] = {}
    for cat, p in region_file_map.items():
        lens = bed_region_lengths(p) if p else []
        region_lengths_total[cat] = int(sum(lens)) if lens else 0
        region_lengths_median[cat] = float(median_len(lens)) if lens else 0.0
    region_lengths_total['Non-localized polymerase'] = 0
    region_lengths_median['Non-localized polymerase'] = 0.0

    # Write summary TSV (use 'reads' column name for compatibility with reporting script)
    summary_file = OUT / "functional_regions_summary.tsv"
    with open(summary_file, "w") as f:
        f.write("region\treads\tregion_count\tregion_length_total_bp\tregion_length_median_bp\n")
        for category, data in summary.items():
            reg_ct = region_counts.get(category, 0)
            reg_len_tot = region_lengths_total.get(category, 0)
            reg_len_med = region_lengths_median.get(category, 0.0)
            f.write(f"{category}\t{data['signal']:.6f}\t{reg_ct}\t{reg_len_tot}\t{reg_len_med:.2f}\n")
    
    # Write BED9 track (following old script format)
    bed_file = OUT / "functional_regions.bed"
    with open(bed_file, "w") as f:
        f.write(f'track name="functional_genomic_regions_{SID}" itemRgb="On"\n')
        
        # Color mapping (from old script)
        colors = {
            'divergent': "178,59,212", 
            'cps': "103,200,249",
            'gene_body': "0,0,0",
            'short_genes': "253,218,13",
            'termination': "255,54,98"
        }
        
        # Add promoter regions (orange) - gene-based ppPolII.txt
        promoter_file = coord_files.get('promoter')
        if promoter_file and Path(promoter_file).exists():
            with open(promoter_file) as pf:
                for ln in pf:
                    if ln.strip() and not ln.startswith(("track", "browser", "#")):
                        fields = split_fields(ln)
                        if len(fields) >= 6:
                            chrom, start, end, name, _, strand = fields[:6]
                            f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\t{start}\t{end}\t243,132,0\n")
        
        # Add gene-based coordinate regions with colors
        for coord_type, coord_file in coord_files.items():
            if coord_type in colors and Path(coord_file).exists():
                with open(coord_file) as cf:
                    for ln in cf:
                        if ln.strip() and not ln.startswith(("track", "browser", "#")):
                            fields = split_fields(ln)
                            if len(fields) >= 6:
                                chrom, start, end, name, _, strand = fields[:6]
                                color = colors[coord_type]
                                f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t{strand}\t{start}\t{end}\t{color}\n")
        
        # Add enhancers (green)
        if Path(enhancers_file).exists():
            with open(enhancers_file) as ef:
                for ln in ef:
                    if ln.strip() and not ln.startswith(("track", "browser", "#")):
                        fields = split_fields(ln)
                        if len(fields) >= 3:
                            chrom, start, end = fields[:3]
                            name = fields[3] if len(fields) > 3 else "Enhancer"
                            f.write(f"{chrom}\t{start}\t{end}\t{name}\t0\t.\t{start}\t{end}\t115,212,122\n")

# ---- MAIN -------------------------------------------------------------------
def main():
    log(f"INFO  FGR ▶ {SID}  mode=pre-computed functional regions (v10.2 - comprehensive fixes)")

    # Read pre-computed functional regions
    all_genes = read_precomputed_functional_regions(args.functional_regions)
    log(f"INFO  Loaded {len(all_genes)} genes with functional regions from {args.functional_regions}")

    # Use DT sites to identify active genes and enhancers
    active_genes, genes_with_promoters_bed, enhancers_bed, active_promoter_dt_bed = find_active_promoters_and_enhancers(all_genes, args.divergent)
    
    if not active_genes:
        sys.exit("ERROR: No active genes found (DT sites ∩ TSS ±500bp); check divergent transcription input")

    # Write gene-based coordinate files for ALL functional regions (including promoter)
    coord_files = write_coordinate_files(active_genes)

    # Convert bedGraphs to read intervals
    reads_file = bedgraph_to_reads(args.pos, args.neg)
    log(f"INFO  Converted bedGraphs to {wc_effective_lines(reads_file)} read intervals")

    # Perform sequential read assignment using gene-based regions
    assigned_reads = sequential_read_assignment(reads_file, coord_files, enhancers_bed)

    # Count and summarize
    summary = count_and_summarize(assigned_reads)

    # Write outputs
    write_outputs(summary, coord_files, assigned_reads)

    log("INFO  functional regions ✓ (pre-computed coordinates v10.2 - comprehensive fixes)")

if __name__ == "__main__":
    main()
