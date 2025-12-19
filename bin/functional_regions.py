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
ap.add_argument("--tw-length",  type=int, default=10_000)

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

def sort_bed(src: str, buffer_size: str = "1G", parallel: int = None):
    """
    Sort a BED file in-place using LC_ALL=C sort with memory limits.
    
    Args:
        src: Path to BED file to sort
        buffer_size: Memory buffer size for sort (default: 1G)
        parallel: Number of parallel threads (default: from THREADS env or 1)
    """
    # Use environment variable for threads if available, or default to 1
    if parallel is None:
        try:
            parallel = int(os.environ.get('THREADS', '1'))
        except (ValueError, TypeError):
            parallel = 1
    
    # Build sort command with buffer size limit and parallel processing
    sort_cmd = (
        f"LC_ALL=C sort -k1,1 -k2,2n -k3,3n "
        f"-S {buffer_size} "
        f"--parallel={parallel} "
        f"'{src}' -o '{src}'"
    )
    
    run(["bash", "-c", sort_cmd], None)

def split_fields(ln: str):  # tolerate spaces/tabs
    return ln.rstrip("\n").split("\t") if "\t" in ln else ln.strip().split()

def wc_effective_lines(p: str) -> int:
    if not (Path(p).exists() and Path(p).stat().st_size>0): return 0
    n=0
    with open(p) as fh:
        for ln in fh:
            if ln.strip() and not ln.startswith(("track","browser","#")): n+=1
    return n

def clamp(a: int, b: int) -> tuple[int,int]:
    """Ensure valid BED coordinates: a <= b and both >= 0, with minimum 1bp length"""
    a, b = (a, b) if a <= b else (b, a)
    a = max(0, a)  # Clamp to chromosome start
    b = max(0, b)
    # Ensure valid interval: end must be > start (minimum 1bp)
    if b <= a:
        b = a + 1
    return (a, b)

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

# ---- coordinate generation (following old script logic) --------------------
def build_coordinate_lists(genes_tsv: str, tss_map: dict, tes_map: dict):
    """
    Build coordinate lists following the old TrackTx.sh logic exactly.
    This replaces the signal-based approach with coordinate list generation.
    """
    req = ["gene_id","gene_name","chr","strand","start","end","tss","tes"]
    with open(genes_tsv) as fh:
        hdr = split_fields(fh.readline())
        cols = [c.strip().lstrip("\ufeff").lower() for c in hdr]
    
    def idx(name: str):
        alts = {
            "gene_id":{"gene_id","id"},
            "gene_name":{"gene_name","name","symbol"},
            "chr":{"chr","chrom","chromosome","seqname"},
            "strand":{"strand","orientation"},
            "start":{"start","gene_start"},
            "end":{"end","gene_end"},
            "tss":{"tss","tx_start","txstart"},
            "tes":{"tes","tx_end","txend"},
        }[name]
        for i,c in enumerate(cols):
            if c in alts: return i
        return None
    
    idxs = {k: idx(k) for k in req}
    miss = [k for k,v in idxs.items() if v is None]
    if miss: sys.exit("ERROR genes.tsv missing: " + ",".join(miss))

    # Store all gene info for coordinate calculation
    all_genes = []
    
    with open(genes_tsv) as fh:
        next(fh)  # skip header
        for ln in fh:
            if not ln.strip(): continue
            f = split_fields(ln)
            try:
                gid = f[idxs["gene_id"]]
                gname = f[idxs["gene_name"]] or gid
                chrom = f[idxs["chr"]]
                strand = f[idxs["strand"]] if f[idxs["strand"]] in {"+","-"} else "+"
                gstart = int(float(f[idxs["start"]]))
                gend = int(float(f[idxs["end"]]))
                tss_d = int(float(f[idxs["tss"]]))
                tes_d = int(float(f[idxs["tes"]]))
            except Exception:
                continue

            # Apply overrides (by name or id)
            tss = tss_map.get(gname, tss_map.get(gid, (None,tss_d,None)))[1]
            tes = tes_map.get(gname, tes_map.get(gid, (None,tes_d,None)))[1]

            # Calculate all coordinates following old script logic EXACTLY
            # CRITICAL: For minus strand, old script swaps TSS/CPS assignments!
            if strand == "+":
                TSS = tss    # For +: TSS = txStart
                CPS = tes    # For +: CPS = txEnd
                DIVs = TSS - args.div_outer    # TSS - 750
                DIVe = TSS - args.div_inner    # TSS - 251  
                PPs = TSS - args.prom_up       # TSS - 250
                PPe = TSS + 249                # TSS + 249 (old script exact value)
                GBs = TSS + args.prom_down     # TSS + 250
                GBe = CPS - 501                # CPS - 501
                CPSs = CPS - 500
                CPSe = CPS + 499
                TWs = CPS + 500
                TWe = CPS + 10499              # CPS + 10499 (old script exact value)
            else:  # strand == "-"
                # CRITICAL FIX: For minus strand, old script uses txEnd as TSS and txStart as CPS
                TSS = tes    # For -: TSS = txEnd (old script: df$TSS <- df$txEnd)
                CPS = tss    # For -: CPS = txStart (old script: df$CPS <- df$txStart)
                DIVs = TSS + args.div_inner    # TSS + 251
                DIVe = TSS + args.div_outer    # TSS + 750
                PPs = TSS - 249                # TSS - 249 (old script exact value)
                PPe = TSS + args.prom_up       # TSS + 250
                GBs = CPS + 501                # CPS + 501
                GBe = TSS - args.prom_down     # TSS - 250
                CPSs = CPS - 499
                CPSe = CPS + 500
                TWs = CPS - 10499              # CPS - 10499 (old script exact value)
                TWe = CPS - 500

            # Store gene info (use transcript length like old script: txEnd - txStart)
            # Old script always uses txEnd - txStart regardless of strand
            tx_length = tes - tss  # transcript length (old script: txEnd - txStart)
            gene_info = {
                'gid': gid, 'gname': gname, 'chrom': chrom, 'strand': strand,
                'gstart': gstart, 'gend': gend, 'TSS': TSS, 'CPS': CPS,
                'DIVs': DIVs, 'DIVe': DIVe, 'PPs': PPs, 'PPe': PPe,
                'GBs': GBs, 'GBe': GBe, 'CPSs': CPSs, 'CPSe': CPSe,
                'TWs': TWs, 'TWe': TWe, 'gene_length': tx_length
            }
            all_genes.append(gene_info)

    return all_genes

# ---- active genes & enhancers -----------------------------------------------
def find_active_promoters_and_enhancers(all_genes: list, dt_bed: str) -> tuple[list, str, str]:
    """
    v9.0: Use DT sites to identify active genes, return gene-based regions.
    
    1. Create promoter regions for all genes (TSS ±500bp for active gene detection)
    2. DT sites overlapping promoters → identify active genes
    3. DT sites NOT overlapping → Enhancers
    4. Return active genes (for creating gene-based regions later)
    """
    # Create promoter regions for ALL genes (TSS ±500bp like old pipeline)
    # This is ONLY for identifying active genes, not for read assignment
    promoter_regions_500bp = OUT / "all_promoter_regions_TSS_pm500.bed"
    with open(promoter_regions_500bp, "w") as f:
        for gene in all_genes:
            # Fixed ±500bp window for active gene detection
            prom_start = gene['TSS'] - 500
            prom_end = gene['TSS'] + 500
            prom_start, prom_end = clamp(prom_start, prom_end)
            f.write(f"{gene['chrom']}\t{prom_start}\t{prom_end}\t{gene['gname']}\t.\t{gene['strand']}\n")
    
    sort_bed(str(promoter_regions_500bp))
    
    # Identify which genes have DT sites at their promoters (TSS ±500bp)
    genes_with_promoters_bed = OUT / "genes_with_active_promoters.bed"
    run([BT, "intersect", "-u", "-wa", "-a", str(promoter_regions_500bp), "-b", dt_bed], 
        str(genes_with_promoters_bed))
    
    # DT sites NOT overlapping gene promoters = Enhancers
    enhancers_bed = OUT / "enhancers.bed"
    run([BT, "intersect", "-v", "-a", dt_bed, "-b", str(promoter_regions_500bp)], 
        str(enhancers_bed))
    
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
    
    log(f"INFO  Found {n_dt_at_promoters} genes with DT sites at promoters (TSS ±500bp)")
    log(f"INFO  Found {n_enhancers} DT sites NOT at promoters (enhancers)")
    log(f"INFO  Active genes: {len(active_genes)} (out of {len(all_genes)} total)")
    
    # Return: active_genes, path_to_genes_with_promoters, path_to_enhancers
    # Note: We don't return active_promoters_bed anymore - we'll use gene-based regions
    return active_genes, str(genes_with_promoters_bed), str(enhancers_bed)

# ---- coordinate file generation ---------------------------------------------
def write_coordinate_files(active_genes: list):
    """
    v10.0: Write ALL gene-based regions including promoter (ppPolII).
    chr, start, end, name, name, strand (BED6-like)
    
    ALL genes are processed uniformly - no special handling for short genes.
    With improved sequencing quality (reads as short as 12bp), we can now
    accurately map reads to small genes.
    """
    log(f"INFO  {len(active_genes)} active genes to process")
    
    # Write promoter regions (ppPolII.txt) - TSS-250 to TSS+249 - all active genes
    pp_file = OUT / "ppPolII.txt"
    with open(pp_file, "w") as f:
        for gene in active_genes:
            pps, ppe = clamp(gene['PPs'], gene['PPe'])
            f.write(f"{gene['chrom']}\t{pps}\t{ppe}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(pp_file))
    
    # Write divergent regions (divTx.txt) - TSS-750 to TSS-251 - all active genes
    div_file = OUT / "divTx.txt"
    with open(div_file, "w") as f:
        for gene in active_genes:
            divs, dive = clamp(gene['DIVs'], gene['DIVe'])
            f.write(f"{gene['chrom']}\t{divs}\t{dive}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(div_file))
    
    # Write CPS regions (CPS.txt) - all genes
    cps_file = OUT / "CPS.txt"
    with open(cps_file, "w") as f:
        for gene in active_genes:
            cpss, cpse = clamp(gene['CPSs'], gene['CPSe'])
            f.write(f"{gene['chrom']}\t{cpss}\t{cpse}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(cps_file))
    
    # Write termination windows (TW.txt) - all genes
    tw_file = OUT / "TW.txt"
    with open(tw_file, "w") as f:
        for gene in active_genes:
            tws, twe = clamp(gene['TWs'], gene['TWe'])
            f.write(f"{gene['chrom']}\t{tws}\t{twe}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(tw_file))
    
    # Write gene body regions (geneBody.txt) - all genes where GBe > GBs + 1
    gb_file = OUT / "geneBody.txt"
    with open(gb_file, "w") as f:
        for gene in active_genes:
            gbs, gbe = clamp(gene['GBs'], gene['GBe'])
            if gbe > gbs + 1:  # Only write if body region exists
                f.write(f"{gene['chrom']}\t{gbs}\t{gbe}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(gb_file))
    
    return {
        'promoter': str(pp_file),
        'divergent': str(div_file), 
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
                                    out.write(f"{chrom}\t{start}\t{end}\t{signal}\t{signal}\t+\n")
                            except (ValueError, IndexError):
                                continue
        
        # Process negative strand (keep original signal values like old script)
        if Path(neg_bg).exists() and Path(neg_bg).stat().st_size > 0:
            with open(neg_bg) as f:
                for ln in f:
                    if ln.strip() and not ln.startswith(("track", "browser", "#")):
                        fields = split_fields(ln)
                        if len(fields) >= 4:
                            try:
                                chrom, start, end, signal = fields[0], int(fields[1]), int(fields[2]), float(fields[3])
                                if abs(signal) > thr and end > start:
                                    # Keep original signal (negative values) like old script
                                    out.write(f"{chrom}\t{start}\t{end}\t{signal}\t{signal}\t-\n")
                            except (ValueError, IndexError):
                                continue
    
    # CRITICAL: This is sorting the large read file - use increased buffer
    log(f"INFO  Sorting large read file (this may take a few minutes)...")
    sort_bed(str(reads_file), buffer_size="2G")  # Use 2G buffer for large file
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
    5. Enhancers:    assign unstranded, remove unstranded
    6. Termination:  assign -s, remove UNSTRANDED
    7. Non-localized: remaining
    
    v10.0: Removed "short genes" step - all genes processed uniformly.
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
    # CRITICAL: Remove ALL overlapping reads (no strand flag) like old script
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['promoter']], str(prom_removed))
    assigned_reads['Promoter'] = str(prom_reads)
    current_reads = str(prom_removed)
    log(f"DEBUG Assigned: {wc_effective_lines(str(prom_reads))}, Remaining: {wc_effective_lines(current_reads)}")
    
    # Step 2: Divergent (gene-based TSS-750..-251, opposite-strand)
    # OLD SCRIPT EXACT LOGIC: assign with -S, remove UNSTRANDED
    log("INFO  Step 2: Assigning divergent reads (opposite-strand from divTx)...")
    log(f"DEBUG Input reads: {wc_effective_lines(current_reads)}, Divergent regions: {wc_effective_lines(coord_files['divergent'])}")
    div_reads = OUT / "PROseq_ppDiv.bed"
    div_removed = OUT / "ppdivRemoved.bed"
    run([BT, "intersect", "-S", "-u", "-wa", "-a", current_reads, "-b", coord_files['divergent']], str(div_reads))
    # CRITICAL: Remove ALL overlapping reads (no strand flag) like old script
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['divergent']], str(div_removed))
    assigned_reads['DivergentTx'] = str(div_reads)
    current_reads = str(div_removed)
    log(f"DEBUG Assigned: {wc_effective_lines(str(div_reads))}, Remaining: {wc_effective_lines(current_reads)}")

    
    # Step 3: CPS (same-strand: -s)
    # OLD SCRIPT: assign with -s, remove UNSTRANDED
    log("INFO  Step 3: Assigning CPS reads...")
    cps_reads = OUT / "PROseq_CPS.bed"
    cps_removed = OUT / "ppdivCPSRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['cps']], str(cps_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['cps']], str(cps_removed))
    assigned_reads['CPS'] = str(cps_reads)
    current_reads = str(cps_removed)
    
    # Step 4: Gene Body (same-strand: -s)
    # OLD SCRIPT: assign with -s, remove UNSTRANDED
    log("INFO  Step 4: Assigning gene body reads...")
    gb_reads = OUT / "PROseq_GB.bed"
    gb_removed = OUT / "ppdivCPSgbRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['gene_body']], str(gb_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['gene_body']], str(gb_removed))
    assigned_reads['Gene body'] = str(gb_reads)
    current_reads = str(gb_removed)
    
    # Step 5: Enhancers (unstranded)
    # OLD SCRIPT: assign unstranded, remove unstranded
    log("INFO  Step 5: Assigning enhancer reads...")
    enh_reads = OUT / "PROseq_enhancers.bed"
    enh_removed = OUT / "ppdivCPSgbEnhRemoved.bed"
    run([BT, "intersect", "-u", "-wa", "-a", current_reads, "-b", enhancers_bed], str(enh_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", enhancers_bed], str(enh_removed))
    assigned_reads['Enhancers'] = str(enh_reads)
    current_reads = str(enh_removed)
    
    # Step 6: Termination Window (same-strand: -s)
    # OLD SCRIPT: assign with -s, remove UNSTRANDED
    log("INFO  Step 6: Assigning termination window reads...")
    tw_reads = OUT / "PROseq_TW.bed"
    tw_removed = OUT / "PROseq_noGene_noEnh.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['termination']], str(tw_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['termination']], str(tw_removed))
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
        
        # DEBUG: Keep original values to see what's happening
        # signal_sum = max(0.0, signal_sum)  # Temporarily disabled to debug
        
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
        'Enhancers': str(enhancers_file),  # DT sites NOT overlapping promoters
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
    log(f"INFO  FGR ▶ {SID}  mode=gene-based assignment (v9.1)")

    # Read gene annotations and build coordinate lists
    tss_map = read_sites(args.tss)
    tes_map = read_sites(args.tes)
    all_genes = build_coordinate_lists(args.genes, tss_map, tes_map)
    log(f"INFO  Loaded {len(all_genes)} genes from {args.genes}")

    # Use DT sites to identify active genes and enhancers
    active_genes, genes_with_promoters_bed, enhancers_bed = find_active_promoters_and_enhancers(all_genes, args.divergent)
    
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

    log("INFO  functional regions ✓ (gene-based v9.1)")

if __name__ == "__main__":
    main()