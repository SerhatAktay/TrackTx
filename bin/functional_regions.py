#!/usr/bin/env python3
# =============================================================================
# functional_regions.py — v7.0 (read-based, active-only, strand-aware)
# -----------------------------------------------------------------------------
# NEW APPROACH: Read-based assignment following the old TrackTx.sh logic
# 
# Key differences from v6.0:
#   - Uses read-based assignment instead of signal-based
#   - Converts bedGraphs to read intervals first
#   - Uses sequential masking with bedtools intersect -v on reads
#   - Follows exact logic from old TrackTx.sh script
#   - Generates coordinate lists like the old script
#
# Default geometry (matches old TrackTx script):
#   Promoter:  TSS -250 .. +249 (for +) / TSS -249 .. +250 (for -)
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

# Geometry
ap.add_argument("--prom-up",    type=int, default=250)
ap.add_argument("--prom-down",  type=int, default=250)
ap.add_argument("--div-inner",  type=int, default=250)
ap.add_argument("--div-outer",  type=int, default=750)
ap.add_argument("--tw-length",  type=int, default=10_000)

# Behavior
ap.add_argument("--min-signal",      type=float, default=0.0)
ap.add_argument("--allow-unstranded",action="store_true")
ap.add_argument("--count-mode",      choices=["signal","event"], default="signal")
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
            subprocess.run(cmd, check=check, text=True, stdout=o, stderr=subprocess.PIPE)
    else:
        subprocess.run(cmd, check=check, text=True, stderr=subprocess.PIPE)

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

def clamp(a: int, b: int) -> tuple[int,int]:
    """Ensure valid BED coordinates: a <= b and both >= 0"""
    a, b = (a, b) if a <= b else (b, a)
    a = max(0, a)  # Clamp to chromosome start
    b = max(0, b)
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
def find_active_genes_and_enhancers(all_genes: list, dt_bed: str) -> tuple[list, str]:
    """
    Following old script logic:
    1. Create TSS±500 regions for all genes
    2. Find enhancers: divergent_transcription MINUS TSS±500
    3. Find active genes: TSS±500 INTERSECT divergent_transcription
    4. Filter all_genes to active genes only
    """
    # Create TSS±500 regions (refGenes_TSSpm500.txt equivalent)
    tss_pm500 = OUT / "refGenes_TSSpm500.bed"
    with open(tss_pm500, "w") as f:
        for gene in all_genes:
            tss_start = gene['TSS'] - 500
            tss_end = gene['TSS'] + 500
            tss_start, tss_end = clamp(tss_start, tss_end)
            f.write(f"{gene['chrom']}\t{tss_start}\t{tss_end}\t{gene['gname']}\n")
    
    sort_bed(str(tss_pm500))
    
    # Find enhancers: divergent_transcription MINUS TSS±500
    enhancers_bed = OUT / "enhancers.bed"
    run([BT, "intersect", "-v", "-a", dt_bed, "-b", str(tss_pm500)], str(enhancers_bed))
    
    # Find active genes: TSS±500 INTERSECT divergent_transcription  
    active_genes_bed = OUT / "activeGenes.bed"
    run([BT, "intersect", "-u", "-wa", "-a", str(tss_pm500), "-b", dt_bed], str(active_genes_bed))
    
    # Read active gene names
    active_names = set()
    if Path(active_genes_bed).exists():
        with open(active_genes_bed) as f:
            for ln in f:
                if ln.strip() and not ln.startswith(("#", "track", "browser")):
                    fields = split_fields(ln)
                    if len(fields) >= 4:
                        active_names.add(fields[3])
    
    # Filter to active genes only
    active_genes = [g for g in all_genes if g['gname'] in active_names]
    
    log(f"INFO  Found {len(active_genes)} active genes out of {len(all_genes)} total genes")
    
    return active_genes, str(enhancers_bed)

# ---- coordinate file generation ---------------------------------------------
def write_coordinate_files(active_genes: list):
    """
    Write coordinate files following old script format:
    chr, start, end, name, name, strand (BED6-like)
    """
    # Separate short and long genes (following old script: <= 750bp is short)
    short_genes = [g for g in active_genes if g['gene_length'] <= 750]
    long_genes = [g for g in active_genes if g['gene_length'] > 750]
    
    log(f"INFO  {len(short_genes)} short genes, {len(long_genes)} long genes")
    
    # Write promoter regions (ppPolII.txt) - all active genes
    prom_file = OUT / "ppPolII.txt"
    with open(prom_file, "w") as f:
        for gene in active_genes:
            pps, ppe = clamp(gene['PPs'], gene['PPe'])
            f.write(f"{gene['chrom']}\t{pps}\t{ppe}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(prom_file))
    
    # Write divergent regions (divTx.txt) - all active genes
    div_file = OUT / "divTx.txt"
    with open(div_file, "w") as f:
        for gene in active_genes:
            divs, dive = clamp(gene['DIVs'], gene['DIVe'])
            f.write(f"{gene['chrom']}\t{divs}\t{dive}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(div_file))
    
    # Write short genes (shortGenes.txt) - full gene coordinates
    sg_file = OUT / "shortGenes.txt"
    with open(sg_file, "w") as f:
        for gene in short_genes:
            gstart, gend = clamp(gene['gstart'], gene['gend'])
            f.write(f"{gene['chrom']}\t{gstart}\t{gend}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(sg_file))
    
    # Write CPS regions (CPS.txt) - long genes only
    cps_file = OUT / "CPS.txt"
    with open(cps_file, "w") as f:
        for gene in long_genes:
            cpss, cpse = clamp(gene['CPSs'], gene['CPSe'])
            f.write(f"{gene['chrom']}\t{cpss}\t{cpse}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(cps_file))
    
    # Write termination windows (TW.txt) - long genes only
    tw_file = OUT / "TW.txt"
    with open(tw_file, "w") as f:
        for gene in long_genes:
            tws, twe = clamp(gene['TWs'], gene['TWe'])
            f.write(f"{gene['chrom']}\t{tws}\t{twe}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(tw_file))
    
    # Write gene body regions (geneBody.txt) - long genes only, where GBe > GBs + 1
    gb_file = OUT / "geneBody.txt"
    with open(gb_file, "w") as f:
        for gene in long_genes:
            gbs, gbe = clamp(gene['GBs'], gene['GBe'])
            if gbe > gbs + 1:  # Following old script logic
                f.write(f"{gene['chrom']}\t{gbs}\t{gbe}\t{gene['gname']}\t{gene['gname']}\t{gene['strand']}\n")
    sort_bed(str(gb_file))
    
    return {
        'promoter': str(prom_file),
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
                                if abs(signal) > args.min_signal and end > start:
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
                                if abs(signal) > args.min_signal and end > start:
                                    # Keep original signal (negative values) like old script
                                    out.write(f"{chrom}\t{start}\t{end}\t{signal}\t{signal}\t-\n")
                            except (ValueError, IndexError):
                                continue
    
    sort_bed(str(reads_file))
    return str(reads_file)

# ---- sequential read assignment (following old script exactly) -------------
def sequential_read_assignment(reads_file: str, coord_files: dict, enhancers_bed: str):
    """
    Perform sequential read assignment following old TrackTx.sh logic exactly.
    Each step removes assigned reads from the pool.
    """
    current_reads = reads_file
    assigned_reads = {}
    
    # Step 1: Promoters (same-strand: -s)
    log("INFO  Step 1: Assigning promoter reads...")
    log(f"DEBUG Input reads: {wc_effective_lines(current_reads)}, Promoter regions: {wc_effective_lines(coord_files['promoter'])}")
    prom_reads = OUT / "PROseq_ppPolII.bed"
    prom_removed = OUT / "ppRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['promoter']], str(prom_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['promoter']], str(prom_removed))
    assigned_reads['Promoter'] = str(prom_reads)
    current_reads = str(prom_removed)
    log(f"DEBUG Assigned: {wc_effective_lines(str(prom_reads))}, Remaining: {wc_effective_lines(current_reads)}")
    
    # Step 2: Divergent (opposite-strand: -S)  
    log("INFO  Step 2: Assigning divergent reads...")
    log(f"DEBUG Input reads: {wc_effective_lines(current_reads)}, Divergent regions: {wc_effective_lines(coord_files['divergent'])}")
    div_reads = OUT / "PROseq_ppDiv.bed"
    div_removed = OUT / "ppdivRemoved.bed"
    run([BT, "intersect", "-S", "-u", "-wa", "-a", current_reads, "-b", coord_files['divergent']], str(div_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['divergent']], str(div_removed))
    assigned_reads['DivergentTx'] = str(div_reads)
    current_reads = str(div_removed)
    log(f"DEBUG Assigned: {wc_effective_lines(str(div_reads))}, Remaining: {wc_effective_lines(current_reads)}")
    
    # Step 3: CPS (same-strand: -s)
    log("INFO  Step 3: Assigning CPS reads...")
    cps_reads = OUT / "PROseq_CPS.bed"
    cps_removed = OUT / "ppdivCPSRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['cps']], str(cps_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['cps']], str(cps_removed))
    assigned_reads['CPS'] = str(cps_reads)
    current_reads = str(cps_removed)
    
    # Step 4: Gene Body (same-strand: -s)
    log("INFO  Step 4: Assigning gene body reads...")
    gb_reads = OUT / "PROseq_GB.bed"
    gb_removed = OUT / "ppdivCPSgbRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['gene_body']], str(gb_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['gene_body']], str(gb_removed))
    assigned_reads['Gene body'] = str(gb_reads)
    current_reads = str(gb_removed)
    
    # Step 5: Short Genes (same-strand: -s)
    log("INFO  Step 5: Assigning short gene reads...")
    sg_reads = OUT / "PROseq_SG.bed"
    sg_removed = OUT / "ppdivCPSgbSGRemoved.bed"
    run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", coord_files['short_genes']], str(sg_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", coord_files['short_genes']], str(sg_removed))
    assigned_reads['Short genes'] = str(sg_reads)
    current_reads = str(sg_removed)
    
    # Step 6: Enhancers (unstranded - no -s or -S flag)
    log("INFO  Step 6: Assigning enhancer reads...")
    enh_reads = OUT / "PROseq_enhancers.bed"
    enh_removed = OUT / "ppdivCPSgbSGEnhRemoved.bed"
    run([BT, "intersect", "-u", "-wa", "-a", current_reads, "-b", enhancers_bed], str(enh_reads))
    run([BT, "intersect", "-v", "-a", current_reads, "-b", enhancers_bed], str(enh_removed))
    assigned_reads['Enhancers'] = str(enh_reads)
    current_reads = str(enh_removed)
    
    # Step 7: Termination Window (same-strand: -s)
    log("INFO  Step 7: Assigning termination window reads...")
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
    """Write final outputs"""
    
    # Write summary TSV (use 'reads' column name for compatibility with reporting script)
    summary_file = OUT / "functional_regions_summary.tsv"
    with open(summary_file, "w") as f:
        f.write("region\treads\tregion_count\n")
        for category, data in summary.items():
            f.write(f"{category}\t{data['signal']:.6f}\t{data['count']}\n")
    
    # Write BED9 track (following old script format)
    bed_file = OUT / "functional_regions.bed"
    with open(bed_file, "w") as f:
        f.write(f'track name="functional_genomic_regions_{SID}" itemRgb="On"\n')
        
        # Color mapping (from old script)
        colors = {
            'promoter': "243,132,0",
            'divergent': "178,59,212", 
            'cps': "103,200,249",
            'gene_body': "0,0,0",
            'short_genes': "253,218,13",
            'termination': "255,54,98"
        }
        
        # Add coordinate regions with colors
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
        
        # Add enhancers
        enhancers_file = OUT / "enhancers.bed"
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
    log(f"INFO  FGR ▶ {SID}  mode=read-based (v7.0)")

    # Read gene annotations and build coordinate lists
    tss_map = read_sites(args.tss)
    tes_map = read_sites(args.tes)
    all_genes = build_coordinate_lists(args.genes, tss_map, tes_map)
    log(f"INFO  Loaded {len(all_genes)} genes from {args.genes}")

    # Find active genes and enhancers
    active_genes, enhancers_bed = find_active_genes_and_enhancers(all_genes, args.divergent)
    
    if not active_genes:
        sys.exit("ERROR: No active genes found (promoter ∩ DT); check divergent transcription input")

    # Write coordinate files
    coord_files = write_coordinate_files(active_genes)

    # Convert bedGraphs to read intervals
    reads_file = bedgraph_to_reads(args.pos, args.neg)
    log(f"INFO  Converted bedGraphs to {wc_effective_lines(reads_file)} read intervals")

    # Perform sequential read assignment
    assigned_reads = sequential_read_assignment(reads_file, coord_files, enhancers_bed)

    # Count and summarize
    summary = count_and_summarize(assigned_reads)

    # Write outputs
    write_outputs(summary, coord_files, assigned_reads)

    log("INFO  functional regions ✓ (read-based)")

if __name__ == "__main__":
    main()
