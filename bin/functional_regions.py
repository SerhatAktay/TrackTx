#!/usr/bin/env python3
# =============================================================================
# functional_regions.py — v6.0 (signal-based, active-only, strand-aware)
# -----------------------------------------------------------------------------
# Default geometry (matches TrackTx figure):
#   Promoter:  TSS -250 .. +250
#   Divergent: TSS -750 .. -250  (opposite strand)
#   CPS:       TES -500 .. +500  (fixed)
#   TW:        from CPS end +10000 in transcription direction
#   Gene body: contiguous between Promoter and CPS
#
# Counting modes
#   signal (default): sum bedGraph values inside regions (bedtools map -o sum)
#   event           : legacy per-base 1bp expansion (kept for parity/testing)
#
# Outputs (under --outdir):
#   functional_regions.bed, functional_regions_summary.tsv
#   + minimal aux under outdir/debug when --debug
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

def ensure_bg_abs(in_bg: str, out_bg: str, take_abs: bool):
    """Normalize bedGraph values; make NEG magnitudes positive when take_abs=True."""
    if not (Path(in_bg).exists() and Path(in_bg).stat().st_size>0):
        Path(out_bg).write_text(""); return
    with open(in_bg) as i, open(out_bg,"w") as o:
        for ln in i:
            if not ln.strip() or ln.startswith(("track","browser","#")): continue
            f = split_fields(ln)
            if len(f) < 4: continue
            try:
                s=int(float(f[1])); e=int(float(f[2])); v=float(f[3])
            except Exception:
                continue
            if e<=s: continue
            if abs(v) <= args.min_signal: continue
            if take_abs: v=abs(v)
            o.write(f"{f[0]}\t{s}\t{e}\t{v}\n")
    sort_bed(out_bg)

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

def clamp(a: int, b: int) -> tuple[int,int]:
    """Ensure valid BED coordinates: a <= b and both >= 0"""
    a, b = (a, b) if a <= b else (b, a)
    a = max(0, a)  # Clamp to chromosome start
    b = max(0, b)
    return (a, b)

# ---- build regions -----------------------------------------------------------
def build_regions(genes_tsv: str, tss_map: dict, tes_map: dict):
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

    prom = OUT/"ppPolII.txt"
    divw = OUT/"divTx.txt"
    cps  = OUT/"CPS.txt"
    tw   = OUT/"TW.txt"
    gb   = OUT/"geneBody.txt"
    sg   = OUT/"shortGenes.txt"

    PR, DV, CP, TW, GB, SG = [], [], [], [], [], []
    with open(genes_tsv) as fh:
        next(fh)
        for ln in fh:
            if not ln.strip(): continue
            f=split_fields(ln)
            try:
                gid=f[idxs["gene_id"]]; gname=f[idxs["gene_name"]] or gid
                chrom=f[idxs["chr"]]; strand=f[idxs["strand"]] if f[idxs["strand"]] in {"+","-"} else "+"
                gstart=int(float(f[idxs["start"]])); gend=int(float(f[idxs["end"]]))
                tss_d=int(float(f[idxs["tss"]])); tes_d=int(float(f[idxs["tes"]]))
            except Exception:
                continue

            # overrides (by name or id)
            tss = tss_map.get(gname, tss_map.get(gid, (None,tss_d,None)))[1]
            tes = tes_map.get(gname, tes_map.get(gid, (None,tes_d,None)))[1]

            # promoter
            if strand=="+": p1,p2 = tss-args.prom_up, tss+args.prom_down
            else:           p1,p2 = tss-args.prom_down, tss+args.prom_up
            p1,p2 = clamp(p1,p2); PR.append((chrom,p1,p2,gname,0,strand))

            # divergent window (opposite strand)
            if strand=="+": d1,d2 = tss-args.div_outer, tss-args.div_inner
            else:           d1,d2 = tss+args.div_inner, tss+args.div_outer
            d1,d2 = clamp(d1,d2); DV.append((chrom,d1,d2,gname,0,"-" if strand=="+" else "+"))

            # CPS (TES ±500)
            c1,c2 = tes-500, tes+500; c1,c2 = clamp(c1,c2); CP.append((chrom,c1,c2,gname,0,strand))

            # TW from CPS end in transcription direction
            if strand=="+": t1,t2 = c2, c2+args.tw_length
            else:           t2,t1 = c1, c1-args.tw_length
            t1,t2 = clamp(t1,t2); TW.append((chrom,t1,t2,gname,0,strand))

            # gene body contiguous between promoter and CPS
            if strand=="+": gbs, gbe = p2, c1
            else:           gbs, gbe = c2, p1
            gbs,gbe = clamp(gbs,gbe)
            if gbe>gbs: GB.append((chrom,gbs,gbe,gname,0,strand))

            # short genes (≤ div_outer)
            if (gend-gstart) <= args.div_outer:
                SG.append((chrom,gstart,gend,gname,0,strand))

    def write_bed(path: Path, rows):
        with open(path,"w") as o:
            for r in rows: o.write("\t".join(map(str,r))+"\n")
        sort_bed(str(path))

    write_bed(prom, PR); write_bed(divw, DV); write_bed(cps, CP)
    write_bed(tw, TW);   write_bed(gb, GB);  write_bed(sg, SG)
    return str(prom), str(divw), str(cps), str(tw), str(gb), str(sg)

# ---- active genes & enhancers -----------------------------------------------
def active_genes(prom_bed: str, dt_bed: str) -> str:
    act = OUT/"activeGenes.bed"
    run([BT,"intersect","-u","-wa","-a",prom_bed,"-b",dt_bed,"-sorted"], str(act), check=False)
    return str(act)

def enhancers_from_dt(dt_bed: str, prom_bed: str) -> str:
    enh = OUT/"enhancers.bed"
    run([BT,"intersect","-v","-a",dt_bed,"-b",prom_bed,"-sorted"], str(enh), check=False)
    return str(enh)

# ---- mapping helpers (signal mode) ------------------------------------------
def split_by_strand(src: str, out_plus: str, out_minus: str):
    with open(src) as i, open(out_plus,"w") as p, open(out_minus,"w") as m:
        for ln in i:
            if not ln.strip() or ln.startswith(("track","browser","#")): continue
            f=split_fields(ln)
            if len(f)>=6 and f[5] in {"+","-"}:
                (p if f[5]=="+" else m).write(ln)

def bg_map_sum(a_bed: str, b_bg: str, outp: str):
    # bedtools map keeps A and writes sum of 4th col from B into new 4th col
    if not (Path(a_bed).exists() and Path(a_bed).stat().st_size>0 and
            Path(b_bg).exists() and Path(b_bg).stat().st_size>0):
        # keep A with zero score
        with open(outp,"w") as o:
            with open(a_bed) as a:
                for ln in a:
                    if not ln.strip() or ln.startswith(("track","browser","#")): continue
                    f=split_fields(ln)
                    # Ensure we have at least 4 fields before setting f[3]
                    while len(f) < 4:
                        f.append(".")
                    f[3]="0"
                    o.write("\t".join(f[:6])+"\n")
        return
    tmp = outp+".tmp"
    run([BT,"map","-a",a_bed,"-b",b_bg,"-c","4","-o","sum","-sorted"], tmp, check=False)
    # normalize empty sums to 0
    with open(tmp) as i, open(outp,"w") as o:
        for ln in i:
            f=split_fields(ln)
            # bedtools map appends the computed value as the last column; preserve A's first 6 fields
            # and overwrite column 4 (name/score slot) with the numeric sum for downstream parsing.
            if len(f) < 7:
                # Unexpected: missing appended value; fallback to 0
                val = "0"
            else:
                val = f[-1]
                if val == '.' or val == '':
                    val = "0"
            # ensure we have at least 4 fields to place the value
            while len(f) < 4:
                f.append(".")
            # Ensure we have exactly 6 fields for BED format
            while len(f) < 6:
                f.append(".")
            f[3] = val
            o.write("\t".join(f[:6])+"\n")
    Path(tmp).unlink(missing_ok=True)

def bg_subtract(a_bg: str, b_bed: str, out_bg: str):
    if not Path(a_bg).exists() or Path(a_bg).stat().st_size==0:
        Path(out_bg).write_text(""); return
    if not Path(b_bed).exists() or Path(b_bed).stat().st_size==0:
        shutil.copyfile(a_bg, out_bg); return
    # Use bedtools subtract WITHOUT -A flag to remove only overlapping portions
    # -A flag removes entire features, which is wrong for sequential masking
    # Default behavior subtracts only the overlapping portions, which is correct
    tmp_out = out_bg + ".tmp"
    run([BT,"subtract","-a",a_bg,"-b",b_bed,"-sorted"], tmp_out, check=False)
    # Sort the output to ensure proper downstream processing
    sort_bed(tmp_out)
    shutil.move(tmp_out, out_bg)

# ---- EVENT mode (legacy) ----------------------------------------------------
def event_mode_notebook():
    sys.exit("ERROR: event mode not implemented in this rewrite. Use --count-mode signal (default).")

# ---- MAIN -------------------------------------------------------------------
def main():
    log(f"INFO  FGR ▶ {SID}  mode={args.count_mode}")

    # Prepare bedGraphs (abs for NEG, threshold, sort)
    tmpdir = Path(tempfile.mkdtemp(prefix=f".fgr_{SID}_", dir=str(OUT)))
    pos_bg = tmpdir/"pos.bg"; neg_bg = tmpdir/"neg.abs.bg"
    ensure_bg_abs(args.pos, str(pos_bg), take_abs=False)
    ensure_bg_abs(args.neg, str(neg_bg), take_abs=True)

    # Regions
    tss_map = read_sites(args.tss); tes_map = read_sites(args.tes)
    prom, divw, cps, tw, gb, sg = build_regions(args.genes, tss_map, tes_map)

    # Active genes (must be non-empty)
    dt_sorted = tmpdir/"dt.sorted.bed"
    shutil.copyfile(args.divergent, dt_sorted)
    sort_bed(str(dt_sorted))
    prom_act = active_genes(prom, str(dt_sorted))
    if wc_effective_lines(prom_act) == 0:
        sys.exit("ERROR: active promoter set empty (promoter ∩ DT); refusing to fallback to all genes.")

    # Restrict all region sets to active gene names
    names=set()
    with open(prom_act) as fh:
        for ln in fh:
            f=split_fields(ln)
            if len(f)>=4: names.add(f[3])
    def keep_active(src: str, dst: str):
        if wc_effective_lines(src)==0: Path(dst).write_text(""); return
        with open(src) as i, open(dst,"w") as o:
            for ln in i:
                f=split_fields(ln)
                if len(f)>=4 and f[3] in names:
                    o.write(ln)
        sort_bed(dst)
    promA = str(OUT/"ppPolII.ACT.bed"); keep_active(prom, promA)
    divA  = str(OUT/"divTx.ACT.bed");   keep_active(divw, divA)
    cpsA  = str(OUT/"CPS.ACT.bed");     keep_active(cps,  cpsA)
    twA   = str(OUT/"TW.ACT.bed");      keep_active(tw,   twA)
    gbA   = str(OUT/"geneBody.ACT.bed");keep_active(gb,   gbA)
    sgA   = str(OUT/"shortGenes.ACT.bed"); keep_active(sg, sgA)

    enh = enhancers_from_dt(str(dt_sorted), promA)  # DT minus promoter (unstranded)

    if args.count_mode == "event":
        event_mode_notebook()

    # ---------- SIGNAL MODE (strand-aware with sequential masking) -------------
    # Split region sets by strand
    promP, promM = str(OUT/"pp.plus.bed"), str(OUT/"pp.minus.bed")
    divP,  divM  = str(OUT/"div.plus.bed"), str(OUT/"div.minus.bed")
    cpsP,  cpsM  = str(OUT/"cps.plus.bed"), str(OUT/"cps.minus.bed")
    gbP,   gbM   = str(OUT/"gb.plus.bed") , str(OUT/"gb.minus.bed")
    sgP,   sgM   = str(OUT/"sg.plus.bed") , str(OUT/"sg.minus.bed")
    twP,   twM   = str(OUT/"tw.plus.bed") , str(OUT/"tw.minus.bed")

    for src,(op,om) in [(promA,(promP,promM)),(divA,(divP,divM)),(cpsA,(cpsP,cpsM)),
                        (gbA,(gbP,gbM)),(sgA,(sgP,sgM)),(twA,(twP,twM))]:
        split_by_strand(src, op, om)

    # Residual bedGraphs: used for sequential masking to assign each signal to one category
    pos_res = tmpdir/"pos.res.bg"; neg_res = tmpdir/"neg.res.bg"
    shutil.copyfile(pos_bg, pos_res); shutil.copyfile(neg_bg, neg_res)

    def sum_and_mask_same(region_plus: str, region_minus: str, label: str):
        # same-strand: POS with +, NEG with -
        outP = OUT/f"{label}.plus.sum.bed"
        outM = OUT/f"{label}.minus.sum.bed"
        bg_map_sum(region_plus, str(pos_res), str(outP))
        bg_map_sum(region_minus, str(neg_res), str(outM))
        # Sequential masking: subtract assigned signal from residuals
        # This ensures each signal is counted in exactly ONE functional category
        bg_subtract(str(pos_res), region_plus, str(pos_res))
        bg_subtract(str(neg_res), region_minus, str(neg_res))
        return str(outP), str(outM)

    def sum_opp(region_plus: str, region_minus: str, label: str):
        # opposite-strand: POS with -, NEG with +
        outP = OUT/f"{label}.plus.sum.bed"
        outM = OUT/f"{label}.minus.sum.bed"
        bg_map_sum(region_plus, str(neg_res), str(outP))
        bg_map_sum(region_minus, str(pos_res), str(outM))
        # no masking needed across strands here (regions are upstream vs GB/CPS/TW)
        return str(outP), str(outM)

    # Assign signal hierarchically with sequential masking:
    # 1. Promoter (same-strand) → mask
    ppP, ppM = sum_and_mask_same(promP, promM, "pp")
    # 2. Divergent (opposite-strand) → no masking needed (different strand)
    dvP, dvM = sum_opp(divP, divM, "div")
    # 3. Gene body (same-strand) → mask
    gbp, gbm = sum_and_mask_same(gbP, gbM, "gb")
    # 4. CPS (same-strand) → mask
    cpP, cpM = sum_and_mask_same(cpsP, cpsM, "cps")
    # 5. Short genes (same-strand) → mask
    sgp, sgm = sum_and_mask_same(sgP, sgM, "sg")
    # 6. Enhancers (unstranded DT - promoter): sum both strands from current residuals
    enh_sum = OUT/"enh.sum.bed"
    # Build enh with strand '.' and carry name
    if wc_effective_lines(enh)==0:
        Path(enh_sum).write_text("")
    else:
        # Map POS residual then NEG residual; merge sums
        pos_enh = OUT/"enh.pos.map.bed"; bg_map_sum(enh, str(pos_res), str(pos_enh))
        neg_enh = OUT/"enh.neg.map.bed"; bg_map_sum(enh, str(neg_res), str(neg_enh))
        with open(pos_enh) as a, open(neg_enh) as b, open(enh_sum,"w") as o:
            for la,lb in zip(a,b):
                fa,fb=split_fields(la),split_fields(lb)
                # fa[3], fb[3] are sums; add
                s = (0.0 if fa[3]=='.' else float(fa[3])) + (0.0 if fb[3]=='.' else float(fb[3]))
                o.write("\t".join([fa[0],fa[1],fa[2],f"{s:.6f}",fa[5] if len(fa)>=6 else "."])+"\n")
        # Sequential masking: enhancers are assigned from residuals
        bg_subtract(str(pos_res), enh, str(pos_res))
        bg_subtract(str(neg_res), enh, str(neg_res))

    # 7. Termination window (same-strand) → mask
    twp, twm = sum_and_mask_same(twP, twM, "tw")
    # 8. Non-localized: total remaining signal in residual bedGraphs after all masking
    rest_sum = 0.0
    for p in (pos_res, neg_res):
        if Path(p).exists() and Path(p).stat().st_size>0:
            with open(p) as fh:
                for ln in fh:
                    if ln.startswith(("track","browser","#")) or not ln.strip(): continue
                    f=split_fields(ln)
                    rest_sum += (float(f[3]) * max(0, int(f[2])-int(f[1])) / max(1, int(f[2])-int(f[1])))
                    # note: bedGraph is already integrated per span; here we just sum values (not length-weighted)

    # ---- Summaries + final BED9 -------------------------------------------
    # Collapse per-class sums and counts
    def total_and_count_from(mapped_plus: str, mapped_minus: str) -> tuple[float, int]:
        s=0.0
        count=0
        for p in (mapped_plus, mapped_minus):
            if not Path(p).exists(): continue
            with open(p) as fh:
                for ln in fh:
                    f=split_fields(ln)
                    if len(f)>=4 and f[3] not in (".",""):
                        try:
                            val = float(f[3])
                            s += val
                            count += 1
                        except ValueError:
                            continue
        return s, count

    # Calculate signal and region counts
    prom_signal, prom_count = total_and_count_from(ppP, ppM)
    div_signal, div_count = total_and_count_from(dvP, dvM)
    cps_signal, cps_count = total_and_count_from(cpP, cpM)
    gb_signal, gb_count = total_and_count_from(gbp, gbm)
    sg_signal, sg_count = total_and_count_from(sgp, sgm)
    tw_signal, tw_count = total_and_count_from(twp, twm)
    
    # Handle enhancers separately
    enh_signal = 0.0
    enh_count = 0
    if wc_effective_lines(enh_sum) > 0:
        with open(enh_sum) as fh:
            for ln in fh:
                f = split_fields(ln)
                if len(f) >= 4 and f[3] not in (".", ""):
                    try:
                        enh_signal += float(f[3])
                        enh_count += 1
                    except ValueError:
                        continue
    
    sums = {
        "Promoter" : prom_signal,
        "DivergentTx": div_signal,
        "CPS"      : cps_signal,
        "Gene body": gb_signal,
        "Short genes": sg_signal,
        "Enhancers": enh_signal,
        "Termination window": tw_signal,
        "Non-localized polymerase": rest_sum,
    }
    
    counts = {
        "Promoter" : prom_count,
        "DivergentTx": div_count,
        "CPS"      : cps_count,
        "Gene body": gb_count,
        "Short genes": sg_count,
        "Enhancers": enh_count,
        "Termination window": tw_count,
        "Non-localized polymerase": 0,  # This is residual signal, not discrete regions
    }

    # Write functional_regions_summary.tsv
    with open(OUT/"functional_regions_summary.tsv","w") as o:
        o.write("region\tsignal\tregion_count\n")
        for k in sums.keys():
            signal = sums[k]
            count = counts[k]
            o.write(f"{k}\t{signal:.6f}\t{count}\n")

    # Build BED9 track by stitching region BEDs with per-feature counts (-c)
    # (Counts not needed here for browser view; we still attach a score = 0)
    header = OUT/"headerLine.txt"
    header.write_text(f'track name="functional_genomic_regions_{SID}" itemRgb="On"\n')
    final = OUT/"functional_regions.bed"
    with open(final,"w") as out:
        with open(header) as h: shutil.copyfileobj(h, out)
        # Concatenate region defs (active) with a dummy score and rgb per class
        rgb = dict(ppPol="243,132,0", ppDiv="178,59,212", enh="115,212,122",
                   GB="0,0,0", CPS="103,200,249", SG="253,218,13", TW="255,54,98")
        def add(path,label,color):
            if wc_effective_lines(path)==0: return
            with open(path) as fh:
                for ln in fh:
                    f=split_fields(ln);  # chr start end name . strand
                    name = f[3] if len(f)>=4 else label
                    strand = f[5] if len(f)>=6 else "."
                    out.write("\t".join([f[0],f[1],f[2],name,"0",strand,f[1],f[2],color])+"\n")
        add(promA,"Promoter",rgb["ppPol"])
        add(divA,"DivergentTx",rgb["ppDiv"])
        add(cpsA,"CPS",rgb["CPS"])
        add(gbA,"Gene body",rgb["GB"])
        add(sgA,"Short genes",rgb["SG"])
        add(enh,"Enhancer",rgb["enh"])
        add(twA,"Termination window",rgb["TW"])

    log("INFO  functional regions ✓")
    shutil.rmtree(tmpdir, ignore_errors=True)

if __name__ == "__main__":
    main()
