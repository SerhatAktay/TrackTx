#!/usr/bin/env python3
# =============================================================================
# normalize_tracks.py — CPM, Spike-in (cRPMsi), and RPKM scaling for bedGraphs
# -----------------------------------------------------------------------------
# Math:
#   CPM      = score * (1e6 / N_sample)
#   cRPMsi   = score * (si_control / si_sample) * (1e6 / N_ctrl_baseline)
#              (baseline picked per condition: tp=0 & rep=1 preferred; else tp=0
#               (min rep); else min(tp) then min(rep). Baseline sample → fac_si=1.0)
#   RPKM     = score * (1e9 / (N_sample * length_bp))   # per-interval
#
# I/O:
#   Inputs:
#     --sample --condition --timepoint --replicate
#     --counts counts_master.tsv     (TSV: sample, condition, timepoint, replicate,
#                                     spike_reads, mapped_reads)
#     --bam    sample.bam            (for N_sample via samtools idxstats)
#     --fasta  genome.fa             (.fai is used or built via samtools faidx)
#     Optional bedGraphs: --pos3 --neg3 --pos5 --neg5  (.bg or .bedgraph; gz OK)
#     --make-bw (optional)           (writes BigWig if bedGraphToBigWig available)
#
#   Outputs (cwd):
#     <SID>_<label>_cpm.bedgraph     (sorted)
#     <SID>_<label>_cpm_si.bedgraph  (sorted; empty if factors unavailable)
#     <SID>_<label>_rpkm.bedgraph    (sorted)
#     *_cpm.bw, *_cpm_si.bw, *_rpkm.bw (optional)
#     normalization_factors.tsv      (key/value summary)
#
# Implementation notes:
#   • All processing is STREAMED; no large Python lists (low memory).
#   • bedGraphs are clipped to chrom sizes and then sorted with GNU sort
#     (tunable via --sort-mem and --tmpdir). Sorting is required for BigWig.
#   • NEG labels (neg3/neg5) are mirrored (value := -value) before scaling.
#   • Requires: samtools; optional: bedGraphToBigWig.
# =============================================================================

from __future__ import annotations
import argparse, csv, os, re, subprocess as sp, sys, shutil, tempfile, gzip, datetime
from typing import Dict, Iterable, Tuple, Optional, Iterator

# ────────────────────────────────────────────────────────────────────────────
# 0) CLI
# ────────────────────────────────────────────────────────────────────────────

def get_args():
    p = argparse.ArgumentParser(description="Normalize bedGraph tracks to CPM, cRPMsi, and RPKM (streaming).")
    p.add_argument("--sample", required=True)
    p.add_argument("--condition", required=True)
    p.add_argument("--timepoint", required=True)
    p.add_argument("--replicate", required=True)
    p.add_argument("--counts", required=True, help="counts_master.tsv")
    p.add_argument("--bam", required=True, help="Aligned BAM (indexed or indexable)")
    p.add_argument("--fasta", required=True, help="Reference FASTA (needs .fai)")
    p.add_argument("--pos3"); p.add_argument("--neg3"); p.add_argument("--pos5"); p.add_argument("--neg5")
    p.add_argument("--make-bw", action="store_true", help="Also write BigWig if bedGraphToBigWig exists")
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--sort-mem", default="25%", help="GNU sort memory (e.g. 2G or 25%)")
    p.add_argument("--tmpdir", default=None, help="Temporary directory for sorting")
    return p.parse_args()

# ────────────────────────────────────────────────────────────────────────────
# 1) Small utilities (tools, run, IO)
# ────────────────────────────────────────────────────────────────────────────

def have(cmd: str) -> bool:
    return shutil.which(cmd) is not None

def run(cmd: list[str], check=True, capture=False, text=True) -> sp.CompletedProcess:
    return sp.run(cmd, check=check, capture_output=capture, text=text)

def open_text(path: str) -> Iterator[str]:
    """Open plain or gz text as an iterator of lines."""
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            for ln in f: yield ln
    else:
        with open(path, "r") as f:
            for ln in f: yield ln

# ────────────────────────────────────────────────────────────────────────────
# 2) Counts/metadata + baseline selection
# ────────────────────────────────────────────────────────────────────────────

def tp_to_minutes(x: str) -> int:
    s = str(x).strip().lower()
    if s in {"0","0m","0min","0h"}: return 0
    m = re.match(r"^\s*(\d+)\s*(m|min)\s*$", s)
    h = re.match(r"^\s*([\d.]+)\s*h\s*$", s)
    if m: return int(m.group(1))
    if h: return int(float(h.group(1))*60)
    try: return int(float(s))
    except Exception: return 10**9

def to_int(x, default=10**6) -> int:
    try: return int(x)
    except Exception:
        m = re.search(r"(\d+)", str(x) if x is not None else "")
        return int(m.group(1)) if m else default

def read_counts_table(path: str) -> list[dict]:
    rows = []
    with open(path, newline="") as f:
        for d in csv.DictReader(f, delimiter="\t"):
            if not d or not d.get("sample"): continue
            d = {k.strip(): v for k, v in d.items()}
            d["condition_l"] = str(d.get("condition","")).lower()
            d["tp_min"] = tp_to_minutes(d.get("timepoint",""))
            d["rep_i"]  = to_int(d.get("replicate",None))
            d["spike"]  = int(d.get("spike_reads",0) or 0)
            d["mapped"] = int(d.get("mapped_reads",0) or 0)
            rows.append(d)
    if not rows:
        sys.exit(f"ERROR: empty/invalid counts table: {path}")
    return rows

def pick_baseline(rows: list[dict], condition: str) -> dict:
    cond_l = str(condition).lower()
    crows = [r for r in rows if r["condition_l"] == cond_l] or rows
    base = [r for r in crows if r["tp_min"]==0 and r["rep_i"]==1] or \
           [r for r in crows if r["tp_min"]==0]
    if not base:
        min_tp = min(r["tp_min"] for r in crows)
        base = [r for r in crows if r["tp_min"]==min_tp]
    return sorted(base, key=lambda r: r["rep_i"])[0]

def compute_spike_factors(rows: list[dict], sid: str, condition: str) -> Tuple[float,int,str]:
    base = pick_baseline(rows, condition)
    baseline_id = base.get("sample","NA")
    si_ctrl = base["spike"]
    n_ctrl  = base["mapped"]
    samp = next((r for r in rows if str(r.get("sample","")) == sid), None)
    si_samp = samp["spike"] if samp else 0
    if not samp or si_ctrl <= 0 or si_samp <= 0:
        fac_si = 0.0                         # disable cRPMsi if missing info
    elif str(samp.get("sample","")) == baseline_id:
        fac_si = 1.0
    else:
        fac_si = float(si_ctrl) / float(si_samp)
    return fac_si, (n_ctrl if n_ctrl>0 else 0), baseline_id

# ────────────────────────────────────────────────────────────────────────────
# 3) BAM & reference helpers
# ────────────────────────────────────────────────────────────────────────────

def ensure_bam_index(bam: str, threads: int=2):
    if os.path.exists(bam + ".bai") or os.path.exists(bam + ".csi"): return
    if not have("samtools"): sys.exit("ERROR: samtools not found")
    hdr = run(["samtools","view","-H", bam], capture=True).stdout
    sq = sum(1 for line in hdr.splitlines() if line.startswith("@SQ"))
    if sq > 512:
        run(["samtools","index","-@", str(threads), "-c", bam])
    else:
        try:
            run(["samtools","index","-@", str(threads), bam])
        except sp.CalledProcessError:
            run(["samtools","index","-@", str(threads), "-c", bam])

def mapped_reads(bam: str) -> int:
    ensure_bam_index(bam)
    out = run(["samtools","idxstats", bam], capture=True).stdout
    m = 0
    for line in out.splitlines():
        if not line or line.startswith("*\t"): continue
        parts = line.split("\t")
        if len(parts) >= 3:
            try: m += int(parts[2])
            except: pass
    return m

def ensure_fai(fasta: str):
    fai = fasta + ".fai"
    if os.path.exists(fai): return
    if not have("samtools"): sys.exit("ERROR: samtools not found to build .fai")
    run(["samtools","faidx", fasta])

def load_chrom_sizes(fasta: str) -> Dict[str,int]:
    ensure_fai(fasta)
    sizes = {}
    with open(fasta + ".fai") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                sizes[p[0]] = int(p[1])
    if not sizes:
        sys.exit("ERROR: failed to read chrom sizes from .fai")
    return sizes

# ────────────────────────────────────────────────────────────────────────────
# 4) bedGraph processing (streaming + external sort)
# ────────────────────────────────────────────────────────────────────────────

def iter_bg_lines(path: str) -> Iterable[Tuple[str,int,int,float]]:
    for ln in open_text(path):
        if not ln.strip() or ln.startswith(("track","browser","#")): continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 4: continue
        try:
            chrom = p[0]; start = int(p[1]); end = int(p[2]); val = float(p[3])
        except Exception:
            continue
        yield chrom, start, end, val

def write_bg_stream(src: Iterable[Tuple[str,int,int,float]],
                    out_path: str,
                    sizes: Dict[str,int],
                    mirror_neg: bool,
                    scale: float,
                    rpkm_const: Optional[float] = None):
    """
    Stream input rows → clip → scale → write UNSORTED bedGraph.
    If rpkm_const is not None, write an RPKM file instead of CPM (out_path is the target).
    """
    with open(out_path, "w") as out:
        for chrom, s, e, v in src:
            mlen = sizes.get(chrom, 0)
            if mlen <= 0: continue
            s2 = 0 if s < 0 else s
            e2 = e if e <= mlen else mlen
            if e2 <= s2: continue
            val = -v if mirror_neg else v
            if rpkm_const is None:
                out.write(f"{chrom}\t{s2}\t{e2}\t{val * scale}\n")
            else:
                length = e2 - s2
                if length > 0:
                    out.write(f"{chrom}\t{s2}\t{e2}\t{val * (rpkm_const / length)}\n")

def sort_bg_inplace(path: str, sort_mem: str = "25%", tmpdir: Optional[str] = None):
    """
    Sort bedGraph by chrom,start using GNU sort (low memory, external).
    """
    if not have("sort"):
        # Fallback in-memory (only if sort is missing)
        rows = []
        with open(path) as f:
            for ln in f:
                p = ln.rstrip("\n").split("\t")
                if len(p) == 4:
                    rows.append((p[0], int(p[1]), int(p[2]), float(p[3])))
        rows.sort(key=lambda x: (x[0], x[1]))
        with open(path, "w") as out:
            for c,s,e,v in rows:
                out.write(f"{c}\t{s}\t{e}\t{v}\n")
        return
    tmp = path + ".sorted.tmp"
    cmd = ["sort", "-k1,1", "-k2,2n", "--stable", f"--buffer-size={sort_mem}"]
    if tmpdir:
        cmd.append(f"--temporary-directory={tmpdir}")
    cmd.append(path)
    with open(tmp, "w") as w:
        sp.run(cmd, check=True, stdout=w)
    os.replace(tmp, path)

def make_bw(sorted_bg: str, sizes: Dict[str,int], out_bw: str):
    """
    Convert sorted bedGraph to BigWig if tool exists. Else emit stub.
    """
    if not have("bedGraphToBigWig"):
        open(out_bw, "w").close()
        return
    with tempfile.NamedTemporaryFile("w", delete=False) as t:
        for k,v in sizes.items():
            t.write(f"{k}\t{v}\n")
        sizes_path = t.name
    try:
        sp.run(["bedGraphToBigWig", sorted_bg, sizes_path, out_bw], check=False)
        if not os.path.exists(out_bw):
            open(out_bw,"w").close()
    finally:
        os.unlink(sizes_path)

# ────────────────────────────────────────────────────────────────────────────
# 5) Main
# ────────────────────────────────────────────────────────────────────────────

def main():
    a = get_args()
    print(f"[normalize_py] start ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)

    # ── 5.1) Factors
    rows = read_counts_table(a.counts)
    fac_si, n_ctrl, baseline_id = compute_spike_factors(rows, a.sample, a.condition)

    n_sample = mapped_reads(a.bam)
    if n_sample <= 0:
        sys.exit("ERROR: no mapped reads in BAM (N_sample==0)")

    fac_cpm_sample = 1_000_000.0 / n_sample
    fac_cpm_ctrl   = (1_000_000.0 / n_ctrl) if n_ctrl > 0 else 0.0
    fac_rpkm_const = 1_000_000_000.0 / n_sample  # 1e9 / N_sample

    sizes = load_chrom_sizes(a.fasta)

    # ── 5.2) Per-track streaming pipeline
    tracks = [
        ("pos3", a.pos3, False),
        ("neg3", a.neg3, True),
        ("pos5", a.pos5, False),
        ("neg5", a.neg5, True),
    ]

    wrote_any = False
    for label, in_bg, is_neg in tracks:
        if not in_bg: continue
        if not os.path.exists(in_bg):
            print(f"WARN: missing {in_bg} ({label})", file=sys.stderr); continue

        # Sources are iterators; we re-open per output to keep streaming & avoid RAM growth.
        # 1) CPM
        out_cpm = f"{a.sample}_{label}_cpm.bedgraph"
        write_bg_stream(iter_bg_lines(in_bg), out_cpm, sizes, is_neg, fac_cpm_sample, rpkm_const=None)
        sort_bg_inplace(out_cpm, a.sort_mem, a.tmpdir)

        # 2) cRPMsi (si*controlCPM) only when info present
        out_si = f"{a.sample}_{label}_cpm_si.bedgraph"
        if fac_si > 0 and n_ctrl > 0:
            write_bg_stream(iter_bg_lines(in_bg), out_si, sizes, is_neg, fac_si * fac_cpm_ctrl, rpkm_const=None)
            sort_bg_inplace(out_si, a.sort_mem, a.tmpdir)
        else:
            open(out_si, "w").close()

        # 3) RPKM
        out_rpkm = f"{a.sample}_{label}_rpkm.bedgraph"
        write_bg_stream(iter_bg_lines(in_bg), out_rpkm, sizes, is_neg, scale=0.0, rpkm_const=fac_rpkm_const)
        sort_bg_inplace(out_rpkm, a.sort_mem, a.tmpdir)

        # Optional BigWigs
        if a.make_bw:
            make_bw(out_cpm,  sizes, f"{a.sample}_{label}_cpm.bw")
            make_bw(out_si,   sizes, f"{a.sample}_{label}_cpm_si.bw")
            make_bw(out_rpkm, sizes, f"{a.sample}_{label}_rpkm.bw")

        wrote_any = True
        print(f"INFO: wrote {out_cpm}, {out_si}, {out_rpkm}", file=sys.stderr)

    if not wrote_any:
        print("WARN: no input tracks provided", file=sys.stderr)

    # ── 5.3) Factors file
    with open("normalization_factors.tsv","w") as f:
        w = f.write
        w("key\tvalue\n")
        w(f"sample\t{a.sample}\n")
        w(f"condition\t{a.condition}\n")
        w(f"timepoint\t{a.timepoint}\n")
        w(f"replicate\t{a.replicate}\n")
        w(f"baseline_sample\t{baseline_id}\n")
        w(f"N_sample\t{n_sample}\n")
        w(f"fac_cpm_sample\t{fac_cpm_sample:.10f}\n")
        w(f"N_ctrl_baseline\t{n_ctrl}\n")
        w(f"fac_cpm_control\t{fac_cpm_ctrl:.10f}\n")
        w(f"fac_si\t{fac_si:.6f}\n")
    print(f"[normalize_py] done ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)

if __name__ == "__main__":
    # Defensive locale for sort/bw tools
    os.environ.setdefault("LC_ALL", "C")
    main()
