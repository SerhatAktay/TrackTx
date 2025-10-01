#!/usr/bin/env python3
# =============================================================================
# pol2_pausing_metrics.py — per-sample Pol II pausing + gene metrics (fast)
# -----------------------------------------------------------------------------
# • Builds TSS ±win and Body windows from GTF (feature types: gene,transcript…)
# • Counts reads with `bedtools coverage -counts` against a prefiltered BAM
# • Emits:
#     - pausing_index.tsv  (lean)
#     - pol2_gene_metrics.tsv (rich; CPM, per-bp densities, normalized PI)
#     - pol2_qc.json
# =============================================================================

from __future__ import annotations
import argparse, csv, json, os, subprocess, tempfile, datetime, sys
from pathlib import Path

def run(cmd: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(cmd, text=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

ap = argparse.ArgumentParser()
ap.add_argument("--bam", required=True)
ap.add_argument("--gtf", required=True)
ap.add_argument("--tss-win", type=int, default=50)
ap.add_argument("--body-offset-min", type=int, default=2000)
ap.add_argument("--body-offset-frac", type=float, default=0.10)
ap.add_argument("--feature-types", default="gene,transcript")
ap.add_argument("--out-pausing", required=True)
ap.add_argument("--out-genes",   required=True)
ap.add_argument("--out-qc",      required=True)
ap.add_argument("--threads", type=int, default=1)
ap.add_argument("--fail-if-empty", default="false")
args = ap.parse_args()
print(f"[pol2_calc] start ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)

# ---- helpers ---------------------------------------------------------------
def parse_attrs(attr: str) -> dict:
    d={}
    for k in ("gene_id","gene_name","ID","Name"):
        import re
        m = re.search(rf'{k}\s+"([^"]+)"', attr) or re.search(rf'{k}=([^;]+)', attr)
        if m: d[k]=m.group(1)
    d.setdefault("gene_id", d.get("ID","NA"))
    d.setdefault("gene_name", d.get("Name", d.get("gene_id","NA")))
    return d

def is_gz(p: str) -> bool:
    try:
        import mimetypes, gzip
        with gzip.open(p, "rb") as g:
            g.read(1)
        return True
    except Exception:
        return str(p).endswith(".gz")

def open_text(path: str):
    import gzip, io
    return io.TextIOWrapper(gzip.open(path, "rb")) if is_gz(path) else open(path, "r", encoding="utf-8", errors="replace")

# ---- total mapped reads (for CPM) -----------------------------------------
def mapped_reads(bam: str) -> int:
    try:
        p = subprocess.run(["samtools","idxstats", bam], text=True, check=True, stdout=subprocess.PIPE)
        tot = 0
        for ln in p.stdout.splitlines():
            parts = ln.strip().split("\t")
            if len(parts) >= 3 and parts[0] != "*":
                try: tot += int(parts[2])
                except: pass
        return tot
    except Exception:
        # Fallback: samtools view -c -F4
        p = subprocess.run(["samtools","view","-c","-F","4", bam], text=True, check=True, stdout=subprocess.PIPE)
        return int(p.stdout.strip())

# ---- build windows from GTF -----------------------------------------------
ft_accept = {t.strip() for t in args.feature_types.split(",") if t.strip()}
tssw = int(args.tss_win); body_min = int(args.body_offset_min); body_fr = float(args.body_offset_frac)

tmpd = Path(tempfile.mkdtemp(prefix=".pol2_", dir=".")).resolve()
genome_file = tmpd/"genome.faidx.tsv"  # contig order/sizes from BAM
tss_bed = tmpd/"tss.bed"
body_bed= tmpd/"body.bed"

# Dictionary to store gene information and aggregate across transcripts
gene_data = {}  # gid -> {gname, chrom, strand, transcripts: [(start, end, tss, body_lo, body_hi, body_len)]}

with open_text(args.gtf) as fh:
    for ln in fh:
        if not ln.strip() or ln.startswith("#"): continue
        p = ln.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] not in ft_accept: continue
        chrom = p[0]; strand = p[6]
        try:
            start = int(p[3]); end = int(p[4])
        except Exception:
            continue
        if end <= start: continue
        attrs = parse_attrs(p[8])
        gid, gname = attrs["gene_id"], attrs["gene_name"]

        tss = start if strand == "+" else end
        glen = max(1, end - start)
        offset = max(body_min, int(glen * body_fr))

        if strand == "+":
            body_lo, body_hi = tss + offset, end
        else:
            body_lo, body_hi = start, tss - offset

        body_lo = max(0, body_lo)
        body_len = max(0, body_hi - body_lo)
        tss_lo, tss_hi = max(0, tss - tssw), tss + tssw

        # Store transcript information
        if gid not in gene_data:
            gene_data[gid] = {
                'gname': gname,
                'chrom': chrom,
                'strand': strand,
                'transcripts': []
            }
        gene_data[gid]['transcripts'].append((tss_lo, tss_hi, body_lo, body_hi, body_len))

# Convert to genes list, aggregating transcript coordinates
genes = []  # (gid,gname,chrom,strand, tss_lo, tss_hi, body_lo, body_hi, body_len)
for gid, data in gene_data.items():
    transcripts = data['transcripts']
    
    # For TSS: take the union of all TSS windows
    tss_lo = min(t[0] for t in transcripts)
    tss_hi = max(t[1] for t in transcripts)
    
    # For body: take the union of all body regions
    body_lo = min(t[2] for t in transcripts)
    body_hi = max(t[3] for t in transcripts)
    body_len = max(0, body_hi - body_lo)
    
    genes.append((gid, data['gname'], data['chrom'], data['strand'], 
                 tss_lo, tss_hi, body_lo, body_hi, body_len))

if not genes:
    if args.fail_if_empty.lower() in ("true","1","yes"):
        raise SystemExit("ERROR: No acceptable loci parsed from GTF.")
    # still emit empty outputs + qc
    Path(args.out_pausing).write_text("gene_id\tchrom\tstrand\ttss_count\tgene_body_count\tpausing_index\tis_truncated\n")
    Path(args.out_genes).write_text("gene_id\tgene_name\tchrom\tstrand\ttss_bp\ttss_lo\ttss_hi\ttss_width\tbody_lo\tbody_hi\tbody_len\ttss_count\ttss_cpm\ttss_density_per_bp\tbody_count\tbody_cpm\tbody_density_per_bp\tpi_raw\tpi_len_norm\tis_truncated\n")
    Path(args.out_qc).write_text(json.dumps(dict(total_mapped=0, genes_seen=0), indent=2))
    raise SystemExit(0)

# write BED6 windows
with open(tss_bed, "w") as o1, open(body_bed, "w") as o2:
    for (gid,gname,chrom,strand, tss_lo,tss_hi, body_lo,body_hi,body_len) in genes:
        o1.write(f"{chrom}\t{tss_lo}\t{tss_hi}\t{gid}\t0\t{strand}\n")
        if body_len > 0:
            o2.write(f"{chrom}\t{body_lo}\t{body_hi}\t{gid}\t0\t{strand}\n")

def write_genome_file_from_bam(bam: str, out_path: Path):
    # Use idxstats for names and lengths in BAM header order
    p = subprocess.run(["samtools","idxstats", bam], text=True, check=True, stdout=subprocess.PIPE)
    lines = []
    for ln in p.stdout.splitlines():
        parts = ln.split("\t")
        if len(parts) >= 2 and parts[0] != "*":
            lines.append(f"{parts[0]}\t{parts[1]}")
    out_path.write_text("\n".join(lines)+"\n")

def sort_bed_with_genome_order(p: Path, genome: Path):
    if not (p.exists() and p.stat().st_size > 0):
        return
    try:
        # bedtools sort with genome file; write stdout back to file
        q = subprocess.run(["bedtools","sort","-g", str(genome), "-i", str(p)],
                           text=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        p.write_text(q.stdout)
    except subprocess.CalledProcessError:
        # Fallback to POSIX sort if bedtools sort is unavailable in this container
        run(["bash","-lc", f"LC_ALL=C sort -k1,1 -k2,2n -k3,3n '{p}' -o '{p}'"])

# Build genome file in BAM order then sort BEDs accordingly
write_genome_file_from_bam(args.bam, genome_file)
sort_bed_with_genome_order(tss_bed, genome_file)
sort_bed_with_genome_order(body_bed, genome_file)

# coverage counts
def coverage_counts(bed: Path, bam: str, genome: Path) -> dict[str,int]:
    if not bed.exists() or bed.stat().st_size == 0:
        return {}
    # Use -sorted with a genome file for consistent contig order
    p = subprocess.run(["bedtools","coverage","-counts","-sorted","-g",str(genome),"-a",str(bed),"-b",bam],
                       text=True, check=True, stdout=subprocess.PIPE)
    out={}
    for ln in p.stdout.splitlines():
        if not ln.strip() or ln.startswith(("track","browser","#")): continue
        f = ln.split("\t")
        if len(f) < 7: continue
        gid = f[3]; cnt = int(f[-1])
        out[gid] = out.get(gid, 0) + cnt
    return out

tss_counts  = coverage_counts(tss_bed,  args.bam, genome_file)
body_counts = coverage_counts(body_bed, args.bam, genome_file)

# totals for CPM
mapped = mapped_reads(args.bam)
denom  = (mapped/1_000_000.0) if mapped>0 else 1e-9

# write outputs
with open(args.out_pausing, "w") as p_out, open(args.out_genes, "w") as g_out:
    pw = csv.writer(p_out, delimiter="\t", lineterminator="\n")
    gw = csv.writer(g_out, delimiter="\t", lineterminator="\n")
    pw.writerow(["gene_id","chrom","strand","tss_count","gene_body_count","pausing_index","is_truncated"])
    gw.writerow([
        "gene_id","gene_name","chrom","strand",
        "tss_lo","tss_hi","tss_width",
        "body_lo","body_hi","body_len",
        "tss_count","tss_cpm","tss_density_per_bp",
        "body_count","body_cpm","body_density_per_bp",
        "pi_raw","pi_len_norm","is_truncated"
    ])
    for (gid,gname,chrom,strand, tss_lo,tss_hi, body_lo,body_hi,body_len) in genes:
        tcnt = int(tss_counts.get(gid, 0))
        bcnt = int(body_counts.get(gid, 0))
        t_w  = max(1, tss_hi - tss_lo)
        pi_raw = (tcnt / bcnt) if bcnt > 0 else float("nan")
        pi_len = ((tcnt / t_w) / (bcnt / body_len)) if (bcnt>0 and body_len>0) else float("nan")
        trunc = int(bcnt==0 or body_len==0)
        pw.writerow([gid, chrom, strand, tcnt, bcnt, pi_raw, trunc])
        gw.writerow([
            gid, gname, chrom, strand,
            tss_lo, tss_hi, t_w,
            body_lo, body_hi, body_len,
            tcnt, (tcnt/denom), (tcnt/max(1,t_w)),
            bcnt, (bcnt/denom), (bcnt/max(1,body_len) if body_len>0 else 0.0),
            pi_raw, pi_len, trunc
        ])

# QC JSON
qc = dict(
    total_mapped=int(mapped),
    genes_seen=len(genes),
    tss_window_bp=int(tssw),
    body_offset_min_bp=int(body_min),
    body_offset_frac=float(body_fr),
)
Path(args.out_qc).write_text(json.dumps(qc, indent=2))
print(f"[pol2_calc] done ts={datetime.datetime.utcnow().isoformat()}Z", file=sys.stderr)

# cleanup tmp
try:
    import shutil; shutil.rmtree(tmpd, ignore_errors=True)
except Exception:
    pass
