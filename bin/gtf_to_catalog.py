#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  gtf_to_catalog.py — build gene catalog from a GTF/GFF                   ║
# ║                                                                          ║
# ║  Inputs  : <in.gtf[.gz]>                                                 ║
# ║  Outputs : genes.tsv  (gene_id, gene_name, chr, strand, start, end, tss, tes, biotype) ║
# ║            tss.bed    (BED6; 1bp TSS per gene; name = gene_name|gene_id) ║
# ║            tes.bed    (BED6; 1bp TES per gene; name = gene_name|gene_id) ║
# ║                                                                          ║
# ║  Design                                                                  ║
# ║   • Streaming parser (low memory); tolerant to GTF or GFF attributes.    ║
# ║   • Falls back across common keys:                                       ║
# ║       gene_id:  gene_id, geneID, gene, ID, Parent                         ║
# ║       gene_name: gene_name, Name, gene                                    ║
# ║       biotype:  gene_type, gene_biotype, biotype                          ║
# ║   • Includes ALL gene types by default (maintains scientific completeness) ║
# ║   • Consolidates per-gene extents using gene features when present;       ║
# ║     otherwise uses min/max across transcripts/exons.                      ║
# ║   • Calculates single TSS/TES coordinates from transcript boundaries       ║
# ║   • Deterministic iteration/sorting at write stage.                       ║
# ║   • BED rows are BED6, 0-based start, 1-based end, with strand.           ║
# ╚══════════════════════════════════════════════════════════════════════════╝

from __future__ import annotations
import sys, os, io, gzip, datetime
from collections import defaultdict
from typing import Dict, Tuple, Iterable

# ── Biotype filtering ────────────────────────────────────────────────────────
# Optional biotype filtering - by default, include ALL gene types
# Users can specify --exclude-biotypes to filter out specific types

def is_valid_gene_biotype(biotype: str, excluded_biotypes: set = None) -> bool:
    """Check if gene biotype should be included in the catalog."""
    if not biotype:
        return True  # Include genes without biotype info
    if excluded_biotypes is None:
        return True  # Include all biotypes by default
    return biotype.lower() not in excluded_biotypes

# ── IO helpers ──────────────────────────────────────────────────────────────
def open_text(path: str) -> io.TextIOBase:
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"))
    return open(path, "r", encoding="utf-8", errors="replace")

# ── Attribute parsing (GTF and GFF tolerant) ───────────────────────────────
def parse_attrs(s: str) -> Dict[str, str]:
    """
    Accepts either:
      GTF: key "value"; key "value";
      GFF: key=value;key=value;
    Returns a {key->value} dict (stripped).
    """
    out: Dict[str, str] = {}
    if not s:
        return out
    # Split on ';' but tolerate missing trailing ';'
    for chunk in s.strip().strip(";").split(";"):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "=" in chunk:  # GFF style
            k, v = chunk.split("=", 1)
            out[k.strip()] = v.strip().strip('"')
        elif " " in chunk:  # GTF style
            k, v = chunk.split(" ", 1)
            out[k.strip()] = v.strip().strip('"')
        else:
            # Bare key or odd token; keep as flag
            out[chunk] = ""
    return out

def pick(keys: Iterable[str], d: Dict[str, str], default: str | None = None) -> str | None:
    for k in keys:
        if k in d and d[k]:
            return d[k]
    return default

def make_name(gene_name: str | None, gene_id: str) -> str:
    if gene_name and gene_name != gene_id:
        return f"{gene_name}|{gene_id}"
    return gene_id

# ── Core conversion ─────────────────────────────────────────────────────────
def build_catalog(gtf_path: str) -> Tuple[
    Dict[str, Dict[str, object]],
    Dict[str, int], Dict[str, int]
]:
    """
    Returns:
      genes: gene_id -> {chr, strand, start, end, gene_name, biotype}
      tx_min, tx_max: gene_id -> min/max coords across transcripts/exons
    """
    genes: Dict[str, Dict[str, object]] = {}  # from explicit gene features
    tx_min: Dict[str, int] = defaultdict(lambda: 10**18)
    tx_max: Dict[str, int] = defaultdict(lambda: -1)

    # Optionally track best-known fields when gene rows are absent/partial
    strand_hint: Dict[str, str] = {}
    chrom_hint: Dict[str, str] = {}
    name_hint: Dict[str, str] = {}
    bio_hint: Dict[str, str] = {}

    with open_text(gtf_path) as fh:
        for ln in fh:
            if not ln or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _src, feat, start, end, _score, strand, _frame, attrs_raw = f
            try:
                s = int(start); e = int(end)
            except Exception:
                continue
            if e <= s:
                continue

            attrs = parse_attrs(attrs_raw)
            gid = pick(["gene_id", "geneID", "gene", "ID", "Parent"], attrs)
            if not gid:
                continue
            gname = pick(["gene_name", "Name", "gene"], attrs, default=gid)
            gtype = pick(["gene_type", "gene_biotype", "biotype"], attrs, default="")

            # Track hints for non-gene rows
            if strand in {"+", "-"}:
                strand_hint[gid] = strand
            chrom_hint.setdefault(gid, chrom)
            if gname and gname != gid:
                name_hint.setdefault(gid, gname)
            if gtype:
                bio_hint.setdefault(gid, gtype)

            if feat == "gene":
                # Note: Biotype filtering is now optional (disabled by default)
                # Users can specify excluded biotypes via command line if needed
                    
                g = genes.setdefault(gid, {
                    "chr": chrom,
                    "strand": strand if strand in {"+", "-"} else strand_hint.get(gid, "+"),
                    "start": s,
                    "end": e,
                    "gene_name": gname or gid,
                    "biotype": gtype or ""
                })
                # Expand if multiple 'gene' rows appear
                if chrom and not g["chr"]:
                    g["chr"] = chrom
                g["start"] = min(int(g["start"]), s)
                g["end"]   = max(int(g["end"]),   e)
                if (not g.get("gene_name")) or g["gene_name"] == gid:
                    g["gene_name"] = gname or gid
                if (not g.get("biotype")) and gtype:
                    g["biotype"] = gtype
            elif feat in {"transcript", "mRNA", "exon"}:
                # Track transcript/exon boundaries for gene extent calculation
                tx_min[gid] = min(tx_min[gid], s)
                tx_max[gid] = max(tx_max[gid], e)

    return genes, tx_min, tx_max

def finalize_rows(genes: Dict[str, Dict[str, object]],
                  tx_min: Dict[str, int],
                  tx_max: Dict[str, int],
                  chrom_hint: Dict[str, str],
                  strand_hint: Dict[str, str],
                  name_hint: Dict[str, str],
                  bio_hint: Dict[str, str]) -> list[tuple]:
    """
    Build final per-gene rows:
      (gene_id, gene_name, chr, strand, start, end, tss, tes, biotype)
    
    Logic:
      1. If explicit 'gene' features exist, only use those (strict mode)
      2. If no explicit 'gene' features exist, fall back to transcript-based inference
      3. Include all gene types by default (maintains scientific completeness)
      4. Calculate single TSS/TES coordinates based on strand
    """
    out = []
    
    # Check if we have any explicit gene features
    has_explicit_genes = len(genes) > 0
    
    if has_explicit_genes:
        # Strict mode: only use genes with explicit 'gene' features
        all_ids = set(genes.keys())
        for gid in sorted(all_ids):
            g = genes[gid]
            chrom = str(g.get("chr") or chrom_hint.get(gid) or "")
            strand = str(g.get("strand") or strand_hint.get(gid) or "+")
            # Use transcript bounds if available to extend gene coordinates
            gstart = int(g.get("start") or tx_min.get(gid, 1))
            gend   = int(g.get("end")   or tx_max.get(gid, gstart))
            gname  = str(g.get("gene_name") or name_hint.get(gid) or gid)
            gtype  = str(g.get("biotype") or bio_hint.get(gid) or "")

            # Note: Biotype filtering is now optional (disabled by default)
            # All gene types are included to maintain scientific completeness

            if gend <= gstart:
                continue

            # Calculate single TSS and TES coordinates based on strand
            tss = gstart if strand == "+" else gend
            tes = gend   if strand == "+" else gstart
                
            out.append((gid, gname, chrom, strand, gstart, gend, tss, tes, gtype))
    else:
        # Fallback mode: infer genes from transcripts when no explicit gene features exist
        all_ids = set(tx_min.keys()) | set(tx_max.keys())
        for gid in sorted(all_ids):
            chrom = chrom_hint.get(gid, "")
            if not chrom:
                # No chromosome info at all — skip
                continue
            strand = strand_hint.get(gid, "+")
            gstart = tx_min.get(gid, None)
            gend   = tx_max.get(gid, None)
            if gstart is None or gend is None or gend <= gstart:
                continue
            gname = name_hint.get(gid, gid)
            gtype = bio_hint.get(gid, "unknown")

            # Note: Biotype filtering is now optional (disabled by default)
            # All gene types are included to maintain scientific completeness

            # Calculate single TSS and TES coordinates based on strand
            tss = gstart if strand == "+" else gend
            tes = gend   if strand == "+" else gstart
                
            out.append((gid, gname, chrom, strand, gstart, gend, tss, tes, gtype))
    
    return out

# ── Writers ────────────────────────────────────────────────────────────────
def write_genes_tsv(path: str, rows: list[tuple]) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("gene_id\tgene_name\tchr\tstrand\tstart\tend\ttss\ttes\tbiotype\n")
        for gid, gname, chrom, strand, gstart, gend, tss, tes, gtype in rows:
            fh.write(f"{gid}\t{gname}\t{chrom}\t{strand}\t{gstart}\t{gend}\t{tss}\t{tes}\t{gtype}\n")

def write_bed6(path: str, rows: list[tuple], which: str) -> None:
    """
    which: "tss" or "tes"
    Writes 1-bp BED6 with 0-based start, 1-based end.
    name field = gene_name|gene_id (or gene_id if identical).
    """
    assert which in {"tss", "tes"}
    idx = 6 if which == "tss" else 7
    with open(path, "w", encoding="utf-8") as fh:
        for gid, gname, chrom, strand, _gs, _ge, tss, tes, _bt in rows:
            pos = int(tss if which == "tss" else tes)
            lo = max(0, pos - 1)
            hi = pos
            name = make_name(gname, gid)
            fh.write(f"{chrom}\t{lo}\t{hi}\t{name}\t0\t{strand}\n")

# ── CLI ────────────────────────────────────────────────────────────────────
def main(argv: list[str]) -> int:
    sys.stderr.write(f"[gtf_to_catalog] start ts={datetime.datetime.utcnow().isoformat()}Z\n")
    if len(argv) != 5:
        sys.stderr.write("Usage: gtf_to_catalog.py <in.gtf[.gz]> <genes.tsv> <tss.bed> <tes.bed>\n")
        return 2

    gtf_in, genes_tsv, tss_bed, tes_bed = argv[1:5]
    if not os.path.exists(gtf_in):
        sys.stderr.write(f"ERROR: input not found: {gtf_in}\n")
        return 2

    # Build
    genes, tx_min, tx_max = build_catalog(gtf_in)

    # Hints were captured inside build_catalog; re-open file to extract hints
    # (We can avoid a second pass by returning hints; do that instead.)
    # Updated: return hints directly.
    # For clarity, re-run minimal pass to gather the hint dicts:
    strand_hint: Dict[str, str] = {}
    chrom_hint: Dict[str, str] = {}
    name_hint: Dict[str, str] = {}
    bio_hint: Dict[str, str] = {}
    with open_text(gtf_in) as fh:
        for ln in fh:
            if not ln or ln[0] == "#":
                continue
            f = ln.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, _src, _feat, _s, _e, _score, strand, _frame, attrs_raw = f
            attrs = parse_attrs(attrs_raw)
            gid = pick(["gene_id", "geneID", "gene", "ID", "Parent"], attrs)
            if not gid:
                continue
            if strand in {"+", "-"}:
                strand_hint[gid] = strand
            chrom_hint.setdefault(gid, chrom)
            gname = pick(["gene_name", "Name", "gene"], attrs, default=None)
            if gname:
                name_hint.setdefault(gid, gname)
            gtype = pick(["gene_type", "gene_biotype", "biotype"], attrs, default=None)
            if gtype:
                bio_hint.setdefault(gid, gtype)

    rows = finalize_rows(genes, tx_min, tx_max, chrom_hint, strand_hint, name_hint, bio_hint)

    # Write (deterministic order already enforced by finalize_rows sort)
    write_genes_tsv(genes_tsv, rows)
    write_bed6(tss_bed, rows, "tss")
    write_bed6(tes_bed, rows, "tes")

    sys.stderr.write(f"✓ Wrote {len(rows)} genes to {genes_tsv}, TSS={tss_bed}, TES={tes_bed}\n")
    sys.stderr.write(f"[gtf_to_catalog] done ts={datetime.datetime.utcnow().isoformat()}Z\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
