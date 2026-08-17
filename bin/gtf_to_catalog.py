#!/usr/bin/env python3
# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  gtf_to_catalog.py — build gene catalog from a GTF/GFF                   ║
# ║                                                                          ║
# ║  Inputs  : <in.gtf[.gz]>                                                 ║
# ║  Outputs : genes.tsv  (gene_id, gene_name, chr, strand, start, end, …)   ║
# ║            tss.bed    (BED6; 1bp TSS per gene; name = gene_name|gene_id) ║
# ║            tes.bed    (BED6; 1bp TES per gene; name = gene_name|gene_id) ║
# ║                                                                          ║
# ║  Options : --exclude-biotypes rRNA,tRNA  (exclude genes by biotype)       ║
# ║            --chr-add-prefix|--chr-remove-prefix (normalize chr names)    ║
# ║                                                                          ║
# ║  Design                                                                  ║
# ║   • Streaming parser (low memory); tolerant to GTF or GFF attributes.    ║
# ║   • Source-agnostic: Ensembl ("1"/"MT"), GENCODE ("chr1"), RefSeq        ║
# ║     (renamed "chr1") — incl. files carrying stray CR / CRLF line         ║
# ║     endings from NCBI assembly-report based renaming.                    ║
# ║   • Falls back across common keys: gene_id, gene_name, biotype           ║
# ║   • Consolidates per-gene extents using gene features when present;     ║
# ║     otherwise uses min/max across transcripts/exons.                      ║
# ║   • Coordinates are pooled ONLY within a single contig. A gene_id that   ║
# ║     appears on more than one contig (e.g. a primary chromosome plus an   ║
# ║     alt-haplotype/patch/random duplicate) is resolved to ONE locus via   ║
# ║     _select_primary_locus() — never merged into one nonsensical window.  ║
# ║   • Deterministic iteration/sorting at write stage.                       ║
# ║   • BED rows are BED6, 0-based start, 1-based end, with strand.           ║
# ║   • Aborts (non-zero exit) if any output gene has an empty/whitespace     ║
# ║     chromosome — this can never silently corrupt downstream steps.        ║
# ╚══════════════════════════════════════════════════════════════════════════╝

from __future__ import annotations
import sys, os, io, gzip, datetime, argparse
from collections import defaultdict
from typing import Dict, Tuple, Iterable, List, Set

# ── IO helpers ──────────────────────────────────────────────────────────────
def open_text(path: str) -> io.TextIOBase:
    """Open plain or gzipped text.

    newline="\\n" disables Python's universal-newline translation, so a stray
    carriage return *inside* a line (e.g. a mis-renamed seqid like "chr1\\r")
    is NOT turned into a line break that would split the record and blank the
    chromosome column. Any remaining CR is stripped per-line in iter_records().
    """
    if path.endswith(".gz"):
        return io.TextIOWrapper(
            gzip.open(path, "rb"), encoding="utf-8", errors="replace", newline="\n"
        )
    return open(path, "r", encoding="utf-8", errors="replace", newline="\n")


def iter_records(fh: io.TextIOBase) -> Iterable[List[str]]:
    """Yield cleaned 9-column GTF/GFF records.

    Robust to:
      • CRLF (Windows) and lone-CR endings, anywhere in the line
      • leading/trailing whitespace around the chromosome name
      • comment ('#') and blank lines
      • rows with fewer than 9 columns (skipped)
    """
    for ln in fh:
        # Strip *all* carriage returns first — GTF/GFF never contain a literal
        # CR as data, so this safely repairs both "...\r\n" and mid-line "chr1\r\t".
        if "\r" in ln:
            ln = ln.replace("\r", "")
        if not ln or ln[0] == "#":
            continue
        ln = ln.rstrip("\n")
        if not ln:
            continue
        f = ln.split("\t")
        if len(f) < 9:
            continue
        # Defensively trim the structural columns most prone to stray whitespace.
        f[0] = f[0].strip()  # chrom
        f[2] = f[2].strip()  # feature type
        f[6] = f[6].strip()  # strand
        yield f


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


# GFF3 commonly prefixes IDs ("gene:ENSG..", "transcript:ENST.."). Strip those
# so a gene's exon/transcript rows resolve to the same gene_id where possible.
_GFF_ID_PREFIXES = ("gene:", "transcript:", "rna:", "mrna:")


def clean_gid(gid: str | None) -> str | None:
    if not gid:
        return gid
    low = gid.lower()
    for p in _GFF_ID_PREFIXES:
        if low.startswith(p):
            return gid[len(p):]
    return gid


# Attribute key search orders, shared by both passes for consistency.
GID_KEYS   = ["gene_id", "geneID", "gene", "ID", "Parent"]
GNAME_KEYS = ["gene_name", "Name", "gene"]
GTYPE_KEYS = ["gene_type", "gene_biotype", "biotype"]

# Feature types that, when no explicit "gene" row exists, bound a gene's extent.
# Kept broad on purpose so ncRNA/primary_transcript/lnc_RNA/etc. all count.
_TX_LIKE = {"transcript", "mrna", "exon", "primary_transcript",
            "lnc_rna", "lncrna", "ncrna", "trna", "rrna", "cds"}


def make_name(gene_name: str | None, gene_id: str) -> str:
    if gene_name and gene_name != gene_id:
        return f"{gene_name}|{gene_id}"
    return gene_id


def normalize_chrom(chrom: str, mode: str | None) -> str:
    """Normalize chromosome name for UCSC/Ensembl compatibility.

    mode: 'add'    = Ensembl -> UCSC  (1 -> chr1,  MT -> chrM)
          'remove' = UCSC -> Ensembl  (chr1 -> 1,  chrM/chrMT -> MT)
    """
    if not chrom:
        return chrom
    if mode == "add":
        if chrom in ("MT", "mt"):
            return "chrM"
        return chrom if chrom.startswith("chr") else f"chr{chrom}"
    if mode == "remove":
        if chrom in ("chrM", "chrMT"):
            return "MT"
        if chrom.startswith("chr"):
            return chrom[3:]
    return chrom


# ── Core conversion ─────────────────────────────────────────────────────────
class Hints:
    """Best-known per-gene fields gathered from any row (gene or otherwise),
    used only as a last-resort fallback for genes with no explicit 'gene'
    feature row at all. Deliberately does NOT track chromosome — chrom is
    always resolved per-contig (see build_catalog/_select_primary_locus) so
    it can never be used to silently merge two different locations."""
    __slots__ = ("strand", "name", "bio")

    def __init__(self) -> None:
        self.strand: Dict[str, str] = {}
        self.name: Dict[str, str] = {}
        self.bio: Dict[str, str] = {}


def _is_alt_contig(chrom: str) -> bool:
    """True for non-primary-assembly sequences: alt haplotypes, patches, and
    unplaced/unlocalized scaffolds.

    Every genome build this pipeline uses (UCSC-style hs1/mm10/dm6/canFam6/
    TAIR10, and GENCODE hg19/hg38) names primary chromosomes as a bare token
    with no underscore (chr1, chr17, chrX, chrM, 2L, Chr1, ...) and marks
    every alt/patch/random/unplaced sequence with an underscore-joined
    suffix (chr17_GL456022_random, chr1_KI270706v1_alt, chrUn_GL456239,
    ...). That's the same convention UCSC/GENCODE use to separate primary
    from non-primary sequences, so this generalizes across organisms
    without hardcoding any specific contig name.

    Not exhaustive for every possible annotation source (a bare Ensembl
    scaffold accession like "KI270706.1" has no underscore and would slip
    through undetected) — but that only affects the tie-break order in
    _select_primary_locus, which always makes *some* deterministic choice
    and logs it; it can never again silently pool coordinates the way the
    bug this replaces did.
    """
    return "_" in chrom


def _select_primary_locus(gid: str, loci: Dict[str, dict], kind: str) -> Tuple[str, dict]:
    """
    Pick ONE chromosome for a gene_id that has rows on more than one contig.

    Background: a gene_id occasionally appears on both a primary chromosome
    and an alt-haplotype/patch/random duplicate of the same region. The bug
    this replaces: the previous version of this script pooled start/end
    coordinates with a plain min()/max() across every row sharing a gene_id,
    regardless of which contig each row was on. When a gene has real loci on
    two different contigs, that produces one nonsensical window spanning
    both locations (real case found: mm10 Hspa1a/Hspa1b landed on a THIRD,
    unrelated alt contig with a ~30 Mb span, purely from pooled coordinates
    — see analysis/scripts/patch_mm10_hspa1.py for the diagnosis this fix
    replaces).

    Rule: prefer a primary (non-alt, see _is_alt_contig) contig if any
    candidate is one. If more than one candidate remains (all primary, or
    all alt/patch — no primary copy at all), fall back to the locus with the
    larger transcript-supported span, then chrom name, for a deterministic
    pick. Always logs which contig was chosen and which were dropped — this
    can never again fail silently the way the pooling bug did.
    """
    if len(loci) == 1:
        chrom, loc = next(iter(loci.items()))
        return chrom, loc

    primary = {c: v for c, v in loci.items() if not _is_alt_contig(c)}
    pool = primary if primary else loci
    chosen = max(pool, key=lambda c: (pool[c]["end"] - pool[c]["start"], c))
    dropped = sorted(c for c in loci if c != chosen)
    sys.stderr.write(
        f"[gtf_to_catalog] WARNING: {kind} {gid!r} has rows on multiple contigs "
        f"{sorted(loci.keys())} -- using {chosen!r}, ignoring {dropped} "
        f"(likely an alt-haplotype/patch/scaffold duplicate; coordinates are "
        f"NOT pooled across contigs)\n"
    )
    return chosen, loci[chosen]


def build_catalog(gtf_path: str) -> Tuple[
    Dict[str, Dict[str, Dict[str, object]]],
    Dict[str, Dict[str, int]], Dict[str, Dict[str, int]],
    Hints,
]:
    """
    Single streaming pass. Returns:
      genes  : gene_id -> chrom -> {strand, start, end, gene_name, biotype}
               (kept PER CHROMOSOME — a gene_id's rows are only ever pooled
               with other rows on the SAME contig; see _select_primary_locus
               for how one locus is chosen when a gene_id spans more than one)
      tx_min : gene_id -> chrom -> min coord across transcript/exon-like rows
      tx_max : gene_id -> chrom -> max coord across transcript/exon-like rows
      hints  : Hints() with strand/name/biotype fallbacks per gene_id, used
               only when a gene has no explicit 'gene' feature row at all
    """
    genes: Dict[str, Dict[str, Dict[str, object]]] = defaultdict(dict)
    tx_min: Dict[str, Dict[str, int]] = defaultdict(dict)
    tx_max: Dict[str, Dict[str, int]] = defaultdict(dict)
    hints = Hints()

    with open_text(gtf_path) as fh:
        for f in iter_records(fh):
            chrom, _src, feat, start, end, _score, strand, _frame, attrs_raw = f[:9]
            try:
                s = int(start); e = int(end)
            except Exception:
                continue
            if e <= s:
                continue

            attrs = parse_attrs(attrs_raw)
            gid = clean_gid(pick(GID_KEYS, attrs))
            if not gid:
                continue
            gname = pick(GNAME_KEYS, attrs, default=gid)
            gtype = pick(GTYPE_KEYS, attrs, default="")

            # Track hints from every row. Chromosome is intentionally NOT
            # tracked here (see Hints docstring) — it is always resolved
            # per-contig below, so two different loci for the same gene_id
            # can never be silently blended into one.
            if strand in {"+", "-"}:
                hints.strand[gid] = strand
            if gname and gname != gid:
                hints.name.setdefault(gid, gname)
            if gtype:
                hints.bio.setdefault(gid, gtype)

            feat_l = feat.lower()
            if feat_l == "gene":
                loc = genes[gid].get(chrom)
                if loc is None:
                    genes[gid][chrom] = {
                        "strand": strand if strand in {"+", "-"} else hints.strand.get(gid, "+"),
                        "start": s,
                        "end": e,
                        "gene_name": gname or gid,
                        "biotype": gtype or "",
                    }
                else:
                    # Multiple 'gene' rows on the SAME contig for this
                    # gene_id (rare but valid, e.g. split annotation
                    # records) — safe to pool since they share a location.
                    loc["start"] = min(loc["start"], s)
                    loc["end"] = max(loc["end"], e)
                    if loc["gene_name"] == gid and gname:
                        loc["gene_name"] = gname
                    if not loc["biotype"] and gtype:
                        loc["biotype"] = gtype
            elif feat_l in _TX_LIKE:
                cur_min = tx_min[gid].get(chrom, 10**18)
                cur_max = tx_max[gid].get(chrom, -1)
                tx_min[gid][chrom] = min(cur_min, s)
                tx_max[gid][chrom] = max(cur_max, e)

    return genes, tx_min, tx_max, hints


def finalize_rows(genes: Dict[str, Dict[str, Dict[str, object]]],
                  tx_min: Dict[str, Dict[str, int]],
                  tx_max: Dict[str, Dict[str, int]],
                  hints: Hints) -> list[tuple]:
    """
    Build final per-gene rows:
      (gene_id, gene_name, chr, strand, start, end, tss, tes, biotype)

    Logic:
      1. If explicit 'gene' features exist anywhere in the GTF, use those
         (per gene_id, resolved to a single contig via _select_primary_locus
         when a gene_id has rows on more than one).
      2. Otherwise infer genes from transcript/exon-like rows, resolved to a
         single contig the same way.
    """
    out = []
    has_explicit_genes = len(genes) > 0

    if has_explicit_genes:
        for gid in sorted(genes.keys()):
            chrom, g = _select_primary_locus(gid, genes[gid], "gene_id")
            strand = str(g.get("strand") or hints.strand.get(gid) or "+")
            gstart = int(g.get("start") or tx_min.get(gid, {}).get(chrom, 1))
            gend   = int(g.get("end")   or tx_max.get(gid, {}).get(chrom, gstart))
            gname  = str(g.get("gene_name") or hints.name.get(gid) or gid)
            gtype  = str(g.get("biotype") or hints.bio.get(gid) or "")
            if gend <= gstart:
                continue
            tss = gstart if strand == "+" else gend
            tes = gend   if strand == "+" else gstart
            out.append((gid, gname, chrom, strand, gstart, gend, tss, tes, gtype))
    else:
        all_ids = set(tx_min.keys()) | set(tx_max.keys())
        for gid in sorted(all_ids):
            chroms = set(tx_min.get(gid, {})) | set(tx_max.get(gid, {}))
            loci = {}
            for c in chroms:
                t_min = tx_min.get(gid, {}).get(c)
                t_max = tx_max.get(gid, {}).get(c)
                if t_min is None or t_max is None or t_max <= t_min:
                    continue
                loci[c] = {"start": t_min, "end": t_max}
            if not loci:
                # No chromosome info, or no valid span, on any contig — skip
                continue
            chrom, loc = _select_primary_locus(gid, loci, "gene_id (inferred from transcripts)")
            strand = hints.strand.get(gid, "+")
            gstart, gend = loc["start"], loc["end"]
            gname = hints.name.get(gid, gid)
            gtype = hints.bio.get(gid, "")
            tss = gstart if strand == "+" else gend
            tes = gend   if strand == "+" else gstart
            out.append((gid, gname, chrom, strand, gstart, gend, tss, tes, gtype))

    return out


# ── Writers ────────────────────────────────────────────────────────────────
def write_genes_tsv(path: str, rows: list[tuple], chr_mode: str | None = None) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("gene_id\tgene_name\tchr\tstrand\tstart\tend\ttss\ttes\tbiotype\n")
        for gid, gname, chrom, strand, gstart, gend, tss, tes, gtype in rows:
            c = normalize_chrom(chrom, chr_mode)
            fh.write(f"{gid}\t{gname}\t{c}\t{strand}\t{gstart}\t{gend}\t{tss}\t{tes}\t{gtype}\n")


def write_bed6(path: str, rows: list[tuple], which: str, chr_mode: str | None = None) -> None:
    """
    which: "tss" or "tes"
    Writes 1-bp BED6 with 0-based start, 1-based end.
    name field = gene_name|gene_id (or gene_id if identical).
    """
    assert which in {"tss", "tes"}
    with open(path, "w", encoding="utf-8") as fh:
        for gid, gname, chrom, strand, _gs, _ge, tss, tes, _bt in rows:
            c = normalize_chrom(chrom, chr_mode)
            pos = int(tss if which == "tss" else tes)
            lo = max(0, pos - 1)
            hi = pos
            name = make_name(gname, gid)
            fh.write(f"{c}\t{lo}\t{hi}\t{name}\t0\t{strand}\n")


# ── Validation ──────────────────────────────────────────────────────────────
def validate_rows(rows: list[tuple]) -> list[str]:
    """Return a list of fatal problems (empty list == OK)."""
    problems: list[str] = []
    if not rows:
        problems.append("no genes were parsed from the GTF")
        return problems
    empty_chrom = sum(1 for r in rows if not r[2])
    ws_chrom    = sum(1 for r in rows if r[2] and r[2] != r[2].strip())
    if empty_chrom:
        problems.append(
            f"{empty_chrom}/{len(rows)} genes have an EMPTY chromosome — the GTF is "
            f"likely corrupt (e.g. a stray CR from seqid renaming split the records)"
        )
    if ws_chrom:
        problems.append(f"{ws_chrom}/{len(rows)} genes have whitespace in the chromosome name")
    return problems


# ── CLI ────────────────────────────────────────────────────────────────────
def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Build gene catalog from GTF/GFF")
    ap.add_argument("gtf_in", help="Input GTF/GFF (plain or .gz)")
    ap.add_argument("genes_tsv", help="Output genes.tsv")
    ap.add_argument("tss_bed", help="Output TSS BED")
    ap.add_argument("tes_bed", help="Output TES BED")
    ap.add_argument("--exclude-biotypes", default="", help="Comma-separated biotypes to exclude (e.g. rRNA,tRNA)")
    ap.add_argument("--chr-add-prefix", action="store_true", help="Add chr prefix to chromosome names (Ensembl->UCSC)")
    ap.add_argument("--chr-remove-prefix", action="store_true", help="Remove chr prefix (UCSC->Ensembl)")
    args = ap.parse_args(argv[1:])

    gtf_in = args.gtf_in
    genes_tsv = args.genes_tsv
    tss_bed = args.tss_bed
    tes_bed = args.tes_bed

    exclude_biotypes: Set[str] = set()
    if args.exclude_biotypes:
        exclude_biotypes = {b.strip().lower() for b in args.exclude_biotypes.split(",") if b.strip()}
        sys.stderr.write(f"[gtf_to_catalog] Excluding biotypes: {exclude_biotypes}\n")

    chr_mode: str | None = None
    if args.chr_add_prefix:
        chr_mode = "add"
        sys.stderr.write("[gtf_to_catalog] Adding chr prefix to chromosome names\n")
    elif args.chr_remove_prefix:
        chr_mode = "remove"
        sys.stderr.write("[gtf_to_catalog] Removing chr prefix from chromosome names\n")

    sys.stderr.write(f"[gtf_to_catalog] start ts={datetime.datetime.utcnow().isoformat()}Z\n")
    if not os.path.exists(gtf_in):
        sys.stderr.write(f"ERROR: input not found: {gtf_in}\n")
        return 2

    # Single streaming pass collects genes + hints together (no fragile re-read).
    genes, tx_min, tx_max, hints = build_catalog(gtf_in)
    rows = finalize_rows(genes, tx_min, tx_max, hints)

    # Filter by biotype if requested
    if exclude_biotypes:
        rows = [r for r in rows if (r[8] or "").strip().lower() not in exclude_biotypes]
        sys.stderr.write(f"[gtf_to_catalog] After biotype filter: {len(rows)} genes\n")

    # Fail loud rather than emit a silently-corrupt catalog.
    problems = validate_rows(rows)
    if problems:
        for p in problems:
            sys.stderr.write(f"ERROR: {p}\n")
        sys.stderr.write("[gtf_to_catalog] aborting without writing outputs\n")
        return 3

    # Write (deterministic order already enforced by finalize_rows sort)
    write_genes_tsv(genes_tsv, rows, chr_mode)
    write_bed6(tss_bed, rows, "tss", chr_mode)
    write_bed6(tes_bed, rows, "tes", chr_mode)

    n_chroms = len({normalize_chrom(r[2], chr_mode) for r in rows})
    sys.stderr.write(
        f"✓ Wrote {len(rows)} genes across {n_chroms} sequences to "
        f"{genes_tsv}, TSS={tss_bed}, TES={tes_bed}\n"
    )
    sys.stderr.write(f"[gtf_to_catalog] done ts={datetime.datetime.utcnow().isoformat()}Z\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
