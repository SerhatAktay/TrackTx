#!/usr/bin/env python3
"""
detect_barcode_umi.py

QC-based detection of inline sample barcodes and UMIs (Unique Molecular
Identifiers) at the 5' and 3' ends of raw sequencing reads, for use by
TrackTx's preprocessing step (03_preprocess_and_quality_filter_reads.nf)
in "auto" and "verify" detect_mode.

WHY THIS EXISTS
----------------
A fixed inline barcode and a random UMI both sit at the very start or end
of a read, and both get trimmed off before alignment -- but the pipeline
needs to know, per dataset: is anything there at all, which read is it on
(R1 or R2), which end (5' or 3'), and how long is it? Published methods
sections are not always explicit about this (see the E. coli MG1655
cohort entry, where the UMI is on R2, not R1, and every prior cohort
entry only ever supported R1). Guessing wrong silently corrupts the
active-site coordinate PRO-seq exists to measure -- so this script gives
the pipeline a way to check a claim against the data, or to propose one
when the user does not know.

THE CORE IDEA (plain terms)
----------------------------
Look at the letter-mix (A/C/G/T proportions) at each position near both
ends of a subsample of raw reads, and compare it to the letter-mix from
a safely-interior "bulk" window of the same reads (real biological
sequence, unambiguously past any adapter/barcode/UMI):

  - If a terminal position's letter-mix is basically the SAME as bulk
    -> nothing artificial there; real sequence starts immediately.
  - If a terminal position's letter-mix is dominated by ONE letter far
    more than bulk is (low entropy, high divergence from bulk)
    -> looks like a FIXED tag (a constant multiplexing barcode).
  - If a terminal position's letter-mix is MORE uniform (higher entropy)
    than this dataset's own bulk letter-mix -> looks like a RANDOM tag (a
    synthesized UMI). This is compared to the dataset's own bulk entropy,
    not to a fixed 25/25/25/25 assumption, because real DNA is almost
    never perfectly uniform at every position but the exact degree of
    bias varies a lot by organism -- a synthesized random tag is reliably
    closer to uniform than the organism's own sequence is, even when
    that organism happens to be fairly GC-balanced.

The longest unbroken run of FIXED-looking (or RANDOM-looking) positions,
counted inward from position 1 of the 5' end or from the last base of
the 3' end, is reported as the candidate barcode (or UMI) length at that
corner. All four corners (R1 5', R1 3', R2 5', R2 3') are scored the
same way in one pass.

LIMITATIONS (read before trusting a result)
--------------------------------------------
  - This is a heuristic, not a certainty. A very short (<=3 nt) tag can
    be statistically invisible against noise; a biologically extreme
    5'/3' base bias (e.g. a strong consensus start site) can look like a
    weak fixed signal. Auto/verify results should be read as "the QC data
    is/isn't consistent with X", not as ground truth on their own.
  - Needs enough reads to get stable position-wise base frequencies;
    the default of 200,000 reads is comfortably enough for any real
    sequencing run, but very small test/toy FASTQs may look noisier.
  - Assumes barcode/UMI, if present, is contiguous from the true read
    edge (true for every published PRO-seq/GRO-seq protocol reviewed for
    this pipeline). A tag with a fixed offset from the edge is not
    something this script looks for.

Dependencies: Python 3 standard library only (no numpy/pandas needed).

Usage:
    # Full profile + best guesses, written as JSON
    python detect_barcode_umi.py --r1 R1.fastq.gz --r2 R2.fastq.gz \\
        --out profile.json

    # Same, but also print bash-sourceable variables (used by the
    # Nextflow module so it doesn't need a JSON parser in bash)
    python detect_barcode_umi.py --r1 R1.fastq.gz --r2 R2.fastq.gz \\
        --out profile.json --emit-shell

    # Check one specific (read, end, length) claim against the profile;
    # exits 0 on MATCH, 1 on MISMATCH (for the pipeline's "verify" mode)
    python detect_barcode_umi.py --r1 R1.fastq.gz --r2 R2.fastq.gz \\
        --check-kind umi --check-read R2 --check-location 5 \\
        --check-length 6 --tolerance 1

    # Run the built-in synthetic self-test (no FASTQ needed)
    python detect_barcode_umi.py --self-test

Author: Serhat Aktay (TrackTx pipeline)
Version: 1.0
"""

import argparse
import gzip
import json
import math
import sys
from collections import Counter

__version__ = "1.0"

BASES = ("A", "C", "G", "T")


# ─────────────────────────────────────────────────────────────────────────
# FASTQ reading
# ─────────────────────────────────────────────────────────────────────────

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def sample_reads(path, n_reads):
    """Yield up to n_reads sequence strings from a FASTQ (gz or plain)."""
    seqs = []
    with open_maybe_gz(path) as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                seqs.append(line.strip())
                if len(seqs) >= n_reads:
                    break
    return seqs


# ─────────────────────────────────────────────────────────────────────────
# Composition / entropy math
# ─────────────────────────────────────────────────────────────────────────

def composition(seqs, pos):
    """ACGT counts at 0-based position `pos` across all seqs long enough."""
    c = Counter()
    for s in seqs:
        if pos < len(s):
            b = s[pos]
            if b in BASES:
                c[b] += 1
    total = sum(c.values())
    if total == 0:
        return None, 0
    frac = {b: c.get(b, 0) / total for b in BASES}
    return frac, total


def entropy_bits(frac):
    """Shannon entropy of a base-fraction dict, in bits (0..2)."""
    h = 0.0
    for b in BASES:
        p = frac.get(b, 0.0)
        if p > 0:
            h -= p * math.log2(p)
    return h


def divergence(frac_a, frac_b):
    """Total variation distance between two base-fraction dicts (0..1)."""
    return 0.5 * sum(abs(frac_a.get(b, 0.0) - frac_b.get(b, 0.0)) for b in BASES)


def bulk_composition(seqs, skip, window):
    """Pooled ACGT composition over an interior window, used as the
    'this is what real biological sequence looks like in this dataset'
    reference. Positions are 0-based, counted from the 5' end."""
    c = Counter()
    for s in seqs:
        end = min(len(s), skip + window)
        if end <= skip:
            continue
        for b in s[skip:end]:
            if b in BASES:
                c[b] += 1
    total = sum(c.values())
    if total == 0:
        return {b: 0.25 for b in BASES}
    return {b: c.get(b, 0) / total for b in BASES}


# ─────────────────────────────────────────────────────────────────────────
# Per-end scan
# ─────────────────────────────────────────────────────────────────────────

def classify_position(frac, bulk, bulk_entropy, fixed_frac_min, div_min, random_entropy_excess_min):
    div = divergence(frac, bulk)
    max_frac = max(frac.values())
    ent = entropy_bits(frac)
    if div >= div_min and max_frac >= fixed_frac_min:
        cls = "FIXED"
    elif ent >= bulk_entropy + random_entropy_excess_min:
        # A synthesized random tag is closer to a perfectly uniform ACGT mix
        # than real biological sequence from THIS dataset ever is -- compare
        # against the dataset's own bulk entropy, not a fixed absolute bar,
        # so this still works for high- or low-GC organisms.
        cls = "RANDOM"
    else:
        cls = "BULK"
    return {
        "max_base_frac": round(max_frac, 4),
        "entropy_bits": round(ent, 4),
        "divergence_from_bulk": round(div, 4),
        "class": cls,
    }


def longest_run(classes, target):
    """Length of the longest unbroken prefix of `classes` equal to `target`."""
    n = 0
    for c in classes:
        if c == target:
            n += 1
        else:
            break
    return n


def scan_end(seqs, window, bulk, bulk_entropy, thresholds, from_three_prime=False):
    """Scan `window` positions inward from one end of the reads.

    Returns per-position stats (position 1 = the very edge base) plus the
    longest contiguous FIXED-looking and RANDOM-looking run starting at
    position 1.
    """
    positions = []
    read_len_min = min((len(s) for s in seqs), default=0)
    usable_window = min(window, max(read_len_min, 0))
    for i in range(usable_window):
        if from_three_prime:
            # position i (0-based from the end) maps to len(s)-1-i for each read;
            # composition() takes a fixed 0-based index, so build the slice manually.
            frac, n = _composition_from_end(seqs, i)
        else:
            frac, n = composition(seqs, i)
        if frac is None:
            break
        stats = classify_position(frac, bulk, bulk_entropy, **thresholds)
        stats["position"] = i + 1
        stats["n_reads"] = n
        positions.append(stats)

    classes = [p["class"] for p in positions]
    return {
        "positions": positions,
        "fixed_run_len": longest_run(classes, "FIXED"),
        "random_run_len": longest_run(classes, "RANDOM"),
    }


def _composition_from_end(seqs, i):
    c = Counter()
    for s in seqs:
        if i < len(s):
            b = s[len(s) - 1 - i]
            if b in BASES:
                c[b] += 1
    total = sum(c.values())
    if total == 0:
        return None, 0
    return {b: c.get(b, 0) / total for b in BASES}, total


# ─────────────────────────────────────────────────────────────────────────
# Full profile for one read (R1 or R2)
# ─────────────────────────────────────────────────────────────────────────

def profile_read(seqs, window, bulk_skip, bulk_window, thresholds):
    bulk = bulk_composition(seqs, bulk_skip, bulk_window)
    bulk_entropy = entropy_bits(bulk)
    five = scan_end(seqs, window, bulk, bulk_entropy, thresholds, from_three_prime=False)
    three = scan_end(seqs, window, bulk, bulk_entropy, thresholds, from_three_prime=True)
    return {
        "n_reads_sampled": len(seqs),
        "read_len_min": min((len(s) for s in seqs), default=0),
        "read_len_max": max((len(s) for s in seqs), default=0),
        "bulk_composition": {b: round(bulk[b], 4) for b in BASES},
        "bulk_entropy_bits": round(bulk_entropy, 4),
        "5prime": five,
        "3prime": three,
    }


def best_guess(profile, kind, exclude=None):
    """kind: 'fixed_run_len' (barcode candidate) or 'random_run_len' (UMI candidate).
    Picks the (read, end) with the longest run > 0 across whatever reads
    are present in `profile` (r1, r2). `exclude` is an optional set of
    (READ, LOCATION) pairs (e.g. {("R1", "5")}) to skip -- used so a second
    barcode slot's auto-detection doesn't just repeat the first slot's pick."""
    exclude = exclude or set()
    best = None
    for read_name in ("r1", "r2"):
        if read_name not in profile or profile[read_name] is None:
            continue
        for end_name, loc in (("5prime", "5"), ("3prime", "3")):
            if (read_name.upper(), loc) in exclude:
                continue
            run_len = profile[read_name][end_name][kind]
            if run_len > 0 and (best is None or run_len > best["length"]):
                best = {
                    "read": read_name.upper(),
                    "location": loc,
                    "length": run_len,
                }
    return best


# ─────────────────────────────────────────────────────────────────────────
# Top-level profile builder
# ─────────────────────────────────────────────────────────────────────────

def build_profile(r1_path, r2_path, n_reads, window, bulk_skip, bulk_window,
                   fixed_frac_min, div_min, random_entropy_excess_min):
    thresholds = dict(
        fixed_frac_min=fixed_frac_min,
        div_min=div_min,
        random_entropy_excess_min=random_entropy_excess_min,
    )
    seqs1 = sample_reads(r1_path, n_reads)
    profile = {"r1": profile_read(seqs1, window, bulk_skip, bulk_window, thresholds)}
    if r2_path:
        seqs2 = sample_reads(r2_path, n_reads)
        profile["r2"] = profile_read(seqs2, window, bulk_skip, bulk_window, thresholds)
    else:
        profile["r2"] = None

    profile["best_guess"] = {
        "barcode": best_guess(profile, "fixed_run_len"),
        "umi": best_guess(profile, "random_run_len"),
    }
    profile["params"] = {
        "n_reads_requested": n_reads,
        "window": window,
        "bulk_skip": bulk_skip,
        "bulk_window": bulk_window,
        "fixed_frac_min": fixed_frac_min,
        "divergence_min": div_min,
        "random_entropy_excess_min": random_entropy_excess_min,
        "version": __version__,
    }
    return profile


def run_len_at(profile, kind, read, location):
    read_key = read.lower()
    end_key = "5prime" if str(location) == "5" else "3prime"
    field = "fixed_run_len" if kind == "barcode" else "random_run_len"
    r = profile.get(read_key)
    if r is None:
        return 0
    return r[end_key][field]


# ─────────────────────────────────────────────────────────────────────────
# Self-test (synthetic reads, no FASTQ files needed)
# ─────────────────────────────────────────────────────────────────────────

def _self_test():
    import random
    random.seed(42)

    def make_read(prefix, bulk_len=60, gc=0.65):
        # Biological "bulk" body: GC-biased (not uniform), simulating a
        # real genome's base composition rather than a synthesized random tag.
        at = (1 - gc) / 2
        gcf = gc / 2
        body = "".join(random.choices(BASES, weights=(at, gcf, gcf, at), k=bulk_len))
        return prefix + body

    n = 5000
    # R1: 6nt FIXED barcode (all 'A's, i.e. constant) at 5', nothing at 3'
    r1 = [make_read("AAAAAA") for _ in range(n)]
    # R2: 6nt RANDOM UMI at 5' (independently random per read), nothing at 3'
    r2 = ["".join(random.choices(BASES, k=6)) + make_read("") for _ in range(n)]

    thresholds = dict(fixed_frac_min=0.75, div_min=0.35, random_entropy_excess_min=0.03)
    p1 = profile_read(r1, window=20, bulk_skip=15, bulk_window=40, thresholds=thresholds)
    p2 = profile_read(r2, window=20, bulk_skip=15, bulk_window=40, thresholds=thresholds)

    ok = True

    def check(label, cond):
        nonlocal ok
        status = "PASS" if cond else "FAIL"
        if not cond:
            ok = False
        print(f"[self-test] {label}: {status}", file=sys.stderr)

    check("R1 5' fixed run == 6 (synthetic constant barcode)",
          p1["5prime"]["fixed_run_len"] == 6)
    check("R1 3' fixed run == 0 (no tag there)",
          p1["3prime"]["fixed_run_len"] == 0)
    check("R1 3' random run == 0 (biological sequence, not random)",
          p1["3prime"]["random_run_len"] == 0)
    check("R2 5' random run == 6 (synthetic random UMI)",
          p2["5prime"]["random_run_len"] == 6)
    check("R2 5' fixed run == 0 (random != fixed)",
          p2["5prime"]["fixed_run_len"] == 0)
    check("R2 3' random run == 0 (no tag there)",
          p2["3prime"]["random_run_len"] == 0)

    # No tag anywhere: pure bulk sequence at both ends of both reads.
    r3 = [make_read("") for _ in range(n)]
    p3 = profile_read(r3, window=20, bulk_skip=15, bulk_window=40, thresholds=thresholds)
    check("No-tag control: R 5' fixed run == 0", p3["5prime"]["fixed_run_len"] == 0)
    check("No-tag control: R 5' random run == 0", p3["5prime"]["random_run_len"] == 0)
    check("No-tag control: R 3' fixed run == 0", p3["3prime"]["fixed_run_len"] == 0)
    check("No-tag control: R 3' random run == 0", p3["3prime"]["random_run_len"] == 0)

    # Low-GC organism (opposite bias) -- classification should be symmetric.
    r4 = ["".join(random.choices(BASES, k=6)) + make_read("", gc=0.30) for _ in range(n)]
    p4 = profile_read(r4, window=20, bulk_skip=15, bulk_window=40, thresholds=thresholds)
    check("Low-GC organism: 5' random run == 6 (still detected)",
          p4["5prime"]["random_run_len"] == 6)

    # Stacked tags: 6nt fixed barcode then 5nt random UMI, both on the 5' end
    # of the same read -- the FIXED run must stop at the boundary, not bleed
    # into the RANDOM region (and vice versa).
    r5 = ["AAAAAA" + "".join(random.choices(BASES, k=5)) + make_read("") for _ in range(n)]
    p5 = profile_read(r5, window=20, bulk_skip=20, bulk_window=40, thresholds=thresholds)
    check("Stacked tags: fixed run == 6 (barcode, not bleeding into UMI)",
          p5["5prime"]["fixed_run_len"] == 6)

    if ok:
        print("[self-test] ALL CHECKS PASSED", file=sys.stderr)
        sys.exit(0)
    else:
        print("[self-test] ONE OR MORE CHECKS FAILED", file=sys.stderr)
        sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────

def parse_args():
    ap = argparse.ArgumentParser(
        description="QC-based detection of inline barcodes/UMIs at read ends",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--r1", help="Path to R1 FASTQ (.fastq or .fastq.gz)")
    ap.add_argument("--r2", help="Path to R2 FASTQ (paired-end only)")
    ap.add_argument("--profile-in", help="Reuse a previously written --out JSON profile instead of "
                     "re-scanning the FASTQs (fast path for a second call in the same pipeline run; "
                     "--r1/--r2 are not needed when this is given)")
    ap.add_argument("--out", help="Write full JSON profile to this path")
    ap.add_argument("--exclude-read", choices=["R1", "R2", "r1", "r2"],
                     help="Exclude this (read, location) corner from best-guess picks "
                          "(pairs with --exclude-location; used for a second barcode slot "
                          "so it doesn't just repeat the first slot's pick)")
    ap.add_argument("--exclude-location", choices=["5", "3"],
                     help="Paired with --exclude-read")
    ap.add_argument("--emit-shell", action="store_true",
                     help="Also print bash-sourceable BEST_BARCODE_*/BEST_UMI_* vars to stdout")
    ap.add_argument("--n-reads", type=int, default=200000,
                     help="Reads to subsample per file (default: 200000)")
    ap.add_argument("--window", type=int, default=20,
                     help="Positions scanned from each end (default: 20)")
    ap.add_argument("--bulk-skip", type=int, default=25,
                     help="Bases skipped from the 5' end before the bulk reference window (default: 25)")
    ap.add_argument("--bulk-window", type=int, default=40,
                     help="Width of the interior bulk reference window (default: 40)")
    ap.add_argument("--fixed-frac-min", type=float, default=0.75,
                     help="Min dominant-base fraction to call a position FIXED (default: 0.75)")
    ap.add_argument("--divergence-min", type=float, default=0.35,
                     help="Min divergence-from-bulk to call a position FIXED or RANDOM (default: 0.35)")
    ap.add_argument("--random-entropy-excess-min", type=float, default=0.03,
                     help="Min entropy ABOVE this dataset's own bulk entropy to call a position RANDOM (default: 0.03 bits)")

    # "verify" mode: check one specific claim against the profile
    ap.add_argument("--check-kind", choices=["barcode", "umi"],
                     help="Check a specific barcode or UMI claim (enables verify mode)")
    ap.add_argument("--check-read", choices=["R1", "R2", "r1", "r2"],
                     help="Read the claim refers to")
    ap.add_argument("--check-location", choices=["5", "3"],
                     help="End the claim refers to")
    ap.add_argument("--check-length", type=int,
                     help="Claimed length in bp")
    ap.add_argument("--tolerance", type=int, default=1,
                     help="Allowed +/- bp mismatch before verify fails (default: 1)")

    ap.add_argument("--self-test", action="store_true",
                     help="Run the built-in synthetic self-test and exit (no FASTQ needed)")
    return ap.parse_args()


def main():
    args = parse_args()

    if args.self_test:
        _self_test()
        return

    if args.profile_in:
        with open(args.profile_in) as fh:
            profile = json.load(fh)
    else:
        if not args.r1:
            print("error: --r1 is required unless --profile-in is given (or use --self-test)",
                  file=sys.stderr)
            sys.exit(2)
        profile = build_profile(
            r1_path=args.r1, r2_path=args.r2, n_reads=args.n_reads, window=args.window,
            bulk_skip=args.bulk_skip, bulk_window=args.bulk_window,
            fixed_frac_min=args.fixed_frac_min, div_min=args.divergence_min,
            random_entropy_excess_min=args.random_entropy_excess_min,
        )

    if args.exclude_read and args.exclude_location:
        exclude = {(args.exclude_read.upper(), args.exclude_location)}
        profile["best_guess"] = {
            "barcode": best_guess(profile, "fixed_run_len", exclude=exclude),
            "umi": best_guess(profile, "random_run_len", exclude=exclude),
        }

    if args.out:
        with open(args.out, "w") as fh:
            json.dump(profile, fh, indent=2)

    if args.emit_shell:
        bc = profile["best_guess"]["barcode"]
        umi = profile["best_guess"]["umi"]
        print(f"BEST_BARCODE_READ={bc['read'] if bc else ''}")
        print(f"BEST_BARCODE_LOCATION={bc['location'] if bc else ''}")
        print(f"BEST_BARCODE_LENGTH={bc['length'] if bc else 0}")
        print(f"BEST_UMI_READ={umi['read'] if umi else ''}")
        print(f"BEST_UMI_LOCATION={umi['location'] if umi else ''}")
        print(f"BEST_UMI_LENGTH={umi['length'] if umi else 0}")

    if args.check_kind:
        if not (args.check_read and args.check_location and args.check_length is not None):
            print("error: --check-kind requires --check-read, --check-location, --check-length",
                  file=sys.stderr)
            sys.exit(2)
        detected_len = run_len_at(profile, args.check_kind, args.check_read, args.check_location)
        diff = abs(detected_len - args.check_length)
        match = diff <= args.tolerance
        best = profile["best_guess"]["barcode" if args.check_kind == "barcode" else "umi"]
        result = {
            "check_kind": args.check_kind,
            "claimed": {"read": args.check_read.upper(), "location": args.check_location,
                        "length": args.check_length},
            "detected_run_len_at_claimed_location": detected_len,
            "tolerance": args.tolerance,
            "match": match,
            "best_guess_elsewhere": best,
        }
        print(json.dumps(result, indent=2))
        sys.exit(0 if match else 1)

    if not args.out and not args.emit_shell and not args.check_kind:
        # default: dump profile to stdout
        print(json.dumps(profile, indent=2))


if __name__ == "__main__":
    main()
