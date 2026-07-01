#!/usr/bin/awk -f
# =============================================================================
# add_nh_tags.awk — deterministic NH:i tagging for bowtie2 -k output
# =============================================================================
#
# Why:
#   bowtie2 does not emit the NH (number-of-hits) tag, and in -k mode it sets
#   MAPQ to 255 ("unavailable"). Downstream "unique read" filtering therefore
#   cannot rely on MAPQ when -k is active. This script counts, per read segment
#   (i.e. per mate for paired-end), how many alignments bowtie2 reported and
#   writes that count as an NH:i tag on every record. A uniquely-mapped read
#   then has NH:i:1; a multimapper has NH:i:>1 on all of its records.
#
# Input requirement:
#   Records MUST be grouped by read name. Pipe through `samtools collate`
#   first (fast, O(n), no full coordinate sort needed). Header lines (@...) are
#   passed through untouched.
#
# Portability:
#   Uses only arithmetic bit tests (no gawk and()/or() builtins), so it runs
#   under mawk, the default awk on Debian/Ubuntu containers.
#
# Usage:
#   samtools collate -@ T -O -u in.bam - \
#     | samtools view -h - \
#     | awk -f add_nh_tags.awk \
#     | samtools view -@ T -b - \
#     | samtools sort -@ T -o out.bam -
# =============================================================================

function bit(f, b) { return int(f / b) % 2 }   # is SAM flag bit `b` set in `f`?

function flush(   i, f, n1, n2, n0, nh) {
  n1 = 0; n2 = 0; n0 = 0
  # First pass: count alignments per segment within this read-name group.
  for (i = 1; i <= cnt; i++) {
    f = flg[i]
    if      (bit(f, 64))  n1++          # first-in-pair  (0x40)
    else if (bit(f, 128)) n2++          # second-in-pair (0x80)
    else                  n0++          # single-end (neither mate bit set)
  }
  # Second pass: emit each record with the NH count for its segment.
  for (i = 1; i <= cnt; i++) {
    f = flg[i]
    if      (bit(f, 64))  nh = n1
    else if (bit(f, 128)) nh = n2
    else                  nh = n0
    print rec[i] "\tNH:i:" nh
  }
  cnt = 0
}

/^@/ { print; next }                    # SAM header → pass through

{
  if (cnt > 0 && $1 != cur) flush()     # read name changed → emit previous group
  cur = $1
  cnt++
  rec[cnt] = $0
  flg[cnt] = $2
}

END { if (cnt > 0) flush() }
