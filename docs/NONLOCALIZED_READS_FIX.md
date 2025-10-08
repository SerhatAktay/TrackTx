# Non-Localized Reads Issue - Root Cause & Fix

## Date: 2025-10-08
## Version: TrackTx v8.1

---

## Overview

This document describes the bug causing excessive "Non-localized polymerase" reads and the fix implemented in `functional_regions.py` v8.1.

---

## The Bug

### Root Cause

In `functional_regions.py` v8.0, the **promoter read assignment** (Step 1 of sequential assignment) was using **unstranded intersection**:

```python
# OLD CODE (v8.0):
run([BT, "intersect", "-u", "-wa", "-a", current_reads, "-b", active_promoters_bed], str(prom_reads))
# No -s or -S flag = unstranded intersection
```

### Why This Caused Mis-assignment

**Divergent transcription sites represent BIDIRECTIONAL transcription:**
- **Same-strand reads** (relative to gene) = Promoter-proximal paused Pol II
- **Opposite-strand reads** = Divergent/upstream antisense transcription

**The bug:**
1. Step 1 (Promoter) took **ALL reads (both strands)** overlapping active promoter DT sites
2. This incorrectly assigned **opposite-strand reads** to "Promoter" category
3. Step 2 (Divergent Tx) then had no reads left to assign to "Divergent Tx"
4. Later steps couldn't properly categorize these reads
5. Result: Many reads ended up as "Non-localized polymerase"

### Visual Example

```
Gene on + strand, TSS at position 1000:

DT site at promoter: chr1:750-1250

Reads:
  Read A: chr1:900-901  strand=+  → Should be "Promoter" ✓
  Read B: chr1:900-901  strand=-  → Should be "Divergent Tx" ✗

OLD BEHAVIOR:
  Step 1 (Promoter): Both Read A and Read B assigned to "Promoter"
  Step 2 (Divergent): No reads left to assign
  Result: Read B incorrectly categorized as "Promoter"

NEW BEHAVIOR:
  Step 1 (Promoter): Only Read A assigned to "Promoter" (same-strand)
  Step 2 (Divergent): Read B assigned to "Divergent Tx" (opposite-strand)
  Result: Both reads correctly categorized ✓
```

---

## The Fix

### Changes Made

**1. Active promoter BED format** (lines 244-267)

Changed from BED4 (no strand) to BED6 (with strand):

```python
# NEW: Use -wa -wb to get DT coords AND gene strand
run([BT, "intersect", "-wa", "-wb", "-a", dt_bed, "-b", str(promoter_regions)], 
    str(active_promoters_tmp))

# Reformat to BED6: chr, start, end, gene_name, signal, strand
# Strand comes from the overlapping gene promoter region
```

**2. Promoter read assignment** (lines 455-456)

Added `-s` flag for same-strand intersection:

```python
# NEW: Use -s flag to only assign same-strand reads
run([BT, "intersect", "-s", "-u", "-wa", "-a", current_reads, "-b", active_promoters_bed], 
    str(prom_reads))
run([BT, "intersect", "-s", "-v", "-a", current_reads, "-b", active_promoters_bed], 
    str(prom_removed))
```

**3. Updated BED9 output** (lines 651-653)

Now handles BED6 format with strand information.

---

## Expected Impact

### Before Fix (v8.0)
- **Promoter**: Inflated (includes both same- and opposite-strand reads)
- **Divergent Tx**: Deflated (opposite-strand reads already removed)
- **Non-localized polymerase**: Inflated (mis-assigned reads end up here)

### After Fix (v8.1)
- **Promoter**: Reduced to only same-strand reads (correct)
- **Divergent Tx**: Increased to include opposite-strand reads (correct)
- **Non-localized polymerase**: Reduced significantly (proper assignment)

### Biological Interpretation

This fix ensures that:
- **Promoter** category truly represents paused Pol II on the gene
- **Divergent Tx** category captures upstream antisense transcription
- **Non-localized polymerase** only contains reads that genuinely don't match any functional category

---

## Testing

### To verify the fix worked:

1. **Check promoter counts decreased:**
   ```bash
   # Compare old vs new functional_regions_summary.tsv
   grep "Promoter" functional_regions_summary.tsv
   ```

2. **Check divergent Tx counts increased:**
   ```bash
   grep "DivergentTx" functional_regions_summary.tsv
   ```

3. **Check non-localized decreased:**
   ```bash
   grep "Non-localized" functional_regions_summary.tsv
   ```

4. **Verify active_promoters.bed has strand column:**
   ```bash
   head active_promoters.bed
   # Should show: chr start end gene_name signal strand
   ```

5. **Check debug logs:**
   ```bash
   grep "Step 1: Assigning promoter reads" functional_regions.log
   # Should say "same-strand only"
   ```

---

## Technical Details

### Sequential Assignment Order (Unchanged)

1. **Promoter** (same-strand, -s) ← FIXED
2. **Divergent Tx** (opposite-strand, -S)
3. **CPS** (same-strand, -s)
4. **Gene body** (same-strand, -s)
5. **Short genes** (same-strand, -s)
6. **Enhancers** (unstranded)
7. **Termination window** (same-strand, -s)
8. **Non-localized polymerase** (remaining)

### Files Modified

- `/Volumes/samsung/TrackTx/bin/functional_regions.py` v8.0 → v8.1

### Backward Compatibility

The fix is backward compatible:
- Old results were incorrect but pipeline will run
- New results are correct
- No changes to input/output file formats (except active_promoters.bed gains strand column)

---

## References

- Issue: "Non-localized reads issue"
- Related: Divergent transcription detection (module 09)
- Related: Functional region assignment logic (module 10)

---

## Author Notes

The key insight is that **divergent transcription sites are bidirectional** - they contain signal from both the gene strand (paused Pol II at promoter) and the opposite strand (upstream divergent/antisense transcription). The assignment logic must respect this strandedness to correctly categorize reads.

