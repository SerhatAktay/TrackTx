# Non-Localized Reads Issue - Complete Root Cause & Fix

## Date: 2025-10-08
## Version: TrackTx v8.2

---

## Executive Summary

**Problem**: ~40% of reads were being assigned to "Non-localized polymerase" instead of proper functional categories.

**Root Cause**: Mismatch between DT site widths (~1-2kb) and gene-based divergent region widths (500bp fixed).

**Fix**: Use actual DT sites for both Promoter AND Divergent Tx assignment, partitioned by strand.

---

## The Journey: Two Bugs, One Solution

### Bug #1: Unstranded Promoter Assignment (v8.0 → v8.1)

**Issue**: Promoter read assignment was unstranded, capturing both same- and opposite-strand reads.

**Impact**: Moderate - some mis-assignment but not the main cause of non-localized reads.

**Fix (v8.1)**:
- Made active_promoters.bed include strand information (BED4 → BED6)
- Added `-s` flag to promoter read assignment (same-strand only)

**Result**: Helped, but non-localized reads still ~40%!

---

### Bug #2: Narrow Divergent Regions vs Wide DT Sites (v8.1 → v8.2)

**The Discovery**:

After implementing v8.1, we found:
- DT sites: Mean width ~2,051 bp, Median ~1,445 bp
- Gene-based divergent regions: Fixed 500 bp (TSS-750 to TSS-250)

**What Was Happening**:

```
Example: Gene on + strand, TSS = 10,000

DT Site (detected from data):
  chr1:9,000-11,500  (2,500 bp wide, bidirectional signal)
  Assigned strand = + (from overlapping gene)

Step 1 (Promoter, v8.1):
  bedtools intersect -s (same-strand)
  → Takes + strand reads from chr1:9,000-11,500
  ✓ Correctly captures promoter-proximal Pol II

Step 2 (Divergent Tx, v8.1):
  bedtools intersect -S (opposite-strand)
  Target: Gene-based region chr1:9,250-9,750 (500 bp)
  → Takes - strand reads from chr1:9,250-9,750 ONLY
  
  ✗ PROBLEM: DT site extends 9,000-11,500
             But divergent region is only 9,250-9,750
             - strand reads from 9,000-9,250 = NOT assigned!
             - strand reads from 9,750-11,500 = NOT assigned!
             → These become "Non-localized"!
```

**Visual**:

```
Position:     9000        9250   9750       10000      11500
              |-----------|-----|-----------|----------|
DT Site:      [=========================================]  (2.5kb)
              
Promoter (v8.1):
(+ strand)    [=========================================]  ✓ CORRECT

Divergent (v8.1):
(- strand)              [=====]                            ✗ TOO NARROW!
              ^^^^^^^            ^^^^^^^^^^^^^^^^^^^^^
              LOST!              LOST!

Divergent (v8.2):
(- strand)    [=========================================]  ✓ CORRECT
```

---

## The Fix (v8.2)

### Change

In `sequential_read_assignment()`, Step 2 now uses **the same active_promoters_bed** as Step 1:

```python
# OLD (v8.1): Used gene-based divergent regions
run([BT, "intersect", "-S", "-u", "-wa", 
     "-a", current_reads, 
     "-b", coord_files['divergent']],  # ← 500bp gene-based regions
     str(div_reads))

# NEW (v8.2): Use actual DT sites
run([BT, "intersect", "-S", "-u", "-wa", 
     "-a", current_reads, 
     "-b", active_promoters_bed],      # ← Actual DT sites (~1-2kb)
     str(div_reads))
```

### Logic

Now the active promoter DT sites are **partitioned by strand**:
- **Same-strand reads** → "Promoter" (paused Pol II on gene)
- **Opposite-strand reads** → "Divergent Tx" (upstream antisense transcription)

This correctly captures the entire bidirectional transcription signal from the DT site.

---

## Expected Results

### Before (v8.1)
```
Promoter:        ~15% (same-strand from DT sites)
Divergent Tx:    ~1%  (opposite-strand from narrow 500bp regions)
Gene body:       ~20%
Short genes:     ~9%
Enhancers:       ~12%
CPS:             ~1%
Termination:     ~2%
Non-localized:   ~40% ← PROBLEM!
```

### After (v8.2)
```
Promoter:        ~10-15% (same-strand from DT sites)
Divergent Tx:    ~10-15% (opposite-strand from SAME DT sites)
Gene body:       ~20%
Short genes:     ~9%
Enhancers:       ~12%
CPS:             ~1%
Termination:     ~2%
Non-localized:   ~5-10% ← FIXED!
```

**Key change**: Divergent Tx should increase ~10-15x, Non-localized should decrease ~4x.

---

## Biological Interpretation

### What DT Sites Represent

Divergent transcription sites are regions of **bidirectional transcription**:
- **Gene-direction signal**: Promoter-proximal paused Pol II (productive)
- **Opposite-direction signal**: Upstream antisense/divergent transcription

### Why This Matters

The DT sites are detected from the data and represent the **actual extent** of bidirectional transcription. Using gene-based fixed windows (500bp) was biologically incorrect because:

1. Not all promoters have the same divergent transcription extent
2. DT site width correlates with transcriptional activity
3. Gene-based windows were arbitrary and too narrow

By using the actual DT sites, we now respect the **biological reality** of each locus.

---

## Testing

### To verify the fix:

1. **Re-run the pipeline**:
   ```bash
   nextflow run main.nf -resume
   ```

2. **Check functional regions summary**:
   ```bash
   cat results/07_functional_regions/*/functional_regions_summary.tsv
   ```

3. **Expected changes**:
   - ✓ Divergent Tx reads: Should increase ~10-15x (from ~10k to ~100-150k)
   - ✓ Non-localized reads: Should decrease ~4x (from ~330k to ~50-80k)
   - ✓ Promoter reads: May decrease slightly (some reads reallocated)

4. **Verify DT site usage**:
   ```bash
   grep "Step 2: Assigning divergent reads" results/07_functional_regions/*/functional_regions.log
   # Should say: "opposite-strand from active promoter DT sites"
   ```

---

## Files Modified

- `bin/functional_regions.py`: v8.0 → v8.1 → v8.2
  - v8.1: Added strand-aware promoter assignment
  - v8.2: Fixed divergent assignment to use actual DT sites

---

## Technical Notes

### Why Not Use Gene-Based Regions at All?

Gene-based regions (CPS, gene body, termination window) are still used because:
- They represent different functional zones ALONG the gene
- They're properly sized for their biological function (e.g., CPS ±500bp around TES)

But for **divergent transcription**, the DT sites ARE the ground truth. Using gene-based windows was conceptually wrong.

### Sequential Assignment Still Correct

The sequential masking order is still appropriate:
1. Promoter/Divergent (from DT sites) - highest priority
2. CPS, Gene body, Short genes - gene structure
3. Enhancers - distal DT sites
4. Termination window - downstream
5. Non-localized - truly unassigned

---

## Conclusion

The non-localized reads issue was caused by a **mismatch in scale**: DT sites (~1-2kb) vs gene-based divergent regions (500bp). By using the actual DT sites for both promoter and divergent assignment and partitioning them by strand, we now correctly capture the full bidirectional transcription signal.

This fix is both **technically correct** and **biologically meaningful**.

