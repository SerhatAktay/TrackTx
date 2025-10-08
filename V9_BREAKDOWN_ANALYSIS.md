# v9.0/v9.1 Complete Breakdown Analysis

## Summary

**Both v9.0 and v9.1 have ~40% non-localized reads due to narrow gene-based regions.**

---

## The Fundamental Problem

**Fixed 500bp gene-based windows are too narrow to capture actual transcription extent.**

### Region Sizes:
- **Promoter** (ppPolII): TSS-250 to TSS+249 = 500bp
- **Divergent** (divTx): TSS-750 to TSS-251 = 500bp  
- **DT sites** (actual data): Mean ~2,051bp, Median ~1,445bp

### Result:
- **DT sites are 3-4x wider than gene-based regions!**
- Most transcription falls outside the narrow 500bp windows
- Those reads become "Non-localized"

---

## v9.1 Step-by-Step Breakdown

### Starting Point
- **Total reads**: 13,601,113

### Step 1: Promoter Assignment
```
bedtools intersect -s -u -wa -a reads -b ppPolII.txt
  → Assigned: 308,317 same-strand reads (2.3%)

bedtools intersect -v -a reads -b ppPolII.txt  
  → Removed: 474,760 total (both strands)
  → Lost: 166,443 opposite-strand reads (removed but not assigned)
  
Remaining: 13,126,353
```

**What happened:**
- Assigned only same-strand reads from narrow 500bp promoter
- Removed ALL reads (both strands) from that region
- Lost 166k opposite-strand reads (antisense, correctly filtered)

### Step 2: Divergent Assignment  
```
bedtools intersect -S -u -wa -a reads -b divTx.txt
  → Assigned: 180,597 opposite-strand reads (1.3%)

bedtools intersect -v -a reads -b divTx.txt
  → Removed: 288,966 total (both strands)
  → Lost: 108,369 same-strand reads (removed but not assigned)
  
Remaining: 12,837,387
```

**What happened:**
- Assigned only opposite-strand reads from narrow 500bp divergent region
- Removed ALL reads (both strands) from that region
- Lost 108k same-strand reads (correctly filtered)

### Total Lost to Unstranded Removal
- **1,000,000 reads** removed across all steps
- These are antisense/non-functional reads in wrong regions
- **This filtering is correct!**

### Final Distribution

| Category | Reads | % | Signal |
|----------|-------|---|--------|
| Promoter | 308,317 | 2.3% | 27,499 |
| Divergent | 180,597 | 1.3% | 13,304 |
| Gene body | 3,351,488 | 24.6% | 195,759 |
| Short genes | 1,701,104 | 12.5% | 89,899 |
| CPS | 147,724 | 1.1% | 8,823 |
| Enhancers | 1,450,662 | 10.7% | 76,548 |
| Termination | 436,473 | 3.2% | 19,860 |
| **Non-localized** | **5,024,748** | **36.9%** | **301,870** |

---

## Why v9.1 is Still Bad

### Promoter: Only 2.3%
- 16,803 genes × 500bp = ~8.4M bp total
- 308,317 reads = **18 reads per promoter**
- Most promoter transcription falls outside TSS-250..+249

### Divergent: Only 1.3%  
- 16,803 genes × 500bp = ~8.4M bp total
- 180,597 reads = **11 reads per divergent region**
- Most divergent transcription falls outside TSS-750..-251

### Non-localized: 36.9%
- These are likely reads from:
  - Promoters outside the TSS±250 window
  - Divergent transcription outside TSS-750..-251
  - Real intergenic transcription

---

## Comparison: All Versions

| Version | Strategy | Promoter % | Divergent % | Non-loc % |
|---------|----------|-----------|-------------|-----------|
| **v8.2** | DT sites (1-2kb) | 14.8% | 5.8% | 38.1% |
| **v9.0** | Gene 500bp, strand removal | 3.3% | 1.6% | 43.3% |
| **v9.1** | Gene 500bp, unstranded removal | 2.3% | 1.3% | 36.9% |

**v8.2 was actually the best!**
- Used actual DT site extents (~1-2kb)
- Captured more real transcription
- Still had 38% non-localized, but better than v9.x

---

## The Dilemma

### Option A: DT Sites (v8.2)
✅ Captures actual transcription extent  
✅ Biologically accurate  
✗ Inconsistent region sizes  
✗ Still 38% non-localized

### Option B: Gene-Based 500bp (v9.1)
✅ Consistent region sizes  
✅ Clean conceptual model  
✗ Misses most transcription  
✗ 37% non-localized

### Option C: Wider Gene-Based Regions
✅ Consistent region sizes  
✅ Might capture more transcription  
? How wide to make them?  
? Just converge back to DT sites?

---

## Key Insights

1. **Fixed 500bp windows don't match biology**
   - Transcription extent varies by gene/condition
   - DT sites show actual extent (~1-2kb average)
   
2. **Unstranded removal is correct**
   - Filters antisense reads in wrong regions
   - ~1M reads correctly removed as non-functional

3. **Non-localized reads are the real mystery**
   - 5M reads (37%) don't fit any category
   - Are they:
     - Real promoter reads outside narrow windows?
     - Intergenic transcription?
     - Noise/artifacts?

4. **The old bash script likely has the same issue**
   - Uses same 500bp windows
   - Would have same low promoter/divergent percentages
   - Unless there's something else we're missing...

---

## Next Steps?

1. **Check what the OLD bash script actually produces**
   - Run it on same data
   - Compare percentages
   
2. **Investigate non-localized reads**
   - Where are they in the genome?
   - Do they overlap DT sites?
   - Are they near genes?

3. **Consider hybrid approach**
   - Use DT sites for promoter/divergent (like v8.2)
   - Use gene-based for gene body/CPS/termination
   - Accept inconsistent sizes as biological reality

4. **Question the assumption**
   - Maybe 35-40% non-localized IS normal?
   - Maybe that's intergenic/enhancer/other transcription?

