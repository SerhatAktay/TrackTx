# v9.0: Return to Gene-Based Assignment Strategy

## Date: 2025-10-08
## Philosophy Change: DT-Driven Discovery, Gene-Based Assignment

---

## Executive Summary

**v9.0 implements the original TrackTx strategy:** Use DT sites to **discover** active genes, then use gene annotations to **define** functional regions for read assignment.

This is a major architectural change from v8.x, which mixed DT sites (variable width) with gene-based regions (fixed width).

---

## The Problem with v8.x

**Mixed coordinate systems created inconsistency:**

```
v8.2 Strategy:
├─ Promoter: Use DT sites (~1-2kb, variable)
├─ Divergent: Use DT sites (~1-2kb, variable)
├─ Gene body: Use gene annotation (fixed)
├─ CPS: Use gene annotation (fixed)
└─ ... inconsistent!
```

**Issues:**
1. Promoter/divergent widths depend on DT detection quality (not biology)
2. Variable region sizes make comparisons difficult
3. DT site quality affects functional assignment (not just discovery)

---

## The v9.0 Solution

**Clean separation of concerns:**

```
v9.0 Strategy:
1. DT sites → DISCOVER which genes are active
2. Gene annotations → DEFINE all functional regions
```

**All regions have consistent, biologically-meaningful sizes:**
- Promoter: TSS-250 to TSS+249 (500bp)
- Divergent: TSS-750 to TSS-251 (500bp, opposite strand)
- Gene body: TSS+250 to CPS-501 (variable, based on gene length)
- CPS: CPS-500 to CPS+499 (1000bp)
- Termination: CPS+500 to CPS+10499 (10kb)

---

## Key Changes

### 1. Active Gene Identification

**OLD (v8.2):** TSS ±250bp
**NEW (v9.0):** TSS ±500bp

```python
# Create promoter regions for ALL genes (TSS ±500bp)
# This is ONLY for identifying active genes
prom_start = gene['TSS'] - 500
prom_end = gene['TSS'] + 500

# Genes with DT sites overlapping this region = active
```

### 2. Promoter Assignment

**OLD (v8.2):** Use DT sites directly (variable width)
**NEW (v9.0):** Use gene-based ppPolII.txt (fixed 500bp)

```python
# Write promoter regions (ppPolII.txt) - all active genes
pp_file = OUT / "ppPolII.txt"
for gene in active_genes:
    pps, ppe = clamp(gene['PPs'], gene['PPe'])  # TSS-250..+249
    f.write(f"{chrom}\t{pps}\t{ppe}\t{gene_name}\t.\t{strand}\n")

# Assignment: same-strand reads from ppPolII.txt
run([BT, "intersect", "-s", "-u", "-wa", 
     "-a", reads, "-b", ppPolII], promoter_reads)
```

### 3. Divergent Assignment

**OLD (v8.2):** Use DT sites (variable width, opposite-strand)
**NEW (v9.0):** Use gene-based divTx.txt (fixed 500bp, opposite-strand)

```python
# Write divergent regions (divTx.txt) - all active genes
div_file = OUT / "divTx.txt"
for gene in active_genes:
    divs, dive = clamp(gene['DIVs'], gene['DIVe'])  # TSS-750..-251
    f.write(f"{chrom}\t{divs}\t{dive}\t{gene_name}\t.\t{strand}\n")

# Assignment: opposite-strand reads from divTx.txt
run([BT, "intersect", "-S", "-u", "-wa", 
     "-a", reads, "-b", divTx], divergent_reads)
```

### 4. Sequential Assignment Order

**Unchanged:**
1. Promoter (ppPolII.txt, same-strand)
2. Divergent (divTx.txt, opposite-strand)
3. CPS (CPS.txt, same-strand)
4. Gene body (geneBody.txt, same-strand)
5. Short genes (shortGenes.txt, same-strand)
6. Enhancers (DT sites not at promoters, unstranded)
7. Termination (TW.txt, same-strand)
8. Non-localized (remaining)

---

## Expected Results

### Region Sizes (Consistent!)

| Region | Width | Notes |
|--------|-------|-------|
| Promoter | 500bp | TSS-250..+249 |
| Divergent | 500bp | TSS-750..-251 |
| CPS | 1000bp | CPS-500..+499 |
| Gene body | Variable | TSS+250 to CPS-501 |
| Short genes | ≤750bp | Full gene if length ≤750bp |
| Termination | 10kb | CPS+500..+10499 |
| Enhancers | Variable | DT sites not at promoters |

### Expected Distribution

Based on original TrackTx pipeline:

| Category | Expected % |
|----------|-----------|
| Promoter | 10-15% |
| Divergent Tx | 5-10% |
| Gene body | 20-30% |
| Short genes | 8-12% |
| CPS | 1-2% |
| Termination | 2-5% |
| Enhancers | 10-15% |
| **Non-localized** | **10-20%** ← Should be MUCH lower! |

---

## Comparison: v8.2 vs v9.0

### v8.2 (DT-site-based)
```
Philosophy: DT sites are the functional elements
- Promoter: DT sites (1-2kb, same-strand)
- Divergent: DT sites (1-2kb, opposite-strand)
- Others: Gene-based

Problem: Inconsistent region definitions
Result: 38% non-localized
```

### v9.0 (Gene-based)
```
Philosophy: DT sites identify active genes; genes define regions
- ALL regions: Gene-based, consistent sizes
- DT sites: Only for discovery + enhancers

Advantage: Consistent, biologically-meaningful
Expected: 10-20% non-localized
```

---

## Implementation Details

### File Structure

```
coord_files = {
    'promoter':    'ppPolII.txt',      # TSS-250..+249 (NEW!)
    'divergent':   'divTx.txt',        # TSS-750..-251
    'gene_body':   'geneBody.txt',     # TSS+250 to CPS-501
    'short_genes': 'shortGenes.txt',   # Full gene if ≤750bp
    'cps':         'CPS.txt',          # CPS±500bp
    'termination': 'TW.txt',           # CPS+500..+10kb
}

enhancers = 'enhancers.bed'  # DT sites NOT at TSS ±500bp
```

### Changes from v8.2

| Aspect | v8.2 | v9.0 |
|--------|------|------|
| Active gene detection | TSS ±250bp | TSS ±500bp |
| Promoter source | DT sites | Gene-based (ppPolII.txt) |
| Divergent source | DT sites | Gene-based (divTx.txt) |
| Region consistency | Mixed | All gene-based (except enhancers) |
| Biological meaning | Detection-dependent | Structure-based |

---

## Testing

### Run the pipeline:

```bash
nextflow run main.nf -resume
```

### Check results:

```bash
# View summary
cat results/07_functional_regions/*/functional_regions_summary.tsv

# Expected improvements:
# - Promoter: ~10-15% (consistent)
# - Divergent: ~5-10% (consistent)
# - Non-localized: ~10-20% (much lower!)
```

### Verify logs:

```bash
grep "mode=gene-based" results/07_functional_regions/*/functional_regions.log
# Should show: "mode=gene-based assignment (v9.0)"
```

---

## Biological Rationale

### Why Gene-Based is Better

1. **Biological Consistency**
   - Promoter size should be based on gene structure, not detection sensitivity
   - All promoters treated equally (500bp)

2. **Comparative Analysis**
   - Consistent regions allow fair comparison across samples
   - Region size doesn't depend on signal strength

3. **Interpretability**
   - "Promoter" means TSS ±250bp (clear biological definition)
   - Not "wherever we detected strong bidirectional signal"

4. **Robustness**
   - Poor DT detection doesn't shrink promoter regions
   - Assignment quality decoupled from detection quality

### DT Sites Still Critical

DT sites remain essential for:
- **Discovery**: Which genes are actually transcribing?
- **Enhancers**: Intergenic transcription sites
- **Validation**: Cross-check with gene annotations

But they no longer define region boundaries for assignment.

---

## Migration Notes

### From v8.2 to v9.0

**What stays the same:**
- DT detection (module 09)
- Gene annotation loading
- Sequential assignment order
- Output format

**What changes:**
- Active gene identification (TSS ±250bp → ±500bp)
- Promoter/divergent use gene-based regions (not DT sites)
- All regions now gene-based (except enhancers)

**Breaking changes:**
- None (output format unchanged)
- Results will differ (more consistent, lower non-localized)

---

## References

- Original TrackTx bash pipeline (feature_coordinates.R, active_genes.R)
- v8.0-v8.2: DT-site-based attempts
- v9.0: Return to gene-based strategy (lessons learned)

---

## Conclusion

v9.0 represents a return to first principles: **use the data to discover, use annotations to define**. This clean separation creates consistent, interpretable, biologically-meaningful functional regions.

The previous v8.x approaches taught us that mixing detection-based and annotation-based coordinates creates more problems than it solves. v9.0 fixes this by keeping them separate.

