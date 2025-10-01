# Functional Regions Classification

## Overview

TrackTx assigns nascent RNA-seq signal to functional genomic regions using a **hierarchical, sequential masking** approach. Each read/signal is assigned to **exactly ONE** functional category to provide biologically meaningful quantification.

## Biological Hierarchy

Signal is assigned in the following order (earlier categories have priority):

```
1. Promoter (ppPolII)
   ↓
2. Divergent Transcription (opposite strand)
   ↓
3. Gene Body
   ↓
4. Cleavage & Polyadenylation Site (CPS)
   ↓
5. Short Genes
   ↓
6. Enhancers (divergent regions not overlapping promoters)
   ↓
7. Termination Window
   ↓
8. Non-localized Polymerase (residual)
```

## Region Definitions

### 1. Promoter (ppPolII - Promoter-Proximal Pol II)
- **Location**: TSS ± window (default: -250 to +250 bp)
- **Strand**: Same as gene
- **Biology**: Paused/initiating Pol II at transcription start sites
- **Expected %**: 10-20% of total signal

### 2. Divergent Transcription
- **Location**: TSS upstream window (default: -750 to -250 bp)
- **Strand**: **Opposite** to gene
- **Biology**: Upstream divergent/antisense transcription
- **Expected %**: 2-5% of total signal

### 3. Gene Body
- **Location**: From promoter end to CPS start
- **Strand**: Same as gene
- **Biology**: Productive elongating Pol II through gene
- **Expected %**: 50-70% of total signal (largest category)

### 4. Cleavage & Polyadenylation Site (CPS)
- **Location**: TES ± 500 bp
- **Strand**: Same as gene
- **Biology**: 3' end processing, termination initiation
- **Expected %**: 5-10% of total signal

### 5. Short Genes
- **Location**: Entire gene if length ≤ div_outer (default: 750 bp)
- **Strand**: Same as gene
- **Biology**: Small genes where regions would overlap
- **Expected %**: 1-3% of total signal

### 6. Enhancers
- **Location**: Divergent transcription regions NOT overlapping promoters
- **Strand**: Unstranded (both strands)
- **Biology**: Active enhancers with bidirectional eRNA production
- **Expected %**: 2-5% of total signal

### 7. Termination Window
- **Location**: Downstream of CPS (default: +10 kb)
- **Strand**: Same as gene
- **Biology**: Terminating/torpedo Pol II readthrough
- **Expected %**: 5-15% of total signal

### 8. Non-localized Polymerase
- **Location**: All remaining genome
- **Strand**: Both
- **Biology**: 
  - Background/non-specific signal
  - Intergenic transcription
  - Unannotated genes/transcripts
  - Technical noise
- **Expected %**: 5-10% of total signal
- **⚠️ WARNING**: If >20%, indicates potential annotation or quality issues

## Sequential Masking Logic

```python
# Start with full signal
residual_signal = total_signal

# 1. Assign promoters (same-strand)
promoter_signal = map_signal(promoter_regions, residual_signal)
residual_signal -= promoter_signal  # MASK

# 2. Assign divergent (opposite-strand, no masking needed)
divergent_signal = map_signal(divergent_regions, residual_signal)

# 3. Assign gene body (same-strand)
genebody_signal = map_signal(genebody_regions, residual_signal)
residual_signal -= genebody_signal  # MASK

# ... continue for all regions ...

# 8. Whatever remains is non-localized
nonlocalized_signal = residual_signal
```

## Quality Control

### Healthy Sample Metrics
- **Non-localized**: 5-10%
- **Gene body**: 50-70% (dominant category)
- **Promoter**: 10-20%
- **Termination window**: 5-15%
- **Divergent + Enhancers**: 5-10% combined

### Warning Signs
- **Non-localized >20%**: Poor annotation, low quality, or wrong genome assembly
- **Promoter >30%**: Excessive pausing (could be biological or technical)
- **Gene body <40%**: Poor elongation or degraded sample
- **Divergent >15%**: Unusual transcription patterns

## Parameters

All parameters are configurable in `params.yaml`:

```yaml
functional_regions:
  prom_up: 250           # Promoter upstream (bp)
  prom_down: 250         # Promoter downstream (bp)
  div_inner: 250         # Divergent inner boundary (bp)
  div_outer: 750         # Divergent outer boundary (bp)
  tw_length: 10000       # Termination window length (bp)
  min_signal: 0.0        # Minimum signal threshold
  count_mode: "signal"   # "signal" (bedGraph sum) or "event" (legacy)
  allow_unstranded: true # Allow unstranded features
```

## Troubleshooting

### High Non-localized Percentage

**If you see >20% non-localized:**

1. **Check genome assembly match**
   - Ensure FASTQ data matches reference genome
   - Example: hg38 data requires hg38 reference (not hg19/GRCh37)

2. **Verify GTF annotation**
   - Use comprehensive annotation (e.g., GENCODE comprehensive, not basic)
   - Check annotation version matches genome assembly

3. **Inspect sample quality**
   - Check duplication rate (<30% expected)
   - Review FastQC/fastp reports
   - Look for adapters or contamination

4. **Confirm library type**
   - PRO-seq, ChIP-seq expected: 5-10% non-localized
   - GRO-seq may have 10-15% (more nascent/unstable transcripts)
   - Total RNA-seq: not appropriate for this pipeline

5. **Review parameter settings**
   - Default windows may not suit all organisms
   - Consider adjusting for very compact genomes

## Implementation Details

### Input
- **Signal tracks**: 3' siCPM bedGraphs (preferred) or CPM
- **Active genes**: Intersection of promoters with divergent transcription
- **Strand-aware**: Plus/minus strand signals matched to gene strand

### Output Files
- `functional_regions.bed` - BED9 track for genome browser
- `functional_regions_summary.tsv` - Per-category signal totals
- `README_functional_regions.txt` - Run parameters

### Algorithm
- Tool: `bedtools map -o sum` for signal integration
- Masking: `bedtools subtract -A` to remove assigned regions
- Sorting: Required for all BED/bedGraph inputs (LC_ALL=C sort)

## References

This hierarchical approach is based on nascent transcription biology:

- Core et al. (2008) "Nascent RNA sequencing reveals widespread pausing..."
- Kwak et al. (2013) "Precise maps of RNA polymerase reveal how promoters direct..."
- Jonkers & Lis (2015) "Getting up to speed with transcription elongation..."

## Version History

- **v6.0**: Sequential masking implemented with clear hierarchy
- **v5.x**: Bug fix - masking was accidentally disabled (led to ~80% non-localized)
- **v4.x**: Signal-based counting (replaced event-based)

