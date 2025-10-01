# Re-run Instructions for Bug Fixes

## Quick Start

```bash
# 1. Delete affected results (keeps alignment and tracks cached)
rm -rf results/07_functional_regions results/10_qc results/11_reports results/14_combined_report

# 2. Re-run pipeline (only affected steps execute)
nextflow run main.nf -profile docker -resume
```

**Estimated time**: 10-20% of original runtime (most work cached)

---

## What Will Happen

### Steps That Will Re-execute ✅
1. **Functional Regions** (step 07) - With corrected sequential masking
2. **QC Module** (step 10) - With real metrics instead of dummy data
3. **Individual Reports** (step 11) - With corrected percentages
4. **Combined Report** (step 14) - With accurate cross-sample comparisons

### Steps That Stay Cached ⚡
- ✅ Read processing and trimming
- ✅ Genome alignment (BAM files)
- ✅ Track generation (bedGraphs, BigWigs)
- ✅ Divergent transcription detection
- ✅ Normalization
- ✅ Pol II metrics calculation

---

## Verification Checklist

After re-running, verify these changes:

### ✅ Functional Regions Fixed
```bash
# Check functional regions summary for one sample
cat results/07_functional_regions/Control_r2/functional_regions_summary.tsv
```

**Expected output:**
```
region                      signal        region_count
Promoter                    150000        7957
DivergentTx                 20000         7957
CPS                         50000         7957
Gene body                   650000        7956
Short genes                 10000         593
Enhancers                   15000         15450
Termination window          100000        7957
Non-localized polymerase    60000         0       ← Should be ~5-10%, not 80%!
```

### ✅ QC Module Implemented
```bash
# Check QC JSON has real data
cat results/10_qc/Control_r2/qc_pol2.json
```

**Expected output** (real metrics, not `{"status": "minimal_test"}`):
```json
{
  "sample_id": "Control_r2",
  "total_reads_raw": 10000000,
  "mapped_reads": 9500000,
  "duplicate_reads": 2000000,
  "mapq_ge_reads_nodup": 7000000,
  ...
}
```

### ✅ Reports Updated
```bash
# View individual sample report
open results/11_reports/samples/Control_r2/Control_r2.report.html

# View combined report
open results/14_combined_report/tracktx_combined_report.html
```

**Expected changes in reports:**
- Non-localized percentage: ~5-10% (not 80-90%)
- Gene body percentage: 50-70% (largest category)
- QC metrics displayed (not "n/a")

---

## Troubleshooting

### "Work directory not found" errors
```bash
# Full re-run (slower but guaranteed to work)
rm -rf work/ results/
nextflow run main.nf -profile docker
```

### Pipeline fails to resume
```bash
# Clean Nextflow cache and try again
nextflow clean -f
rm -rf results/07_functional_regions results/10_qc results/11_reports results/14_combined_report
nextflow run main.nf -profile docker -resume
```

### Results still look wrong
1. Verify you're using the latest code:
   ```bash
   git status
   git log --oneline -5
   ```

2. Check the bug was actually fixed:
   ```bash
   grep -A 2 "Sequential masking" bin/functional_regions.py
   # Should show uncommented bg_subtract calls
   ```

3. Review BUGFIX_SUMMARY.md for detailed explanation

---

## Alternative: Full Clean Re-run

If you want to re-run everything from scratch:

```bash
# WARNING: This deletes ALL results and work files
rm -rf results/ work/ .nextflow*

# Fresh run (will take full pipeline time)
nextflow run main.nf -profile docker
```

---

## Expected Timeline

| System | Affected Steps | Estimated Time |
|--------|---------------|----------------|
| Laptop (4 cores) | 30-60 minutes | ~10% of original |
| Workstation (12 cores) | 15-30 minutes | ~10% of original |
| Server (32+ cores) | 5-15 minutes | ~10% of original |

---

## Questions?

- **Detailed biology**: See `docs/FUNCTIONAL_REGIONS.md`
- **Bug details**: See `BUGFIX_SUMMARY.md`
- **General help**: See `README.md`

---

**Ready?** Run the commands at the top to get started! 🚀
