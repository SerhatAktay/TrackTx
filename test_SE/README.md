# TrackTx test_SE — Single-End PRO-seq

Test with 1 sample: 10% subset of SRR8669162 (GSE127844 K562, Vihervaara lab).

- **Spike-in:** Drosophila S2 (dm6)
- **Adapter trimming:** TruSeq small RNA

**Setup (run once):**
```bash
./scripts/download_and_subset_test_data.sh SE
# Or with Docker (same image as pipeline):
./scripts/download_and_subset_test_data.sh --docker SE
```

**Run pipeline:**
```bash
./run_pipeline.sh --samplesheet test_SE/samplesheet_SE.csv --params-file test_SE/params_SE.yaml --output_dir ./results_test_SE
```
