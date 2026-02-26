# TrackTx test_PE — Paired-End PRO-seq

Test with 1 sample: 10% subset of SRR12267707 (GSE154746 K562, Vihervaara lab).

- **Layout:** Paired-end
- **Spike-in:** Drosophila S2 (dm6)
- **Adapter trimming:** TruSeq small RNA
- **UMI:** 7 nt at 5'

**Setup (run once):**
```bash
./scripts/download_and_subset_test_data.sh PE
# Or with Docker (same image as pipeline):
./scripts/download_and_subset_test_data.sh --docker PE
```

**Run pipeline:**
```bash
./run_pipeline.sh --samplesheet test_PE/samplesheet_PE.csv --params-file test_PE/params_PE.yaml --output_dir ./results_test_PE
```
