# Monitoring a Run

[← Back to README](../README.md)

```bash
cd /path/to/tracktx
python3 nfmon.py --from-start --tail 80
```

- **Trace-aware progress**: nfmon reads the latest `results*/trace/trace.txt` to compute a snapshot-based progress bar (% of tasks that are done or cached out of all known tasks).
- **Task runtimes**: per-task durations are taken from the Nextflow trace when available; running tasks show time since they started, not since nfmon was launched.
- **Logs pane**: for TrackTx modules, nfmon tails the module logs (`preprocess_reads.log`, `align_reads.log`, `tracks.log`, etc.) and surfaces `PREP | ...`, `ALIGN | ...`, `TRACKS | ...` stage lines. When only Nextflow wrapper logs exist, it shows a small `<waiting for module logs>` placeholder instead of noisy `+ set -e` output.

**Custom output dir?** nfmon auto-detects trace files in `results/` and `results_*/` (e.g. `results_test_PE3/trace/trace.txt`). To force a specific trace path:

```bash
python3 nfmon.py --trace results_test_PE3/trace/trace.txt
```

**Install Rich (optional, for enhanced UI):** The monitor works without it (basic curses fallback), but for a better layout, colors, and live updates:
```bash
pip install rich
# or with conda:
conda install -c conda-forge rich
```

Features:
- Real-time task progress and resource usage (based on Nextflow trace)
- Live log tailing for active tasks (TrackTx module logs preferred over `.command.*`)
- CPU/memory/load monitoring
- Per-task performance metrics and slow-task detection

**Header fields in nfmon:**

- `Run` / `Session`: extracted from `.nextflow.log`. If nfmon cannot find that log, these may show as `n/a`.
- `Exec`: Nextflow executor (often `local`, even when using Docker).
- `Container`: best-effort detection of container engine (e.g. `docker`, `singularity`) based on the Nextflow log.
- `Mode`:
  - `full`: both `.nextflow.log` and `trace.txt` available; all metrics accurate.
  - `log_only`: only `.nextflow.log` found; progress/runtimes are approximate.
  - `trace_only`: only trace file available; task metrics are good, but run metadata may be limited.
  - `fs_only`: only `work/` is visible; nfmon falls back to filesystem heuristics.

**Monitor options:**
```bash
python3 nfmon.py --help                          # All options
python3 nfmon.py --filter "alignment"            # Watch specific tasks
python3 nfmon.py --all-logs                      # See all task logs
python3 nfmon.py --oneshot                       # Quick snapshot
python3 nfmon.py --oneshot --json status.json    # Export JSON
```
