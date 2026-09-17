# Troubleshooting

[← Back to README](../README.md)

## Reading error messages

When a process fails, Nextflow prints the captured output. Look for the TRACKTX ERROR block; it names the problem and the fix:

```
═══════════════════════════════════════════════════════════════════════
TRACKTX ERROR
═══════════════════════════════════════════════════════════════════════
Module:  detect_divergent_transcription
Problem: Missing Python dependencies (numpy, pandas, scikit-learn, scipy)
Fix:     pip install numpy pandas scikit-learn scipy | Or use: -profile conda | -profile docker
═══════════════════════════════════════════════════════════════════════
```

**Note:** Progress output is written to log files only. On failure, Nextflow shows stderr (errors and the TRACKTX ERROR block). Full output is in the work dir (see `Work dir:` in the error message).

- **Quick find:** `grep -A 6 "TRACKTX ERROR"` in the output
- **Full log:** Check the `.log` file in the work dir (shown at the end of the error)
- **Resume:** Add `-resume` to continue after fixing the issue

## Common issues

**Docker not running:**
```bash
# Install Docker Desktop: https://docs.docker.com/get-docker/
# Make sure it's running before starting pipeline
```

**Out of memory (exit 137):**

OOM means RAM exhaustion, not disk space. Removing `work` and `.nextflow` frees disk space but does not fix OOM.

Common causes when you "have enough space" (disk):

1. **Docker memory limit:** Docker Desktop has its own RAM limit (Settings → Resources → Memory). The pipeline detects *host* RAM and allocates per-task memory accordingly, but containers only see Docker's limit. If Docker has 8 GB and the pipeline assumes 64 GB, multiple tasks can exceed available RAM.
   - **Fix:** `./run_pipeline.sh` auto-detects Docker memory when using the docker profile. If OOM persists: `export NXF_HOST_MEM=8` (match Docker limit), or increase Docker memory in Settings.
2. **WSL2:** WSL reports host RAM, not its own memory limit. Set `NXF_HOST_MEM` to your WSL memory limit (e.g. in `.wslconfig`).
3. **Parallelism:** Several tasks run at once; total RAM ≈ per-task × forks. If detection is wrong, total can exceed actual RAM.
4. **Large samples:** Many unaligned reads (e.g. 20M+) need more memory for spike-in alignment.

```bash
# Tell Nextflow the actual available RAM (e.g. Docker or WSL limit)
export NXF_HOST_MEM=8   # Use 8 if Docker/WSL has 8 GB
export NXF_HOST_CPUS=4   # Reduce parallelism
./run_pipeline.sh
```

**Conda environment fails:**
```bash
# Use Docker instead (more reliable)
./run_pipeline.sh -profile docker

# Or clean conda cache
conda clean --all --yes
```

**"Missing Python dependencies" (divergent transcription step):**
```bash
# Use Docker (recommended; has all deps pre-installed)
./run_pipeline.sh -profile docker

# Or use conda profile (creates env from envs/tracktx.yaml)
./run_pipeline.sh -profile conda

# Or install manually: pip install -r envs/requirements-divergent.txt
```

**Pipeline seems slow:**
- First run downloads reference genomes (~10-30 min)
- Use SSD storage for better performance
- Monitor with `python3 nfmon.py` to see bottlenecks

**Low unique-read rate (e.g. &lt;30%):**
- QC reports the unique-read rate via `uniqueness_method` in `qc_pol.json`: `NH==1` when `align.multimap_k > 1` (the default), or `MAPQ≥threshold` in single-best mode
- PRO-seq often has 30-60% uniquely-mapped reads; subset data or repetitive genomes can be lower. Multimappers are still retained in the **allMap** tracks even when excluded from unique-read metrics and gene quantification
- Single-best mode only: to use a lower MAPQ threshold add `qc: { mapq: 5 }` (or `pol: { mapq: 5 }`) to params.yaml. In `-k` mode uniqueness is exact (`NH==1`), so MAPQ thresholds don't apply
- A high multimapper fraction in repetitive regions is expected for some datasets and is not an error

**"Failed to publish file [link]"** (external drive / exFAT):

Hard links don't work on exFAT (common on USB drives). Fix:

```bash
./run_pipeline.sh --external-drive        # Cache/work on local; fixes publish + file-lock errors
```

**OverlappingFileLockException** (e.g. `preprocess_and_quality_filter_reads`, `download_genome_annotations`):

Java file-lock conflict. Common causes and fixes:

1. **Multiple runs from same directory:** Only one Nextflow run per directory at a time. Stop other runs or use a separate project copy.
2. **Stale lock from previous run:** If you used Ctrl+Z or killed the process uncleanly:
   ```bash
   rm -rf work .nextflow
   ./run_pipeline.sh
   ```
3. **USB drive with exFAT/FAT32 (macOS) / OverlappingFileLockException:** exFAT and FAT32 do not support file locking. Fix: use `./run_pipeline.sh --external-drive`, which puts cache, temp, and work (~10-50 GB) on local (`~/tmp/tracktx_cache`, `~/tmp/tracktx_work`); results stay on your project. Ensure ~20-50 GB free on internal drive. Or reformat the USB to APFS or Mac OS Extended.
4. **NFS / network / cloud-synced storage:** File locking is unreliable on NFS, SMB, iCloud, Dropbox. Set work dir to internal disk or a USB drive formatted as APFS/ext4: `export NXF_WORK=/tmp/nextflow-work` or `export NXF_WORK=/Volumes/MySSD/nextflow-work` (macOS, SSD must be APFS/HFS+).
5. **Conda profile:** Multiple tasks can contend on the conda cache. Try `./run_pipeline.sh -profile docker`, or set `export NXF_CONDA_CACHEDIR=/tmp/conda-$USER-$$` before running.
6. **Upgrade Nextflow:** Pipeline requires ≥26.04.0; older versions have locking issues.

**"matplotlib is building a font cache" seems stuck:**
- Matplotlib scans system fonts on first import (30s-2min). `umi_tools` (preprocess, coverage) and report/aggregate tasks use it.
- **Docker:** The image pre-builds the cache; pull the latest and rebuild if needed.
- **Conda:** `MPLCONFIGDIR` is set to `$TMPDIR` in affected modules. Ensure `TMPDIR` points to local disk (not NFS).

**Spike-in alignment fails (sample-specific):**
- Samples with many unaligned reads (e.g. 20M+) need more memory for spike-in alignment
- Increase Docker memory (Settings → Resources) or system RAM
- Check `bowtie2_spikein.log` in the failed task's work dir for details

**Finished tasks re-run from sample 1 (even with -resume):**
- Nextflow's cache depends on input file path, size, and timestamp. NFS/network storage can give inconsistent timestamps; add `preprocess_reads_lenient_cache: true` to params.yaml or run with `--preprocess_reads_lenient_cache`.
- Docker `:latest` changes when the image is updated; use a fixed tag (e.g. `tracktx:1.4.0`) for stable caching.
- Debug: `nextflow run ... -resume -dump-hashes 2>&1 | grep "cache hash"` and compare between runs.

**preprocess_and_quality_filter_reads re-runs after stop/restart:**
- When you Ctrl+C, running/queued tasks are cancelled and not cached
- Only completed tasks are reused with `-resume`
- Let the pipeline finish, or stop when no preprocess_and_quality_filter_reads tasks are active

## Getting help

1. **Check logs**: `.nextflow.log` in the working directory
2. **Review trace**: `{output_dir}/trace/report.html` for resource issues
3. **Monitor live**: `python3 nfmon.py` to see what's happening
4. **GitHub Issues**: [Report bugs](https://github.com/serhataktay/tracktx/issues)
