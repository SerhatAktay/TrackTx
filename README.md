# TrackTx — Nascent RNA Analysis Pipeline

<div align="center">

**End-to-end Nextflow pipeline for PRO-seq and nascent RNA-seq data analysis**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-supported-0db7ed.svg)](https://www.docker.com/)
[![Conda](https://img.shields.io/badge/conda-supported-green.svg)](https://conda.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-1.5.0-blue.svg)](https://github.com/serhataktay/tracktx/releases/tag/v1.5.0)

</div>

TrackTx takes raw FASTQ files (or SRA accessions) and produces publication-ready coverage tracks, divergent transcription calls, RNA Polymerase II occupancy metrics, and interactive HTML reports — all in a single run. It is designed for PRO-seq and related nascent RNA-seq protocols and runs on anything from a laptop to an HPC cluster, configured through a point-and-click generator and a single YAML file.

---

## Contents

- [⚡ Quick Start](#-quick-start)
- [🧪 Testing the Pipeline](#-testing-the-pipeline)
- [📊 What TrackTx Does](#-what-tracktx-does)
- [🔬 Pipeline Modules](#-pipeline-modules)
- [🔧 Installation](#-installation)
- [📁 Input Files](#-input-files)
- [📂 Outputs](#-outputs)
- [🎯 Common Use Cases](#-common-use-cases)
- [⚙️ Execution Profiles](#-execution-profiles)
- [⚡ Performance & Storage](#-performance--storage)
- [🛠 Troubleshooting](#-troubleshooting)
- [📖 Resources](#-resources)
- [🧬 Citation](#-citation)
- [📜 License](#-license)

---

## ⚡ Quick Start

### 1️⃣ Generate your configuration (recommended)

Open **`TrackTx_config_generator.html`** in your browser — double-click it, or:

```bash
open TrackTx_config_generator.html   # macOS
```

Fill in your study details, then click **Download ZIP** to get `params.yaml` and `samplesheet.csv`. Place both files in the same directory as `main.nf`.

For local FASTQ files, enter the full path to each file (right-click in your file manager → Copy Path).

### 2️⃣ Run the pipeline

```bash
./run_pipeline.sh
```

The launcher auto-detects Docker, Conda, or a local environment, loads your `params.yaml` and `samplesheet.csv`, and starts the pipeline. Common options:

```bash
./run_pipeline.sh --resume               # Resume a previous run
./run_pipeline.sh -profile docker        # Force a specific profile
./run_pipeline.sh --output_dir my_run    # Custom output directory
./run_pipeline.sh --external-drive       # Running from USB/NAS (see Performance)
./run_pipeline.sh --help                 # All options
```

> **All files default to the project directory.** Use `--external-drive` only when the project lives on exFAT/USB — cache, temp, and work (~10–50 GB) are then redirected to local storage while results stay on the project drive.

### 3️⃣ Monitor progress

```bash
python3 nfmon.py --from-start --tail 80
```

The live monitor shows task progress (drawn from the Nextflow trace), per-task runtimes, CPU/memory usage, and tailed log output from running modules. Install `rich` for an enhanced layout:

```bash
pip install rich          # or: conda install -c conda-forge rich
```

Useful options:

```bash
python3 nfmon.py --filter "align"        # Watch specific tasks
python3 nfmon.py --oneshot               # Single snapshot and exit
python3 nfmon.py --oneshot --json s.json # Export snapshot as JSON
python3 nfmon.py --trace path/trace.txt  # Point to a specific trace file
python3 nfmon.py --help                  # All options
```

---

## 🧪 Testing the Pipeline

A ready-made test setup is included. It downloads a small public PRO-seq dataset (~10% of reads, ~100–200 MB) so you can verify the pipeline works before running your own data.

**Step 1: Download the test data**

```bash
# Standard (requires curl)
./scripts/download_and_subset_test_data.sh

# With Docker (no local tools needed)
./scripts/download_and_subset_test_data.sh --docker
```

The script fetches paired-end data from ENA, subsets to 10%, and writes to `test_PE/test_data/`.

**Step 2: Run the pipeline on the test data**

```bash
./run_pipeline.sh \
  --samplesheet test_PE/samplesheet_PE.csv \
  --params-file  test_PE/params_PE.yaml \
  --output_dir   ./results_test_PE
```

See `test_PE/README.md` for dataset details and expected outputs.

---

## 📊 What TrackTx Does

TrackTx analyses nascent RNA sequencing data to characterise **real-time transcription**: where RNA Polymerase is engaged across the genome, how it pauses at gene promoters, and where divergent (upstream antisense) transcription is occurring.

```mermaid
graph LR
    A[📁 Raw FASTQ / SRA] --> B[✂️ QC & Trimming]
    B --> C[🧭 Alignment]
    C --> D[📈 Coverage Tracks]
    D --> E[⚖️ Normalization]
    E --> F[🔬 Divergent Detection]
    F --> G[📊 Pol-II Metrics]
    G --> H[📋 Reports]
```

**Key capabilities:**

- **Automated** — raw reads to publication-ready outputs in one command
- **Flexible** — handles SE and PE data, UMIs, barcodes, and spike-in normalization
- **Statistical** — divergent transcription detection with Gaussian Mixture Models and FDR control, no manual thresholding required
- **Quantitative** — CPM and spike-in CPM (siCPM) normalization for cross-sample comparisons
- **Comprehensive** — per-sample HTML reports plus a cohort-wide dashboard
- **Scalable** — runs identically on a laptop, workstation, or HPC cluster

---

## 🔬 Pipeline Modules

TrackTx runs 15 modules in sequence. Here is what each one does.

**01 — Download Genome Annotations**
Downloads GTF gene annotation files from Ensembl, RefSeq, or GENCODE for the chosen reference genome. Annotations drive gene-body boundary calls, divergent transcription pairing, and functional region assignment throughout the pipeline.

**02 — Download SRA Samples**
Retrieves FASTQ files from NCBI SRA or the European Nucleotide Archive (ENA) when `sample_source: srr` is set. Uses `fasterq-dump` with multi-threaded conversion and falls back to parallel ENA downloads automatically. Raw FASTQs are stored via Nextflow's `storeDir` at `{output_dir}/00_sra_cache/{SRR}/`, so re-downloads are skipped as long as those files exist — even after the `work/` directory has been deleted between runs.

**03 — Preprocess and Quality Filter Reads**
Trims adapter sequences, removes barcodes, and applies a minimum length filter in a single `cutadapt` pass. UMI extraction with `umi_tools` follows if enabled. FastQC runs on both raw and cleaned reads to confirm quality before alignment.

**04 — Download Genome and Build Alignment Index**
Downloads the reference (and spike-in) genome FASTA, then builds a `bowtie2` alignment index. Pre-built indexes are cached and reused across runs.

**05 — Align Reads to Genome**
Aligns cleaned reads to the reference genome with `bowtie2` (end-to-end mode). Spike-in reads are aligned separately to derive normalization factors. Outputs include coordinate-sorted, indexed BAM files for primary, all-mapping, and spike-in alignments. BAMs and alignment artefacts are stored via `storeDir` at `{output_dir}/02_alignments/{sample_id}/`, so alignment is skipped for any sample whose BAM already exists there — even across runs where `work/` has been cleared. To force realignment for a sample, delete its folder: `rm -rf {output_dir}/02_alignments/{sample_id}`.

**06 — Generate Coverage Tracks**
Produces strand-specific 3′-end and 5′-end coverage tracks in bedGraph format using `bedtools genomecov`. Four independent coverage jobs (3p positive, 3p negative, 5p positive, 5p negative) run in parallel. Each bedGraph is pre-sorted inline and converted to BigWig format for use in genome browsers.

**07 — Quantify Reads Per Gene**
Counts uniquely aligned reads overlapping each gene using `samtools idxstats` and gene coordinate intervals. Produces a counts master file that drives CPM and siCPM normalization factor calculation in the next module.

**08 — Normalize Coverage Tracks**
Scales raw bedGraph signal to CPM (counts per million mapped reads) and siCPM (spike-in CPM) using pre-computed factors. Positive and negative strand tracks are normalized in parallel. Outputs both bedGraph and BigWig formats for all track sets (3p, 5p, main, and allMap).

**09 — Detect Divergent Transcription**
The statistical core of the pipeline. Pairs upstream antisense peaks with downstream sense peaks, computes a suite of features (signal balance, local enrichment, strand specificity), and fits a two-component Gaussian Mixture Model to separate signal from noise. Divergent regions passing the FDR threshold are written as a BED file with confidence scores. No manual thresholds are required — set `divergent_threshold: auto` and the calibration percentile handles it.

**10 — Assign Signal to Functional Regions**
Assigns normalized coverage to a hierarchical set of genomic functional regions: active promoters, gene bodies, cleavage and polyadenylation sites, enhancers, termination windows, and non-localized signal. Each position is assigned to exactly one region by sequential masking, so the categories are mutually exclusive.

**11 — Calculate Polymerase Occupancy Metrics**
Computes two complementary views of Pol II occupancy. The density metrics approach sums normalized bedGraph signal over each functional region. The gene metrics approach operates on the filtered BAM directly, computing per-gene TSS-window and gene-body coverage from which pausing indices (PI = TSS density / body density) are derived. Both approaches run in parallel so neither waits on the other.

**12 — Summarize Polymerase Metrics**
Aggregates per-sample Pol II metrics across the cohort into summary TSVs — pausing index distributions, region density tables, and normalization factor comparisons — for use in the cohort report.

**13 — Quality Control Aligned Reads**
Calculates per-sample alignment QC: total and mapped read counts, duplicate rates, MAPQ pass rates, strand balance (critical for PRO-seq validation), fragment length distribution (PE only), and mean genome coverage depth. Results feed the per-sample HTML reports and cohort outlier detection.

**14 — Generate Per-Sample Reports**
Produces an interactive HTML report for each sample, summarising QC metrics, coverage distributions, divergent transcription statistics, and Pol II pausing results, with inline visualizations.

**15 — Combine Reports into Cohort**
Merges all per-sample outputs into a single cohort-level HTML dashboard. Includes by-condition QC comparisons, divergent transcription patterns, Pol II pausing distributions, functional region composition, normalization factor validation, replicate consistency (coefficient of variation), and an interactive sample metrics table.

---

## 🔧 Installation

### Prerequisites

You need **Nextflow** plus one container or environment manager:

| Tool | Purpose | Install |
|------|---------|---------|
| **Nextflow ≥ 24.04.0** | Runs the pipeline | See below |
| **Docker Desktop** | Easiest — all tools pre-packaged | [docker.com](https://www.docker.com/products/docker-desktop/) |
| **Miniconda** | Alternative when Docker isn't available | [docs.conda.io](https://docs.conda.io/en/latest/miniconda.html) |

**Install Nextflow:**

```bash
# Via Conda (recommended if you already use Conda)
conda install -c bioconda nextflow

# Standalone (works without Conda)
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
```

**Verify:**

```bash
nextflow -version   # must be ≥ 24.04.0
docker --version    # if using Docker
conda --version     # if using Conda
```

---

### Option 1: Docker (Recommended)

Docker packages all required bioinformatics tools — nothing else to install.

```bash
git clone https://github.com/serhataktay/tracktx.git
cd tracktx
./run_pipeline.sh
```

The first run downloads the container image (~2–5 min). After a `git pull`, the launcher automatically pulls the updated image. To skip the pull (e.g. offline): `TRACKTX_SKIP_PULL=1 ./run_pipeline.sh`.

---

### Option 2: Conda

Use Conda when Docker is unavailable (e.g. restricted HPC environments).

```bash
git clone https://github.com/serhataktay/tracktx.git
cd tracktx
./run_pipeline.sh
```

The launcher auto-detects Conda and creates the pipeline environment on first run (~10–20 min).

---

### Option 3: Manual / HPC

If you already have the required tools in your PATH:

```bash
git clone https://github.com/serhataktay/tracktx.git
cd tracktx
nextflow run main.nf -entry TrackTx -profile local \
  --samplesheet samplesheet.csv -params-file params.yaml
```

---

### Windows (WSL2)

On Windows, run TrackTx inside **WSL2** (Windows Subsystem for Linux) with Ubuntu.

**Step 1: Install WSL with Ubuntu**

```powershell
wsl --install -d Ubuntu
```

Restart if prompted, then complete the Ubuntu first-boot setup (username and password).

**Step 2: Install Nextflow inside WSL**

```bash
sudo apt update && sudo apt install -y openjdk-17-jdk curl
curl -s https://get.nextflow.io | bash
chmod +x nextflow && sudo mv nextflow /usr/local/bin/
```

**Step 3: Install Docker Desktop**

Download from [docker.com](https://www.docker.com/products/docker-desktop/), enable **"Use the WSL 2 based engine"** during setup, and confirm Docker shows "Running".

**Step 4: Clone and run**

```bash
cd /mnt/c/Users/YourUsername   # or use ~ for the Linux home directory
git clone https://github.com/serhataktay/tracktx.git
cd tracktx
./run_pipeline.sh
```

> **Tip:** Store your data on the Linux filesystem (`~` or `/home/you`) rather than `/mnt/c` for significantly better I/O performance.

**Verify the setup:**

```bash
nextflow -version && docker --version && docker run --rm hello-world
```

---

### System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| **OS** | Linux, macOS, Windows (WSL2) | Linux or macOS |
| **CPU** | 2+ cores | 8+ cores |
| **RAM** | 8 GB | 32+ GB |
| **Storage** | 50 GB | 200+ GB SSD |

The first run downloads reference genomes (1–5 GB depending on species). Ensure sufficient free disk space before starting.

---

## 📁 Input Files

### Sample Sheet (`samplesheet.csv`)

Use the [config generator](TrackTx_config_generator.html) to create this, or format manually:

```csv
sample,condition,timepoint,replicate,file1,file2
ctrl_rep1,control,0,1,data/ctrl_R1.fastq,data/ctrl_R2.fastq
heat_rep1,treatment,30,1,data/heat_R1.fastq,data/heat_R2.fastq
```

- **Single-end data:** use `file1` only; leave `file2` empty.
- **SRA accessions:** put the SRR ID in `file1` and set `sample_source: srr` in `params.yaml`.
- **Local files:** use absolute paths or paths relative to the project directory.

### Parameters (`params.yaml`)

Generate this with the interactive [`TrackTx_config_generator.html`](TrackTx_config_generator.html) — it validates all parameters and explains each option inline. To write one manually:

```yaml
# ── Core settings ────────────────────────────────────────────────────────────
samplesheet:        "samplesheet.csv"
sample_source:      "local"          # "local" or "srr"
reference_genome:   "hg38"           # hg38, hs1, mm39, mm10, dm6, …
paired_end:         false
output_dir:         "./results"

# ── Spike-in normalization (optional) ────────────────────────────────────────
spikein_genome:     "dm6"

# ── Annotation ────────────────────────────────────────────────────────────────
annotation_source:          refseq   # refseq | ensembl | gencode
annotation_exclude_biotypes: ""      # e.g. "rRNA,tRNA"
annotation_chr_naming:      none     # none | add | remove

# ── Adapter trimming ──────────────────────────────────────────────────────────
adapter_trimming:
  enabled: true
  preset:  "illumina"
  adapter1: "TGGAATTCTCGGGTGCCAAGG"

# ── UMI processing (if your library has UMIs) ─────────────────────────────────
umi:
  enabled:  false
  length:   8
  location: 5                        # 5 or 3

# ── Divergent transcription detection ────────────────────────────────────────
advanced:
  divergent_threshold:               auto   # auto or float
  divergent_sum_thr:                 auto
  divergent_fdr:                     0.08
  divergent_calibration_percentile:  65
  divergent_merge_gap:               150
  divergent_nt_window:               1000
  divergent_balance:                 0.0
  divergent_qc:                      true

# ── Functional region windows ─────────────────────────────────────────────────
functional_regions:
  tss_active_pm: 500
  prom_up:       500
  prom_down:     250
```

**Disk space tips:**

| Flag | Savings |
|------|---------|
| `publish_trimmed_fastq: false` | ~920 MB per sample |
| `publish_references_gtf: false` | ~790 MB (one-time) |
| `output.raw_tracks: false` | ~150 MB per sample |
| `output.bedgraph: false` | ~100 MB per sample |
| `norm.emit_allmap: false` | ~50 MB per sample |

> **Note:** `publish_alignments` no longer suppresses BAM output. Module 05 uses `storeDir`, so BAMs are always written to `{output_dir}/02_alignments/` and serve as the persistent cache. Delete individual sample folders manually if you need to reclaim space.

---

## 📂 Outputs

The pipeline writes results to `{output_dir}/` (default: `./results`).

```
results/
├── 00_sra_cache/                  # storeDir — raw FASTQs (SRA runs only)
│   └── <SRR>/
│       ├── <SRR>_R1.fastq[.gz]   # Persisted to survive work/ deletion
│       ├── <SRR>_R2.fastq[.gz]   # (PE only)
│       └── README_fastq.txt
│
├── 02_alignments/                 # storeDir — BAMs and alignment artefacts
│   └── <sample>/
│       ├── <sample>.bam / .bai
│       ├── <sample>_allMap.bam
│       ├── <sample>_spikein.bam / .bai
│       ├── aligner_summary.tsv
│       └── bowtie2_*.log
│
├── 05_normalized_tracks/          # ⭐ Load in IGV or UCSC Genome Browser
│   └── <sample>/
│       ├── 3p/*.cpm.bw            # CPM BigWig tracks (3' end)
│       ├── 3p/*.sicpm.bw          # Spike-in CPM BigWig tracks
│       └── 5p/*.cpm.bw            # 5' end tracks (paired-end)
│
├── 06_divergent_tx/               # Divergent transcription regions
│   └── <sample>/
│       ├── divergent_transcription.bed     # High-confidence regions (BED5)
│       └── divergent_transcription_qc.txt  # Statistical QC report
│
├── 07_functional_regions/         # Genomic region annotations
│   └── <sample>/
│       ├── functional_regions.bed
│       └── functional_regions_summary.tsv
│
├── 08_pol_metrics/                # Pol II pausing and density
│   └── <sample>/
│       ├── pausing_index.tsv
│       ├── pol_gene_metrics.tsv
│       └── pol_density.tsv
│
├── 10_qc/                         # Per-sample alignment QC
│   └── <sample>/qc_pol.json
│
└── 11_reports/                    # ⭐ Interactive HTML reports
    ├── cohort/
    │   └── global_summary.html    # Cohort dashboard (start here)
    └── samples/
        └── <sample>/<sample>.report.html
```

Intermediate outputs (genome references, trimmed FASTQs, raw alignments, raw tracks, read counts, and aggregated metrics) are also produced in their numbered folders. Nextflow trace files live in `{output_dir}/trace/`.

**Where to start:**

1. **`11_reports/cohort/global_summary.html`** — the cohort dashboard. Covers QC assessment and outlier detection, divergent transcription patterns by condition, Pol II pausing index distributions, functional region composition, normalization factor validation, replicate consistency, and an interactive sample table.

2. **`11_reports/samples/<sample>/<sample>.report.html`** — detailed per-sample reports with inline visualizations.

3. **`05_normalized_tracks/<sample>/3p/*.cpm.bw`** — load directly into IGV or the UCSC Genome Browser for visual inspection.

4. **`{output_dir}/trace/report.html`** — Nextflow's pipeline performance report: task runtimes, CPU/memory usage, and resource efficiency.

---

## 🎯 Common Use Cases

### Time-course experiment

```bash
cat > timecourse.csv << EOF
sample,condition,timepoint,replicate,file1,file2
t0_r1,control,0,1,data/t0_r1.fastq,
t0_r2,control,0,2,data/t0_r2.fastq,
t30_r1,heat,30,1,data/t30_r1.fastq,
t30_r2,heat,30,2,data/t30_r2.fastq,
EOF

./run_pipeline.sh --samplesheet timecourse.csv
```

### Drug treatment with spike-in normalization

```yaml
# params.yaml
reference_genome: "mm39"
spikein_genome:   "dm6"
paired_end:       false
```

```bash
./run_pipeline.sh --params-file params.yaml
```

### Downloading from SRA

```csv
sample,condition,timepoint,replicate,file1,file2
ctrl,control,0,1,SRR4454567,
treated,treatment,30,1,SRR4454568,
```

```yaml
# params.yaml
sample_source:    srr
reference_genome: "hg38"
output_dir:       "./results"
```

```bash
./run_pipeline.sh --samplesheet sra_samples.csv --params-file params.yaml
```

---

## ⚙️ Execution Profiles

The launcher **auto-detects** your environment. Force a specific profile with `-profile`:

| Profile | When to use |
|---------|-------------|
| `docker` | Docker is available (recommended) |
| `conda` | Docker not available |
| `conda_server` | Conda on NFS / network storage |
| `singularity` | HPC with Singularity |
| `slurm` | SLURM scheduler (combine with a container profile) |
| `performance` | External drive — auto-appended by `--external-drive` |
| `local` | Tools already in PATH |

```bash
./run_pipeline.sh                             # Auto-detect
./run_pipeline.sh -profile docker             # Force Docker
./run_pipeline.sh -profile slurm,singularity  # HPC + containers
./run_pipeline.sh --external-drive            # USB / NAS
```

---

## ⚡ Performance & Storage

### Running from external storage (USB / NAS / exFAT)

If your project lives on an external drive or network storage, use:

```bash
./run_pipeline.sh --external-drive
```

This automatically:

- Keeps **results** on your project directory
- Moves cache, temp, and work (~10–50 GB) to local disk (`~/tmp/tracktx_*`)
- Fixes publish errors (`copy` instead of hard links — required on exFAT)
- Fixes `OverlappingFileLockException` (exFAT has no file-locking support)
- Increases task parallelism for better I/O utilisation

Ensure ~20–50 GB free on your internal drive before using this mode.

### Storage footprint

Typical sizes for a single paired-end sample (full run):

| Location | Size | Contents |
|----------|------|----------|
| `results/` | ~2.3 GB | Published outputs + storeDir folders |
| `results/00_sra_cache/` | ~100 GB+ | Raw SRA FASTQs (SRA runs only, storeDir) |
| `results/02_alignments/` | ~270 MB/sample | BAMs (storeDir) |
| `work/` | ~1.7 GB | Nextflow intermediates (safe to delete) |
| `.cache/genomes` | ~7 GB | Reference genomes (shared across runs) |
| Input data | ~325 MB | Raw FASTQ (test subset) |

Largest result folders: `01_trimmed_fastq` (~920 MB), `00_references` (~800 MB GTF), `02_alignments` (~270 MB per sample), `05_normalized_tracks` (~165 MB).

The `work/` directory can be safely deleted between runs — `storeDir` folders (`00_sra_cache/` and `02_alignments/`) persist independently and will be used as cache on the next run.

### Expected runtimes

| Samples | Internal SSD (optimised) | External SSD |
|---------|--------------------------|--------------|
| 2 | 30–60 min | 2–3 hours ⚠️ |
| 4 | 1–2 hours | 4–6 hours ⚠️ |
| 8 | 2–4 hours | 8–12 hours ⚠️ |

### Memory tuning

If you hit out-of-memory errors (exit code 137), the pipeline may be allocating based on host RAM while Docker has a lower memory limit. Fix:

```bash
export NXF_HOST_MEM=8   # match your Docker or WSL memory limit (GB)
export NXF_HOST_CPUS=4  # reduce parallelism if needed
./run_pipeline.sh
```

In Docker Desktop, check Settings → Resources → Memory. In WSL2, set a limit in `.wslconfig`.

### Clean up after a run

```bash
nextflow clean -f   # removes work/ directories for completed runs
```

> **`nextflow clean` does not touch `storeDir` folders.** Raw SRA FASTQs (`00_sra_cache/`) and BAMs (`02_alignments/`) will remain on disk and continue to act as cache for future runs. Delete them manually only when you are certain they are no longer needed:
> ```bash
> rm -rf results/00_sra_cache/         # remove all cached SRA FASTQs
> rm -rf results/02_alignments/SAMPLE  # remove one sample's BAMs to force realignment
> ```

---

## 🛠 Troubleshooting

### Reading error messages

When a module fails, look for the **TRACKTX ERROR block** in the output:

```
═══════════════════════════════════════════════════════════════════════
TRACKTX ERROR
═══════════════════════════════════════════════════════════════════════
Module:  detect_divergent_transcription
Problem: Missing Python dependencies (numpy, pandas, scikit-learn, scipy)
Fix:     pip install numpy pandas scikit-learn scipy
         Or use: -profile conda | -profile docker
═══════════════════════════════════════════════════════════════════════
```

- **Quick find:** `grep -A 6 "TRACKTX ERROR"` in the Nextflow output
- **Full log:** check the `.log` file in the work directory shown at the bottom of the error
- **Resume:** add `--resume` once you have fixed the issue

---

### Environment issues

**Docker not running:**
Start Docker Desktop and wait until it shows "Running" before launching the pipeline.

**Conda environment fails:**
```bash
# Use Docker (more reliable across platforms)
./run_pipeline.sh -profile docker

# Or clean the Conda cache and retry
conda clean --all --yes
./run_pipeline.sh -profile conda
```

**"Missing Python dependencies" (divergent transcription step):**
```bash
./run_pipeline.sh -profile docker   # recommended — all deps pre-installed
# Or: ./run_pipeline.sh -profile conda
# Or: pip install -r envs/requirements-divergent.txt
```

**"matplotlib is building a font cache" appears stuck:**
Matplotlib scans system fonts on first import (30 s–2 min). Docker pre-builds the cache. Under Conda, ensure `TMPDIR` points to local disk (not NFS), as `MPLCONFIGDIR` is set to `$TMPDIR` in affected modules.

---

### File-locking and publish errors

**"Failed to publish file [link]" errors:**
Hard links do not work on exFAT (common on USB drives).
```bash
./run_pipeline.sh --external-drive
```

**`OverlappingFileLockException`** (`preprocess_and_quality_filter_reads`, `download_genome_annotations`):

1. **Two runs from the same directory:** only one Nextflow run per directory at a time. Stop other runs first.
2. **Stale lock from a killed run:** `rm -rf work .nextflow && ./run_pipeline.sh`
3. **USB with exFAT/FAT32:** use `--external-drive` (moves work to local disk).
4. **NFS / cloud-synced storage (iCloud, Dropbox):** file locking is unreliable. Set the work directory to internal disk: `export NXF_WORK=/tmp/nextflow-work`
5. **Conda profile contention:** try Docker, or `export NXF_CONDA_CACHEDIR=/tmp/conda-$USER-$$`
6. **Old Nextflow version:** pipeline requires ≥ 24.04.0.

---

### Memory and performance

**Out of memory (exit 137):**
OOM means RAM exhaustion, not disk space. Common causes:

- Docker memory limit is lower than host RAM — adjust in Settings → Resources → Memory, then `export NXF_HOST_MEM=<limit_in_gb>`
- WSL2 reports host RAM, not its own limit — set `NXF_HOST_MEM` to match your `.wslconfig` limit
- Multiple parallel tasks exceed total available RAM — reduce with `export NXF_HOST_CPUS=4`

**Pipeline running slowly:**
- First run downloads reference genomes (10–30 min) — this is normal
- Use SSD storage; external/network drives can add 3–4× overhead
- Monitor with `python3 nfmon.py` to identify bottlenecks
- Use `--external-drive` if your project is on USB/NAS

**Spike-in alignment fails (sample-specific):**
Samples with many unaligned reads (20 M+) need more memory. Increase Docker memory or check `bowtie2_spikein.log` in the failed task's work directory.

---

### Quality and biological issues

**Low MAPQ pass rate (< 30%):**
PRO-seq data typically shows 30–60% MAPQ ≥ 10. Lower rates are expected for datasets with repetitive regions or low-complexity genomes. To lower the threshold:
```yaml
qc:  { mapq: 5 }
pol: { mapq: 5 }
```

**Strand bias looks wrong:**
PRO-seq should be ~45–55% per strand. A severe imbalance (> 70% one strand) usually indicates a library preparation issue or incorrect `paired_end` setting. Check the strand bias table in the per-sample HTML report.

---

### Caching and resume issues

**Understanding the two caching mechanisms:**
TrackTx uses two complementary caching strategies:

- **`storeDir`** (modules 02, 04, 05 and genome annotation/index steps) — caches by **file existence**. If the output files are present at the store path, the process is skipped entirely, regardless of whether `work/` exists. Survives `nextflow clean` and deleting `work/`.
- **`-resume`** (all other modules) — caches by **input checksum**. Requires `work/` to be intact. Skips tasks whose inputs have not changed since the last run.

**`-resume` re-runs tasks that already completed:**
Nextflow's resume cache depends on input file path, size, and timestamp. NFS/network storage can give inconsistent timestamps. Fix:
```yaml
preprocess_reads_lenient_cache: true
```
or run with `--preprocess_reads_lenient_cache`. Debug with:
```bash
nextflow run ... -resume -dump-hashes 2>&1 | grep "cache hash"
```

**SRA download or alignment re-runs even though outputs exist:**
These modules use `storeDir`. If a rerun triggers them unexpectedly, check that the expected output files are actually present and non-empty at:
- `{output_dir}/00_sra_cache/{SRR}/{SRR}_R1.fastq`
- `{output_dir}/02_alignments/{sample_id}/{sample_id}.bam`

If outputs are missing or zero-size, delete the folder and rerun.

**Forcing a re-download or realignment:**
Delete the relevant `storeDir` subfolder — `-resume` alone is not enough, as `storeDir` bypasses the checksum cache:
```bash
rm -rf results/00_sra_cache/SRR123456      # re-download one accession
rm -rf results/02_alignments/my_sample     # realign one sample
```

**`preprocess_and_quality_filter_reads` re-runs after a Ctrl+C:**
Only completed tasks are cached. Tasks that were running or queued when you interrupted are not. Let the pipeline finish naturally, or stop it between preprocessing and alignment.

---

### Getting help

1. **Check logs:** `.nextflow.log` in the working directory
2. **Review trace:** `{output_dir}/trace/report.html` for resource usage and task timings
3. **Monitor live:** `python3 nfmon.py` to see which task is running and what it is doing
4. **GitHub Issues:** [github.com/serhataktay/tracktx/issues](https://github.com/serhataktay/tracktx/issues)

---

## 📖 Resources

| Resource | Description |
|----------|-------------|
| [TrackTx_config_generator.html](TrackTx_config_generator.html) | Interactive parameter and samplesheet generator |
| [GitHub Issues](https://github.com/serhataktay/tracktx/issues) | Bug reports and feature requests |
| [Nextflow documentation](https://www.nextflow.io/docs/latest/) | Nextflow executor, profile, and resume reference |

---

## 🧬 Citation

If TrackTx is useful for your research, please cite:
[https://github.com/serhataktay/tracktx](https://github.com/serhataktay/tracktx)

### Key references

- **PRO-seq:** Kwak et al. (2013). *Precise maps of RNA polymerase reveal how promoters direct initiation and pausing.* Science. [doi:10.1126/science.1229386](https://doi.org/10.1126/science.1229386)
- **Nascent RNA-seq:** Core et al. (2008). *Nascent RNA sequencing reveals widespread pausing and divergent initiation at human promoters.* Science. [doi:10.1126/science.1162228](https://doi.org/10.1126/science.1162228)
- **Nextflow:** Di Tommaso et al. (2017). *Nextflow enables reproducible computational workflows.* Nature Biotechnology. [doi:10.1038/nbt.3820](https://doi.org/10.1038/nbt.3820)

---

## 📜 License

TrackTx is released under the [MIT License](LICENSE).

---

<div align="center">

**⭐ Star this repo if TrackTx is useful for your research!**

[Issues](https://github.com/serhataktay/tracktx/issues) · [Releases](https://github.com/serhataktay/tracktx/releases)

</div>
