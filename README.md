# TrackTx — Nascent RNA Analysis Pipeline

<div align="center">

**A powerful Nextflow pipeline for PRO-seq and nascent RNA-seq analysis**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.0-23aa62.svg)](https://www.nextflow.io/)
[![Docker](https://img.shields.io/badge/docker-supported-0db7ed.svg)](https://www.docker.com/)
[![Conda](https://img.shields.io/badge/conda-supported-green.svg)](https://conda.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-3.0-blue.svg)](https://github.com/serhataktay/tracktx/releases)

</div>

---

## 📚 Contents

- [🎉 What's New in v3.0](#-whats-new-in-v30)
- [⚡ Quick Start](#-quick-start)
- [🧪 Testing the Pipeline](#-testing-the-pipeline)
- [📊 What Does TrackTx Do?](#-what-does-tracktx-do)
- [🔧 Installation](#-installation)
- [📁 Input Files](#-input-files)
- [📊 Outputs](#-outputs)
- [🎯 Common Use Cases](#-common-use-cases)
- [⚙️ Execution Profiles](#-execution-profiles)
- [⚡ Performance Optimization](#-performance-optimization)
- [🔧 Troubleshooting](#-troubleshooting)
- [📖 Documentation](#-documentation)
- [🧬 Citation](#-citation)
- [📜 License](#-license)

---

## 🎉 What's New in v3.0

### 🔬 Statistical Divergent Transcription Detection
- **Gaussian Mixture Models (GMM)** with FDR control for high-confidence region calling
- **Auto-calibration**: No manual threshold tuning needed - set to `auto` and go!
- **Beta-Binomial models** for robust background estimation
- Replaces threshold-based detection with statistically-driven approach
- Comprehensive QC reports with score distributions

### 📊 Comprehensive Cohort Reporting
- **8 detailed analysis sections** with extensive explanations:
  - Quality Control Analysis with outlier detection
  - Divergent Transcription patterns across conditions
  - Pol II Pausing Index distributions and comparisons
  - Functional Region composition analysis
  - Normalization factor validation
  - Sample-level interactive metrics table
- **Beautiful modern UI** with smooth navigation and responsive design
- **Interactive visualizations** with histograms, distributions, and comparisons
- **By-condition analysis** for all major metrics
- **Replicate consistency checks** (coefficient of variation)

### 🎨 Enhanced Configuration Generator
- **Extensive explanatory text** for every setting
- **More spacious layout** for better readability
- **Logical organization** with related settings grouped together
- **Bold section headers** for easier navigation
- **Inline help** explaining when/why to use each option
- Updated for the latest parameters and best practices

### 🐛 Pipeline Improvements
- Fixed caching issues with divergent detection parameters
- Better handling of UMI and barcode processing
- Improved parameter validation and error messages
- Enhanced documentation throughout

---

## ⚡ Quick Start

Get started in **3 simple steps**:

### 1️⃣ Generate Configuration (Interactive, Recommended)

Open the **TrackTx configuration generator** in your browser:

```bash
open TrackTx_config_generator.html   # macOS
# Or double-click the file on any system
```

- Fill in your study details
- For local FASTQ files: enter the **full path** to each file (copy from your file manager → right-click → Copy path)
- Click **Download ZIP** to save `params.yaml` and `samplesheet.csv`
- **Put both files in the same folder as the pipeline** (the tracktx directory where `main.nf` lives)

### 2️⃣ Run Pipeline (Automatic)

Use the **smart launcher** that auto-detects your environment:

```bash
./run_pipeline.sh
```

That's it! The script will:
- ✅ Auto-detect Docker, Conda, or local environment
- ✅ Load your `params.yaml` and `samplesheet.csv`
- ✅ Optimize resource allocation for your system
- ✅ Start the pipeline

**Advanced options:**
```bash
./run_pipeline.sh --help                 # See all options
./run_pipeline.sh                       # All files in project dir
./run_pipeline.sh -profile docker        # Force specific profile
./run_pipeline.sh --resume               # Resume previous run
./run_pipeline.sh --output_dir my_run    # Custom output directory
```

**💡 All files in project dir by default.** Use `--external-drive` only when project is on exFAT/USB — then cache, temp, and work (~10–50 GB) go to local; results stay on project.

### 3️⃣ Monitor Progress (Real-time)

Watch your pipeline in action with the **live monitor**:

```bash
python3 nfmon.py
```

**Install Rich (optional, for enhanced UI):** The monitor works without it (basic curses fallback), but for a better layout, colors, and live updates:
```bash
pip install rich
# or with conda:
conda install -c conda-forge rich
```

Features:
- 📊 Real-time task progress and resource usage
- 🎯 Live log tailing for active tasks
- 💻 CPU/memory/load monitoring
- ⚡ Per-task performance metrics

**Monitor options:**
```bash
python3 nfmon.py --help                          # All options
python3 nfmon.py --filter "alignment"            # Watch specific tasks
python3 nfmon.py --all-logs                      # See all task logs
python3 nfmon.py --oneshot                       # Quick snapshot
python3 nfmon.py --oneshot --json status.json    # Export JSON
```

---

## 🧪 Testing the Pipeline

Want to verify the pipeline works before running your own data? Use the bundled test setup with readymade samplesheets, params, and a script that downloads small test datasets.

### Step 1: Download test data

Run the download script to fetch and subset public PRO-seq data (~10% of reads, ~100–200 MB per test):

```bash
# Single-end (1 sample, ~2 min)
./scripts/download_and_subset_test_data.sh SE

# Paired-end (1 sample, ~3 min)
./scripts/download_and_subset_test_data.sh PE

# Both SE and PE
./scripts/download_and_subset_test_data.sh all
```

**With Docker** (same image as the pipeline; no local curl/gzip needed):

```bash
./scripts/download_and_subset_test_data.sh --docker SE
./scripts/download_and_subset_test_data.sh --docker PE
./scripts/download_and_subset_test_data.sh --docker all
```

The script downloads from ENA, subsets to 10%, and removes the full files. Outputs go to `test_SE/test_data/` and `test_PE/test_data/`.

### Step 2: Run the pipeline with the readymade configs

**Single-end test:**
```bash
./run_pipeline.sh \
  --samplesheet test_SE/samplesheet_SE.csv \
  --params-file test_SE/params_SE.yaml \
  --output_dir ./results_test_SE
```

**Paired-end test:**
```bash
./run_pipeline.sh \
  --samplesheet test_PE/samplesheet_PE.csv \
  --params-file test_PE/params_PE.yaml \
  --output_dir ./results_test_PE
```

Both configs use `sample_source: "local"` and point to the subset FASTQs. See `test_SE/README.md` and `test_PE/README.md` for dataset details.

---

## 📊 What Does TrackTx Do?

TrackTx analyzes nascent RNA sequencing data (PRO-seq, GRO-seq, etc.) to understand **real-time transcription**:

```mermaid
graph LR
    A[📁 FASTQ Files] --> B[✂️ QC & Trimming]
    B --> C[🧭 Alignment]
    C --> D[📈 Track Generation]
    D --> E[⚖️ Normalization]
    E --> F[🔬 Statistical Detection]
    F --> G[📊 Comprehensive Reports]
```

**Key Features:**
- 🎯 **Automated**: From raw reads to publication-ready figures
- 🚀 **Fast**: Optimized for any system (laptop → HPC)
- 📊 **Comprehensive**: Alignment, tracks, divergent transcription, Pol-II metrics
- ⚖️ **Quantitative**: Spike-in normalization for cross-sample comparisons
- 🎨 **Beautiful**: Interactive HTML reports with extensive visualizations
- 🔬 **Statistical**: GMM-based divergent detection with FDR control
- 📈 **Insightful**: Cohort-wide analysis with by-condition comparisons
- 💡 **Guided**: Extensive explanations and quality assessment throughout

---

## 🔧 Installation

### Prerequisites

You need **Nextflow** (the workflow engine) plus **one** of Docker or Conda (for the tools):

| Requirement | Purpose | Install |
|-------------|---------|---------|
| **Nextflow** (≥24.04.0) | Runs the pipeline | See below |
| **Docker Desktop** | Easiest—packages all tools | [Get Docker](https://docs.docker.com/get-docker/) |
| **Miniconda** | Alternative if Docker unavailable | [Get Miniconda](https://docs.conda.io/en/latest/miniconda.html) |

**Install Nextflow** (choose one):

```bash
# Option A: Conda (recommended if you use Conda)
conda install -c bioconda nextflow

# Option B: Standalone (works without Conda)
curl -s https://get.nextflow.io | bash
# Moves nextflow to your PATH, e.g.:
sudo mv nextflow /usr/local/bin/   # Linux/macOS
```

**Verify:**
```bash
nextflow -version   # Must show 24.04.0 or higher
docker --version   # If using Docker
conda --version    # If using Conda
```

---

### Windows (WSL)

On Windows, use **WSL2** (Windows Subsystem for Linux) with Ubuntu. This gives you a Linux environment where the pipeline runs natively.

**Step 1: Install WSL with Ubuntu**
```powershell
wsl --install -d Ubuntu
```
Restart if prompted. After reboot, Ubuntu will open; complete the initial setup (username, password).

**Step 2: Install dependencies and Nextflow** (run inside WSL/Ubuntu)
```bash
sudo apt update
sudo apt install -y openjdk-17-jdk curl

cd ~
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/
```

**Step 3: Install Docker Desktop**
- Download from [docker.com](https://www.docker.com/products/docker-desktop/)
- During setup, enable **“Use the WSL 2 based engine”**
- Start Docker Desktop and ensure it shows “Running”

**Step 4: Clone and run**
```bash
cd /mnt/c/Users/YourUsername   # Replace with your Windows username; or use ~ for home
git clone https://github.com/serhataktay/tracktx.git
cd tracktx
./run_pipeline.sh
```
*Tip:* In WSL, `C:\Users\YourName` is `/mnt/c/Users/YourName`. Store data on the Linux filesystem (`~` or `/home/you`) for better performance than `/mnt/c`.

**Verify everything works:**
```bash
nextflow -version
docker --version
docker run --rm hello-world
java -version
```

---

### Option 1: Docker (Recommended for Novices)

Docker packages everything needed—no manual tool installation.

**Step 1: Install Docker Desktop**
- **macOS/Windows:** Download from [docker.com/get-started](https://www.docker.com/products/docker-desktop/)
- **Linux:** `curl -fsSL https://get.docker.com | sh` (or use your package manager)
- Start Docker Desktop and wait until it shows "Running"

**Step 2: Install Git** (if not already installed)
- **macOS:** `xcode-select --install` or install [Xcode Command Line Tools](https://developer.apple.com/xcode/)
- **Windows:** Install [Git for Windows](https://git-scm.com/download/win)
- **Linux:** `sudo apt install git` (Ubuntu/Debian) or equivalent

**Step 3: Clone and run**
```bash
git clone https://github.com/serhataktay/tracktx.git
cd tracktx
./run_pipeline.sh
```

The script auto-detects Docker and runs the pipeline. First run will download the container image (~2–5 min).

**Updating after git pull:** When you `git pull` and run again, `run_pipeline.sh` automatically pulls the Docker image (tag `tracktx:3.0` by default) so pipeline and container stay in sync. To skip the pull (e.g. offline): `TRACKTX_SKIP_PULL=1 ./run_pipeline.sh`

---

### Option 2: Conda

Use Conda if Docker is not available (e.g. restricted HPC, no admin rights).

**Step 1: Install Miniconda**
- Download the installer for your OS: [docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
- Run the installer and follow prompts (accept license, choose install location)
- Restart your terminal, then run `conda --version` to verify

**Step 2: Clone and run**
```bash
git clone https://github.com/serhataktay/tracktx.git
cd tracktx
./run_pipeline.sh
```

The script auto-detects Conda and creates the pipeline environment on first run (~10–20 min).

---

### Option 3: Manual (Advanced)

If you already have Nextflow and the required tools installed:

```bash
git clone https://github.com/serhataktay/tracktx.git
cd tracktx
nextflow run main.nf -entry TrackTx -profile local --samplesheet samplesheet.csv -params-file params.yaml
```

---

### System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| **OS** | Linux, macOS, Windows (WSL2) | Linux or macOS |
| **CPU** | 2+ cores | 8+ cores |
| **RAM** | 8+ GB | 32+ GB |
| **Storage** | 50+ GB | 200+ GB (SSD) |

**Note:** First run downloads reference genomes (~1–5 GB depending on species). Ensure enough free disk space.

---

## 📁 Input Files

### Sample Sheet (`samplesheet.csv`)

Use the [config generator](TrackTx_config_generator.html) to create this, or format manually:

```csv
sample,condition,timepoint,replicate,file1,file2
ctrl_rep1,control,0,1,data/ctrl_R1.fastq,data/ctrl_R2.fastq
heat_rep1,treatment,30,1,data/heat_R1.fastq,data/heat_R2.fastq
```

**Single-end:** Use `file1` only; `file2` can be empty.

**For SRR downloads:** Put the SRR accession in `file1`:
```csv
sample,condition,timepoint,replicate,file1,file2
sample1,control,0,1,SRR123456,
```

**Local files:** Use absolute paths (e.g. `/Users/you/data/sample_R1.fastq.gz`) or paths relative to the project directory.

### Parameters (`params.yaml`)

**Recommended:** Generate this using the interactive generator [`TrackTx_config_generator.html`](TrackTx_config_generator.html), which exports a validated `params.yaml` for you.  
You can also create it manually if you prefer:

```yaml
# Basic Settings
samplesheet: "samplesheet.csv" # Sample sheet CSV (or override via --samplesheet)
sample_source: "local"        # "local" for FASTQ paths, "srr" for SRA downloads
reference_genome: "hg38"       # Human (hg38, hs1), Mouse (mm39, mm10), Fly (dm6), etc.
paired_end: false              # true for paired-end data
output_dir: "./results"        # Where to save results

# Optional: Spike-in Normalization
spikein_genome: "dm6"          # Drosophila spike-in

# Annotation options (affects divergent read assignment)
annotation_source: refseq      # refseq, ensembl, gencode
annotation_exclude_biotypes: ""  # e.g. rRNA,tRNA
annotation_chr_naming: none    # none, add (Ensembl->UCSC), remove (UCSC->Ensembl)

# Optional: Adapter Trimming
adapter_trimming:
  enabled: true
  preset: "illumina"           # illumina, nextera, or custom
  adapter1: "TGGAATTCTCGGGTGCCAAGG"

# Optional: UMI Processing (if your library has UMIs)
umi:
  enabled: false               # Enable if using UMIs
  length: 8                    # Typical: 6-10 bases
  location: 5                  # 5' or 3'

# Advanced: Divergent Transcription (Statistical)
advanced:
  divergent_threshold: auto    # auto = 75th percentile, or specify float
  divergent_sum_thr: auto     # auto = 3x threshold, or specify float
  divergent_fdr: 0.05         # False discovery rate (0.01-0.10)
  divergent_calibration_percentile: 75   # More permissive (targets 50-100K sites)
  divergent_calibration_sum_multiplier: 3
  divergent_merge_gap: 500    # Merge overlapping regions (bp); 0=disabled
  divergent_nt_window: 1000   # Max pairing distance (bp)
  divergent_balance: 0.0      # Balance ratio (0 = max sensitivity)
  divergent_qc: true          # Generate QC reports

# Functional regions (affects divergent read assignment)
functional_regions:
  tss_active_pm: 500          # TSS ± bp for active gene detection
  prom_up: 500
  prom_down: 250
```

**💡 Pro Tip:** Use the interactive config generator (`TrackTx_config_generator.html`) for guided parameter selection with detailed explanations!

---

## 📊 Outputs

```
results/
├── 📈 05_normalized_tracks/        # Load in IGV/UCSC Browser
│   └── <sample>/3p/*.cpm.bw        # CPM and siCPM normalized BigWigs
├── 🔬 06_divergent_tx/             # Divergent transcription (statistical)
│   └── <sample>/
│       ├── divergent_transcription.bed   # High-confidence regions (BED5)
│       └── divergent_transcription_qc.txt # Statistical QC report
├── 🧬 07_functional_regions/       # Genomic region annotations
│   └── <sample>/
│       ├── functional_regions.bed
│       └── functional_regions_summary.tsv
├── 📊 08_pol_metrics/             # Pol-II pausing & density
│   └── <sample>/
│       ├── pausing_index.tsv
│       └── pol_density.tsv
├── 🔍 10_qc/                       # Quality control metrics
│   └── <sample>/qc_pol.json
├── 📋 11_reports/                  # Interactive HTML reports
│   ├── cohort/
│   │   └── global_summary.html     # ⭐ Comprehensive cohort dashboard
│   └── samples/
│       └── <sample>/<sample>.report.html
└── 🔍 trace/                       # Pipeline performance (in output_dir)
    ├── report.html
    ├── timeline.html
    └── trace.txt
```

**Note:** Intermediate outputs (00_references, 01_trimmed_fastq, 02_alignments, 03_genome_tracks, 04_counts, 09_pol_aggregate) are also produced. Trace files live in `{output_dir}/trace/`.

**🎯 Start Here:**
1. **`11_reports/cohort/global_summary.html`** - Comprehensive cohort analysis with:
   - Quality control assessment and outlier detection
   - Divergent transcription patterns across conditions
   - Pol II pausing index distributions
   - Functional region composition analysis
   - Normalization factor validation
   - Interactive sample metrics table
   
2. **`11_reports/samples/<sample>/<sample>.report.html`** - Detailed per-sample reports

3. **`05_normalized_tracks/<sample>/3p/*.cpm.bw`** - Load directly in IGV/UCSC genome browsers

4. **`trace/report.html`** - Pipeline performance and resource usage

---

## 🎯 Common Use Cases

### Time-Course Heat Shock

```bash
# Create sample sheet (or use config generator)
cat > timecourse.csv << EOF
sample,condition,timepoint,replicate,file1,file2
t0_r1,control,0,1,data/t0_r1.fastq,
t0_r2,control,0,2,data/t0_r2.fastq,
t30_r1,heat,30,1,data/t30_r1.fastq,
t30_r2,heat,30,2,data/t30_r2.fastq,
EOF

./run_pipeline.sh --samplesheet timecourse.csv
```

### Drug Treatment with Spike-in

```yaml
# params.yaml
reference_genome: "mm39"
spikein_genome: "dm6"
paired_end: false
```

```bash
./run_pipeline.sh --params-file params.yaml
```

### Download from SRA

**Samplesheet** (SRR accessions in `file1`):
```csv
sample,condition,timepoint,replicate,file1,file2
sample1,control,0,1,SRR4454567,
sample2,treatment,30,1,SRR4454568,
```

**params.yaml** must include `sample_source: srr`:
```yaml
sample_source: srr
reference_genome: "hg38"
output_dir: "./results"
```

```bash
./run_pipeline.sh --samplesheet sra_samples.csv --params-file params.yaml
```

---

## ⚙️ Execution Profiles

The pipeline **auto-detects** your environment, but you can force a specific profile:

| Profile | Description | Use When |
|---------|-------------|----------|
| **docker** | Everything included (recommended) | Docker is available |
| **conda** | Auto environment setup | Docker not available |
| **conda_server** | For NFS/network storage | Conda has filesystem issues |
| **singularity** | HPC containers | On HPC clusters |
| **slurm** | SLURM scheduler | Combine with containers |
| **performance** | External drive optimizations | Use with docker/conda via `--external-drive` |
| **local** | System tools | Tools already installed |

**Examples:**
```bash
./run_pipeline.sh                          # Auto-detect
./run_pipeline.sh -profile docker          # Force Docker
./run_pipeline.sh -profile conda           # Force Conda
./run_pipeline.sh -profile slurm,singularity  # HPC with Singularity
./run_pipeline.sh --external-drive         # USB/exFAT: appends performance profile
```

---

## ⚡ Performance Optimization

### Pipeline Too Slow? Try These Fixes

**Problem**: Pipeline running from USB/exFAT/NAS, or getting "Failed to publish file [link]" errors?

**Solution**: Use **external drive mode** for correct behavior:

```bash
./run_pipeline.sh --external-drive
```

This automatically:
- ✅ Keeps **results** on your project directory (no local space needed for outputs)
- ✅ Fixes publish errors (uses copy instead of hard links)
- ✅ Fixes OverlappingFileLockException (cache, temp, work on local — exFAT lacks file locking)
- ✅ Disables scratch space (reduces file copying on slow storage)
- ✅ Increases task parallelism for better I/O utilization

**Note:** `--external-drive` puts cache, temp, and work (~10–50 GB) on local disk (`~/tmp/tracktx_*`); only results stay on the project dir. Ensure ~20–50 GB free on your internal drive.

### Expected Performance

| Sample Count | Optimized (Internal SSD) | Default (External SSD) |
|--------------|-------------------------|------------------------|
| 2 samples    | 30-60 min               | 2-3 hours ⚠️           |
| 4 samples    | 1-2 hours               | 4-6 hours ⚠️           |
| 8 samples    | 2-4 hours               | 8-12 hours ⚠️          |

## 🔧 Troubleshooting

### Reading Error Messages

When a process fails, Nextflow prints verbose output. **Look for the TRACKTX ERROR block** — it summarizes the problem and fix:

```
═══════════════════════════════════════════════════════════════════════
TRACKTX ERROR
═══════════════════════════════════════════════════════════════════════
Module:  detect_divergent_transcription
Problem: Missing Python dependencies (numpy, pandas, scikit-learn, scipy)
Fix:     pip install numpy pandas scikit-learn scipy | Or use: -profile conda | -profile docker
═══════════════════════════════════════════════════════════════════════
```

- **Quick find:** `grep -A 6 "TRACKTX ERROR"` in the output
- **Full log:** Check the `.log` file in the work dir (shown at the end of the error)
- **Resume:** Add `-resume` to continue after fixing the issue

### Common Issues

**Docker not running:**
```bash
# Install Docker Desktop: https://docs.docker.com/get-docker/
# Make sure it's running before starting pipeline
```

**Out of memory (exit 137):**

OOM = **RAM exhaustion**, not disk space. Removing `work` and `.nextflow` frees disk space but does not fix OOM.

Common causes when you "have enough space" (disk):

1. **Docker memory limit:** Docker Desktop has its own RAM limit (Settings → Resources → Memory). The pipeline detects *host* RAM and allocates per-task memory accordingly—but containers only see Docker’s limit. If Docker has 8 GB and the pipeline assumes 64 GB, multiple tasks can exceed available RAM.
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
3. **USB drive with exFAT/FAT32 (macOS) / OverlappingFileLockException:** exFAT and FAT32 do not support file locking. **Fix:** Use `./run_pipeline.sh --external-drive` — it puts cache, temp, and work (~10–50 GB) on local (`~/tmp/tracktx_cache`, `~/tmp/tracktx_work`); results stay on your project. Ensure ~20–50 GB free on internal drive. Or reformat the USB to **APFS** or **Mac OS Extended**.
4. **NFS / network / cloud-synced storage:** File locking is unreliable on NFS, SMB, iCloud, Dropbox. Set work dir to internal disk or a USB drive formatted as APFS/ext4: `export NXF_WORK=/tmp/nextflow-work` or `export NXF_WORK=/Volumes/MySSD/nextflow-work` (macOS, SSD must be APFS/HFS+).
5. **Conda profile:** Multiple tasks can contend on the conda cache. Try `./run_pipeline.sh -profile docker`, or set `export NXF_CONDA_CACHEDIR=/tmp/conda-$USER-$$` before running.
6. **Upgrade Nextflow:** Pipeline requires ≥24.04.0; older versions have locking issues.

**"matplotlib is building a font cache" seems stuck:**
- Matplotlib scans system fonts on first import (30s–2min). **umi_tools** (preprocess, coverage) and report/aggregate tasks use it.
- **Docker:** The image pre-builds the cache; pull the latest and rebuild if needed.
- **Conda:** `MPLCONFIGDIR` is set to `$TMPDIR` in affected modules. Ensure `TMPDIR` points to local disk (not NFS).
**Spike-in alignment fails (sample-specific):**
- Samples with many unaligned reads (e.g. 20M+) need more memory for spike-in alignment
- Increase Docker memory (Settings → Resources) or system RAM
- Check `bowtie2_spikein.log` in the failed task's work dir for details

**Finished tasks re-run from sample 1 (even with -resume):**
- Nextflow’s cache depends on input file path, size, and timestamp. NFS/network storage can give inconsistent timestamps → add `preprocess_reads_lenient_cache: true` to params.yaml or run with `--preprocess_reads_lenient_cache`.
- Docker `:latest` changes when the image is updated → use a fixed tag (e.g. `tracktx:3.0`) for stable caching.
- Debug: `nextflow run ... -resume -dump-hashes 2>&1 | grep "cache hash"` and compare between runs.

**preprocess_and_quality_filter_reads re-runs after stop/restart:**
- When you Ctrl+C, running/queued tasks are cancelled and not cached
- Only completed tasks are reused with `-resume`
- Let the pipeline finish, or stop when no preprocess_and_quality_filter_reads tasks are active

### Getting Help

1. **Check logs**: `.nextflow.log` in the working directory
2. **Review trace**: `{output_dir}/trace/report.html` for resource issues
3. **Monitor live**: `python3 nfmon.py` to see what's happening
4. **GitHub Issues**: [Report bugs](https://github.com/serhataktay/tracktx/issues)

---

## 📖 Documentation

| Document | Description |
|----------|-------------|
| [TrackTx_config_generator.html](TrackTx_config_generator.html) | Interactive config and samplesheet generator |

---

## 🧬 Citation

If TrackTx is useful for your research, please cite: [https://github.com/serhataktay/tracktx](https://github.com/serhataktay/tracktx)

### Key References
- **PRO-seq**: [Kwak et al., 2013](https://doi.org/10.1126/science.1229386) - Precision run-on sequencing
- **Nascent RNA-seq**: [Core et al., 2008](https://doi.org/10.1126/science.1162228) - Global nascent transcription
- **Nextflow**: [Di Tommaso et al., 2017](https://doi.org/10.1038/nbt.3820) - Scalable workflows

---

## 📜 License

TrackTx is released under the [MIT License](LICENSE).

---

<div align="center">

**⭐ Star this repo if TrackTx is useful for your research!**

[Issues](https://github.com/serhataktay/tracktx/issues) • [Releases](https://github.com/serhataktay/tracktx/releases)

</div>
