# Installation

[← Back to README](../README.md)

> **Why one container for every step?** TrackTx ships a single Docker image / conda
> environment (`envs/Dockerfile`, `envs/tracktx.yaml`) covering all 17 pipeline stages,
> rather than nf-core's convention of one container per tool. This is deliberate: every
> run of a given TrackTx version uses the exact same toolchain end to end, so citing one
> image tag (or conda lockfile) in a paper's Methods section fully specifies every tool
> version used, with no per-process container matrix to reconcile. The tradeoff is that
> adding or upgrading one tool rebuilds the whole image, which is fine for a pipeline
> with a fixed, curated toolchain rather than one that composes many independently
> versioned community modules.

## Prerequisites

You need **Nextflow** (the workflow engine) plus **one** of Docker or Conda (for the tools):

| Requirement | Purpose | Install |
|-------------|---------|---------|
| **Nextflow** (≥26.04.0) | Runs the pipeline | See below |
| **Docker Desktop** | Easiest, packages all tools | [Get Docker](https://docs.docker.com/get-docker/) |
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
nextflow -version   # Must show 26.04.0 or higher
docker --version   # If using Docker
conda --version    # If using Conda
```

---

## Windows (WSL)

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
- During setup, enable "Use the WSL 2 based engine"
- Start Docker Desktop and ensure it shows "Running"

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

## Option 1: Docker (recommended for novices)

Docker packages everything needed, so there's no manual tool installation.

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

The script auto-detects Docker and runs the pipeline. First run will download the container image (~2-5 min).

**Updating after git pull:** When you `git pull` and run again, `run_pipeline.sh` automatically pulls the Docker image (tag `tracktx:1.4.0` by default) so pipeline and container stay in sync. To skip the pull (e.g. offline): `TRACKTX_SKIP_PULL=1 ./run_pipeline.sh`

---

## Option 2: Conda

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

The script auto-detects Conda and creates the pipeline environment on first run (~10-20 min).

---

## Option 3: Manual (advanced)

If you already have Nextflow and the required tools installed:

```bash
git clone https://github.com/serhataktay/tracktx.git
cd tracktx
nextflow run main.nf -entry TrackTx -profile local --samplesheet samplesheet.csv -params-file params.yaml
```

---

## System Requirements

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| **OS** | Linux, macOS, Windows (WSL2) | Linux or macOS |
| **CPU** | 2+ cores | 8+ cores |
| **RAM** | 8+ GB | 32+ GB |
| **Storage** | 50+ GB | 200+ GB (SSD) |

**Note:** First run downloads reference genomes (~1-5 GB depending on species). Ensure enough free disk space.
