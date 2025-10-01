# TrackTx Deployment Guide

This guide covers deploying the TrackTx pipeline on different computing environments, from local development to HPC clusters.

## Table of Contents

1. [Local Development](#local-development)
2. [Server Deployment](#server-deployment)
3. [HPC Cluster Deployment](#hpc-cluster-deployment)
4. [Configuration Profiles](#configuration-profiles)
5. [Storage Considerations](#storage-considerations)
6. [Troubleshooting](#troubleshooting)

## Local Development

### Prerequisites

- **Nextflow**: Version 24.04.0 or higher
- **Java**: OpenJDK 11 or higher
- **Container Runtime**: Docker, Singularity, or Conda
- **Resources**: 16GB RAM, 8 CPU cores minimum

### Quick Start

```bash
# Clone the repository
git clone https://github.com/serhataktay/tracktx-nf.git
cd tracktx-nf

# Run with Docker (recommended)
nextflow run main.nf -entry TrackTx -profile docker -params-file params.yaml

# Run with Conda
nextflow run main.nf -entry TrackTx -profile conda -params-file params.yaml

# Run with Singularity
nextflow run main.nf -entry TrackTx -profile singularity -params-file params.yaml
```

### Local Configuration

Create a `params.local.yaml` file:

```yaml
# Local development parameters
output_dir: './results'
work_dir: './work'

# Resource settings
host_cpus: 8
host_mem: 32

# Container settings
container:
  enabled: true
  engine: 'docker'
  cache_dir: '~/.docker'
```

## Server Deployment

### Automated Setup

Use the provided setup script for automated server configuration:

```bash
# Download and run the setup script
curl -O https://raw.githubusercontent.com/serhataktay/tracktx-nf/main/scripts/setup_server.sh
chmod +x setup_server.sh
./setup_server.sh
```

### Manual Setup

#### 1. Install Dependencies

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install -y curl wget git build-essential ca-certificates gnupg lsb-release software-properties-common
```

**CentOS/RHEL:**
```bash
sudo yum update -y
sudo yum install -y curl wget git gcc gcc-c++ make ca-certificates gnupg2
```

#### 2. Install Nextflow

```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
sudo chmod +x /usr/local/bin/nextflow
```

#### 3. Install Container Runtime

**Docker:**
```bash
# Ubuntu/Debian
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin

# Add user to docker group
sudo usermod -aG docker $USER
```

**Singularity:**
```bash
# Install from conda-forge
conda install -y -c conda-forge singularity
```

#### 4. Configure Storage

```bash
# Create storage directories
sudo mkdir -p /opt/tracktx/{assets,genomes,cache,projects}
sudo mkdir -p /var/tmp/nextflow-work
sudo chown -R $USER:$USER /opt/tracktx
sudo chmod -R 755 /opt/tracktx
```

#### 5. Server Configuration

Create `params.server.yaml`:

```yaml
# Server-specific parameters
output_dir: '/opt/tracktx/projects/results'
work_dir: '/var/tmp/nextflow-work'

# Storage configuration
storage:
  shared_data: '/opt/tracktx/data'
  shared_scratch: '/var/tmp'
  project_dir: '/opt/tracktx/projects'
  cache_root: '/opt/tracktx/cache'
  assets_root: '/opt/tracktx/assets'

# Asset and genome caching
assets_dir: '/opt/tracktx/assets'
genome_cache: '/opt/tracktx/genomes'

# Server resource settings
host_cpus: 16
host_mem: 64

# Container settings
container:
  enabled: true
  engine: 'docker'
  cache_dir: '/opt/tracktx/containers'
```

### Running on Server

```bash
# Run with server profile
nextflow run main.nf -entry TrackTx -profile server,docker -params-file params.server.yaml --samplesheet samplesheet.csv
```

## HPC Cluster Deployment

### Supported Schedulers

- **SLURM**: Most common in academic HPC centers
- **SGE/UGE**: Common in some academic and commercial environments
- **PBS/Torque**: Used in some HPC environments

### Automated Setup

Use the provided HPC setup script:

```bash
# Download and run the setup script
curl -O https://raw.githubusercontent.com/serhataktay/tracktx-nf/main/scripts/setup_hpc.sh
chmod +x setup_hpc.sh
./setup_hpc.sh
```

### Manual HPC Setup

#### 1. Load Required Modules

```bash
# Load environment modules
module load singularity
module load nextflow  # if available

# Or install Nextflow locally
curl -s https://get.nextflow.io | bash
export PATH=$PWD:$PATH
```

#### 2. Configure Environment

```bash
# Set up environment variables
export NXF_HOME="$HOME/.nextflow"
export SINGULARITY_CACHEDIR="$HOME/.singularity"
export NXF_TEMP="/tmp/nextflow-work"

# Create directories
mkdir -p "$NXF_HOME" "$SINGULARITY_CACHEDIR" "$NXF_TEMP"
```

#### 3. HPC Configuration

Create `params.hpc.yaml`:

```yaml
# HPC-specific parameters
output_dir: './results'
work_dir: '/tmp/nextflow-work'

# HPC storage configuration
storage:
  shared_data: '/shared/data'
  shared_scratch: '/scratch'
  project_dir: '/shared/projects/tracktx'
  cache_root: '/shared/cache'
  assets_root: '/shared/assets'

# Scheduler settings
scheduler:
  type: 'slurm'  # or 'sge', 'pbs'
  partition: 'compute'
  account: 'default'
  time_limit: '48:00:00'
  max_jobs: 50

# Container settings
container:
  enabled: true
  engine: 'singularity'
  cache_dir: '/shared/containers'
```

### Job Submission

#### SLURM

```bash
# Submit job
sbatch scripts/run_tracktx.slurm --samplesheet samplesheet.csv

# Monitor jobs
squeue -u $USER

# Check job details
scontrol show job <job_id>
```

#### SGE

```bash
# Submit job
qsub scripts/run_tracktx.sge --samplesheet samplesheet.csv

# Monitor jobs
qstat -u $USER

# Check job details
qstat -j <job_id>
```

#### PBS

```bash
# Submit job
qsub scripts/run_tracktx.pbs --samplesheet samplesheet.csv

# Monitor jobs
qstat -u $USER

# Check job details
qstat -f <job_id>
```

## Configuration Profiles

### Available Profiles

| Profile | Description | Use Case |
|---------|-------------|----------|
| `conda` | Local execution with Conda | Development, small datasets |
| `docker` | Local execution with Docker | Development, reproducible environments |
| `server` | Server-optimized settings | Dedicated server deployment |
| `slurm` | SLURM scheduler | HPC clusters with SLURM |
| `sge` | SGE scheduler | HPC clusters with SGE |
| `pbs` | PBS scheduler | HPC clusters with PBS |
| `singularity` | Singularity containers | HPC without Docker |
| `apptainer` | Apptainer containers | Modern Singularity fork |

### Profile Combinations

```bash
# Server with Docker
nextflow run main.nf -profile server,docker

# HPC with SLURM and Singularity
nextflow run main.nf -profile slurm,singularity

# HPC with SGE and Apptainer
nextflow run main.nf -profile sge,apptainer
```

### Custom Profiles

Create custom profiles in `nextflow.config`:

```groovy
profiles {
  custom {
    process.executor = 'local'
    process.cpus = 16
    process.memory = '64 GB'
    process.container = 'your-custom-image:latest'
  }
}
```

## Storage Considerations

### Local Storage

- **Work Directory**: Fast SSD recommended for temporary files
- **Output Directory**: Sufficient space for results (100GB+)
- **Cache Directory**: Persistent storage for genome indices

### HPC Storage

- **Shared Storage**: Network-attached storage for persistent data
- **Scratch Space**: High-performance temporary storage
- **User Home**: Limited space, use for configuration only

### Storage Optimization

```yaml
# Optimize for HPC storage
storage:
  shared_data: '/shared/data'      # Persistent, shared access
  shared_scratch: '/scratch'       # High-performance, temporary
  user_scratch: '/scratch/$USER'   # User-specific scratch
  cache_root: '/shared/cache'      # Persistent cache
```

## Troubleshooting

### Common Issues

#### 1. Container Issues

**Problem**: Container not found or permission denied
```bash
# Solution: Check container availability
singularity exec docker://ghcr.io/serhataktay/tracktx:latest echo "Container works"

# Or use local image
docker pull ghcr.io/serhataktay/tracktx:latest
```

#### 2. Resource Issues

**Problem**: Out of memory or CPU limits
```bash
# Solution: Adjust resource settings in params file
host_cpus: 8
host_mem: 32
```

#### 3. Storage Issues

**Problem**: No space left on device
```bash
# Solution: Clean up old work directories
nextflow clean -f
rm -rf work/
```

#### 4. HPC Scheduler Issues

**Problem**: Jobs not starting or failing
```bash
# Solution: Check scheduler configuration
squeue -u $USER  # SLURM
qstat -u $USER   # SGE/PBS

# Check logs
tail -f logs/tracktx-*.out
```

### Performance Optimization

#### 1. Resource Tuning

```yaml
# Optimize for your hardware
process {
  withName: 'run_alignment' {
    cpus = 12
    memory = '48 GB'
  }
}
```

#### 2. Storage Optimization

```yaml
# Use fast storage for work directory
work_dir: '/fast/scratch/nextflow-work'
```

#### 3. Container Optimization

```yaml
# Use local container cache
singularity {
  cacheDir = '/shared/containers'
}
```

### Getting Help

1. **Check Logs**: Review Nextflow and process logs
2. **Resource Usage**: Monitor CPU, memory, and disk usage
3. **Community Support**: GitHub issues and discussions
4. **Documentation**: Check this guide and inline help

### Support Contacts

- **GitHub Issues**: [TrackTx Issues](https://github.com/serhataktay/tracktx-nf/issues)
- **Email**: serhat.aktay@scilifelab.se
- **Documentation**: [TrackTx Docs](https://github.com/serhataktay/tracktx-nf/docs)

## Best Practices

### 1. Resource Planning

- **CPU**: Plan for 8-16 cores per sample
- **Memory**: 32-64GB RAM recommended
- **Storage**: 100GB+ for results and temporary files
- **Time**: 24-48 hours for typical datasets

### 2. Data Management

- **Backup**: Regular backups of important results
- **Version Control**: Track pipeline versions and parameters
- **Documentation**: Document custom configurations

### 3. Monitoring

- **Job Monitoring**: Use provided monitoring scripts
- **Resource Usage**: Monitor CPU, memory, and disk usage
- **Log Analysis**: Regular review of pipeline logs

### 4. Security

- **Access Control**: Proper file permissions
- **Data Privacy**: Secure handling of sensitive data
- **Container Security**: Regular updates of container images
