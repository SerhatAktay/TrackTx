#!/bin/bash
# ════════════════════════════════════════════════════════════════════════════
# TrackTx HPC Setup Script
# Configures an HPC cluster environment for running the TrackTx pipeline
# ════════════════════════════════════════════════════════════════════════════

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SETUP_USER="${SETUP_USER:-$(whoami)}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Logging functions
log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[SUCCESS]${NC} $1"
}

log_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Detect HPC scheduler
detect_scheduler() {
    log_info "Detecting HPC scheduler..."
    
    if command -v squeue >/dev/null 2>&1; then
        SCHEDULER="slurm"
        log_success "Detected SLURM scheduler"
    elif command -v qstat >/dev/null 2>&1; then
        if qstat -help 2>&1 | grep -q "SGE"; then
            SCHEDULER="sge"
            log_success "Detected SGE/UGE scheduler"
        else
            SCHEDULER="pbs"
            log_success "Detected PBS scheduler"
        fi
    else
        log_error "No supported HPC scheduler detected (SLURM, SGE, PBS)"
        log_info "Please ensure you're running this on an HPC cluster"
        exit 1
    fi
}

# Check HPC requirements
check_hpc_requirements() {
    log_info "Checking HPC cluster requirements..."
    
    # Check if we're on a compute node or login node
    if [[ -n "${SLURM_JOB_ID:-}" ]] || [[ -n "${PBS_JOBID:-}" ]] || [[ -n "${JOB_ID:-}" ]]; then
        log_info "Running on compute node - good for testing"
    else
        log_warning "Running on login node - ensure you have sufficient resources"
    fi
    
    # Check available modules
    if command -v module >/dev/null 2>&1; then
        log_success "Environment modules available"
        
        # Check for common modules
        local modules=("singularity" "apptainer" "docker" "nextflow" "java")
        for module in "${modules[@]}"; do
            if module avail "$module" 2>&1 | grep -q "$module"; then
                log_success "Module available: $module"
            else
                log_warning "Module not available: $module"
            fi
        done
    else
        log_warning "Environment modules not available"
    fi
    
    # Check scratch space
    if [[ -d "/scratch" ]]; then
        local scratch_gb=$(df -BG /scratch 2>/dev/null | awk 'NR==2{print $4}' | sed 's/G//' || echo "unknown")
        log_success "Scratch space available: ${scratch_gb}GB"
    else
        log_warning "No /scratch directory found"
    fi
}

# Setup user environment
setup_user_environment() {
    log_info "Setting up user environment..."
    
    # Create user directories
    mkdir -p "$HOME/.nextflow"
    mkdir -p "$HOME/.singularity"
    mkdir -p "$HOME/.config/nextflow"
    
    # Create project workspace
    mkdir -p "$PROJECT_DIR/workspace"
    mkdir -p "$PROJECT_DIR/logs"
    mkdir -p "$PROJECT_DIR/tmp"
    
    # Set up environment variables
    cat >> "$HOME/.bashrc" << EOF

# TrackTx HPC Environment
export NXF_HOME="$HOME/.nextflow"
export SINGULARITY_CACHEDIR="$HOME/.singularity"
export NXF_TEMP="/tmp/nextflow-work"
export NXF_OPTS="-Xmx4g -XX:+UseG1GC"
EOF
    
    # Create Nextflow config
    cat > "$HOME/.nextflow/config" << EOF
// TrackTx HPC Configuration
workDir = '/tmp/nextflow-work'
cacheDir = '$HOME/.nextflow/cache'

// HPC-specific settings
executor {
  queueSize = 100
  submitRateLimit = '10/1min'
}

// Container settings
singularity {
  enabled = true
  cacheDir = '$HOME/.singularity'
  autoMounts = true
}

// Storage settings
params {
  scratch_dir = '/tmp'
  shared_storage = '/shared'
}
EOF
    
    log_success "User environment configured"
}

# Create HPC-specific configuration
create_hpc_config() {
    log_info "Creating HPC-specific configuration..."
    
    # Create HPC params file
    cat > "$PROJECT_DIR/params.hpc.yaml" << EOF
# HPC-specific parameters
output_dir: './results'
work_dir: '/tmp/nextflow-work'

# HPC storage configuration
storage:
  shared_data: '/shared/data'
  shared_scratch: '/scratch'
  shared_tmp: '/tmp'
  project_dir: '/shared/projects/tracktx'
  user_scratch: '/scratch/\$USER'
  user_tmp: '/tmp/\$USER'
  cache_root: '/shared/cache'
  assets_root: '/shared/assets'

# Asset and genome caching
assets_dir: '/shared/assets/tracktx'
genome_cache: '/shared/cache/genomes'

# HPC resource settings (will be auto-detected)
host_cpus: \$(nproc)
host_mem: \$(free -g | awk '/^Mem:/{print \$2}')

# Scheduler-specific settings
scheduler:
  type: '$SCHEDULER'
  partition: 'compute'
  account: 'default'
  time_limit: '48:00:00'
  max_jobs: 50
  max_cpus_per_job: 32

# Container settings
container:
  enabled: true
  engine: 'singularity'
  cache_dir: '/shared/containers'
  bind_paths:
    - '/scratch:/scratch'
    - '/tmp:/tmp'
    - '/shared:/shared'
EOF
    
    log_success "HPC configuration created"
}

# Create job submission scripts
create_job_scripts() {
    log_info "Creating job submission scripts..."
    
    # Create SLURM job script
    if [[ "$SCHEDULER" == "slurm" ]]; then
        cat > "$PROJECT_DIR/scripts/run_tracktx.slurm" << 'EOF'
#!/bin/bash
#SBATCH --job-name=tracktx
#SBATCH --output=logs/tracktx-%j.out
#SBATCH --error=logs/tracktx-%j.err
#SBATCH --time=48:00:00
#SBATCH --partition=compute
#SBATCH --account=default
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1

# Load required modules
module load singularity
module load nextflow 2>/dev/null || true

# Set up environment
export NXF_HOME="$HOME/.nextflow"
export SINGULARITY_CACHEDIR="$HOME/.singularity"
export NXF_TEMP="/tmp/nextflow-work"

# Create work directory
mkdir -p "$NXF_TEMP"

# Run TrackTx
nextflow run main.nf \
    -profile slurm,singularity \
    -params-file params.hpc.yaml \
    -work-dir "$NXF_TEMP" \
    -with-trace \
    -with-report \
    -with-timeline \
    -with-dag \
    "$@"
EOF
        chmod +x "$PROJECT_DIR/scripts/run_tracktx.slurm"
        log_success "SLURM job script created"
    fi
    
    # Create SGE job script
    if [[ "$SCHEDULER" == "sge" ]]; then
        cat > "$PROJECT_DIR/scripts/run_tracktx.sge" << 'EOF'
#!/bin/bash
#$ -N tracktx
#$ -o logs/tracktx-$JOB_ID.out
#$ -e logs/tracktx-$JOB_ID.err
#$ -l h_rt=48:00:00
#$ -l slots=8
#$ -l h_vmem=32G
#$ -V
#$ -j y
#$ -cwd

# Load required modules
module load singularity
module load nextflow 2>/dev/null || true

# Set up environment
export NXF_HOME="$HOME/.nextflow"
export SINGULARITY_CACHEDIR="$HOME/.singularity"
export NXF_TEMP="/tmp/nextflow-work"

# Create work directory
mkdir -p "$NXF_TEMP"

# Run TrackTx
nextflow run main.nf \
    -profile sge,singularity \
    -params-file params.hpc.yaml \
    -work-dir "$NXF_TEMP" \
    -with-trace \
    -with-report \
    -with-timeline \
    -with-dag \
    "$@"
EOF
        chmod +x "$PROJECT_DIR/scripts/run_tracktx.sge"
        log_success "SGE job script created"
    fi
    
    # Create PBS job script
    if [[ "$SCHEDULER" == "pbs" ]]; then
        cat > "$PROJECT_DIR/scripts/run_tracktx.pbs" << 'EOF'
#!/bin/bash
#PBS -N tracktx
#PBS -o logs/tracktx-$PBS_JOBID.out
#PBS -e logs/tracktx-$PBS_JOBID.err
#PBS -l walltime=48:00:00
#PBS -l nodes=1:ppn=8
#PBS -l mem=32gb
#PBS -V
#PBS -j oe

# Load required modules
module load singularity
module load nextflow 2>/dev/null || true

# Set up environment
export NXF_HOME="$HOME/.nextflow"
export SINGULARITY_CACHEDIR="$HOME/.singularity"
export NXF_TEMP="/tmp/nextflow-work"

# Create work directory
mkdir -p "$NXF_TEMP"

# Run TrackTx
nextflow run main.nf \
    -profile pbs,singularity \
    -params-file params.hpc.yaml \
    -work-dir "$NXF_TEMP" \
    -with-trace \
    -with-report \
    -with-timeline \
    -with-dag \
    "$@"
EOF
        chmod +x "$PROJECT_DIR/scripts/run_tracktx.pbs"
        log_success "PBS job script created"
    fi
}

# Create monitoring script
create_monitoring_script() {
    log_info "Creating monitoring script..."
    
    cat > "$PROJECT_DIR/scripts/monitor_hpc.sh" << 'EOF'
#!/bin/bash
# TrackTx HPC Monitoring Script

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
LOG_FILE="$PROJECT_DIR/logs/monitor.log"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

log() {
    echo -e "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

log_info() {
    log "${GREEN}[INFO]${NC} $1"
}

log_warning() {
    log "${YELLOW}[WARNING]${NC} $1"
}

log_error() {
    log "${RED}[ERROR]${NC} $1"
}

# Monitor job status
monitor_jobs() {
    log_info "Monitoring TrackTx jobs..."
    
    # Check for running jobs
    if command -v squeue >/dev/null 2>&1; then
        local running_jobs=$(squeue -u "$USER" --name=tracktx --format="%.18i %.9P %.8j %.8u %.2t %.10M %.6D %R" --noheader | wc -l)
        log_info "SLURM jobs running: $running_jobs"
    elif command -v qstat >/dev/null 2>&1; then
        local running_jobs=$(qstat -u "$USER" | grep tracktx | wc -l)
        log_info "SGE/PBS jobs running: $running_jobs"
    fi
    
    # Check disk usage
    local disk_usage=$(df -h . | awk 'NR==2{print $5}' | sed 's/%//')
    if [[ $disk_usage -gt 80 ]]; then
        log_warning "Disk usage is high: ${disk_usage}%"
    else
        log_info "Disk usage: ${disk_usage}%"
    fi
    
    # Check scratch space
    if [[ -d "/scratch" ]]; then
        local scratch_usage=$(df -h /scratch | awk 'NR==2{print $5}' | sed 's/%//')
        if [[ $scratch_usage -gt 90 ]]; then
            log_warning "Scratch space usage is high: ${scratch_usage}%"
        else
            log_info "Scratch space usage: ${scratch_usage}%"
        fi
    fi
}

# Clean up old files
cleanup() {
    log_info "Cleaning up old files..."
    
    # Clean up old logs (older than 30 days)
    find "$PROJECT_DIR/logs" -name "*.out" -mtime +30 -delete 2>/dev/null || true
    find "$PROJECT_DIR/logs" -name "*.err" -mtime +30 -delete 2>/dev/null || true
    
    # Clean up old Nextflow work directories
    find "/tmp/nextflow-work" -type d -mtime +7 -exec rm -rf {} + 2>/dev/null || true
    
    log_info "Cleanup completed"
}

# Main monitoring loop
main() {
    log_info "Starting TrackTx HPC monitoring..."
    
    while true; do
        monitor_jobs
        cleanup
        sleep 300  # Check every 5 minutes
    done
}

# Run main function
main "$@"
EOF
    
    chmod +x "$PROJECT_DIR/scripts/monitor_hpc.sh"
    log_success "Monitoring script created"
}

# Main setup function
main() {
    log_info "Starting TrackTx HPC setup..."
    
    detect_scheduler
    check_hpc_requirements
    setup_user_environment
    create_hpc_config
    create_job_scripts
    create_monitoring_script
    
    log_success "TrackTx HPC setup completed!"
    log_info "Next steps:"
    log_info "1. Load required modules: module load singularity nextflow"
    log_info "2. Test the installation: nextflow run $PROJECT_DIR/main.nf -profile ${SCHEDULER},singularity -help"
    log_info "3. Submit a job: $PROJECT_DIR/scripts/run_tracktx.${SCHEDULER} --samplesheet samplesheet.csv"
    
    log_warning "Remember to:"
    log_warning "- Configure your cluster-specific settings in params.hpc.yaml"
    log_warning "- Set appropriate partition/queue names"
    log_warning "- Adjust resource limits based on your cluster's policies"
}

# Run main function
main "$@"
