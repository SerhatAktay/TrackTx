#!/bin/bash
# ═══════════════════════════════════════════════════════════════════════════
# TrackTx Pipeline Runner - Universal Launcher
# ═══════════════════════════════════════════════════════════════════════════
#
# PURPOSE:
#   Smart launcher that automatically detects your environment and runs
#   the pipeline with optimal settings. Works on:
#   • Laptops (Mac/Linux)
#   • Workstations
#   • Servers (standard or problematic)
#   • HPC clusters
#
# USAGE:
#   ./run_pipeline.sh                    # Auto-detect everything
#   ./run_pipeline.sh --fast             # Performance mode (external drives)
#   ./run_pipeline.sh -profile docker    # Force specific profile
#   ./run_pipeline.sh --help             # Show full help
#
# FEATURES:
#   • Divergent transcription: edge-to-edge distance algorithm
#   • Automatic resource optimization for any system size
#   • Smart profile detection (Docker, Conda, Singularity, HPC)
#   • Network storage compatibility detection
#
# ═══════════════════════════════════════════════════════════════════════════

set -euo pipefail

# ═══════════════════════════════════════════════════════════════════════════
# COLORS & FORMATTING
# ═══════════════════════════════════════════════════════════════════════════

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
BOLD='\033[1m'
NC='\033[0m'

info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[✓]${NC} $1"; }
warning() { echo -e "${YELLOW}[⚠]${NC} $1"; }
error()   { echo -e "${RED}[✗]${NC} $1"; }
header()  { echo -e "${BOLD}${CYAN}$1${NC}"; }

separator() {
    echo "═══════════════════════════════════════════════════════════════════════════"
}

# ═══════════════════════════════════════════════════════════════════════════
# SYSTEM DETECTION UTILITIES
# ═══════════════════════════════════════════════════════════════════════════

# Check if command exists
has_command() { command -v "$1" >/dev/null 2>&1; }

# Detect number of CPUs
detect_cpus() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        sysctl -n hw.ncpu 2>/dev/null || echo 4
    else
        nproc 2>/dev/null || grep -c ^processor /proc/cpuinfo 2>/dev/null || echo 4
    fi
}

# Detect total memory in GB
detect_memory_gb() {
    if [[ "$OSTYPE" == "darwin"* ]]; then
        echo $(($(sysctl -n hw.memsize 2>/dev/null || echo 8589934592) / 1073741824))
    else
        echo $(($(grep MemTotal /proc/meminfo 2>/dev/null | awk '{print $2}' || echo 8388608) / 1048576))
    fi
}

# Categorize system size
categorize_system() {
    local cpus=$(detect_cpus)
    local mem_gb=$(detect_memory_gb)
    
    if [[ $cpus -le 4 && $mem_gb -le 8 ]]; then
        echo "laptop"
    elif [[ $cpus -le 8 && $mem_gb -le 16 ]]; then
        echo "desktop"
    elif [[ $cpus -le 16 && $mem_gb -le 32 ]]; then
        echo "workstation"
    elif [[ $cpus -le 32 && $mem_gb -le 64 ]]; then
        echo "server"
    else
        echo "hpc"
    fi
}

# Detect if running on HPC scheduler
detect_hpc() {
    if has_command sbatch || has_command sinfo; then
        echo "slurm"
    elif has_command qsub || has_command qstat; then
        echo "pbs"
    elif has_command bsub || has_command bjobs; then
        echo "lsf"
    else
        echo "none"
    fi
}

# Check if path is on network storage (NFS, SMB, etc.)
is_network_storage() {
    local path="${1:-.}"
    
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        # Check mount type
        local mount_point=$(df "$path" 2>/dev/null | tail -1 | awk '{print $1}')
        local fs_type=$(df -T "$path" 2>/dev/null | tail -1 | awk '{print $2}')
        
        case "$fs_type" in
            nfs|nfs4|cifs|smb|smbfs|afs|gpfs|lustre|beegfs)
                return 0
                ;;
        esac
    fi
    
    return 1
}

# ═══════════════════════════════════════════════════════════════════════════
# CONTAINER DETECTION
# ═══════════════════════════════════════════════════════════════════════════

# Check if Docker daemon is running (silent check)
docker_daemon_running() {
    if ! has_command docker; then
        return 1
    fi
    timeout 5s docker info >/dev/null 2>&1
}

# Attempt to start Docker (with user interaction)
docker_start() {
    # Already running
    if docker_daemon_running; then
        return 0
    fi
    
    warning "Docker is installed but not running"
    
    # Non-interactive mode
    if [[ ! -t 0 ]] || [[ ${NO_DOCKER_PROMPT:-0} -eq 1 ]]; then
        info "Non-interactive mode - please start Docker manually"
        if [[ "$OSTYPE" == "darwin"* ]] && [[ -d "/Applications/Docker.app" ]]; then
            info "  Command: open -a Docker"
        elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
            info "  Command: sudo systemctl start docker"
        fi
        return 1
    fi
    
    # macOS Docker Desktop
    if [[ "$OSTYPE" == "darwin"* ]] && [[ -d "/Applications/Docker.app" ]]; then
        echo -n "Start Docker Desktop? [Y/n]: "
        read -r -t 60 ANSWER || ANSWER="n"
        
        case "${ANSWER:-y}" in
            [Yy]*)
                info "Starting Docker Desktop..."
                open -a Docker
                
                # Wait for daemon (max 60s)
                for i in {1..30}; do
                    sleep 2
                    if docker_daemon_running; then
                        success "Docker started!"
                        return 0
                    fi
                    echo -ne "Waiting... ${i}/30\r"
                done
                echo ""
                error "Docker failed to start within 60 seconds"
                return 1
                ;;
            *)
                info "Please start Docker Desktop manually"
                return 1
                ;;
        esac
    fi
    
    # Linux
    if [[ "$OSTYPE" == "linux-gnu"* ]]; then
        info "To start Docker daemon:"
        info "  sudo systemctl start docker"
    fi
    
    return 1
}

# Check Singularity/Apptainer availability
singularity_available() {
    has_command singularity || has_command apptainer
}

# Check Podman availability
podman_available() {
    has_command podman && podman info >/dev/null 2>&1
}

# ═══════════════════════════════════════════════════════════════════════════
# CONDA DETECTION
# ═══════════════════════════════════════════════════════════════════════════

# Detect conda variant (conda, mamba, micromamba)
detect_conda_variant() {
    if has_command micromamba; then
        echo "micromamba"
    elif has_command mamba; then
        echo "mamba"
    elif has_command conda; then
        echo "conda"
    else
        echo "none"
    fi
}

# Check if conda is available (any variant)
conda_available() {
    [[ "$(detect_conda_variant)" != "none" ]]
}

# Check if conda environment needs special handling (server mode)
needs_conda_server_profile() {
    # Check for network storage
    if is_network_storage "."; then
        return 0
    fi
    
    # Check for common HPC indicators
    if [[ -d "/gpfs" ]] || [[ -d "/lustre" ]] || [[ -d "/beegfs" ]]; then
        return 0
    fi
    
    return 1
}

# ═══════════════════════════════════════════════════════════════════════════
# SMART PROFILE DETECTION
# ═══════════════════════════════════════════════════════════════════════════

detect_best_profile() {
    local hpc_type=$(detect_hpc)
    local system_type=$(categorize_system)
    
    # HPC Environment
    if [[ "$hpc_type" != "none" ]]; then
        if singularity_available; then
            echo "${hpc_type},singularity"
            return
        elif conda_available; then
            echo "${hpc_type},conda_server"
            return
        fi
    fi
    
    # Desktop/Workstation/Server
    if docker_daemon_running; then
        echo "docker"
    elif podman_available; then
        echo "podman"
    elif conda_available; then
        if needs_conda_server_profile; then
            echo "conda_server"
        else
            echo "conda"
        fi
    elif singularity_available; then
        echo "singularity"
    else
        echo "local"
    fi
}

# ═══════════════════════════════════════════════════════════════════════════
# SYSTEM RESOURCE RECOMMENDATIONS
# ═══════════════════════════════════════════════════════════════════════════

show_system_info() {
    local cpus=$(detect_cpus)
    local mem_gb=$(detect_memory_gb)
    local system_type=$(categorize_system)
    
    separator
    header "System Information"
    separator
    echo "  CPUs:        ${cpus} cores"
    echo "  Memory:      ${mem_gb} GB"
    echo "  Category:    ${system_type}"
    
    # System-specific recommendations
    case "$system_type" in
        laptop)
            echo "  Expected:    Sequential processing (1 sample at a time)"
            echo "  Runtime:     ~1-2 hours per sample"
            ;;
        desktop)
            echo "  Expected:    2-3 samples in parallel"
            echo "  Runtime:     ~30-60 minutes per batch"
            ;;
        workstation)
            echo "  Expected:    3-5 samples in parallel"
            echo "  Runtime:     ~20-40 minutes per batch"
            ;;
        server)
            echo "  Expected:    5-10 samples in parallel"
            echo "  Runtime:     ~15-30 minutes per batch"
            ;;
        hpc)
            echo "  Expected:    10+ samples in parallel"
            echo "  Runtime:     ~10-20 minutes per batch"
            ;;
    esac
    
    # Storage recommendations
    echo ""
    echo -e "  Storage:     Ensure ${BOLD}10-20 GB free${NC} per sample"
    
    if is_network_storage "."; then
        warning "Detected network storage - using conda_server profile for compatibility"
    fi
    
    separator
}

# ═══════════════════════════════════════════════════════════════════════════
# PROFILE VALIDATION & DESCRIPTION
# ═══════════════════════════════════════════════════════════════════════════

describe_profile() {
    local profile="$1"
    
    case "$profile" in
        docker)
            echo -e "🐳 ${BOLD}Docker${NC} - Containerized (recommended for desktops)"
            echo "   Requirements: Docker installed and running"
            echo "   Benefits: Reproducible, no dependency management"
            echo "   Best for: Mac, Linux desktop/workstation"
            ;;
        conda)
            local variant=$(detect_conda_variant)
            echo -e "🐍 ${BOLD}Conda${NC} - Environment management (${variant})"
            echo "   Requirements: conda/mamba/micromamba installed"
            echo "   Benefits: Works on most Linux systems, no containers"
            echo "   Best for: Linux workstation, standard servers"
            ;;
        conda_server)
            echo -e "🐍 ${BOLD}Conda Server${NC} - Network-storage-safe conda"
            echo "   Requirements: conda/mamba/micromamba installed"
            echo "   Benefits: Handles NFS, network storage, HPC filesystems"
            echo "   Best for: HPC servers, problematic storage systems"
            ;;
        singularity|apptainer)
            echo -e "📦 ${BOLD}Singularity/Apptainer${NC} - HPC containers"
            echo "   Requirements: Singularity or Apptainer installed"
            echo "   Benefits: Container reproducibility on HPC"
            echo "   Best for: HPC clusters without Docker"
            ;;
        slurm,*)
            echo -e "🖥️  ${BOLD}Slurm + ${profile#slurm,}${NC} - HPC scheduler + containers"
            echo "   Requirements: Slurm access + ${profile#slurm,}"
            echo "   Benefits: Automatic job submission, resource management"
            echo "   Best for: University/research HPC clusters"
            ;;
        podman)
            echo -e "🐋 ${BOLD}Podman${NC} - Rootless containers"
            echo "   Requirements: Podman installed"
            echo "   Benefits: Like Docker but rootless"
            echo "   Best for: Systems where Docker requires root"
            ;;
        local)
            echo -e "🖥️  ${BOLD}Local${NC} - System-installed tools"
            echo "   Requirements: All tools manually installed (bowtie2, samtools, etc.)"
            echo "   Benefits: No overhead, direct tool control"
            echo "   Best for: Development, testing, custom setups"
            ;;
        *)
            echo "Custom profile: $profile"
            ;;
    esac
}

validate_profile() {
    local profile="$1"
    
    # Handle composite profiles (e.g., slurm,singularity)
    if [[ "$profile" == *","* ]]; then
        local scheduler="${profile%%,*}"
        local container="${profile#*,}"
        
        # Validate scheduler
        case "$scheduler" in
            slurm)
                if ! has_command sbatch; then
                    error "Slurm requested but sbatch not found"
                    return 1
                fi
                ;;
            pbs)
                if ! has_command qsub; then
                    error "PBS requested but qsub not found"
                    return 1
                fi
                ;;
        esac
        
        # Validate container
        return $(validate_profile "$container")
    fi
    
    # Validate single profiles
    case "$profile" in
        docker)
            if ! docker_start; then
                error "Docker profile requested but Docker unavailable"
                info "Install: https://docs.docker.com/get-docker/"
                return 1
            fi
            ;;
        conda|conda_server)
            if ! conda_available; then
                error "Conda profile requested but conda/mamba unavailable"
                info "Install: https://docs.conda.io/en/latest/miniconda.html"
                return 1
            fi
            ;;
        singularity|apptainer)
            if ! singularity_available; then
                error "Singularity/Apptainer profile requested but unavailable"
                info "Install: https://apptainer.org/docs/user/latest/quick_start.html"
                return 1
            fi
            ;;
        podman)
            if ! podman_available; then
                error "Podman profile requested but unavailable"
                info "Install: https://podman.io/getting-started/installation"
                return 1
            fi
            ;;
        local)
            warning "Local profile - ensure tools installed manually"
            ;;
        *)
            error "Unknown profile: $profile"
            return 1
            ;;
    esac
    
    return 0
}

# ═══════════════════════════════════════════════════════════════════════════
# PRE-FLIGHT CHECKS
# ═══════════════════════════════════════════════════════════════════════════

check_nextflow() {
    if ! has_command nextflow; then
        error "Nextflow not found"
        info "Install: curl -s https://get.nextflow.io | bash"
        info "Or: conda install -c bioconda nextflow"
        return 1
    fi
    
    local nf_version=$(nextflow -version 2>&1 | grep version | awk '{print $3}')
    success "Nextflow ${nf_version}"
    return 0
}

check_disk_space() {
    local required_gb=20  # Minimum recommended
    
    if [[ "$OSTYPE" == "darwin"* ]]; then
        # macOS: df output in 512-byte blocks by default, convert to GB
        local avail_blocks=$(df . | tail -1 | awk '{print $4}')
        local avail_gb=$((avail_blocks / 2097152))  # 512-byte blocks to GB
    else
        # Linux: use -BG for GB output
        local avail_gb=$(df -BG . | tail -1 | awk '{print $4}' | sed 's/G//')
    fi
    
    if [[ $avail_gb -lt $required_gb ]]; then
        warning "Low disk space: ${avail_gb} GB available (recommend ${required_gb}+ GB)"
        warning "Pipeline may fail if processing many samples"
    else
        success "Disk space: ${avail_gb} GB available"
    fi
}

check_write_permissions() {
    local test_file=".tracktx_write_test_$$"
    
    if touch "$test_file" 2>/dev/null; then
        rm -f "$test_file"
        success "Write permissions OK"
        return 0
    else
        error "Cannot write to current directory"
        return 1
    fi
}

check_input_files() {
    local samplesheet="$1"
    local params_file="$2"
    
    if [[ ! -f "$samplesheet" ]]; then
        error "Sample sheet not found: $samplesheet"
        info "Create template:"
        info "  Use TrackTx_config_generator.html or: echo 'sample,condition,timepoint,replicate,file1,file2' > samplesheet.csv"
        return 1
    fi
    success "Sample sheet: $samplesheet"
    
    # Count samples (excluding header)
    local n_samples=$(($(wc -l < "$samplesheet") - 1))
    if [[ $n_samples -le 0 ]]; then
        error "Sample sheet is empty (no samples after header)"
        return 1
    fi
    info "  Samples: ${n_samples}"
    
    if [[ -f "$params_file" ]]; then
        success "Parameters: $params_file"
    else
        info "No custom parameters (using pipeline defaults)"
    fi
    
    return 0
}

# ═══════════════════════════════════════════════════════════════════════════
# HELP & USAGE
# ═══════════════════════════════════════════════════════════════════════════

show_help() {
    cat << 'EOF'
═══════════════════════════════════════════════════════════════════════════
TrackTx Pipeline Runner
═══════════════════════════════════════════════════════════════════════════

QUICK START:
    ./run_pipeline.sh                    # Auto-detect & run
    ./run_pipeline.sh -profile docker    # Force Docker
    ./run_pipeline.sh --resume           # Resume previous run

USAGE:
    ./run_pipeline.sh [OPTIONS]

OPTIONS:
    -h, --help                      Show this help message
    -profile PROFILE               Force execution profile
    --samplesheet FILE             Sample sheet CSV (default: samplesheet.csv)
    --params-file FILE             Parameters YAML (default: params.yaml)
    --output_dir DIR               Override output directory (from params.yaml)
    --resume                       Resume previous run
    --fast, --performance          ⚡ Enable performance mode (for external drives)
    --no-auto-resume               Disable auto-resume detection
    --no-clear                     Keep terminal history visible
    --clear-delay SEC              Seconds before starting (default: 30, press Enter to skip)
    --skip-countdown               Start immediately without countdown
    --no-docker-prompt             Skip Docker auto-start prompt
    --dry-run                      Show command without executing
    --show-system-info             Display system resources and exit

PERFORMANCE MODE (--fast):
    Optimizes pipeline for external storage (USB SSD, NAS):
    • Disables scratch space (reduces I/O)
    • Increases parallel task limit
    • Uses fast work directory (~/tmp/tracktx_work)
    ⚠️  Recommended when running from external drives
    📊 Expected speedup: 2-3x faster execution

PROFILES:
    docker         🐳 Containers via Docker (recommended for desktops)
    conda          🐍 Conda/Mamba environment (standard servers)
    conda_server   🐍 Conda for network storage (HPC/problematic servers)
    singularity    📦 Singularity/Apptainer containers (HPC)
    slurm,*        🖥️  Slurm scheduler + container (HPC)
    podman         🐋 Rootless containers
    local          🖥️  System-installed tools

PROFILE AUTO-DETECTION:
    Laptop/Mac     → docker (if available)
    Workstation    → docker or conda
    Standard Server→ conda
    HPC Server     → slurm,singularity or slurm,conda_server
    Network Storage→ conda_server (automatic)

FEATURES:
    • Divergent transcription: edge-to-edge distance algorithm
    • Primary/unique mapper tracks for clean detection
    • 500bp default pairing window (optimized)
    • Automatic resource scaling (4 CPUs to 128+ CPUs)
    • Universal compatibility (laptop → HPC)

EXAMPLES:
    # Auto-detect environment (recommended)
    ./run_pipeline.sh

    # Desktop with Docker
    ./run_pipeline.sh -profile docker

    # Workstation with Conda
    ./run_pipeline.sh -profile conda

    # HPC with Slurm + Singularity
    ./run_pipeline.sh -profile slurm,singularity

    # Server with problematic network storage
    ./run_pipeline.sh -profile conda_server

    # Custom parameters
    ./run_pipeline.sh --params-file my_config.yaml

    # Resume previous run
    ./run_pipeline.sh --resume

    # Dry run to preview command
    ./run_pipeline.sh --dry-run

TROUBLESHOOTING:
    No Docker?     Use: -profile conda
    NFS issues?    Use: -profile conda_server
    HPC cluster?   Use: -profile slurm,singularity
    Permissions?   Check: ls -la . && touch test.txt

MORE HELP:
    GitHub: https://github.com/SerhatAktay/TrackTx
    Docs:   https://tracktx.readthedocs.io

═══════════════════════════════════════════════════════════════════════════
EOF
}

# ═══════════════════════════════════════════════════════════════════════════
# MAIN EXECUTION
# ═══════════════════════════════════════════════════════════════════════════

main() {
    # Parse arguments
    local PROFILE=""
    local SAMPLESHEET="samplesheet.csv"
    local PARAMS_FILE="params.yaml"
    local RESUME=""
    local DRY_RUN=""
    local NO_AUTO_RESUME=0
    local NO_CLEAR=0
    local CLEAR_DELAY=30
    local NO_DOCKER_PROMPT=0
    local WANT_PROMPT_RESUME=1
    local SHOW_SYSTEM_INFO_ONLY=0
    local PERFORMANCE_MODE=0
    local EXTRA_ARGS=()
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            -profile)
                PROFILE="$2"
                shift 2
                ;;
            --samplesheet)
                SAMPLESHEET="$2"
                shift 2
                ;;
            --params-file)
                PARAMS_FILE="$2"
                shift 2
                ;;
            -resume|--resume)
                RESUME="-resume"
                shift
                ;;
            --dry-run)
                DRY_RUN="-preview"
                shift
                ;;
            --no-auto-resume)
                NO_AUTO_RESUME=1
                shift
                ;;
            --no-prompt)
                WANT_PROMPT_RESUME=0
                shift
                ;;
            --no-clear)
                NO_CLEAR=1
                shift
                ;;
            --clear-delay)
                CLEAR_DELAY="$2"
                shift 2
                ;;
            --skip-countdown)
                NO_CLEAR=1
                shift
                ;;
            --no-docker-prompt)
                NO_DOCKER_PROMPT=1
                shift
                ;;
            --show-system-info)
                SHOW_SYSTEM_INFO_ONLY=1
                shift
                ;;
            --fast|--performance)
                PERFORMANCE_MODE=1
                shift
                ;;
            *)
                EXTRA_ARGS+=("$1")
                shift
                ;;
        esac
    done
    
    # ═══════════════════════════════════════════════════════════════════════
    # BANNER
    # ═══════════════════════════════════════════════════════════════════════
    
    echo ""
    separator
    header "  TrackTx Pipeline — Nascent RNA Analysis"
    separator
    echo ""
    
    # Show system info if requested
    if [[ $SHOW_SYSTEM_INFO_ONLY -eq 1 ]]; then
        show_system_info
        exit 0
    fi
    
    # ═══════════════════════════════════════════════════════════════════════
    # PRE-FLIGHT CHECKS
    # ═══════════════════════════════════════════════════════════════════════
    
    header "Pre-Flight Checks"
    separator
    
    check_nextflow || exit 1
    check_write_permissions || exit 1
    check_disk_space
    check_input_files "$SAMPLESHEET" "$PARAMS_FILE" || exit 1
    
    echo ""
    
    # ═══════════════════════════════════════════════════════════════════════
    # SYSTEM INFORMATION
    # ═══════════════════════════════════════════════════════════════════════
    
    show_system_info
    echo ""
    
    # ═══════════════════════════════════════════════════════════════════════
    # PROFILE DETECTION & VALIDATION
    # ═══════════════════════════════════════════════════════════════════════
    
    header "Profile Selection"
    separator
    
    if [[ -z "$PROFILE" ]]; then
        PROFILE=$(detect_best_profile)
        info "Auto-detected profile: ${BOLD}$PROFILE${NC}"
    else
        info "User-specified profile: ${BOLD}$PROFILE${NC}"
    fi
    
    echo ""
    describe_profile "$PROFILE"
    echo ""
    
    if ! validate_profile "$PROFILE"; then
        exit 1
    fi
    
    success "Profile validated: $PROFILE"
    echo ""
    
    # ═══════════════════════════════════════════════════════════════════════
    # RESUME DETECTION (Updated for UNBOUND variable safety)
    # ═══════════════════════════════════════════════════════════════════════

    # RESUME=""  <-- Removed to preserve -resume flag if set earlier
    new_extra_args=()

    # Only proceed with filtering if the array is not empty
    if ((${#EXTRA_ARGS[@]} > 0)); then
        # Loop through the existing arguments to safely find and remove "-resume"
        for arg in "${EXTRA_ARGS[@]}"; do
            if [[ "$arg" == "-resume" ]]; then
                RESUME="-resume"
            else
                new_extra_args+=("$arg")
            fi
        done
        
        # Update the main array with the filtered arguments
        EXTRA_ARGS=("${new_extra_args[@]}")
    fi

    #---

    # Auto-detect resume
    if [[ $NO_AUTO_RESUME -eq 0 && -z "$RESUME" ]]; then
        if [[ -d .nextflow || -f .nextflow.log || -f results/trace/trace.txt ]]; then
            if [[ $WANT_PROMPT_RESUME -eq 1 && -t 0 ]]; then
                echo -n "Found previous run. Resume? [Y/n]: "
                read -r -t 30 ANSWER || ANSWER="y"
                case "${ANSWER:-y}" in
                    [Yy]*) RESUME="-resume" ;;
                esac
            else
                info "Auto-detected previous run → enabling resume"
                RESUME="-resume"
            fi
        fi
    fi
    
    # ═══════════════════════════════════════════════════════════════════════
    # PERFORMANCE MODE SETUP
    # ═══════════════════════════════════════════════════════════════════════
    
    FAST_WORK_DIR=""
    if [[ $PERFORMANCE_MODE -eq 1 ]]; then
        # Create fast work directory on internal storage
        FAST_WORK_DIR="$HOME/tmp/tracktx_work"
        mkdir -p "$FAST_WORK_DIR"
        
        success "Performance mode enabled!"
        info "Using fast work directory: $FAST_WORK_DIR"
        info "This can make the pipeline 2-3x faster on external drives"
    fi
    
    # ═══════════════════════════════════════════════════════════════════════
    # BUILD COMMAND
    # ═══════════════════════════════════════════════════════════════════════
    
    CMD=(
        nextflow run main.nf
        -entry TrackTx
        -profile "$PROFILE"
        --samplesheet "$SAMPLESHEET"
    )
    
    # Note: output_dir is now controlled by params.yaml (not hardcoded!)
    # You can override it with: ./run_pipeline.sh --output_dir my_custom_results
    # Execution reports (trace, timeline, DAG) will be in ${output_dir}/trace/
    
    [[ -f "$PARAMS_FILE" ]] && CMD+=(-params-file "$PARAMS_FILE")
    [[ -n "$RESUME" ]] && CMD+=("$RESUME")
    [[ -n "$DRY_RUN" ]] && CMD+=("$DRY_RUN")
    
    # Add performance optimizations
    if [[ $PERFORMANCE_MODE -eq 1 ]]; then
        CMD+=(-work-dir "$FAST_WORK_DIR")
        [[ -f "performance.config" ]] && CMD+=(-c performance.config)
    fi

    
    [[ ${#EXTRA_ARGS[@]} -gt 0 ]] && CMD+=("${EXTRA_ARGS[@]}")
    
    # ═══════════════════════════════════════════════════════════════════════
    # EXECUTION SUMMARY
    # ═══════════════════════════════════════════════════════════════════════
    
    separator
    header "Execution Plan"
    separator
    echo -e "  Profile:      ${BOLD}$PROFILE${NC}"
    echo "  Sample sheet: $SAMPLESHEET"
    if [[ -f "$PARAMS_FILE" ]]; then
        echo "  Parameters:   $PARAMS_FILE"
        # Extract output_dir from params.yaml if present
        local OUTPUT_FROM_PARAMS=$(grep '^output_dir:' "$PARAMS_FILE" 2>/dev/null | sed 's/output_dir: *//;s/["'\'']//g' || echo "")
        if [[ -n "$OUTPUT_FROM_PARAMS" ]]; then
            echo "  Output:       ${OUTPUT_FROM_PARAMS}"
        else
            echo "  Output:       ./results/ (default)"
        fi
    else
        echo "  Output:       ./results/ (default)"
    fi
    [[ -n "$RESUME" ]] && echo -e "  Mode:         ${BOLD}Resume${NC}"
    if [[ $PERFORMANCE_MODE -eq 1 ]]; then
        echo -e "  Performance:  ${BOLD}${GREEN}⚡ FAST MODE${NC} (work: $FAST_WORK_DIR)"
    fi
    separator
    
    if [[ -n "$DRY_RUN" ]]; then
        echo ""
        info "DRY RUN MODE - Command preview:"
        echo "  ${CMD[*]}"
        exit 0
    fi
    
    # ═══════════════════════════════════════════════════════════════════════
    # EXECUTE
    # ═══════════════════════════════════════════════════════════════════════
    
    echo ""
    info "Full command:"
    echo "  ${CMD[*]}"
    echo ""
    
    # Optional terminal clear with countdown
    if [[ $NO_CLEAR -eq 0 ]]; then
        echo ""
        warning "Starting in ${CLEAR_DELAY}s (Press Enter to start now, or Ctrl+C to abort)"
        
        # Countdown with Enter to skip
        for ((i=CLEAR_DELAY; i>0; i--)); do
            # Check if Enter was pressed
            read -t 1 -n 1 -s key && break
            printf "\r  Starting in %2ds... (Press Enter to skip)" "$i"
        done
        printf "\r\033[K"  # Clear the countdown line
        echo ""
        clear 2>/dev/null || true
    fi
    
    separator
    success "Starting TrackTx Pipeline"
    separator
    echo ""
    
    exec "${CMD[@]}"
}

# ═══════════════════════════════════════════════════════════════════════════
# SIGNAL HANDLING
# ═══════════════════════════════════════════════════════════════════════════

trap 'echo ""; warning "Pipeline interrupted by user"; exit 130' INT TERM

# ═══════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════

main "$@"