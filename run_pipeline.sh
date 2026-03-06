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
#   ./run_pipeline.sh --external-drive   # When project is on exFAT/USB
#   ./run_pipeline.sh -profile docker    # Force specific profile
#   ./run_pipeline.sh --setup-completion # Install tab completion
#   ./run_pipeline.sh --help             # Show full help
#
# FEATURES:
#   • Divergent transcription: edge-to-edge distance algorithm
#   • Automatic resource optimization for any system size
#   • Smart profile detection (Docker, Conda, Singularity, HPC)
#   • Network storage compatibility detection
#   • Tab completion support (run --setup-completion once)
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
BOLD='\033[1m'
NC='\033[0m'

info()    { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[✓]${NC} $1"; }
warning() { echo -e "${YELLOW}[⚠]${NC} $1"; }
error()   { echo -e "${RED}[✗]${NC} $1"; }
header()  { echo -e "${BOLD}${CYAN}$1${NC}"; }

# Uniform 71-char separator for all banners
SEP="═══════════════════════════════════════════════════════════════════════════"
separator() { echo "$SEP"; }

# Section header: separator + title + separator (uniform block style)
section_header() {
    separator
    header "$1"
    separator
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
    local cpus
    local mem_gb
    cpus=$(detect_cpus)
    mem_gb=$(detect_memory_gb)

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
# Uses timeout to avoid hanging on stale/slow NFS mounts
is_network_storage() {
    local path="${1:-.}"

    if [[ "$OSTYPE" != "linux-gnu"* ]]; then
        return 1
    fi

    local df_cmd=""
    if has_command timeout; then
        df_cmd="timeout 3s df"
    elif has_command gtimeout; then
        df_cmd="gtimeout 3s df"
    else
        df_cmd="df"
    fi

    local fs_type
    fs_type=$($df_cmd -T "$path" 2>/dev/null | tail -1 | awk '{print $2}')
    case "$fs_type" in
        nfs|nfs4|cifs|smb|smbfs|afs|gpfs|lustre|beegfs)
            return 0
            ;;
    esac
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
    docker info >/dev/null 2>&1
}

# Check if Docker can actually run containers (catches file-sharing, permission issues)
# Uses 60s timeout to avoid hanging on slow pulls or exFAT/external drives
docker_can_run_containers() {
    local timeout_cmd=""
    if has_command timeout; then
        timeout_cmd="timeout 60s"
    elif has_command gtimeout; then
        timeout_cmd="gtimeout 60s"
    fi
    if [[ -n "$timeout_cmd" ]]; then
        $timeout_cmd docker run --rm hello-world >/dev/null 2>&1
    else
        # No timeout available (e.g. macOS without coreutils) - run in background with kill
        docker run --rm hello-world >/dev/null 2>&1 &
        local pid=$!
        local i=0
        while kill -0 "$pid" 2>/dev/null && [[ $i -lt 60 ]]; do
            sleep 1
            i=$((i + 1))
        done
        if kill -0 "$pid" 2>/dev/null; then
            kill "$pid" 2>/dev/null
            wait "$pid" 2>/dev/null
            return 1
        fi
        wait "$pid" 2>/dev/null
        return $?
    fi
}

# Detect container runtime memory limit (Docker/Podman) in GB
# Used to set NXF_HOST_MEM so Nextflow allocates correctly inside containers
detect_container_memory_gb() {
    local mem_bytes=""
    if docker_daemon_running; then
        mem_bytes=$(docker info --format '{{.MemTotal}}' 2>/dev/null)
    elif has_command podman && podman info >/dev/null 2>&1; then
        mem_bytes=$(podman info --format '{{.Host.MemTotal}}' 2>/dev/null)
    fi
    if [[ -n "$mem_bytes" && "$mem_bytes" =~ ^[0-9]+$ && "$mem_bytes" -gt 0 ]]; then
        echo $((mem_bytes / 1073741824))
    else
        echo ""
    fi
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
    local hpc_type
    hpc_type=$(detect_hpc)

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

    # Linux server: prefer conda_server over Docker since servers commonly run on NFS
    # and conda_server is more reliable in that context
    if [[ "$OSTYPE" == "linux-gnu"* ]] && conda_available && needs_conda_server_profile; then
        echo "conda_server"
        return
    fi

    # Desktop/Workstation/Server — Docker first if available
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
    local cpus
    local mem_gb
    local system_type
    cpus=$(detect_cpus)
    mem_gb=$(detect_memory_gb)
    system_type=$(categorize_system)

    section_header "System Information"
    echo "  CPUs:        ${cpus} cores"
    echo "  Memory:      ${mem_gb} GB"
    echo "  Category:    ${system_type}"
    echo ""

    # System-specific recommendations
    # NOTE: First-run times include genome download (~5 min) + index build (~20-40 min for hg38).
    # Subsequent runs skip these steps and go straight to alignment.
    case "$system_type" in
        laptop)
            echo "  Parallelism: Sequential (1 sample at a time)"
            echo "  Runtime:     First run:  ~3-5h (includes genome download + index build)"
            echo "               Subsequent: ~60-90 min per sample"
            ;;
        desktop)
            echo "  Parallelism: 2-3 samples in parallel"
            echo "  Runtime:     First run:  ~2-3h (includes genome download + index build)"
            echo "               Subsequent: ~30-60 min per batch"
            ;;
        workstation)
            echo "  Parallelism: 3-5 samples in parallel"
            echo "  Runtime:     First run:  ~1.5-2.5h (includes genome download + index build)"
            echo "               Subsequent: ~20-40 min per batch"
            ;;
        server)
            echo "  Parallelism: 5-10 samples in parallel"
            echo "  Runtime:     First run:  ~1-2h (includes genome download + index build)"
            echo "               Subsequent: ~15-30 min per batch"
            ;;
        hpc)
            echo "  Parallelism: 10+ samples in parallel"
            echo "  Runtime:     First run:  ~45-90 min (includes genome download + index build)"
            echo "               Subsequent: ~10-20 min per batch"
            ;;
    esac

    echo ""
    echo -e "  Storage:     Ensure ${BOLD}50 GB free${NC} (index ~5 GB + ~10-20 GB per sample)"

    if is_network_storage "."; then
        warning "Detected network storage — using conda_server profile for compatibility"
    fi

    separator
}

# ═══════════════════════════════════════════════════════════════════════════
# PROFILE VALIDATION & DESCRIPTION
# ═══════════════════════════════════════════════════════════════════════════

describe_profile() {
    local profile="$1"

    # Composite with performance (e.g., docker,performance from --external-drive)
    if [[ "$profile" == *",performance" ]]; then
        local base="${profile%,performance}"
        describe_profile "$base"
        echo "     + performance profile (scratch off, publish copy, increased parallelism)"
        return
    fi

    case "$profile" in
        docker)
            echo -e "  🐳 ${BOLD}Docker${NC} - Containerized (recommended for desktops)"
            echo "     Requirements: Docker installed and running"
            echo "     Benefits: Reproducible, no dependency management"
            echo "     Best for: Mac, Linux desktop/workstation"
            ;;
        conda)
            local variant
            variant=$(detect_conda_variant)
            echo -e "  🐍 ${BOLD}Conda${NC} - Environment management (${variant})"
            echo "     Requirements: conda/mamba/micromamba installed"
            echo "     Benefits: Works on most Linux systems, no containers"
            echo "     Best for: Linux workstation, standard servers (local filesystem)"
            ;;
        conda_server)
            echo -e "  🐍 ${BOLD}Conda Server${NC} - Network-storage-safe conda"
            echo "     Requirements: conda/mamba/micromamba installed"
            echo "     Benefits: Handles NFS, network storage, HPC filesystems"
            echo "     Best for: HPC servers, NFS/shared storage (your typical server setup)"
            ;;
        singularity|apptainer)
            echo -e "  📦 ${BOLD}Singularity/Apptainer${NC} - HPC containers"
            echo "     Requirements: Singularity or Apptainer installed"
            echo "     Benefits: Container reproducibility on HPC"
            echo "     Best for: HPC clusters without Docker"
            ;;
        slurm,*)
            echo -e "  🖥️  ${BOLD}Slurm + ${profile#slurm,}${NC} - HPC scheduler + containers"
            echo "     Requirements: Slurm access + ${profile#slurm,}"
            echo "     Benefits: Automatic job submission, resource management"
            echo "     Best for: University/research HPC clusters"
            ;;
        podman)
            echo -e "  🐋 ${BOLD}Podman${NC} - Rootless containers"
            echo "     Requirements: Podman installed"
            echo "     Benefits: Like Docker but rootless"
            echo "     Best for: Systems where Docker requires root"
            ;;
        local)
            echo -e "  🖥️  ${BOLD}Local${NC} - System-installed tools"
            echo "     Requirements: All tools manually installed (bowtie2, samtools, etc.)"
            echo "     Benefits: No overhead, direct tool control"
            echo "     Best for: Development, testing, custom setups"
            ;;
        *)
            echo "  Custom profile: $profile"
            ;;
    esac
}

validate_profile() {
    local profile="$1"

    # Handle composite profiles (e.g., slurm,singularity or docker,performance)
    if [[ "$profile" == *","* ]]; then
        local first="${profile%%,*}"
        local second="${profile#*,}"

        # performance is additive (used with docker/conda via --external-drive); validate base only
        if [[ "$second" == "performance" ]]; then
            validate_profile "$first"
            return $?
        fi

        # Scheduler,container pattern (e.g., slurm,singularity)
        case "$first" in
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

        validate_profile "$second"
        return $?
    fi

    # Validate single profiles
    case "$profile" in
        docker)
            if ! docker_start; then
                error "Docker profile requested but Docker unavailable"
                info "Install: https://docs.docker.com/get-docker/"
                return 1
            fi
            # Verify Docker can run containers (catches file-sharing issues on Mac)
            # Skipped by default: hello-world test can hang 60s+ on external drives
            if [[ ${SKIP_DOCKER_RUN_TEST:-1} -eq 0 ]]; then
                info "Verifying Docker can run containers (may take up to 60s)..."
                if ! docker_can_run_containers; then
                    error "Docker daemon responds but cannot run containers"
                    info "Common causes:"
                    info "  • Mac: Docker Desktop → Settings → Resources → File Sharing (add your project path)"
                    info "  • External drive: Try running from internal drive, or use: -profile conda"
                    info "  • Permissions: Ensure your user is in the 'docker' group (Linux)"
                    info "  • Skip this check: SKIP_DOCKER_RUN_TEST=1 (default: skipped)"
                    return 1
                fi
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

    local nf_version
    nf_version=$(nextflow -version 2>&1 | grep version | awk '{print $3}')
    success "Nextflow ${nf_version}"

    # Validate minimum version (requires >=24.04.0)
    local required_major=24 required_minor=4
    local actual_major actual_minor
    actual_major=$(echo "$nf_version" | cut -d. -f1)
    actual_minor=$(echo "$nf_version" | cut -d. -f2)
    if [[ -n "$actual_major" && -n "$actual_minor" ]]; then
        if [[ "$actual_major" -lt "$required_major" ]] || \
           [[ "$actual_major" -eq "$required_major" && "$actual_minor" -lt "$required_minor" ]]; then
            warning "Nextflow ${nf_version} is below recommended >=24.04.0 — pipeline may not work correctly"
            info "Upgrade: nextflow self-update"
        fi
    fi

    return 0
}

check_disk_space() {
    local n_samples="${1:-1}"
    local required_gb=$(( 50 + n_samples * 15 ))  # 50 GB base (index etc.) + 15 GB per sample
    local df_out=""

    # Use timeout to avoid hang on external/slow drives (df can block indefinitely)
    local df_tmp
    df_tmp=$(mktemp 2>/dev/null || echo "/tmp/tracktx_df_$$")
    local df_args="."
    [[ "$OSTYPE" != "darwin"* ]] && df_args="-BG ."

    if has_command timeout; then
        df_out=$(timeout 5s df $df_args 2>/dev/null | tail -1)
    elif has_command gtimeout; then
        df_out=$(gtimeout 5s df $df_args 2>/dev/null | tail -1)
    else
        # macOS: no timeout; run in background and kill if slow
        df $df_args 2>/dev/null > "$df_tmp" &
        local pid=$!
        local i=0
        while kill -0 "$pid" 2>/dev/null && [[ $i -lt 5 ]]; do
            sleep 1
            i=$((i + 1))
        done
        if kill -0 "$pid" 2>/dev/null; then
            kill "$pid" 2>/dev/null
            wait "$pid" 2>/dev/null
            df_out=""
        else
            wait "$pid" 2>/dev/null
            df_out=$(tail -1 "$df_tmp" 2>/dev/null)
        fi
        rm -f "$df_tmp" 2>/dev/null
    fi

    local avail_gb
    if [[ "$OSTYPE" == "darwin"* ]]; then
        local avail_blocks
        avail_blocks=$(echo "$df_out" | awk '{print $4}')
        avail_gb=$((avail_blocks / 2097152))  # 512-byte blocks to GB
    else
        avail_gb=$(echo "$df_out" | awk '{print $4}' | sed 's/G//')
    fi

    # Fallback if df failed or timed out
    avail_gb=${avail_gb:-0}
    if [[ ! "$avail_gb" =~ ^[0-9]+$ ]]; then
        avail_gb=0
    fi

    if [[ $avail_gb -lt $required_gb ]]; then
        warning "Low disk space: ${avail_gb} GB available (recommend ${required_gb}+ GB for ${n_samples} sample(s))"
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
    local n_samples
    n_samples=$(($(wc -l < "$samplesheet") - 1))
    if [[ $n_samples -le 0 ]]; then
        error "Sample sheet is empty (no samples after header)"
        return 1
    fi
    info "Samples: ${n_samples}"

    if [[ -f "$params_file" ]]; then
        success "Parameters: $params_file"
    else
        info "No custom parameters (using pipeline defaults)"
    fi

    # Return sample count for use in disk check
    echo "$n_samples"
    return 0
}

# ═══════════════════════════════════════════════════════════════════════════
# TAB COMPLETION SETUP
# ═══════════════════════════════════════════════════════════════════════════

setup_completion() {
    local completion_dir="${HOME}/.bash_completion.d"
    local completion_file="${completion_dir}/tracktx"

    mkdir -p "$completion_dir"

    cat > "$completion_file" << 'COMP'
_tracktx_complete() {
    local cur prev
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    local opts="--samplesheet --params-file --params_file --output_dir --output-dir
                -profile --resume --external-drive --skip-countdown
                --clear-delay --no-clear --no-auto-resume --no-resume-prompt
                --no-docker-prompt --dry-run --show-system-info
                --setup-completion --help"
    local profiles="docker conda conda_server singularity apptainer podman local
                    slurm,singularity slurm,conda_server slurm,apptainer"

    case "$prev" in
        -profile)
            COMPREPLY=( $(compgen -W "$profiles" -- "$cur") )
            return 0
            ;;
        --samplesheet)
            COMPREPLY=( $(compgen -f -X '!*.csv' -- "$cur") )
            return 0
            ;;
        --params-file|--params_file)
            COMPREPLY=( $(compgen -f -X '!*.yaml' -- "$cur") )
            COMPREPLY+=( $(compgen -f -X '!*.yml' -- "$cur") )
            return 0
            ;;
        --output_dir|--output-dir)
            COMPREPLY=( $(compgen -d -- "$cur") )
            return 0
            ;;
        --clear-delay)
            COMPREPLY=( $(compgen -W "5 10 15 30" -- "$cur") )
            return 0
            ;;
    esac

    COMPREPLY=( $(compgen -W "$opts" -- "$cur") )
}
complete -F _tracktx_complete run_pipeline.sh
COMP

    # Source from .bashrc if not already present
    local bashrc="${HOME}/.bashrc"
    local source_line="source ${completion_file}"
    if ! grep -qF "$source_line" "$bashrc" 2>/dev/null; then
        echo "" >> "$bashrc"
        echo "# TrackTx tab completion" >> "$bashrc"
        echo "$source_line" >> "$bashrc"
    fi

    # Activate immediately in current shell
    # shellcheck disable=SC1090
    source "$completion_file" 2>/dev/null || true

    success "Tab completion installed to: ${completion_file}"
    info "Restart your shell or run: source ${completion_file}"
    info "Usage: ./run_pipeline.sh --<TAB> or -profile <TAB>"
}

# ═══════════════════════════════════════════════════════════════════════════
# HELP & USAGE
# ═══════════════════════════════════════════════════════════════════════════

show_help() {
    cat << 'EOF'
═══════════════════════════════════════════════════════════════════════════
TrackTx Pipeline Runner
═══════════════════════════════════════════════════════════════════════════

GETTING STARTED (First time?):
    1. Generate config:  Open TrackTx_config_generator.html in your browser
    2. Download:         Save params.yaml and samplesheet.csv to this folder
    3. Run:              ./run_pipeline.sh
    That's it! The script auto-detects Docker/Conda and runs the pipeline.

QUICK START:
    ./run_pipeline.sh                    # Auto-detect & run
    ./run_pipeline.sh -profile docker    # Force Docker
    ./run_pipeline.sh --resume           # Resume previous run
    ./run_pipeline.sh --setup-completion # Enable tab completion (run once)

USAGE:
    ./run_pipeline.sh [OPTIONS]

OPTIONS:
  Input & output:
    --samplesheet FILE             Sample sheet CSV (default: samplesheet.csv)
    --params-file FILE             Parameters YAML (default: params.yaml)
    --params_file FILE             Same as --params-file
    --output_dir DIR               Override output directory (from params.yaml)
    --output-dir DIR               Same as --output_dir

  Execution:
    -profile PROFILE               Force profile: docker, conda, conda_server, singularity,
                                   slurm,singularity, podman, local
    --resume                       Resume from last run (reuse completed tasks)
    --external-drive               Use when project is on exFAT/USB/NAS (fixes publish errors)

  Prompts & timing:
    --skip-countdown               Start immediately (no countdown)
    --clear-delay SEC              Countdown seconds before start (default: 5)
    --no-clear                     Keep terminal history (no clear screen)
    --no-auto-resume               Don't auto-detect or prompt for resume
    --no-resume-prompt             Don't prompt to resume, auto-resume if detected
    --no-docker-prompt             Don't prompt to start Docker if not running

  Other:
    -h, --help                     Show this help message
    --dry-run                      Preview command without executing
    --show-system-info             Show CPU, memory, and exit
    --setup-completion             Install bash tab completion (run once)

  Extra args after options are passed to Nextflow (e.g. --reference_genome hg38)

EXTERNAL DRIVE MODE (--external-drive):
    For USB, exFAT, or network drives:
    • Fixes "Failed to publish file [link]" error (uses copy instead of hard links)
    • Fixes OverlappingFileLockException (cache, temp, work on local ~/tmp/tracktx_*)
    • Results stay on project dir; work dir ~10–50 GB on local (exFAT lacks file locking)
    • Disables scratch space, increases parallelism for slow I/O

PROFILES:
    docker         🐳 Containers via Docker (recommended for desktops/Mac)
    conda          🐍 Conda/Mamba environment (local filesystem servers)
    conda_server   🐍 Conda for network storage (NFS/HPC — recommended for servers)
    singularity    📦 Singularity/Apptainer containers (HPC)
    slurm,*        🖥️  Slurm scheduler + container (HPC)
    podman         🐋 Rootless containers
    local          🖥️  System-installed tools

PROFILE AUTO-DETECTION:
    Mac/Desktop    → docker (if available)
    Linux Server   → conda_server (if on NFS), else conda or docker
    HPC            → slurm,singularity or slurm,conda_server

RUNTIME ESTIMATES:
    These include genome download + index build on first run (~30-45 min for hg38).
    Subsequent runs reuse the cached index and are significantly faster.

    Laptop         First run: ~3-5h  |  Subsequent: ~60-90 min/sample
    Desktop        First run: ~2-3h  |  Subsequent: ~30-60 min/batch
    Workstation    First run: ~1.5h  |  Subsequent: ~20-40 min/batch
    Server         First run: ~1-2h  |  Subsequent: ~15-30 min/batch
    HPC            First run: ~45min |  Subsequent: ~10-20 min/batch

FEATURES:
    • Divergent transcription: edge-to-edge distance algorithm
    • Primary/unique mapper tracks for clean detection
    • 500bp default pairing window (optimized)
    • Automatic resource scaling (4 CPUs to 128+ CPUs)
    • Universal compatibility (laptop → HPC)

EXAMPLES:
    # Basic usage
    ./run_pipeline.sh
    ./run_pipeline.sh -profile docker
    ./run_pipeline.sh -profile conda_server

    # Custom input files
    ./run_pipeline.sh --samplesheet my_samples.csv
    ./run_pipeline.sh --params-file my_params.yaml
    ./run_pipeline.sh --samplesheet data.csv --params-file config.yaml

    # Output and resume
    ./run_pipeline.sh --output_dir my_results
    ./run_pipeline.sh --resume
    ./run_pipeline.sh --output_dir run2 --resume

    # USB/exFAT/external drive
    ./run_pipeline.sh --external-drive

    # Skip prompts and countdown
    ./run_pipeline.sh --skip-countdown
    ./run_pipeline.sh --no-auto-resume

    # HPC and servers
    ./run_pipeline.sh -profile slurm,singularity
    ./run_pipeline.sh -profile conda_server

    # Preview and info
    ./run_pipeline.sh --dry-run
    ./run_pipeline.sh --show-system-info

    # Enable tab completion (once)
    ./run_pipeline.sh --setup-completion

ENVIRONMENT VARIABLES:
    TRACKTX_SKIP_PULL=1           Skip Docker image pull (offline, slow network)
    NXF_HOST_MEM=N                Tell Nextflow available RAM in GB (e.g. Docker limit)
    NXF_HOST_CPUS=N               Tell Nextflow available CPUs
    NXF_WORK=DIR                  Override work directory

TROUBLESHOOTING:
    No Docker?     Use: -profile conda_server
    NFS issues?    Use: -profile conda_server
    HPC cluster?   Use: -profile slurm,singularity
    Permissions?   Check: ls -la . && touch test.txt
    Mac M1/M2/M3?  Use multi-arch image (see envs/tracktx.yaml build instructions)
    External drive? Use: --external-drive or -profile conda
    Skip Docker test?  SKIP_DOCKER_RUN_TEST=1 ./run_pipeline.sh

MORE HELP:
    GitHub: https://github.com/serhataktay/tracktx
    Issues: https://github.com/serhataktay/tracktx/issues

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
    local CLEAR_DELAY=5
    local NO_DOCKER_PROMPT=0
    local WANT_PROMPT_RESUME=1
    local SHOW_SYSTEM_INFO_ONLY=0
    local EXTERNAL_DRIVE_MODE=0
    local OUTPUT_DIR_OVERRIDE=""
    local EXTRA_ARGS=()

    # Auto-disable countdown in non-interactive mode (e.g. cron, CI)
    if [[ ! -t 1 ]]; then
        NO_CLEAR=1
    fi

    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                show_help
                exit 0
                ;;
            --setup-completion)
                setup_completion
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
            --params-file|--params_file)
                PARAMS_FILE="$2"
                shift 2
                ;;
            -resume|--resume)
                RESUME="-resume"
                shift
                ;;
            --dry-run)
                # Maps to Nextflow's -preview flag
                DRY_RUN="-preview"
                shift
                ;;
            --no-auto-resume)
                NO_AUTO_RESUME=1
                shift
                ;;
            --no-resume-prompt)
                # Auto-resume if detected, but don't prompt
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
            --output_dir|--output-dir)
                OUTPUT_DIR_OVERRIDE="$2"
                shift 2
                ;;
            --show-system-info)
                SHOW_SYSTEM_INFO_ONLY=1
                shift
                ;;
            --external-drive|--fast|--performance|--exfat)
                EXTERNAL_DRIVE_MODE=1
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
    section_header "TrackTx Pipeline — PRO-seq Nascent RNA Analysis"
    echo ""

    # Show system info if requested
    if [[ $SHOW_SYSTEM_INFO_ONLY -eq 1 ]]; then
        show_system_info
        exit 0
    fi

    # ═══════════════════════════════════════════════════════════════════════
    # PRE-FLIGHT CHECKS
    # ═══════════════════════════════════════════════════════════════════════

    section_header "Pre-Flight Checks"

    check_nextflow || exit 1
    check_write_permissions || exit 1

    # check_input_files prints sample count on stdout; capture it
    local N_SAMPLES
    N_SAMPLES=$(check_input_files "$SAMPLESHEET" "$PARAMS_FILE") || exit 1

    check_disk_space "$N_SAMPLES"

    echo ""

    # ═══════════════════════════════════════════════════════════════════════
    # SYSTEM INFORMATION
    # ═══════════════════════════════════════════════════════════════════════

    show_system_info
    echo ""

    # ═══════════════════════════════════════════════════════════════════════
    # PROFILE DETECTION & VALIDATION
    # ═══════════════════════════════════════════════════════════════════════

    section_header "Profile Selection"

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
    # NXF_TEMP: Use project dir for Nextflow staging (avoids filling boot drive
    # when project/work are on USB or external drive; Nextflow defaults to /var/folders)
    # NOTE: In --external-drive mode this is overridden to local (exFAT lacks file locking)
    # ═══════════════════════════════════════════════════════════════════════
    if [[ -z "${NXF_TEMP:-}" ]]; then
        export NXF_TEMP="$(pwd)/.nxf_temp"
        mkdir -p "$NXF_TEMP"
    fi

    # ═══════════════════════════════════════════════════════════════════════
    # DOCKER/PODMAN: Auto-set NXF_HOST_MEM from container limit
    # Nextflow detects host RAM, but containers see Docker's limit—causing OOM.
    # ═══════════════════════════════════════════════════════════════════════
    if [[ -z "${NXF_HOST_MEM:-}" ]] && [[ "$PROFILE" == *docker* || "$PROFILE" == *podman* ]]; then
        local container_mem
        container_mem=$(detect_container_memory_gb)
        if [[ -n "$container_mem" && "$container_mem" -ge 2 ]]; then
            export NXF_HOST_MEM=$container_mem
            info "Docker/Podman memory limit: ${container_mem} GB → NXF_HOST_MEM=${container_mem}"
        fi
    fi

    # ═══════════════════════════════════════════════════════════════════════
    # DOCKER/PODMAN: Pull latest image so git pull → fresh image
    # Ensures users get the updated container after pulling pipeline changes.
    # Set TRACKTX_SKIP_PULL=1 to skip (e.g. offline, slow network).
    # Image version is read from nextflow.config to avoid hardcoding here.
    # ═══════════════════════════════════════════════════════════════════════
    local NF_VERSION
    NF_VERSION=$(grep "^\s*version\s*=" nextflow.config 2>/dev/null | head -1 | grep -oP "[\d.]+" || echo "3.0")
    local TRACKTX_IMAGE="ghcr.io/serhataktay/tracktx:${NF_VERSION}"

    if [[ "${TRACKTX_SKIP_PULL:-0}" -eq 0 ]]; then
        if [[ "$PROFILE" == *docker* ]] && has_command docker; then
            section_header "Docker Image"
            echo -e "  Image:   ${BOLD}${TRACKTX_IMAGE}${NC}"
            echo ""
            info "Pulling image (may take a few minutes on first run)..."
            if docker pull "$TRACKTX_IMAGE"; then
                success "Image ready"
            else
                warning "Pull failed — using cached image (pipeline will still run if cached)"
            fi
            separator
            echo ""
        elif [[ "$PROFILE" == *podman* ]] && has_command podman; then
            section_header "Podman Image"
            echo -e "  Image:   ${BOLD}${TRACKTX_IMAGE}${NC}"
            echo ""
            info "Pulling image (may take a few minutes on first run)..."
            if podman pull "$TRACKTX_IMAGE"; then
                success "Image ready"
            else
                warning "Pull failed — using cached image (pipeline will still run if cached)"
            fi
            separator
            echo ""
        fi
    fi

    # ═══════════════════════════════════════════════════════════════════════
    # RESUME DETECTION
    # ═══════════════════════════════════════════════════════════════════════

    local new_extra_args=()

    # Extract -resume from EXTRA_ARGS if passed through (normalize to our RESUME var)
    if ((${#EXTRA_ARGS[@]} > 0)); then
        for arg in "${EXTRA_ARGS[@]}"; do
            if [[ "$arg" == "-resume" ]]; then
                RESUME="-resume"
            else
                new_extra_args+=("$arg")
            fi
        done
        EXTRA_ARGS=("${new_extra_args[@]}")
    fi

    # Auto-detect resume: check for previous run artifacts
    # Only triggers resume if at least one task completed successfully
    if [[ $NO_AUTO_RESUME -eq 0 && -z "$RESUME" ]]; then
        local trace_path="results/trace/trace.txt"
        if [[ -n "$OUTPUT_DIR_OVERRIDE" ]]; then
            trace_path="${OUTPUT_DIR_OVERRIDE}/trace/trace.txt"
        elif [[ -f "$PARAMS_FILE" ]]; then
            local out_from_params
            out_from_params=$(awk -F': ' '/^output_dir:/{gsub(/["\047]/, "", $2); print $2}' "$PARAMS_FILE" 2>/dev/null || echo "")
            [[ -n "$out_from_params" ]] && trace_path="${out_from_params}/trace/trace.txt"
        fi

        # Check both that a previous run exists AND that at least one task succeeded
        local has_prior_run=0
        if [[ -f "$trace_path" ]] && grep -q "COMPLETED" "$trace_path" 2>/dev/null; then
            has_prior_run=1
        elif [[ -d .nextflow ]] && find work -name ".exitcode" -exec grep -lx "0" {} \; 2>/dev/null | grep -q .; then
            has_prior_run=1
        fi

        if [[ $has_prior_run -eq 1 ]]; then
            if [[ $WANT_PROMPT_RESUME -eq 1 && -t 0 ]]; then
                echo -n "Found previous run with completed tasks. Resume? [Y/n]: "
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
    # EXTERNAL DRIVE MODE SETUP
    # exFAT/USB/NFS lack file locking. Nextflow cache, temp, and work dir all need it.
    # Redirect ALL lock-using dirs to local. Results stay on project dir (external).
    # ═══════════════════════════════════════════════════════════════════════

    if [[ $EXTERNAL_DRIVE_MODE -eq 1 ]]; then
        local TRACKTX_CACHE="${HOME}/tmp/tracktx_cache"
        local proj_name
        proj_name=$(basename "$(pwd)" 2>/dev/null | sed 's/[^a-zA-Z0-9_.-]/_/g')
        proj_name=${proj_name:-default}
        local TRACKTX_WORK="${HOME}/tmp/tracktx_work/${proj_name}"
        mkdir -p "$TRACKTX_CACHE" "$TRACKTX_WORK"
        export NXF_CACHE_DIR="$TRACKTX_CACHE"
        export NXF_TEMP="${TRACKTX_CACHE}/.nxf_temp"
        mkdir -p "$NXF_TEMP"
        export NXF_WORK="$TRACKTX_WORK"
        success "External drive mode: cache, temp, work on local (~10–50 GB); results on project dir"
    fi

    # ═══════════════════════════════════════════════════════════════════════
    # BUILD COMMAND
    # ═══════════════════════════════════════════════════════════════════════

    # External drive mode: append performance profile (scratch off, publish_mode copy)
    [[ $EXTERNAL_DRIVE_MODE -eq 1 ]] && PROFILE="${PROFILE},performance"

    local CMD=(
        nextflow run main.nf
        -entry TrackTx
        -profile "$PROFILE"
        --samplesheet "$SAMPLESHEET"
    )

    # Note: output_dir is controlled by params.yaml (not hardcoded here)
    # Override with: ./run_pipeline.sh --output_dir my_custom_results
    # Execution reports (trace, timeline, DAG) will be in ${output_dir}/trace/

    [[ -f "$PARAMS_FILE" ]] && CMD+=(-params-file "$PARAMS_FILE")
    [[ -n "$OUTPUT_DIR_OVERRIDE" ]] && CMD+=(--output_dir "$OUTPUT_DIR_OVERRIDE")
    [[ -n "$RESUME" ]] && CMD+=("$RESUME")
    [[ -n "$DRY_RUN" ]] && CMD+=("$DRY_RUN")
    [[ ${#EXTRA_ARGS[@]} -gt 0 ]] && CMD+=("${EXTRA_ARGS[@]}")

    # ═══════════════════════════════════════════════════════════════════════
    # EXECUTION SUMMARY
    # ═══════════════════════════════════════════════════════════════════════

    section_header "Execution Plan"
    echo -e "  Profile:      ${BOLD}$PROFILE${NC}"
    echo "  Sample sheet: $SAMPLESHEET (${N_SAMPLES} sample(s))"
    [[ -f "$PARAMS_FILE" ]] && echo "  Parameters:   $PARAMS_FILE"

    if [[ -n "$OUTPUT_DIR_OVERRIDE" ]]; then
        echo "  Output:       ${OUTPUT_DIR_OVERRIDE} (from --output_dir)"
    elif [[ -f "$PARAMS_FILE" ]]; then
        local OUTPUT_FROM_PARAMS
        OUTPUT_FROM_PARAMS=$(awk -F': ' '/^output_dir:/{gsub(/["\047]/, "", $2); print $2}' "$PARAMS_FILE" 2>/dev/null || echo "")
        if [[ -n "$OUTPUT_FROM_PARAMS" ]]; then
            echo "  Output:       ${OUTPUT_FROM_PARAMS}"
        else
            echo "  Output:       ./results/ (default)"
        fi
    else
        echo "  Output:       ./results/ (default)"
    fi

    [[ -n "$RESUME" ]] && echo -e "  Mode:         ${BOLD}Resume${NC}"
    if [[ $EXTERNAL_DRIVE_MODE -eq 1 ]]; then
        echo -e "  External:     ${BOLD}${GREEN}enabled${NC} (cache/temp/work: local; results: project dir)"
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

    # Optional countdown before starting (skipped in non-interactive mode)
    if [[ $NO_CLEAR -eq 0 ]]; then
        echo ""
        warning "Starting in ${CLEAR_DELAY}s (Press Enter to start now, or Ctrl+C to abort)"

        for ((i=CLEAR_DELAY; i>0; i--)); do
            read -t 1 -n 1 -s && break || true
            printf "\r  Starting in %2ds... (Press Enter to skip)" "$i"
        done
        printf "\r\033[K"
        echo ""
        # NOTE: we intentionally do NOT clear the screen here — the execution
        # plan and command printed above are useful context to keep visible.
    fi

    separator
    success "Starting TrackTx Pipeline"
    separator
    echo ""

    # Record start time so elapsed time is visible in terminal scrollback
    info "Started at: $(date '+%Y-%m-%d %H:%M:%S')"
    echo ""

    "${CMD[@]}"
    local exit_code=$?

    echo ""
    if [[ $exit_code -eq 0 ]]; then
        success "Pipeline completed at: $(date '+%Y-%m-%d %H:%M:%S')"
    else
        error "Pipeline exited with code ${exit_code} at: $(date '+%Y-%m-%d %H:%M:%S')"
        info "Check .nextflow.log for details"
    fi

    return $exit_code
}

# ═══════════════════════════════════════════════════════════════════════════
# SIGNAL HANDLING
# ═══════════════════════════════════════════════════════════════════════════

trap 'echo ""; warning "Pipeline interrupted by user"; exit 130' INT TERM

# ═══════════════════════════════════════════════════════════════════════════
# ENTRY POINT
# ═══════════════════════════════════════════════════════════════════════════

main "$@"
