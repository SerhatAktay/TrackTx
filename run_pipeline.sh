#!/bin/bash
# TrackTx Pipeline Runner - Simple & Smart
# Automatically detects environment and runs with optimal settings

set -euo pipefail

# Colors for pretty output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

info() { echo -e "${BLUE}[INFO]${NC} $1"; }
success() { echo -e "${GREEN}[SUCCESS]${NC} $1"; }
warning() { echo -e "${YELLOW}[WARNING]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; }

# Check if command exists
has_command() { command -v "$1" >/dev/null 2>&1; }

# Check if Docker daemon is running (silent check for detection)
docker_daemon_running() {
    # If Docker command doesn't exist, return false
    if ! has_command docker; then
        return 1
    fi
    
    # Test Docker daemon with timeout to avoid hanging
    timeout 10s docker info >/dev/null 2>&1
}

# Handle Docker startup with user interaction
docker_works() {
    # First check if Docker is already running
    if docker_daemon_running; then
        return 0
    fi
    
    # Docker installed but not running - handle based on environment
    echo ""
    warning "Docker is installed but not currently running"
    
    # Check if we're in an interactive terminal or Docker prompt is disabled
    if [[ ! -t 0 ]] || [[ $NO_DOCKER_PROMPT -eq 1 ]]; then
        info "Non-interactive mode detected. Please start Docker manually:"
        if [[ "$OSTYPE" == "darwin"* ]] && [[ -d "/Applications/Docker.app" ]]; then
            info "  open -a Docker"
        elif [[ "$OSTYPE" == "linux-gnu"* ]]; then
            info "  sudo systemctl start docker"
        fi
        return 1
    fi
    
    # Check if we're on macOS with Docker Desktop
    if [[ "$OSTYPE" == "darwin"* ]] && [[ -d "/Applications/Docker.app" ]]; then
        # Use a longer timeout for the prompt (60 seconds)
        echo -n "Would you like to start Docker Desktop? [Y/n]: "
        if ! timeout 60s read -r ANSWER; then
            echo ""
            info "Prompt timeout - assuming 'no'"
            ANSWER="n"
        fi
        
        case "${ANSWER:-y}" in
            [Yy]|[Yy][Ee][Ss])
                info "Starting Docker Desktop..."
                open -a Docker
                
                # Wait for Docker daemon to start (max 60 seconds)
                info "Waiting for Docker daemon to start (this may take up to 60 seconds)..."
                for i in {1..30}; do
                    sleep 2
                    if timeout 5s docker info >/dev/null 2>&1; then
                        success "Docker started successfully!"
                        return 0
                    fi
                    echo -ne "Waiting... ${i}/30\r"
                done
                echo ""
                error "Docker failed to start within 60 seconds"
                info "Please start Docker Desktop manually and try again"
                return 1
                ;;
            *)
                info "Please start Docker Desktop manually and try again"
                return 1
                ;;
        esac
    else
        # Linux or other OS - provide instructions
        info "Please start the Docker daemon:"
        if [[ "$OSTYPE" == "linux-gnu"* ]]; then
            info "  sudo systemctl start docker"
        fi
        info "Then run this script again"
        return 1
    fi
}

# Detect best execution profile
detect_profile() {
    if docker_daemon_running; then
        echo "docker"
    # Use conda profile (works with conda/mamba/micromamba)
    elif has_command conda || has_command mamba || has_command micromamba; then
        echo "conda"
    else
        echo "local"
    fi
}

# Show usage
show_help() {
    cat << EOF
TrackTx Pipeline Runner

USAGE:
    $0 [OPTIONS]

EXAMPLES:
    $0                              # Auto-detect environment and run
    $0 -profile docker              # Force Docker profile
    $0 --samplesheet my_samples.csv # Use custom sample sheet
    $0 --params-file my_params.yaml # Use custom parameters

OPTIONS:
    -h, --help                      Show this help
    -profile PROFILE               Force profile: docker, conda, or local
    --samplesheet FILE             Sample sheet CSV (default: samplesheet.csv)
    --params-file FILE             Parameters YAML (default: params.yaml) 
    --resume                       Resume previous run
    --no-auto-resume               Disable auto-detection of resume
    --no-clear                     Do not clear the terminal before starting
    --clear-delay SEC              Seconds to wait before clearing (default: 6)
    --no-docker-prompt             Skip Docker auto-start prompt (for automated environments)
    --dry-run                      Show what would be executed
    
PROFILES:
    docker       🐳 Everything included, fastest (recommended)
    conda        🐍 Automatic environment setup (works with conda/mamba/micromamba)
    conda_server 🐍 Conda with network-storage-safe settings (for HPC/servers)
    singularity  📦 Container via Singularity/Apptainer (HPC)
    local        🖥️  Use system-installed tools

For more help: https://github.com/SerhatAktay/TrackTx/
EOF
}

# Parse command line arguments
PROFILE=""
SAMPLESHEET="samplesheet.csv"
PARAMS_FILE="params.yaml"
RESUME=""
DRY_RUN=""
NO_AUTO_RESUME=0
NO_CLEAR=0
CLEAR_DELAY=6
NO_DOCKER_PROMPT=0
EXTRA_ARGS=()
WANT_PROMPT_RESUME=1
    
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
        -resume)
            # Accept Nextflow-style single-dash resume for convenience
            RESUME="-resume"
            shift
            ;;
        --samplesheet)
            SAMPLESHEET="$2"
            shift 2
            ;;
        --params-file)
            PARAMS_FILE="$2"
                shift 2
                ;;
            --resume)
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
            # Disable interactive prompts (always auto-resume if detected)
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
        --no-docker-prompt)
            NO_DOCKER_PROMPT=1
            shift
            ;;
        *)
            EXTRA_ARGS+=("$1")
                shift
                ;;
        esac
    done

# Main execution
main() {
    info "TrackTx Pipeline - Nascent RNA Analysis"
        echo ""
        
    # Check prerequisites
    info "Checking prerequisites..."
    if ! has_command nextflow; then
        error "Nextflow not found. Please install: https://nextflow.io/docs/latest/getstarted.html"
                exit 1
    fi
    success "Prerequisites OK"
    
    # Check input files
    info "Checking input files..."
    if [[ ! -f "$SAMPLESHEET" ]]; then
        error "Sample sheet not found: $SAMPLESHEET"
        info "Create it with: echo 'sample_id,condition,timepoint,replicate,file' > $SAMPLESHEET"
            exit 1
        fi
    success "Found sample sheet: $SAMPLESHEET"
    
    if [[ -f "$PARAMS_FILE" ]]; then
        success "Found parameters: $PARAMS_FILE"
    else
        info "No custom parameters found (using defaults)"
    fi
    
    # Detect or validate profile
    if [[ -z "$PROFILE" ]]; then
        PROFILE=$(detect_profile)
        info "Auto-detected profile: $PROFILE"
    else
        info "Using specified profile: $PROFILE"
    fi
    
    # Validate profile
    case $PROFILE in
            docker)
            if ! docker_works; then
                error "Docker profile requested but Docker not available"
                info "Install Docker: https://docs.docker.com/get-docker/"
                    exit 1
                fi
            success "Docker profile validated"
                ;;
            conda|conda_server)
            if ! has_command conda && ! has_command mamba && ! has_command micromamba; then
                error "Conda profile requested but conda/mamba/micromamba not available"
                info "Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
                    exit 1
                fi
            if [[ "$PROFILE" == "conda_server" ]]; then
                success "Conda server profile validated (network-storage-safe)"
            else
                success "Conda profile validated"
            fi
                ;;
            singularity|apptainer)
            if ! has_command singularity && ! has_command apptainer; then
                error "Singularity/Apptainer profile requested but not available"
                info "Install Apptainer: https://apptainer.org/docs/user/latest/quick_start.html"
                    exit 1
                fi
            success "Singularity/Apptainer profile validated"
                ;;
            local)
            warning "Local profile - ensure all tools are installed manually"
                ;;
            *)
            error "Invalid profile: $PROFILE (use: docker, conda, conda_server, singularity, or local)"
                exit 1
                ;;
        esac
    
    # Normalize/avoid duplicate -resume if user passed it via EXTRA_ARGS
    if [[ ${#EXTRA_ARGS[@]} -gt 0 ]] && [[ " ${EXTRA_ARGS[*]} " == *" -resume "* ]]; then
        RESUME="-resume"
        # Remove any duplicate -resume in EXTRA_ARGS to avoid passing twice
        FILTERED_ARGS=()
        for a in "${EXTRA_ARGS[@]}"; do
            [[ "$a" == "-resume" ]] && continue
            FILTERED_ARGS+=("$a")
        done
        EXTRA_ARGS=("${FILTERED_ARGS[@]}")
    fi

    # Conda activation hardening (for conda/conda_server profiles)
    if [[ "$PROFILE" == "conda" || "$PROFILE" == "conda_server" ]]; then
        # Prefer user-provided prebuilt env if available
        if [[ -n "${NXF_CONDA_ENV_PATH:-}" && -d "${NXF_CONDA_ENV_PATH}/bin" ]]; then
            export PATH="${NXF_CONDA_ENV_PATH}/bin:${PATH}"
            export CONDA_PREFIX="${NXF_CONDA_ENV_PATH}"
        fi
        # Help Nextflow find conda frontend
        if [[ -z "${CONDA_EXE:-}" ]] && has_command conda; then
            export CONDA_EXE="$(command -v conda)"
        fi
        # Derive a reliable activate path
        if [[ -z "${NXF_CONDA_ACTIVATE:-}" ]]; then
            if [[ -n "${CONDA_EXE:-}" ]]; then
                _conda_bin_dir="$(dirname "${CONDA_EXE}")"
                if [[ -f "${_conda_bin_dir}/activate" ]]; then
                    export NXF_CONDA_ACTIVATE="${_conda_bin_dir}/activate"
                fi
            fi
            for a in \
                "${HOME}/miniconda3/bin/activate" \
                "${HOME}/anaconda3/bin/activate" \
                "/opt/conda/bin/activate" \
                "/usr/local/miniconda3/bin/activate" \
                "/usr/local/anaconda3/bin/activate"; do
                if [[ -z "${NXF_CONDA_ACTIVATE:-}" && -f "$a" ]]; then
                    export NXF_CONDA_ACTIVATE="$a"
                fi
            done
        fi
        export NXF_CONDA_SHELL=bash
    fi

    # Auto-detect resume (unless disabled or explicitly provided)
    if [[ $NO_AUTO_RESUME -eq 0 && -z "$RESUME" ]]; then
        if [[ -d .nextflow || -f .nextflow.log || -f results/trace/trace.txt ]]; then
            if [[ $WANT_PROMPT_RESUME -eq 1 && -t 0 ]]; then
                # Interactive prompt with sane default (Yes)
                echo -n "Found previous run artifacts. Resume? [Y/n]: "
                read -r ANSWER || ANSWER="y"
                case "${ANSWER:-y}" in
                    [Yy]|[Yy][Ee][Ss]) RESUME="-resume" ;;
                    *) RESUME="" ;;
                esac
            else
                info "Auto-detected previous run artifacts → enabling -resume"
                RESUME="-resume"
            fi
        fi
    fi

    # Build command
    CMD=(
        nextflow run main.nf
        -entry TrackTx
        -profile "$PROFILE"
        --samplesheet "$SAMPLESHEET"
        -with-report       "results/nf_report.html"
        -with-timeline     "results/nf_timeline.html"
        -with-dag          "results/nf_dag.png"
        -with-trace        "results/trace/trace.txt"
    )
    
    # Add optional parameters
    [[ -f "$PARAMS_FILE" ]] && CMD+=(-params-file "$PARAMS_FILE")
    [[ -n "$RESUME" ]] && CMD+=("$RESUME")
    [[ -n "$DRY_RUN" ]] && CMD+=("$DRY_RUN")
    [[ ${#EXTRA_ARGS[@]} -gt 0 ]] && CMD+=("${EXTRA_ARGS[@]}")
    
    # Show what we're doing
        echo ""
    info "Running TrackTx pipeline..."
    info "Profile: $PROFILE"
    info "Sample sheet: $SAMPLESHEET"
    [[ -f "$PARAMS_FILE" ]] && info "Parameters: $PARAMS_FILE"
    [[ -n "$RESUME" ]] && info "Resuming previous run"
    
    if [[ -n "$DRY_RUN" ]]; then
        info "DRY RUN - Command preview:"
        echo "${CMD[*]}"
        exit 0
    fi
    
    # Execute
    echo ""
    info "Command: ${CMD[*]}"
    echo ""

    # Controlled clear of terminal to reduce noise, unless opted out
    if [[ $NO_CLEAR -eq 0 ]]; then
        warning "Clearing terminal in ${CLEAR_DELAY}s to avoid clutter. Press Ctrl+C to abort."
        for ((i=CLEAR_DELAY;i>0;i--)); do
            echo -ne "Clearing in ${i}s...\r"; sleep 1;
        done
        echo ""; clear || true
    fi
    
    exec "${CMD[@]}"
}

# Handle interrupts gracefully
trap 'echo ""; warning "Pipeline interrupted by user"; exit 130' INT TERM

# Run main function
main "$@"
