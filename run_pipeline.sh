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

# Check if Docker is running
docker_works() {
    has_command docker && docker info >/dev/null 2>&1
}

# Detect best execution profile
detect_profile() {
    if docker_works; then
        echo "docker"
    elif has_command conda || has_command mamba; then
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
    --dry-run                      Show what would be executed
    
PROFILES:
    docker     🐳 Everything included, fastest (recommended)
    conda      🐍 Automatic environment setup
    local      🖥️  Use system-installed tools

For more help: https://github.com/your-repo/TrackTx/
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
            conda)
            if ! has_command conda && ! has_command mamba; then
                error "Conda profile requested but conda/mamba not available"
                info "Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
                    exit 1
                fi
            success "Conda profile validated"
                ;;
            local)
            warning "Local profile - ensure all tools are installed manually"
                ;;
            *)
            error "Invalid profile: $PROFILE (use: docker, conda, or local)"
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
        --output_dir "./results"
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
