#!/bin/bash
# ════════════════════════════════════════════════════════════════════════════
# TrackTx Server Setup Script
# Optimizes a server environment for running the TrackTx pipeline
# ════════════════════════════════════════════════════════════════════════════

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SETUP_USER="${SETUP_USER:-$(whoami)}"
INSTALL_DIR="${INSTALL_DIR:-/opt/tracktx}"

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

# Check if running as root
check_root() {
    if [[ $EUID -eq 0 ]]; then
        log_error "This script should not be run as root for security reasons"
        log_info "Please run as a regular user with sudo privileges"
        exit 1
    fi
}

# Check system requirements
check_requirements() {
    log_info "Checking system requirements..."
    
    # Check available memory
    local mem_gb=$(free -g | awk '/^Mem:/{print $2}')
    if [[ $mem_gb -lt 32 ]]; then
        log_warning "System has ${mem_gb}GB RAM. TrackTx recommends at least 32GB for optimal performance"
    else
        log_success "System has ${mem_gb}GB RAM - sufficient for TrackTx"
    fi
    
    # Check available disk space
    local disk_gb=$(df -BG . | awk 'NR==2{print $4}' | sed 's/G//')
    if [[ $disk_gb -lt 100 ]]; then
        log_warning "Available disk space: ${disk_gb}GB. TrackTx recommends at least 100GB"
    else
        log_success "Available disk space: ${disk_gb}GB - sufficient for TrackTx"
    fi
    
    # Check CPU cores
    local cores=$(nproc)
    if [[ $cores -lt 8 ]]; then
        log_warning "System has ${cores} CPU cores. TrackTx recommends at least 8 cores for optimal performance"
    else
        log_success "System has ${cores} CPU cores - sufficient for TrackTx"
    fi
}

# Install system dependencies
install_dependencies() {
    log_info "Installing system dependencies..."
    
    # Detect package manager
    if command -v apt-get >/dev/null 2>&1; then
        PKG_MGR="apt"
        sudo apt-get update
        sudo apt-get install -y \
            curl \
            wget \
            git \
            build-essential \
            ca-certificates \
            gnupg \
            lsb-release \
            software-properties-common \
            unzip \
            bzip2 \
            gzip \
            tar \
            rsync
    elif command -v yum >/dev/null 2>&1; then
        PKG_MGR="yum"
        sudo yum update -y
        sudo yum install -y \
            curl \
            wget \
            git \
            gcc \
            gcc-c++ \
            make \
            ca-certificates \
            gnupg2 \
            unzip \
            bzip2 \
            gzip \
            tar \
            rsync
    elif command -v dnf >/dev/null 2>&1; then
        PKG_MGR="dnf"
        sudo dnf update -y
        sudo dnf install -y \
            curl \
            wget \
            git \
            gcc \
            gcc-c++ \
            make \
            ca-certificates \
            gnupg2 \
            unzip \
            bzip2 \
            gzip \
            tar \
            rsync
    else
        log_error "Unsupported package manager. Please install dependencies manually"
        exit 1
    fi
    
    log_success "System dependencies installed"
}

# Install Nextflow
install_nextflow() {
    log_info "Installing Nextflow..."
    
    # Create installation directory
    sudo mkdir -p "$INSTALL_DIR"
    sudo chown "$SETUP_USER:$SETUP_USER" "$INSTALL_DIR"
    
    # Download and install Nextflow
    cd "$INSTALL_DIR"
    curl -s https://get.nextflow.io | bash
    sudo chmod +x nextflow
    sudo ln -sf "$INSTALL_DIR/nextflow" /usr/local/bin/nextflow
    
    log_success "Nextflow installed to $INSTALL_DIR"
}

# Install container runtime
install_containers() {
    log_info "Installing container runtime..."
    
    # Install Docker
    if ! command -v docker >/dev/null 2>&1; then
        log_info "Installing Docker..."
        
        if [[ "$PKG_MGR" == "apt" ]]; then
            # Add Docker's official GPG key
            sudo mkdir -p /etc/apt/keyrings
            curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
            
            # Set up the repository
            echo "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
            
            sudo apt-get update
            sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
        elif [[ "$PKG_MGR" == "yum" || "$PKG_MGR" == "dnf" ]]; then
            sudo yum install -y yum-utils
            sudo yum-config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
            sudo yum install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
        fi
        
        # Add user to docker group
        sudo usermod -aG docker "$SETUP_USER"
        sudo systemctl enable docker
        sudo systemctl start docker
        
        log_success "Docker installed"
    else
        log_success "Docker already installed"
    fi
    
    # Install Singularity/Apptainer
    if ! command -v singularity >/dev/null 2>&1 && ! command -v apptainer >/dev/null 2>&1; then
        log_info "Installing Singularity..."
        
        # Install Singularity from conda-forge (easier than building from source)
        if command -v conda >/dev/null 2>&1; then
            conda install -y -c conda-forge singularity
        else
            log_warning "Singularity installation requires conda. Please install manually if needed."
        fi
        
        log_success "Singularity installation attempted"
    else
        log_success "Container runtime already available"
    fi
}

# Setup storage directories
setup_storage() {
    log_info "Setting up storage directories..."
    
    # Create project directories
    sudo mkdir -p /opt/tracktx/{assets,genomes,cache}
    sudo mkdir -p /var/tmp/nextflow-work
    sudo mkdir -p /tmp/nextflow-work
    
    # Set permissions
    sudo chown -R "$SETUP_USER:$SETUP_USER" /opt/tracktx
    sudo chmod -R 755 /opt/tracktx
    sudo chmod -R 1777 /var/tmp/nextflow-work
    sudo chmod -R 1777 /tmp/nextflow-work
    
    # Create user directories
    mkdir -p "$HOME/.nextflow"
    mkdir -p "$HOME/.singularity"
    
    log_success "Storage directories created"
}

# Create configuration files
create_configs() {
    log_info "Creating server configuration..."
    
    # Create server-specific params file
    cat > "$PROJECT_DIR/params.server.yaml" << EOF
# Server-specific parameters
output_dir: './results'
work_dir: '/var/tmp/nextflow-work'

# Storage configuration
storage:
  shared_data: '/opt/tracktx/data'
  shared_scratch: '/var/tmp'
  shared_tmp: '/tmp'
  project_dir: '/opt/tracktx/projects'
  user_scratch: '/var/tmp/\$USER'
  user_tmp: '/tmp/\$USER'
  cache_root: '/opt/tracktx/cache'
  assets_root: '/opt/tracktx/assets'

# Asset and genome caching
assets_dir: '/opt/tracktx/assets'
genome_cache: '/opt/tracktx/genomes'

# Server resource settings
host_cpus: $(nproc)
host_mem: $(free -g | awk '/^Mem:/{print $2}')

# Container settings
container:
  enabled: true
  engine: 'docker'
  cache_dir: '/opt/tracktx/containers'
EOF
    
    log_success "Server configuration created"
}

# Setup monitoring
setup_monitoring() {
    log_info "Setting up monitoring..."
    
    # Create log directory
    mkdir -p "$PROJECT_DIR/logs"
    
    # Create systemd service for pipeline monitoring (optional)
    cat > /tmp/tracktx-monitor.service << EOF
[Unit]
Description=TrackTx Pipeline Monitor
After=network.target

[Service]
Type=simple
User=$SETUP_USER
WorkingDirectory=$PROJECT_DIR
ExecStart=$PROJECT_DIR/scripts/monitor.sh
Restart=always
RestartSec=300

[Install]
WantedBy=multi-user.target
EOF
    
    log_success "Monitoring setup completed"
}

# Main installation function
main() {
    log_info "Starting TrackTx server setup..."
    
    check_root
    check_requirements
    install_dependencies
    install_nextflow
    install_containers
    setup_storage
    create_configs
    setup_monitoring
    
    log_success "TrackTx server setup completed!"
    log_info "Next steps:"
    log_info "1. Log out and back in to apply Docker group changes"
    log_info "2. Test the installation: nextflow run $PROJECT_DIR/main.nf -profile server,docker -help"
    log_info "3. Configure your samplesheet and run the pipeline"
    
    log_warning "Remember to configure your firewall and security settings as needed"
}

# Run main function
main "$@"
