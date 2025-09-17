#!/bin/bash

set -euo pipefail

# Mechanosynthesis Pipeline Installer
# Installs DFTB+, Slakos parameters, and validates installation

readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly INSTALL_DIR="$HOME/opt/dftb+"
readonly SLAKOS_DIR="$INSTALL_DIR/slakos"
readonly BENCHMARK_DIR="$SCRIPT_DIR/benchmarks"

# Colors for output
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly NC='\033[0m' # No Color

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
    echo -e "${RED}[ERROR]${NC} $1" >&2
}

die() {
    log_error "$1"
    exit 1
}

# Check system requirements
check_system() {
    log_info "Checking system requirements..."
    
    # Check OS
    case "$(uname -s)" in
        Darwin*)
            log_info "Detected macOS"
            ;;
        Linux*)
            log_info "Detected Linux"
            ;;
        *)
            die "Unsupported operating system: $(uname -s)"
            ;;
    esac
    
    # Check conda/mamba
    if command -v mamba &> /dev/null; then
        readonly CONDA_CMD="mamba"
        log_info "Using mamba for package management"
    elif command -v conda &> /dev/null; then
        readonly CONDA_CMD="conda"
        log_info "Using conda for package management"
    else
        die "conda or mamba not found. Please install miniconda/anaconda first."
    fi
    
    # Check internet connectivity
    if ! ping -c 1 google.com &> /dev/null; then
        die "No internet connection. Required for downloading packages."
    fi
    
    log_success "System requirements check passed"
}

# Install DFTB+ via conda
install_dftb() {
    log_info "Installing DFTB+ via conda-forge..."
    
    if command -v dftb+ &> /dev/null; then
        log_warning "DFTB+ already installed at $(which dftb+)"
        return 0
    fi
    
    $CONDA_CMD install -c conda-forge dftb-plus -y || die "Failed to install DFTB+"
    
    # Verify installation
    if ! command -v dftb+ &> /dev/null; then
        die "DFTB+ installation failed - binary not found in PATH"
    fi
    
    log_success "DFTB+ installed successfully: $(which dftb+)"
}

# Download and setup Slakos parameters
setup_slakos() {
    log_info "Setting up Slakos parameters..."
    
    # Create directories
    mkdir -p "$SLAKOS_DIR"
    
    # Use existing parameters from 3rdparty if available
    if [[ -f "$SCRIPT_DIR/3rdparty/ptbp.zip" ]]; then
        log_info "Using existing PTBP parameters from 3rdparty..."
        if [[ ! -d "$SLAKOS_DIR/ptbp" ]]; then
            cd "$SLAKOS_DIR"
            unzip -q "$SCRIPT_DIR/3rdparty/ptbp.zip" || die "Failed to extract PTBP parameters"
            log_success "PTBP parameters installed from 3rdparty"
        else
            log_info "PTBP parameters already present"
        fi
    else
        # Try downloading from official source
        log_info "Downloading PTBP Slakos parameters..."
        local ptbp_url="https://www.dftb.org/fileadmin/DFTB/public/slako/ptbp/ptbp.tar.xz"
        local ptbp_archive="$SLAKOS_DIR/ptbp.tar.xz"
        
        if [[ ! -d "$SLAKOS_DIR/ptbp" ]]; then
            if curl -L "$ptbp_url" -o "$ptbp_archive" 2>/dev/null; then
                cd "$SLAKOS_DIR"
                tar -xf "$ptbp_archive" && rm "$ptbp_archive"
                log_success "PTBP parameters downloaded and installed"
            else
                log_warning "Failed to download PTBP parameters from official source"
                log_info "Please manually download and extract to: $SLAKOS_DIR/ptbp"
            fi
        fi
    fi
    
    # Use existing prior parameters if available
    if [[ -f "$SCRIPT_DIR/3rdparty/prior.zip" ]]; then
        log_info "Using existing prior parameters from 3rdparty..."
        if [[ ! -d "$SLAKOS_DIR/prior" ]]; then
            cd "$SLAKOS_DIR"
            unzip -q "$SCRIPT_DIR/3rdparty/prior.zip" || die "Failed to extract prior parameters"
            log_success "Prior parameters installed from 3rdparty"
        else
            log_info "Prior parameters already present"
        fi
    else
        log_info "Prior parameters not found in 3rdparty (optional)"
    fi
    
    # Set default to PTBP if available, otherwise try prior
    if [[ -d "$SLAKOS_DIR/ptbp" ]]; then
        export SKDIR="$SLAKOS_DIR/ptbp"
        log_info "Using PTBP parameters as default"
    elif [[ -d "$SLAKOS_DIR/prior" ]]; then
        export SKDIR="$SLAKOS_DIR/prior"
        log_info "Using prior parameters as default"
    else
        die "No Slakos parameters found. Please ensure parameters are available."
    fi
    
    # Add to shell profile
    setup_environment
    
    log_success "Slakos parameters setup complete"
}

# Setup environment variables
setup_environment() {
    log_info "Setting up environment variables..."
    
    local shell_profile=""
    if [[ -n "${ZSH_VERSION:-}" ]]; then
        shell_profile="$HOME/.zshrc"
    elif [[ -n "${BASH_VERSION:-}" ]]; then
        shell_profile="$HOME/.bashrc"
    else
        shell_profile="$HOME/.profile"
    fi
    
    local env_setup="
# DFTB+ Slakos parameters
export SKDIR=\"$SLAKOS_DIR/ptbp\"
export PATH=\"$INSTALL_DIR/bin:\$PATH\"
"
    
    if ! grep -q "SKDIR" "$shell_profile" 2>/dev/null; then
        echo "$env_setup" >> "$shell_profile"
        log_success "Environment variables added to $shell_profile"
    else
        log_info "Environment variables already configured"
    fi
    
    # Export for current session
    export SKDIR="$SLAKOS_DIR/ptbp"
    export PATH="$INSTALL_DIR/bin:$PATH"
}

# Check for ORCA installation
check_orca() {
    log_info "Checking for ORCA installation..."
    
    if command -v orca &> /dev/null; then
        local orca_version=$(orca --version 2>&1 | head -1 || echo "unknown")
        log_success "ORCA found: $(which orca) ($orca_version)"
        return 0
    fi
    
    # Check common installation paths
    local orca_paths=(
        "/Applications/orca*/orca"
        "/opt/orca*/orca"
        "/usr/local/orca*/orca"
        "$HOME/orca*/orca"
    )
    
    for path in "${orca_paths[@]}"; do
        if ls $path &> /dev/null; then
            log_warning "ORCA found but not in PATH: $path"
            log_info "Add to PATH: export PATH=\"$(dirname $path):\$PATH\""
            return 0
        fi
    done
    
    log_warning "ORCA not found"
    log_info "To install ORCA:"
    log_info "1. Register at https://orcaforum.kofo.mpg.de"
    log_info "2. Download ORCA for your platform"
    log_info "3. Extract and add to PATH"
}

# Download benchmark datasets
setup_benchmarks() {
    log_info "Setting up benchmark datasets..."
    
    if [[ ! -d "$BENCHMARK_DIR" ]]; then
        mkdir -p "$BENCHMARK_DIR"
    fi
    
    # Check if benchmarks already exist
    if [[ -f "$BENCHMARK_DIR/run_benchmarks.py" ]]; then
        log_info "Benchmark infrastructure already present"
        return 0
    fi
    
    log_warning "Benchmark datasets not found in $BENCHMARK_DIR"
    log_info "Benchmarks should include:"
    log_info "- CMR adsorption dataset"
    log_info "- NHTBH38 reaction barriers"
    log_info "Please ensure benchmark data is available for validation"
}

# Validate installation
validate_installation() {
    log_info "Validating installation..."
    
    # Test DFTB+
    if ! command -v dftb+ &> /dev/null; then
        die "DFTB+ not found in PATH"
    fi
    
    # Test Slakos parameters
    if [[ ! -d "$SKDIR" ]]; then
        die "SKDIR not set or directory doesn't exist: $SKDIR"
    fi
    
    # Check for essential Slakos files
    local essential_files=("C-C.skf" "H-H.skf" "C-H.skf" "H-C.skf")
    for file in "${essential_files[@]}"; do
        if [[ ! -f "$SKDIR/$file" ]]; then
            die "Essential Slakos file missing: $SKDIR/$file"
        fi
    done
    
    # Test DFTB+ with simple molecule
    local test_dir=$(mktemp -d)
    cat > "$test_dir/dftb_in.hsd" << 'EOF'
Geometry = GenFormat {
2 C
H
1 1 0.0 0.0 0.0
2 1 0.74 0.0 0.0
}

Driver = GeometryOptimization {}

Hamiltonian = DFTB {
  SCC = Yes
  SlaterKosterFiles = Type2FileNames {
    Prefix = ""
    Separator = "-"
    Suffix = ".skf"
  }
  MaxAngularMomentum {
    H = "s"
  }
}

Options {
  WriteResultsTag = Yes
}

Analysis {
  CalculateForces = Yes
}
EOF
    
    cd "$test_dir"
    if dftb+ > dftb.log 2>&1; then
        log_success "DFTB+ validation test passed"
    else
        log_error "DFTB+ validation test failed. Check log: $test_dir/dftb.log"
        cat "$test_dir/dftb.log"
        die "DFTB+ validation failed"
    fi
    
    rm -rf "$test_dir"
}

# Main installation function
main() {
    log_info "Starting Mechanosynthesis Pipeline installation..."
    
    check_system
    install_dftb
    setup_slakos
    check_orca
    setup_benchmarks
    validate_installation
    
    log_success "Installation completed successfully!"
    log_info ""
    log_info "Environment setup:"
    log_info "  SKDIR: $SKDIR"
    log_info "  DFTB+: $(which dftb+)"
    log_info ""
    log_info "To use in new shell sessions, restart your terminal or run:"
    log_info "  source ~/.zshrc  # or ~/.bashrc"
    log_info ""
    log_info "To switch Slakos parameter sets:"
    log_info "  export SKDIR=\"$SLAKOS_DIR/ptbp\"   # for PTBP (default)"
    log_info "  export SKDIR=\"$SLAKOS_DIR/prior\"  # for prior parameters"
}

# Run main function
if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    main "$@"
fi