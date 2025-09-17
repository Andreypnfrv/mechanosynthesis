# Mechanosynthesis Pipeline Installation

Automated installer for DFTB+, Slakos parameters, and computational chemistry tools.

## Quick Start

```bash
# Clone repository
git clone <repository-url>
cd mechanosynthesis

# Run installer
./install.sh
```

The installer will:
1. ✅ Check system requirements (macOS/Linux, conda)
2. ✅ Install DFTB+ via conda-forge
3. ✅ Setup Slakos parameter sets (PTBP + prior)
4. ✅ Configure environment variables
5. ✅ Validate installation with test calculations
6. ⚠️ Check for ORCA (manual installation required)

## Requirements

- **macOS** or **Linux**
- **conda** or **mamba** package manager
- **Internet connection** for downloads

## Manual Installation Steps

If the automated installer fails, follow these manual steps:

### 1. Install DFTB+
```bash
conda install -c conda-forge dftb-plus
```

### 2. Setup Slakos Parameters

Create directory structure:
```bash
mkdir -p ~/opt/dftb+/slakos
```

Download parameter sets:
- **PTBP**: For organic molecules
- **Prior**: Alternative parameter set

### 3. Environment Variables

Add to your shell profile (`~/.zshrc` or `~/.bashrc`):
```bash
export SKDIR="$HOME/opt/dftb+/slakos/ptbp"
export PATH="$HOME/opt/dftb+/bin:$PATH"
```

### 4. Install ORCA (Optional)

1. Register at https://orcaforum.kofo.mpg.de
2. Download ORCA for your platform
3. Extract and add to PATH

## Parameter Sets

### PTBP (Default)
- Optimized for organic molecules
- Best for mechanosynthesis research
- Location: `~/opt/dftb+/slakos/ptbp`

### Prior
- Alternative parameter set
- Broader element coverage
- Location: `~/opt/dftb+/slakos/prior`

Switch between parameter sets:
```bash
export SKDIR="$HOME/opt/dftb+/slakos/ptbp"   # PTBP
export SKDIR="$HOME/opt/dftb+/slakos/prior"  # Prior
```

## Validation

Test your installation:
```bash
# Activate environment
source .pyenv/bin/activate

# Run simple calculation
python py/main.py relax test_h2 --backend dftb

# Run benchmarks
python py/main.py test unified_dftb
```

## Benchmark Datasets

For validation, ensure these datasets are available:
- **CMR adsorption**: Surface interaction benchmarks
- **NHTBH38**: Reaction barrier benchmarks

## Troubleshooting

### DFTB+ not found
```bash
which dftb+
conda list dftb-plus
```

### Slakos parameters missing
```bash
echo $SKDIR
ls -la $SKDIR
```

### Permission errors
```bash
chmod +x install.sh
sudo chown -R $USER ~/opt/dftb+
```

### Environment not persistent
Add exports to shell profile:
```bash
echo 'export SKDIR="$HOME/opt/dftb+/slakos/ptbp"' >> ~/.zshrc
source ~/.zshrc
```

## Architecture

```
~/opt/dftb+/
├── bin/dftb+              # DFTB+ binary
└── slakos/
    ├── ptbp/              # PTBP parameters (default)
    │   ├── C-C.skf
    │   ├── C-H.skf
    │   └── ...
    └── prior/             # Prior parameters
        ├── C-C.skf
        └── ...
```

## Support

For issues:
1. Check this troubleshooting section
2. Verify all requirements are met
3. Run installer with verbose logging
4. Check individual component installations

The installer is idempotent - safe to run multiple times.