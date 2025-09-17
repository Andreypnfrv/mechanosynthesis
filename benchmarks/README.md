# Quantum Chemistry Benchmarks for Mechanosynthesis

This directory contains curated benchmark datasets for evaluating quantum chemistry methods in mechanosynthesis applications.

## 🚀 **NEW COMPREHENSIVE BENCHMARK SUITE**

**One command tests your method on ALL mechanosynthesis-relevant benchmarks:**
- ✅ **Reaction barriers** (NHTBH38) - organic bond formation  
- ✅ **Surface chemistry** (CMR) - STM tips, platforms, electrodes
- 🔄 **Future benchmarks** (BH76, SBH17) - added automatically

```bash
# Test single method comprehensively
PYTHONPATH=../py python run_benchmarks.py --backends dftb

# Compare multiple methods with cross-backend summary
PYTHONPATH=../py python run_benchmarks.py --backends dftb xtb orca-simple
```

## Target Elements

Our mechanosynthesis focus covers:
- **H, O**: passivation
- **C**: diamondoids, substrates
- **W, Au**: STM tips and platforms ⭐ **NEW**
- **Cu, Pt, Ag**: electrodes and alternative tips
- **B**: acceptor/dopant
- **Ge, Ti, Al**: platforms and contact metallurgy
- **P**: substrate
- **N**: vacancies
- **Si**: platforms


## Available Benchmark Datasets

| Dataset | Size | Elements | Method | Accuracy | Status | Auto-Run |
|---------|------|----------|---------|----------|---------|----------|
| **NHTBH38** | 38 reactions | C,H,O,N,F,Cl | CCSD(T) | ±0.5 kcal/mol | ✅ **Implemented** | ✅ **Default** |
| **CMR Adsorption** | 25 metals × 8 adsorbates | H,Au,W,Cu,Pt,Ag + 20 more | PBE/LDA | ±0.1 eV | ✅ **Implemented** | ✅ **Default** |
| BH76  | 76 reactions | C,H,O,N,F,Cl | CCSD(T) | ±0.3 kcal/mol | 🔄 **Next** | 🔄 Planned |
| SBH17 | 17 surface reactions | H + metals | SRP-DFT | ±0.2 eV | 🔄 **Next** | 🔄 Planned |
| HTBH38  | 38 H-transfer | C,H | CCSD(T) | ±0.5 kcal/mol | 🔄 Planned | 🔄 Optional |

*mechanosynthesis-relevant subset

### Dataset Descriptions

#### NHTBH38 - Non-Hydrogen Transfer Barriers
**Source**: Minnesota Database Collection (Truhlar Group)  
**Content**: 38 gas-phase reaction barriers covering:
- **Heavy atom transfer**: C-C, C-O, C-N bond formation reactions
- **Nucleophilic substitution**: SN2 reactions at carbon centers  
- **Association reactions**: Radical-radical coupling mechanisms
- **Unimolecular reactions**: Cyclization and rearrangement processes

**Reference Method**: CCSD(T)/aug-cc-pVTZ single points on B3LYP/6-31G(d) geometries  
**Mechanosynthesis Relevance**: Critical for carbon-carbon bond formation in diamond mechanosynthesis

#### CMR Adsorption - Transition Metal Surface Chemistry  
**Source**: DTU Computational Materials Repository  
**Content**: Adsorption energies for multiple adsorbates on 25 metal surfaces:
- **Metals**: Au, W, Cu, Pt, Ag, Pd, Ir, Os, Re, Ta, Hf, Mo, Nb, Zr, Y, Sc, Ti, V, Cr, Mn, Fe, Co, Ni, Ru, Rh
- **Adsorbates**: H, OH, CH, NO, CO, N₂, N, O (H adsorption currently implemented)
- **Surfaces**: (111) facets for FCC metals, (0001) for HCP, (110) for BCC

**Reference Methods**: LDA, PBE, RPBE, BEEF-vdW functionals  
**Mechanosynthesis Relevance**: Essential for STM tip interactions (W), platform chemistry (Au), and electrode design (Cu, Ag)

## Recent Results Summary

### NHTBH38 Reaction Barriers (5 key reactions)

| Method | MAE (kcal/mol) | RMSE (kcal/mol) | Best For |
|--------|----------------|-----------------|----------|
| **DFTB+** | 2.19 | 2.19 | ✅ **Overall accuracy** |
| xTB | 3.69 | 3.69 | Speed |

### CMR H Adsorption on Metals (vs PBE reference) ⭐ **NEW**
| Metal | DFTB+ Error (eV) | xTB Error (eV) | Relevance | Best Method |
|-------|------------------|----------------|-----------|-------------|
| **W (STM tips)** | 0.157 | **0.075** | HIGH | ✅ **xTB** |
| **Au (platforms)** | 0.594 | **0.456** | HIGH | ✅ **xTB** |
| Cu (electrodes) | 0.589 | **0.472** | MEDIUM | ✅ **xTB** |
| Pt (alt tips) | **0.574** | 0.678 | MEDIUM | ✅ **DFTB+** |
| Ag (electrodes) | 0.929 | **0.860** | MEDIUM | ✅ **xTB** |

**Overall CMR Summary:**
- **xTB**: MAE = 0.508 eV (better for surface chemistry)
- **DFTB+**: MAE = 0.569 eV (better for organic reactions)

## Usage

### 🎯 **Recommended Usage** (New Multi-Backend Architecture)

```bash
# Test single method comprehensively (NHTBH38 + CMR + future benchmarks)
python run_benchmarks.py --backends dftb
python run_benchmarks.py --backends xtb

# Compare multiple methods with cross-backend summary table
python run_benchmarks.py --backends dftb xtb
python run_benchmarks.py --backends dftb xtb orca-simple

# Get method recommendations for mechanosynthesis applications
python run_benchmarks.py --backends dftb xtb  # → Auto-generates comparison summary
```

### 📊 **Information Commands**
```bash
# List available datasets and backends
python run_benchmarks.py --list-datasets
python run_benchmarks.py --list-backends
```

### 🔧 **Advanced Usage**
```bash
# Run specific dataset with multiple backends
python run_benchmarks.py --dataset cmr_adsorption --backends dftb xtb

# Legacy support (comma-separated, still works)
python run_benchmarks.py --compare dftb,xtb
python run_benchmarks.py --backend dftb
```

## Results Storage
- **Raw results**: `benchmarks/results/`
- **CMR data**: `benchmarks/cmr_adsorption/raw/`
- **Summary reports**: Generated on demand

## 🏆 **Auto-Generated Cross-Backend Comparison**

When you run multiple backends, you get an automatic comparison summary:

```
🏆 CROSS-BACKEND COMPARISON SUMMARY
============================================================

📊 REACTION BARRIERS (NHTBH38)
Method       MAE (kcal/mol)  RMSE (kcal/mol)  Status
------------------------------------------------------------
dftb         2.19           2.19              ✅ Good
xtb          3.69           3.69              ⚠️ Fair

🔬 SURFACE ADSORPTION (CMR)  
Method       MAE vs PBE (eV)  MAE vs LDA (eV)  Status
-------------------------------------------------------------
dftb         0.569           0.523             ⚠️ Fair
xtb          0.508           0.461             ⚠️ Fair

🎯 MECHANOSYNTHESIS RECOMMENDATIONS
==================================================
• **Organic reactions**: DFTB+ best for C-H bond formation
• **STM tips (W/Au)**: xTB best for surface interactions
• **Mixed systems**: DFTB+ most robust across applications
• **High accuracy**: Orca when precision is critical
```

## Data Sources
- **NHTBH38**: Curated subset from original 38 reaction barriers database
- **CMR**: DTU Computational Materials Repository (https://cmrdb.fysik.dtu.dk/adsorption/)
- All datasets include full attribution in respective metadata.yaml files

## 🚀 **Quick Start**

```bash
# Complete method evaluation (one command!)
cd mechanosynthesis/benchmarks
PYTHONPATH=../py python run_benchmarks.py --backends dftb

# Compare DFTB+ vs xTB (with auto cross-backend summary!)
PYTHONPATH=../py python run_benchmarks.py --backends dftb xtb
```

**Output**: Comprehensive method evaluation covering:
- ✅ **Reaction barriers** (NHTBH38) - organic mechanosynthesis reactions  
- ✅ **Surface chemistry** (CMR) - STM tips, platforms, electrodes
- 🔄 **Future benchmarks** (BH76, SBH17) - auto-added when implemented
- 🏆 **Cross-method comparison** - automatic when multiple backends specified

## 🗺️ **Roadmap** 
| Priority | Dataset | Description | ETA |
|----------|---------|-------------|-----|
| 🔥 **HIGH** | **BH76** | 76 reaction barriers (broader NHTBH38) | Next |
| 🔥 **HIGH** | **SBH17** | Surface reaction barriers (chemisorption) | Next |  
| 🔶 **MED** | **W4-11** | High-accuracy thermochemistry reference | Later |
| 🔶 **MED** | **Transition states** | Exact barrier calculations | Later |
| 🔷 **LOW** | **ML predictions** | Property estimation models | Future |
