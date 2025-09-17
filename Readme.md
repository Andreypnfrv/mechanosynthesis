# Requirements

dftb+, orca, python venv

# How To Use

## Environment Setup
```bash
source .pyenv/bin/activate
export SKDIR=/Users/andreypanferov/opt/dftb+/slakos/
export PATH=/Users/andreypanferov/opt/dftb+/bin:$PATH
```

## Main Commands

### Testing
```bash
python py/main.py test          # Run all tests
```

### Parts Pipeline (Individual Molecules)
```bash
# Build molecule from SMILES or name
python py/main.py build molecule H2O --project test_h2o

# Relax structure with different backends
python py/main.py relax PROJECT_NAME --backend dftb     # DFTB+ (main workhorse)
python py/main.py relax PROJECT_NAME --backend xtb      # xTB native
python py/main.py relax PROJECT_NAME --backend orca     # ORCA DFT

# Calculate frequencies/Hessian for validation
python py/main.py hessian PROJECT_NAME --backend dftb           # DFTB+ analytical (fast)
python py/main.py hessian PROJECT_NAME --backend xtb            # xTB numerical
python py/main.py hessian PROJECT_NAME --backend orca           # ORCA analytical
python py/main.py hessian PROJECT_NAME --hessian-orca-simple    # ORCA PBE/def2-SVP (stable)

# Thermochemistry analysis
python py/main.py thermo PROJECT_NAME --plot --report --temperature-range 298,1500
python py/main.py thermo PROJECT_NAME --animate --frames 30 --fps 5
```

### Sequences Pipeline (Multi-Component Scenes)
```bash
# Generate sequence frames
python py/main.py sequence SCENE_NAME --generate --frame 0
python py/main.py sequence SCENE_NAME --generate --frame 1
```

### Benchmarking
```bash
# Run benchmarks on datasets
python benchmarks/run_benchmarks.py --dataset nhtbh38 --backend dftb
python benchmarks/run_benchmarks.py --dataset cmr_adsorption --backend xtb
python benchmarks/run_benchmarks.py --compare dftb,xtb
```

## Backend Recommendations
- **Small molecules (<15 atoms)**: `--backend xtb` for speed or `--backend dftb` for accuracy
- **Medium molecules (15-50 atoms)**: `--backend dftb` or `--hessian-orca-simple`
- **Large molecules (>50 atoms)**: `--backend dftb` (fastest analytical methods)


# Log

#### 2025.09.17
Have time to continue. To reproduce CBN's chain I need to add thermo (so added this one, however no time to test yet, will test tomorrow)
Cleaned up the repo, added comprehensive Readme.txt,
Added installation for xTB, DFTB+, ORCA.
Upgrading my benchmarks.

### 2025.09.13
Implemented:
- Metrics for comparing structure before and after the relaxation
- First time tried to assemble the scene
- Benchmarks to compare different backends (taken from different libraries). Actually quite interesting, for mechanosynthesis it looks like that DFTB+ is more accurate then XTB

### 2025.09.12
So, i want to compare how different backends work for the different atoms, that take part in mechanosynthesis. It's important for me to compare, performance, accuracy,stability. The atoms that we'll be definetly using for mechanosynthesis:
- H, O for passivation
- C for diamondoids, as a substrate
- W, Pt, Ir for STM tip
- B as aceptor / dopant
- Ge (Si-compatible platform), Ti (adhesion layers), Al (contact metallurgy)
- P as a substrate
- N for vacancies
- Au, Cu, Ag as a platform
- Si as a platform
- Si as a substrate
Ideally i need to make a small set of molecules with well-known geometries that include these key atoms, to test the method accuracy
Reporting template must include:
-	System, backend+level (func/basis/PP/relativistic/dispersion), settings (cutoff, k-mesh, grids), result(s), Δ vs ref, SCF iters/step, wall time, fail notes
- Ideally ofc i need to find some professional benchmarks, but let's start at least with this
List of benchmarks:
- W(CH3)6

### 2025.08-2025.09
Just started this repo and slowly diving into mechanosynthesis and it's problems.
Researched ORCA, evolution of DFTB+, XTB, G-XTB.
Implemented the initial repo
