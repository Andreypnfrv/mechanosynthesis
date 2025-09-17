# Requirements

dftb+, orca, python venv
TODO: yeah, maybe i should add some kind of installation script

# How To Use

`python3 py/main.py test` - proof that functions do work

# Log

#### 2025.09.17
Have time to continue. To reproduce CBN's chain i need to add thermo.
Also actually i want to deploy the repo

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
