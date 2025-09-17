# SUMMARY
I'm researching mechanosynthesis and building a pipeline for it. This is my learning repository. The goal is to have ability to build some initial molecular parts - tips, platforms, nanoparticles, molecular tools, etc. then calculate them with some different qchem methods , and then build scenese out of them and calculate reactions. The most important reactions and scenes must be calculated using methods with highest prediction accuracy DFT.

## METHODS
- DFTB+ - main workhorse
- XTB - native
- gXTB - in docker
- ORCA


## IMPLEMENTATIONS
- `py/` - Python implementation (current working version). use python from .pyenv to run it `source .pyenv/bin/activate`
- `rust/` - Rust implementation (experimental/future)

# CONSTRAINTS
- We can't use PBC or freeze atoms as we'll need to calc real nanoparticles

# ENV
source .pyenv/bin/activate

# DIRECTORY STRUCTURE
Keep the root folder clean.

**Architecture:**
- **Parts Module** (`py/kernel/parts/`) - Individual molecular pieces
- **Sequences Module** (`py/kernel/sequences/`) - Scene composition and animation
- **Backends Module** (`py/kernel/backends/`) - Calculation engine interfaces
- **Settings** (`py/settings.py`) - Configuration using `.env` file

```
mechanosynthesis/
├── py/                  # Python implementation
│   ├── main.py          # Command runner
│   ├── settings.py      # Configuration and paths (uses .env)
│   ├── kernel/          # Core computational modules
│   │   ├── parts/       # Parts Module - individual molecular pieces
│   │   │   ├── build.py # Molecular building functions
│   │   │   ├── relax.py # Relaxation methods
│   │   │   └── hessian.py # Frequency calculations
│   │   ├── sequences/   # Sequences Module - scene composition
│   │   │   ├── scene.py # Scene assembly and positioning
│   │   │   ├── artifacts.py # Fragment detection and extraction
│   │   │   └── animation.py # Keyframe animation
│   │   ├── backends/    # Calculation Backends
│   │   │   ├── base.py  # Abstract interfaces
│   │   │   ├── xtb_backend.py    # xTB implementations
│   │   │   ├── dftb_backend.py   # DFTB+ implementations
│   │   │   └── orca_backend.py   # Orca implementations
│   │   ├── build/       # Legacy build modules
│   │   └── utils/       # Shared utilities
├── rust/                # Rust implementation (future)
├── tests/               # Testing framework
├── data/                # All project data and scenes
│   ├── PROJECT_NAME/    # Individual project workflows
│   └── scenes/          # Scene definitions (YAML)
├── docs/                # Documentation and lessons
├── 3rdparty/            # Third-party tools (XTB variants, etc.)
└── .env                 # Environment configuration
```

# PIPELINE WORKFLOWS

## Parts Pipeline (Individual Molecules)
1. Create project folder inside `/data/PROJECT_NAME/` (via `settings.py`)
2. Build molecular part with `py/kernel/parts/build.py`
3. Relax part with `py/kernel/parts/relax.py` (xTB, DFTB+, Orca DFT)
4. Validate with `py/kernel/parts/hessian.py` (frequencies/stability)

## Sequences Pipeline (Multi-Component Scenes)
1. Define scene in `data/scenes/SCENE_NAME.yaml`
2. Assemble pieces with `py/kernel/sequences/scene.py`
3. Calculate scene reactions using backends
4. Extract artifacts with `py/kernel/sequences/artifacts.py`
5. Create animations with `py/kernel/sequences/animation.py`

## Backend Integration
- All calculations use `py/kernel/backends/` interfaces
- Unified API for xTB, DFTB+, Orca methods
- Easy to add new calculation engines
- Settings managed via `.env` and `py/settings.py`

## RELAXATION METHODS
- `--relax`: DFTB+ with ASE/BFGS optimizer

## HESSIAN/FREQUENCY ANALYSIS METHODS
- `--hessian`: DFTB+ native analytical Hessian (fast and accurate)
- `--hessian-orca`: Orca DFT Hessian calculation (analytical or numerical)
- `--hessian-orca-simple`: Orca PBE/def2-SVP Hessian (more stable)
- `--hessian-xtb`: xTB Hessian calculation (fast, numerical derivatives)
- `--hessian-gxtb`: g-xTB Hessian via Docker ⚠️ VERY SLOW (numerical, 3N+1 evaluations)

### Hessian Method Recommendations
- **Small molecules (<15 atoms)**: `--hessian-xtb` for fast xTB quality or `--hessian` for DFTB+ accuracy
- **Medium molecules (15-50 atoms)**: `--hessian` (DFTB+ native) or `--hessian-orca-simple` (PBE/def2-SVP)
- **Large molecules (>50 atoms)**: `--hessian` (DFTB+ native, fastest analytical method)

## OUTPUT FORMATS (for visualization compatibility)
DFTB+ relaxation generates multiple trajectory formats:
- `.xyz`: OVITO-compatible multi-frame XYZ with energy information
- `.extxyz`: Extended XYZ with metadata
- `.pdb`: PDB trajectory for VMD/PyMOL
- `.traj`: ASE trajectory format
- `.mp4/.gif`: Optimization video animations

# SCENE SYSTEM (NEW)

## Scene Definitions (YAML)
```yaml
name: "stm_approach"
pieces:
  - name: "tip"
    source: "data/w_stm_tip/outputs/xtb/tip_relaxed.xyz"
    position: [0, 0, 5]  # 5Å above surface
  - name: "surface"
    source: "data/diamond/outputs/diamond_relaxed.xyz"
    position: [0, 0, 0]

calculation:
  method: "xtb"
  backend: "fixed"

extract:
  fixed:
    - name: "modified_tip"
      atoms: [1-50]
  artifacts:
    - method: "connectivity"
      min_size: 1
      displacement_threshold: 2.0
```

## Animation Keyframes
```yaml
animation:
  keyframes:
    - frame: 0
      description: "initial_approach"
      tip_position: [0, 0, 10]
    - frame: 1
      description: "contact"
      tip_position: [0, 0, 3]
    - frame: 2
      description: "reaction"
      calculate: true  # Run optimization
```


# STACK
> which dftb+
/Users/andreypanferov/opt/dftb+/
> slakos files
export SKDIR=/Users/andreypanferov/opt/dftb+/slakos/


## In general
If i ask you to calculate with specific backend - don't calculate with the other backend!
If you need to create some tmp files - put them in tmp folder
Never rewrite tests if i haven't asked you about it directly

# TESTS
All the tests must be in `py/tests`, done in `pytest` and summoned via `python py/main.py`
We don't use unit tests, we only write e2e tests that launch the system with some params with CLI
