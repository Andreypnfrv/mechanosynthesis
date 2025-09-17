"""
Calculation Backends Module

This module provides interfaces to various quantum chemistry calculation engines:
- xTB (semiempirical tight-binding)
- g-xTB (machine learning enhanced tight-binding)  
- DFTB+ (density functional tight-binding)
- Orca (ab initio DFT)
- ASE (Atomic Simulation Environment calculators)

Each backend implements standard interfaces for:
- Single point energy calculations
- Geometry optimization 
- Frequency/Hessian calculations
- Molecular dynamics
"""

from . import xtb_backend
from . import dftb_backend  
from . import orca_backend
from . import ase_backend

__all__ = ['xtb_backend', 'dftb_backend', 'orca_backend', 'ase_backend']