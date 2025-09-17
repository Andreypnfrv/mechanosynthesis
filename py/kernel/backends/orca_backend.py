"""
Orca Backend Implementation

Handles Orca ab initio DFT calculations
"""
import os
from typing import Dict, Any, Optional
from ase import Atoms
from .base import RelaxationBackend, HessianBackend, BackendResult


class OrcaRelaxationBackend(RelaxationBackend):
    """Orca relaxation backend"""
    
    def __init__(self, method: str = 'B3LYP', basis: str = 'def2-SVP', 
                 variant: str = 'native', **kwargs):
        super().__init__('orca_relax', **kwargs)
        self.method = method
        self.basis = basis
        self.variant = variant  # 'native', 'ase'
    
    def is_available(self) -> bool:
        """Check if Orca is available"""
        orca_path = '/Applications/orca-6.1.0/orca'
        return os.path.exists(orca_path)
    
    def single_point(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Single point Orca calculation"""
        return {'energy': 0.0, 'forces': [], 'success': False}
    
    def optimize(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Orca geometry optimization"""
        return {'atoms': atoms, 'converged': False, 'success': False}
    
    def frequencies(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Orca frequency calculation"""
        return {'frequencies': [], 'success': False}
    
    def relax_structure(self, project_name: str, atoms: Optional[Atoms] = None, **kwargs) -> bool:
        """Relax structure using Orca - delegates to existing implementation"""
        try:
            from ..parts.relax import relax_orca_native, relax_orca_dft
            
            nprocs = kwargs.get('nprocs', 4)
            
            if self.variant == 'native':
                return relax_orca_native(project_name, self.method, self.basis, nprocs, False, 'phonopy')
            elif self.variant == 'ase':
                return relax_orca_dft(project_name, self.method, self.basis, nprocs)
            else:
                return False
        except ImportError:
            return False


class OrcaHessianBackend(HessianBackend):
    """Orca Hessian calculation backend"""
    
    def __init__(self, method: str = 'B3LYP', basis: str = 'def2-SVP', **kwargs):
        super().__init__('orca_hessian', **kwargs)
        self.method = method
        self.basis = basis
    
    def is_available(self) -> bool:
        """Check if Orca is available"""
        return OrcaRelaxationBackend().is_available()
    
    def single_point(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Not used for Hessian backend"""
        return {'success': False}
    
    def optimize(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Not used for Hessian backend"""
        return {'success': False}
    
    def frequencies(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Orca frequency calculation"""
        return {'frequencies': [], 'success': False}
    
    def calculate_hessian(self, project_name: str, atoms: Optional[Atoms] = None, **kwargs) -> bool:
        """Calculate Hessian using Orca - delegates to existing implementation"""
        try:
            from ..parts.hessian import calculate_hessian_orca_standalone
            
            nprocs = kwargs.get('nprocs', 4)
            return calculate_hessian_orca_standalone(project_name, self.method, self.basis, nprocs)
        except ImportError:
            return False