"""
xTB Backend Implementation

Handles xTB semiempirical tight-binding calculations with multiple backend variants:
- native: System-installed xTB
- fixed: Our patched xTB version (recommended)
- gxtb: g-xTB via Docker (high quality)
"""
import os
import subprocess
from typing import Dict, Any, Optional
from ase import Atoms
from .base import RelaxationBackend, HessianBackend, BackendResult


class XTBRelaxationBackend(RelaxationBackend):
    """xTB relaxation backend"""
    
    def __init__(self, method: str = 'GFN2-xTB', backend: str = 'fixed', **kwargs):
        super().__init__('xtb_relax', **kwargs)
        self.method = method
        self.backend = backend
    
    def is_available(self) -> bool:
        """Check if xTB is available"""
        if self.backend == 'fixed':
            fixed_path = '/Users/andreypanferov/Documents/mechanosynthesis/3rdparty/xtb/fork/build_cmake/xtb'
            return os.path.exists(fixed_path)
        elif self.backend == 'native':
            try:
                result = subprocess.run(['xtb', '--version'], capture_output=True)
                return result.returncode == 0
            except FileNotFoundError:
                return False
        elif self.backend == 'gxtb':
            # Check for Docker and g-xTB container
            try:
                result = subprocess.run(['docker', 'images', 'gxtb'], capture_output=True)
                return 'gxtb' in result.stdout.decode()
            except FileNotFoundError:
                return False
        return False
    
    def single_point(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Single point xTB calculation"""
        # Implementation would call appropriate xTB backend
        return {'energy': 0.0, 'forces': [], 'success': False}
    
    def optimize(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """xTB geometry optimization"""
        # Implementation would call appropriate xTB backend
        return {'atoms': atoms, 'converged': False, 'success': False}
    
    def frequencies(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """xTB frequency calculation"""
        # Implementation would call appropriate xTB backend
        return {'frequencies': [], 'success': False}
    
    def relax_structure(self, project_name: str, atoms: Optional[Atoms] = None, **kwargs) -> bool:
        """Relax structure using xTB - delegates to existing implementation"""
        # Import the existing function to avoid code duplication
        try:
            from ..parts.relax import relax_xtb
            return relax_xtb(project_name, self.method, kwargs.get('nprocs', 1), self.backend)
        except ImportError:
            return False


class XTBHessianBackend(HessianBackend):
    """xTB Hessian calculation backend"""
    
    def __init__(self, method: str = 'GFN2-xTB', backend: str = 'fixed', **kwargs):
        super().__init__('xtb_hessian', **kwargs)
        self.method = method
        self.backend = backend
    
    def is_available(self) -> bool:
        """Check if xTB is available"""
        return XTBRelaxationBackend(self.method, self.backend).is_available()
    
    def single_point(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Not used for Hessian backend"""
        return {'success': False}
    
    def optimize(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Not used for Hessian backend"""
        return {'success': False}
    
    def frequencies(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """xTB frequency calculation"""
        # Implementation would call appropriate xTB backend
        return {'frequencies': [], 'success': False}
    
    def calculate_hessian(self, project_name: str, atoms: Optional[Atoms] = None, **kwargs) -> bool:
        """Calculate Hessian using xTB - delegates to existing implementation"""
        try:
            from ..parts.hessian import calculate_hessian_xtb_standalone
            return calculate_hessian_xtb_standalone(project_name, self.method, self.backend)
        except ImportError:
            return False