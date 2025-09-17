"""
Base classes for calculation backends

Defines standard interfaces that all backends must implement
"""
from abc import ABC, abstractmethod
from typing import Dict, Any, Optional, List
from ase import Atoms


class CalculationBackend(ABC):
    """Base class for all calculation backends"""
    
    def __init__(self, name: str, **kwargs):
        self.name = name
        self.config = kwargs
    
    @abstractmethod
    def single_point(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Perform single point energy calculation"""
        pass
    
    @abstractmethod
    def optimize(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Perform geometry optimization"""
        pass
    
    @abstractmethod
    def frequencies(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Calculate vibrational frequencies/Hessian"""
        pass
    
    @abstractmethod
    def is_available(self) -> bool:
        """Check if backend is properly installed and configured"""
        pass


class RelaxationBackend(CalculationBackend):
    """Base class for relaxation/optimization backends"""
    
    @abstractmethod
    def relax_structure(self, project_name: str, atoms: Optional[Atoms] = None, **kwargs) -> bool:
        """Relax structure in project directory"""
        pass


class HessianBackend(CalculationBackend):
    """Base class for Hessian/frequency calculation backends"""
    
    @abstractmethod
    def calculate_hessian(self, project_name: str, atoms: Optional[Atoms] = None, **kwargs) -> bool:
        """Calculate Hessian matrix and frequencies"""
        pass


class BackendResult:
    """Standard result container for backend calculations"""
    
    def __init__(self, success: bool, atoms: Optional[Atoms] = None, 
                 energy: Optional[float] = None, forces: Optional[List] = None,
                 trajectory: Optional[List[Atoms]] = None, **kwargs):
        self.success = success
        self.atoms = atoms
        self.energy = energy
        self.forces = forces
        self.trajectory = trajectory
        self.metadata = kwargs
    
    @property
    def converged(self) -> bool:
        """Check if calculation converged"""
        return self.success and self.metadata.get('converged', False)