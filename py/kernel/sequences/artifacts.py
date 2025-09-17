"""
Artifact detection and extraction module

Handles detecting separated molecules and fragments during calculations
"""
import numpy as np
from typing import List, Dict, Any, Tuple
from ase import Atoms
from ase.neighborlist import NeighborList


def analyze_connectivity(atoms: Atoms, cutoff_factor: float = 1.2) -> List[List[int]]:
    """
    Analyze atomic connectivity and return disconnected fragments
    
    Args:
        atoms: ASE Atoms object
        cutoff_factor: Multiplier for covalent radii to determine bonds
        
    Returns:
        List of lists, each containing atom indices for connected fragments
    """
    from ase.data import covalent_radii
    
    # Build neighbor list
    cutoffs = []
    for symbol in atoms.get_chemical_symbols():
        atomic_num = atoms.numbers[atoms.symbols.index(symbol)]
        cutoffs.append(covalent_radii[atomic_num] * cutoff_factor)
    
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    
    # Find connected components using DFS
    visited = [False] * len(atoms)
    fragments = []
    
    for i in range(len(atoms)):
        if not visited[i]:
            fragment = []
            stack = [i]
            
            while stack:
                current = stack.pop()
                if not visited[current]:
                    visited[current] = True
                    fragment.append(current)
                    
                    # Add neighbors to stack
                    indices, offsets = nl.get_neighbors(current)
                    for neighbor in indices:
                        if not visited[neighbor]:
                            stack.append(neighbor)
            
            fragments.append(fragment)
    
    return fragments


def detect_artifacts(initial_atoms: Atoms, final_atoms: Atoms, 
                    displacement_threshold: float = 3.0,
                    min_artifact_size: int = 1,
                    max_artifact_size: int = 20) -> List[Dict[str, Any]]:
    """
    Detect artifacts (separated molecules/fragments) after calculation
    
    Args:
        initial_atoms: Initial atomic configuration
        final_atoms: Final atomic configuration
        displacement_threshold: Minimum displacement to consider an artifact (Å)
        min_artifact_size: Minimum number of atoms for an artifact
        max_artifact_size: Maximum number of atoms for an artifact
        
    Returns:
        List of artifact dictionaries with metadata
    """
    artifacts = []
    
    # Connectivity-based detection
    fragments = analyze_connectivity(final_atoms)
    
    # Filter fragments by size
    potential_artifacts = [
        frag for frag in fragments 
        if min_artifact_size <= len(frag) <= max_artifact_size
    ]
    
    # Analyze displacement for each fragment
    initial_positions = initial_atoms.get_positions()
    final_positions = final_atoms.get_positions()
    
    for i, fragment in enumerate(potential_artifacts):
        fragment_positions = final_positions[fragment]
        initial_fragment_positions = initial_positions[fragment]
        
        # Calculate center of mass displacement
        com_initial = np.mean(initial_fragment_positions, axis=0)
        com_final = np.mean(fragment_positions, axis=0)
        displacement = np.linalg.norm(com_final - com_initial)
        
        if displacement > displacement_threshold:
            symbols = [final_atoms.get_chemical_symbols()[idx] for idx in fragment]
            
            artifact = {
                'id': i,
                'atom_indices': fragment,
                'symbols': symbols,
                'formula': Atoms(symbols=symbols).get_chemical_formula(),
                'size': len(fragment),
                'displacement': displacement,
                'com_position': com_final,
                'atoms': final_atoms[fragment]
            }
            artifacts.append(artifact)
    
    return artifacts


def extract_fragment(atoms: Atoms, indices: List[int]) -> Atoms:
    """Extract a fragment as a new Atoms object"""
    fragment_atoms = atoms[indices]
    
    # Center the fragment and set reasonable cell
    fragment_atoms.center(vacuum=5.0)
    
    return fragment_atoms


def classify_artifacts(artifacts: List[Dict[str, Any]]) -> Dict[str, List[Dict[str, Any]]]:
    """Classify artifacts by type (molecules, atoms, etc.)"""
    classification = {
        'atoms': [],
        'small_molecules': [],
        'large_fragments': []
    }
    
    for artifact in artifacts:
        if artifact['size'] == 1:
            classification['atoms'].append(artifact)
        elif artifact['size'] <= 10:
            classification['small_molecules'].append(artifact)
        else:
            classification['large_fragments'].append(artifact)
    
    return classification


def save_artifacts(artifacts: List[Dict[str, Any]], output_dir: str):
    """Save detected artifacts as separate files"""
    from pathlib import Path
    from ase.io import write
    
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    for artifact in artifacts:
        filename = f"artifact_{artifact['id']}_{artifact['formula']}.xyz"
        filepath = output_path / filename
        
        write(str(filepath), artifact['atoms'])
        print(f"Saved artifact: {filename}")
    
    return len(artifacts)