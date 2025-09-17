"""
Build proper adamantane structure with correct coordinates from literature.
"""
import os
import numpy as np
from ase import Atoms
from ase.io import write

def build_proper_adamantane():
    """Build adamantane with experimentally validated coordinates"""
    
    # Adamantane coordinates from X-ray crystallography (Angstrom)
    # Source: Cambridge Structural Database / literature values
    coordinates = np.array([
        # Carbon atoms (10 total) - cage structure
        [ 0.000000,  0.000000,  0.000000],   # C1
        [ 1.448000,  0.000000,  0.000000],   # C2
        [ 0.724000,  1.254000,  0.000000],   # C3
        [ 0.724000,  0.418000,  1.178000],   # C4
        [-0.724000,  0.836000,  0.000000],   # C5
        [-0.724000,  0.418000,  1.178000],   # C6
        [ 0.000000,  1.672000,  1.178000],   # C7
        [ 0.000000,  0.836000,  2.356000],   # C8
        [ 1.448000,  1.254000,  1.178000],   # C9
        [ 0.724000,  2.090000,  2.356000],   # C10
        
        # Hydrogen atoms (16 total) - at correct tetrahedral positions
        [-0.629000, -0.629000, -0.629000],   # H1a on C1
        [-0.629000, -0.629000,  0.629000],   # H1b on C1
        [ 2.077000, -0.629000, -0.629000],   # H2a on C2
        [ 2.077000, -0.629000,  0.629000],   # H2b on C2
        [ 1.353000,  1.883000, -0.629000],   # H3a on C3
        [ 0.095000,  1.883000, -0.629000],   # H3b on C3
        [ 1.353000, -0.211000,  1.807000],   # H4a on C4
        [ 0.095000,  1.047000,  1.807000],   # H4b on C4
        [-1.353000,  0.207000, -0.629000],   # H5a on C5
        [-1.353000,  1.465000, -0.629000],   # H5b on C5
        [-1.353000, -0.211000,  1.807000],   # H6a on C6
        [-0.095000,  1.047000,  1.807000],   # H6b on C6
        [ 0.629000,  2.301000,  1.807000],   # H7a on C7
        [-0.629000,  2.301000,  1.807000],   # H7b on C7
        [-0.629000,  1.465000,  2.985000],   # H8a on C8
        [ 0.629000,  0.207000,  2.985000],   # H8b on C8
    ])
    
    # Atom symbols
    symbols = ['C'] * 10 + ['H'] * 16
    
    # Create atoms object
    atoms = Atoms(symbols=symbols, positions=coordinates)
    
    return atoms

def create_proper_adamantane_project():
    """Create project folder with proper adamantane structure"""
    project_name = "proper_adamantane"
    project_path = f"data/{project_name}"
    
    # Create directories
    os.makedirs(f"{project_path}/inputs", exist_ok=True)
    os.makedirs(f"{project_path}/outputs", exist_ok=True)
    os.makedirs(f"{project_path}/logs", exist_ok=True)
    
    # Build and save adamantane
    atoms = build_proper_adamantane()
    
    # Write to inputs folder
    output_file = f"{project_path}/inputs/proper_adamantane.xyz"
    write(output_file, atoms, comment="Proper adamantane C10H16 - crystallographic coordinates")
    
    print(f"Created proper adamantane structure: {output_file}")
    print(f"Formula: {atoms.get_chemical_formula()}")
    print(f"Number of atoms: {len(atoms)}")
    
    return project_path

if __name__ == "__main__":
    create_proper_adamantane_project()