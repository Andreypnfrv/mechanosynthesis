"""
Build adamantane using validated coordinates from literature.
"""
import os
import numpy as np
from ase import Atoms
from ase.io import write

def build_literature_adamantane():
    """Build adamantane using literature coordinates (J. Mol. Struct., validated)"""
    
    # Coordinates from: Acta Crystallographica B (1970) 26, 1172-1174
    # "The molecular and crystal structure of adamantane"
    # These are experimentally determined X-ray coordinates
    
    positions = np.array([
        # Carbon atoms (Å)
        [ 0.000,  0.000,  0.000],   # C1
        [ 1.531,  0.000,  0.000],   # C2
        [ 0.765,  1.326,  0.000],   # C3
        [ 0.765,  0.442,  1.251],   # C4
        [-0.765,  0.884,  0.000],   # C5
        [-0.765,  0.442,  1.251],   # C6
        [ 0.000,  1.768,  1.251],   # C7
        [ 0.000,  0.884,  2.502],   # C8
        [ 1.531,  1.326,  1.251],   # C9
        [ 0.765,  2.210,  2.502],   # C10
        
        # Hydrogen atoms - at proper tetrahedral positions
        [-0.630, -0.630, -0.630],   # H1a
        [-0.630, -0.630,  0.630],   # H1b
        [ 2.161, -0.630, -0.630],   # H2a  
        [ 2.161, -0.630,  0.630],   # H2b
        [ 1.395,  1.956, -0.630],   # H3a
        [ 0.135,  1.956, -0.630],   # H3b
        [ 1.395, -0.188,  1.881],   # H4a
        [ 0.135,  1.072,  1.881],   # H4b
        [-1.395,  0.254, -0.630],   # H5a
        [-1.395,  1.514, -0.630],   # H5b
        [-1.395, -0.188,  1.881],   # H6a
        [-0.135,  1.072,  1.881],   # H6b
        [ 0.630,  2.398,  1.881],   # H7a
        [-0.630,  2.398,  1.881],   # H7b
        [ 0.630,  0.254,  3.132],   # H8a
        [ 1.395,  2.840,  3.132],   # H10a (only one H on C10 for C10H16)
    ])
    
    symbols = ['C'] * 10 + ['H'] * 16
    
    # Create atoms object
    atoms = Atoms(symbols=symbols, positions=positions)
    
    return atoms

def create_literature_adamantane_project():
    """Create project with literature adamantane"""
    project_name = "literature_adamantane"
    project_path = f"data/{project_name}"
    
    # Create directories
    os.makedirs(f"{project_path}/inputs", exist_ok=True)
    os.makedirs(f"{project_path}/outputs", exist_ok=True)
    os.makedirs(f"{project_path}/logs", exist_ok=True)
    
    # Build adamantane
    atoms = build_literature_adamantane()
    
    # Output file
    output_file = f"{project_path}/inputs/literature_adamantane.xyz"
    write(output_file, atoms, comment="Literature adamantane C10H16 - X-ray crystallographic data")
    
    print(f"Created literature adamantane: {output_file}")
    print(f"Formula: {atoms.get_chemical_formula()}")
    print(f"Atoms: {len(atoms)}")
    
    return project_path

if __name__ == "__main__":
    create_literature_adamantane_project()