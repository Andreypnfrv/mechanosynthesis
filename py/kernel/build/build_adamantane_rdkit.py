"""
Build proper adamantane structure using RDKit for mechanosynthesis research.
"""
import os
import numpy as np
from ase import Atoms
from ase.io import write

def build_adamantane_rdkit():
    """Build proper adamantane (C10H16) structure using RDKit"""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdMolDescriptors
    except ImportError:
        print("RDKit not available. Installing...")
        os.system("conda install -c conda-forge rdkit -y")
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from rdkit.Chem import rdMolDescriptors
    
    # Create adamantane molecule from SMILES
    smiles = "C12CC3CC(C1)CC(C3)C2"  # Adamantane SMILES
    mol = Chem.MolFromSmiles(smiles)
    
    # Add hydrogens
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.OptimizeMolecule(mol, maxIters=1000)
    
    # Extract coordinates
    conf = mol.GetConformer()
    positions = []
    symbols = []
    
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        positions.append([pos.x, pos.y, pos.z])
        symbols.append(atom.GetSymbol())
    
    # Create ASE atoms object
    atoms = Atoms(symbols=symbols, positions=positions)
    
    return atoms

def create_rdkit_adamantane_project():
    """Create project folder with RDKit-generated adamantane structure"""
    project_name = "rdkit_adamantane" 
    project_path = f"data/{project_name}"
    
    # Create directories
    os.makedirs(f"{project_path}/inputs", exist_ok=True)
    os.makedirs(f"{project_path}/outputs", exist_ok=True)
    os.makedirs(f"{project_path}/logs", exist_ok=True)
    
    # Build and save adamantane
    atoms = build_adamantane_rdkit()
    
    # Write to inputs folder
    output_file = f"{project_path}/inputs/rdkit_adamantane.xyz"
    write(output_file, atoms, comment="RDKit-generated adamantane C10H16 - properly optimized structure")
    
    print(f"Created RDKit adamantane structure: {output_file}")
    print(f"Formula: {atoms.get_chemical_formula()}")
    print(f"Number of atoms: {len(atoms)}")
    
    return project_path

if __name__ == "__main__":
    create_rdkit_adamantane_project()