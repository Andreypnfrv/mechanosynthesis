from ase.build import molecule
from kernel.project_structure import ase_to_xyz_input
from datetime import datetime


def build_molecule(molecule_name, project_name=None, vacuum=6.0):
    """
    Build a molecule using ASE and save it to project structure
    
    Args:
        molecule_name (str): Name of the molecule (e.g., "C6H6", "H2O", "CH4")
        project_name (str): Optional project name, defaults to current date + molecule name
        vacuum (float): Vacuum space around molecule in Angstroms
    
    Returns:
        ase.Atoms: The built molecule
    """
    try:
        mol = molecule(molecule_name, vacuum=vacuum)
    except KeyError:
        raise ValueError(f"ASE doesn't have this structure: {molecule_name}")
    
    if project_name is None:
        date_str = "2025-08-23"
        project_name = f"{date_str}_{molecule_name.lower()}"
    
    ase_to_xyz_input(mol, project_name, molecule_name.lower())
    
    return mol


if __name__ == "__main__":
    # Example usage - build benzene
    benzene = build_molecule("C6H6", "2025-08-22_benzene", 6.0)