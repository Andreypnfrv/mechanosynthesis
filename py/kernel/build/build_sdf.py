from ase.io import read
from kernel.project_structure import create_project_structure, ase_to_xyz_input
from datetime import datetime
import os
import shutil


def build_from_sdf(sdf_file_path, project_name=None):
    """
    Build a project from an SDF file and convert to XYZ format
    
    Args:
        sdf_file_path (str): Path to the SDF file
        project_name (str): Optional project name, defaults to SDF filename without extension
    
    Returns:
        ase.Atoms: The built molecule
    """
    if not os.path.exists(sdf_file_path):
        raise FileNotFoundError(f"SDF file not found: {sdf_file_path}")
    
    # Read the SDF file using ASE
    try:
        mol = read(sdf_file_path)
    except Exception as e:
        raise ValueError(f"Failed to read SDF file {sdf_file_path}: {e}")
    
    # Generate project name if not provided
    if project_name is None:
        base_name = os.path.splitext(os.path.basename(sdf_file_path))[0]
        project_name = base_name.lower()
    
    # Create project structure
    create_project_structure(project_name)
    
    # Copy original SDF file to src directory
    src_dir = f"data/{project_name}/src"
    sdf_filename = os.path.basename(sdf_file_path)
    src_sdf_path = os.path.join(src_dir, sdf_filename)
    shutil.copy2(sdf_file_path, src_sdf_path)
    print(f"Copied SDF file to: {src_sdf_path}")
    
    # Convert to XYZ and save to project inputs
    xyz_filename = f"{os.path.splitext(os.path.basename(sdf_file_path))[0]}.xyz"
    ase_to_xyz_input(mol, project_name, xyz_filename)
    
    print(f"Created project '{project_name}' from SDF file: {sdf_file_path}")
    
    return mol


if __name__ == "__main__":
    # Example usage
    sdf_file = "Conformer3D_COMPOUND_CID_9238.sdf"
    if os.path.exists(sdf_file):
        mol = build_from_sdf(sdf_file)