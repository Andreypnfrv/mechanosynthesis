#!/usr/bin/env python3
import os
from pathlib import Path
from datetime import datetime
from ase.io import write
import sys
sys.path.append(str(Path(__file__).parent.parent))
import settings


def create_project_structure(name):
    """Create a complete project structure for the given name"""
    
    project_paths = settings.ensure_data_structure(name)
    
    # Additional directories for compatibility
    additional_dirs = ["misc", "src"]
    
    for dir_name in additional_dirs:
        additional_path = project_paths['project'] / dir_name
        additional_path.mkdir(exist_ok=True)
        print(f"Created directory: {dir_name}")
        
    print(f"Project structure created at: {project_paths['project']}")
    return str(project_paths['project']) + "/"



def ase_to_xyz_input(atoms, project_folder, filename=None):
    """Convert ASE Atoms object to .xyz file in input directory
    
    Args:
        atoms: ASE Atoms object (or compatible structure)
        project_name: Name of the project (used to determine input directory)
        filename: Optional filename for the .xyz file (default: structure.xyz)
    
    Returns:
        str: Path to the created .xyz file
    """
    if filename is None:
        filename = "structure.xyz"
    
    # Ensure filename has .xyz extension
    if not filename.endswith('.xyz'):
        filename += '.xyz'
    
    # Get project paths using settings
    project_path = settings.get_project_data_path(project_folder)
    input_dir = project_path / "inputs"
    input_dir.mkdir(parents=True, exist_ok=True)
    
    # Full path for the output file
    output_path = input_dir / filename
    
    # Write the structure to .xyz file
    write(str(output_path), atoms, format='xyz')
    
    print(f"Wrote structure to: {output_path}")
    return str(output_path)
