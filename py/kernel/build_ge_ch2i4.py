#!/usr/bin/env python3
"""
Build tetrakis(iodomethyl)germane Ge(CH2I)4 molecule.
Molecular formula: C4H8GeI4
Structure: Tetrahedral Ge center with 4 CH2I groups
"""

import numpy as np
from ase import Atoms
from ase.io import write


def build_ge_ch2i4():
    """Build Ge(CH2I)4 molecule with tetrahedral geometry around Ge center."""
    
    # Bond lengths (Angstroms) - typical values for organo-germanium compounds
    ge_c_bond = 1.95   # Ge-C bond length
    c_h_bond = 1.09    # C-H bond length  
    c_i_bond = 2.14    # C-I bond length
    
    # Tetrahedral angles
    tet_angle = np.radians(109.47)  # tetrahedral angle
    
    # Start with Ge at origin
    positions = []
    symbols = []
    
    # Central Ge atom
    positions.append([0.0, 0.0, 0.0])
    symbols.append('Ge')
    
    # Tetrahedral directions from center
    # Standard tetrahedral vertices
    directions = np.array([
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1], 
        [-1, -1, 1]
    ])
    
    # Normalize to unit vectors
    directions = directions / np.linalg.norm(directions, axis=1)[:, np.newaxis]
    
    # Build 4 CH2I groups
    for i, direction in enumerate(directions):
        # Carbon position
        c_pos = direction * ge_c_bond
        positions.append(c_pos)
        symbols.append('C')
        
        # For each CH2I group, we need 2 H atoms and 1 I atom
        # Create local coordinate system for this CH2I group
        # Use direction as one axis
        z_local = direction
        
        # Create perpendicular vectors for H and I placement
        # Find a vector not parallel to z_local
        if abs(z_local[2]) < 0.9:
            x_temp = np.array([0, 0, 1])
        else:
            x_temp = np.array([1, 0, 0])
            
        # Create perpendicular x_local
        x_local = np.cross(z_local, x_temp)
        x_local = x_local / np.linalg.norm(x_local)
        
        # Create y_local
        y_local = np.cross(z_local, x_local)
        
        # Place H atoms in tetrahedral arrangement around C
        # H-C-Ge angle ~109.47°, H-C-I angle ~109.47°
        
        # Two H atoms positioned symmetrically
        h_angle = np.radians(107)  # H-C-H angle in CH2I
        
        # H1 position - slightly towards Ge direction but offset
        h1_dir = -0.3 * z_local + 0.7 * x_local
        h1_dir = h1_dir / np.linalg.norm(h1_dir)
        h1_pos = c_pos + h1_dir * c_h_bond
        positions.append(h1_pos)
        symbols.append('H')
        
        # H2 position - opposite side
        h2_dir = -0.3 * z_local - 0.7 * x_local  
        h2_dir = h2_dir / np.linalg.norm(h2_dir)
        h2_pos = c_pos + h2_dir * c_h_bond
        positions.append(h2_pos)
        symbols.append('H')
        
        # Iodine position - opposite to Ge
        i_dir = -0.6 * z_local + 0.4 * y_local
        i_dir = i_dir / np.linalg.norm(i_dir)
        i_pos = c_pos + i_dir * c_i_bond
        positions.append(i_pos)
        symbols.append('I')
    
    # Create ASE atoms object
    molecule = Atoms(symbols=symbols, positions=positions)
    
    return molecule


def main():
    """Main function to build and save Ge(CH2I)4 molecule."""
    # Build molecule
    molecule = build_ge_ch2i4()
    
    # Print info
    print(f"Built Ge(CH2I)4 molecule:")
    print(f"Formula: {molecule.get_chemical_formula()}")
    print(f"Number of atoms: {len(molecule)}")
    print(f"Symbols: {molecule.get_chemical_symbols()}")
    
    return molecule


if __name__ == "__main__":
    molecule = main()
    
    # Save to XYZ format
    output_file = "ge_ch2i4.xyz"
    write(output_file, molecule)
    print(f"Saved to {output_file}")