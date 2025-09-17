#!/usr/bin/env python3
"""
Build hexamethyl tungsten W(CH3)6 structure
Octahedral geometry with 6 methyl groups around tungsten center
"""

import numpy as np
from ase import Atoms
from ase.io import write

def build_hexamethyl_tungsten():
    """Build W(CH3)6 with octahedral geometry"""
    
    # W-C bond length (approximate, will be optimized)
    w_c_distance = 2.2  # Angstrom
    
    # C-H bond length in methyl groups
    c_h_distance = 1.09  # Angstrom
    
    # Octahedral directions (±x, ±y, ±z)
    directions = np.array([
        [1, 0, 0],   # +x
        [-1, 0, 0],  # -x
        [0, 1, 0],   # +y
        [0, -1, 0],  # -y
        [0, 0, 1],   # +z
        [0, 0, -1]   # -z
    ])
    
    # Start with tungsten at origin
    positions = [[0, 0, 0]]
    symbols = ['W']
    
    # Add methyl groups in octahedral positions
    for i, direction in enumerate(directions):
        # Carbon position
        c_pos = direction * w_c_distance
        positions.append(c_pos)
        symbols.append('C')
        
        # Three hydrogen atoms around each carbon
        # Create tetrahedral arrangement for CH3
        if abs(direction[2]) > 0.9:  # z-direction methyl groups
            # Use x and y directions for H atoms
            h_dirs = np.array([
                [1, 0, 0],
                [-0.5, np.sqrt(3)/2, 0],
                [-0.5, -np.sqrt(3)/2, 0]
            ])
        else:
            # For x and y direction methyls, use z and perpendicular directions
            if abs(direction[0]) > 0.9:  # x-direction
                h_dirs = np.array([
                    [0, 1, 0],
                    [0, -0.5, np.sqrt(3)/2],
                    [0, -0.5, -np.sqrt(3)/2]
                ])
            else:  # y-direction
                h_dirs = np.array([
                    [1, 0, 0],
                    [-0.5, 0, np.sqrt(3)/2],
                    [-0.5, 0, -np.sqrt(3)/2]
                ])
        
        # Add hydrogen atoms
        for h_dir in h_dirs:
            h_pos = c_pos + h_dir * c_h_distance
            positions.append(h_pos)
            symbols.append('H')
    
    # Create ASE atoms object
    atoms = Atoms(symbols=symbols, positions=positions)
    
    return atoms

if __name__ == "__main__":
    # Build the molecule
    molecule = build_hexamethyl_tungsten()
    
    # Save in multiple formats
    write('hexamethyl_tungsten.xyz', molecule, format='xyz')
    write('hexamethyl_tungsten.pdb', molecule, format='proteindatabank')
    
    print(f"Built hexamethyl tungsten W(CH3)6 with {len(molecule)} atoms")
    print(f"Formula: {molecule.get_chemical_formula()}")
    print("Files created:")
    print("- hexamethyl_tungsten.xyz")
    print("- hexamethyl_tungsten.pdb")