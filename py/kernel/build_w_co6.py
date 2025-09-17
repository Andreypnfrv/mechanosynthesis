#!/usr/bin/env python3
"""
Build tungsten hexacarbonyl W(CO)6 with experimental geometry
Based on NIST B3LYP/GENECP optimized structure
"""

import os
import sys
import numpy as np
from ase import Atoms
from ase.io import write

def build_w_co6_from_sdf(sdf_file):
    """Convert NIST SDF coordinates to ASE Atoms object"""
    
    # Read SDF file and extract coordinates
    with open(sdf_file, 'r') as f:
        lines = f.readlines()
    
    # Find coordinate block (starts after counts line)
    coord_start = None
    for i, line in enumerate(lines):
        if ' V2000' in line:
            coord_start = i + 1
            break
    
    if coord_start is None:
        raise ValueError("Could not find V2000 block in SDF file")
    
    # Parse atoms from coordinate block
    symbols = []
    positions = []
    
    # Count line tells us how many atoms
    count_line = lines[coord_start - 1].strip()
    n_atoms = int(count_line.split()[0])
    
    for i in range(n_atoms):
        line = lines[coord_start + i].strip().split()
        x, y, z = float(line[0]), float(line[1]), float(line[2])
        symbol = line[3]
        
        symbols.append(symbol)
        positions.append([x, y, z])
    
    # Create ASE atoms object
    atoms = Atoms(symbols=symbols, positions=positions)
    
    return atoms

def build_w_co6_experimental():
    """Build W(CO)6 with experimental bond lengths from literature"""
    
    # Experimental values: W-C = 2.058 Å, C-O = 1.148 Å
    w_c_dist = 2.058  # Angstrom
    c_o_dist = 1.148  # Angstrom
    
    symbols = ['W']
    positions = [[0.0, 0.0, 0.0]]  # W at origin
    
    # Octahedral CO positions (+/-x, +/-y, +/-z)
    co_directions = [
        [1, 0, 0], [-1, 0, 0],  # x-axis
        [0, 1, 0], [0, -1, 0],  # y-axis  
        [0, 0, 1], [0, 0, -1]   # z-axis
    ]
    
    for direction in co_directions:
        # Carbon position
        c_pos = [w_c_dist * direction[i] for i in range(3)]
        symbols.append('C')
        positions.append(c_pos)
        
        # Oxygen position (further along same direction)
        total_dist = w_c_dist + c_o_dist
        o_pos = [total_dist * direction[i] for i in range(3)]
        symbols.append('O')
        positions.append(o_pos)
    
    atoms = Atoms(symbols=symbols, positions=positions)
    return atoms

def main():
    """Build W(CO)6 molecules and save in multiple formats"""
    
    output_dir = "data/w_co6_validation/inputs"
    os.makedirs(output_dir, exist_ok=True)
    
    print("Building W(CO)6 molecules...")
    
    # Method 1: From NIST SDF file
    sdf_file = "/tmp/w_co6.sdf"
    if os.path.exists(sdf_file):
        try:
            atoms_nist = build_w_co6_from_sdf(sdf_file)
            
            # Save NIST structure
            write(f"{output_dir}/w_co6_nist.xyz", atoms_nist)
            write(f"{output_dir}/w_co6_nist.gen", atoms_nist)
            print(f"✓ NIST B3LYP/GENECP structure saved")
            print(f"  Atoms: {len(atoms_nist)}")
            
        except Exception as e:
            print(f"✗ Failed to parse NIST SDF: {e}")
    
    # Method 2: Experimental ideal octahedral geometry
    atoms_exp = build_w_co6_experimental()
    
    write(f"{output_dir}/w_co6_experimental.xyz", atoms_exp)
    write(f"{output_dir}/w_co6_experimental.gen", atoms_exp)
    print(f"✓ Experimental geometry structure saved")
    print(f"  W-C: 2.058 Å, C-O: 1.148 Å")
    print(f"  Atoms: {len(atoms_exp)}")
    
    # Analysis
    print("\n=== Structure Comparison ===")
    if 'atoms_nist' in locals():
        nist_wc_bonds = []
        exp_wc_bonds = []
        
        for atoms, name in [(atoms_nist, "NIST"), (atoms_exp, "Experimental")]:
            w_idx = 0  # W is first atom
            c_indices = [i for i, s in enumerate(atoms.get_chemical_symbols()) if s == 'C']
            
            bonds = []
            for c_idx in c_indices:
                dist = atoms.get_distance(w_idx, c_idx)
                bonds.append(dist)
            
            avg_bond = np.mean(bonds)
            std_bond = np.std(bonds)
            
            print(f"{name:12s}: W-C = {avg_bond:.3f} ± {std_bond:.3f} Å")
            
            if name == "NIST":
                nist_wc_bonds = bonds
            else:
                exp_wc_bonds = bonds
    
    print(f"\nFiles saved to: {output_dir}/")
    print("Ready for method validation calculations!")

if __name__ == "__main__":
    main()