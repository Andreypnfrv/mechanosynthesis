#!/usr/bin/env python3
"""
Analyze W(CO)6 calculation results and compare with experimental data
"""

import numpy as np
from ase.io import read
import os

def analyze_structure(filename, label):
    """Analyze bond lengths in a structure"""
    atoms = read(filename)
    symbols = atoms.get_chemical_symbols()
    
    # Find W atom (should be first)
    w_idx = None
    for i, symbol in enumerate(symbols):
        if symbol == 'W':
            w_idx = i
            break
    
    if w_idx is None:
        print(f"Warning: No W atom found in {filename}")
        return None
    
    # Find C and O atoms
    c_indices = [i for i, s in enumerate(symbols) if s == 'C']
    o_indices = [i for i, s in enumerate(symbols) if s == 'O']
    
    # Calculate W-C bond lengths
    w_c_distances = []
    for c_idx in c_indices:
        dist = atoms.get_distance(w_idx, c_idx)
        w_c_distances.append(dist)
    
    # Calculate C-O bond lengths (each C should have one O)
    c_o_distances = []
    for c_idx in c_indices:
        c_pos = atoms.positions[c_idx]
        min_dist = float('inf')
        closest_o = None
        
        for o_idx in o_indices:
            o_pos = atoms.positions[o_idx]
            dist = np.linalg.norm(c_pos - o_pos)
            if dist < min_dist:
                min_dist = dist
                closest_o = o_idx
        
        if closest_o is not None:
            c_o_distances.append(min_dist)
    
    # Statistics
    w_c_mean = np.mean(w_c_distances)
    w_c_std = np.std(w_c_distances)
    c_o_mean = np.mean(c_o_distances)
    c_o_std = np.std(c_o_distances)
    
    print(f"\n{label}:")
    print(f"  W-C bonds: {w_c_mean:.3f} ± {w_c_std:.3f} Å")
    print(f"  C-O bonds: {c_o_mean:.3f} ± {c_o_std:.3f} Å")
    print(f"  Individual W-C distances: {[f'{d:.3f}' for d in w_c_distances]}")
    print(f"  Individual C-O distances: {[f'{d:.3f}' for d in c_o_distances]}")
    
    return {
        'w_c_mean': w_c_mean,
        'w_c_std': w_c_std,
        'c_o_mean': c_o_mean,
        'c_o_std': c_o_std,
        'w_c_distances': w_c_distances,
        'c_o_distances': c_o_distances
    }

def main():
    print("=== W(CO)6 Structure Analysis ===")
    
    # Experimental reference
    exp_w_c = 2.058  # Å
    exp_c_o = 1.148  # Å
    
    print(f"\nExperimental reference:")
    print(f"  W-C: {exp_w_c:.3f} Å")  
    print(f"  C-O: {exp_c_o:.3f} Å")
    
    # Analyze structures
    structures = []
    
    # Initial experimental geometry
    initial_exp = "inputs/w_co6_experimental.xyz"
    if os.path.exists(initial_exp):
        result = analyze_structure(initial_exp, "Initial experimental geometry")
        if result:
            structures.append(("Initial experimental", result))
    
    # Initial NIST geometry  
    initial_nist = "inputs/w_co6_nist.xyz"
    if os.path.exists(initial_nist):
        result = analyze_structure(initial_nist, "Initial NIST B3LYP/GENECP geometry")
        if result:
            structures.append(("Initial NIST", result))
    
    # xTB optimized structures
    xtb_exp = "outputs/xtb/w_co6_experimental_relaxed.xyz"
    if os.path.exists(xtb_exp):
        result = analyze_structure(xtb_exp, "xTB-optimized (from experimental)")
        if result:
            structures.append(("xTB (from exp)", result))
    
    xtb_nist = "outputs/xtb/w_co6_nist_relaxed.xyz"
    if os.path.exists(xtb_nist):
        result = analyze_structure(xtb_nist, "xTB-optimized (from NIST)")
        if result:
            structures.append(("xTB (from NIST)", result))
    
    # Comparison with experimental
    print(f"\n=== Comparison with Experimental Values ===")
    print(f"{'Method':<20} {'W-C Error':<12} {'C-O Error':<12} {'W-C %':<10} {'C-O %':<10}")
    print("-" * 64)
    
    for name, result in structures:
        w_c_error = result['w_c_mean'] - exp_w_c
        c_o_error = result['c_o_mean'] - exp_c_o
        w_c_percent = (w_c_error / exp_w_c) * 100
        c_o_percent = (c_o_error / exp_c_o) * 100
        
        print(f"{name:<20} {w_c_error:>+8.3f} Å   {c_o_error:>+8.3f} Å   {w_c_percent:>+6.1f}%   {c_o_percent:>+6.1f}%")
    
    # Check convergence
    print(f"\n=== Method Validation ===")
    for name, result in structures:
        w_c_std = result['w_c_std']
        symmetry_quality = "Perfect" if w_c_std < 0.001 else "Good" if w_c_std < 0.01 else "Poor"
        print(f"{name:<20}: W-C std = {w_c_std:.4f} Å ({symmetry_quality} octahedral symmetry)")

if __name__ == "__main__":
    os.chdir('/Users/andreypanferov/Documents/mechanosynthesis/data/w_co6')
    main()