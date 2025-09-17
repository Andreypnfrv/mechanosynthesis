from ase.build import diamond111
from ase.neighborlist import NeighborList
import numpy as np


def find_atoms_with_hanging_bonds(atoms, cutoff=1.8, expected_coordination=None):
    """
    Find atoms with hanging bonds (fewer neighbors than expected).
    
    Args:
        atoms: ASE Atoms object
        cutoff: Distance cutoff for neighbor detection
        expected_coordination: Dict mapping atomic numbers to expected coordination numbers.
                             If None, uses default values for common elements.
    """
    if expected_coordination is None:
        # Default coordination numbers for common elements
        expected_coordination = {
            1: 1,   # H
            6: 4,   # C
            14: 4,  # Si
            32: 4,  # Ge
            7: 3,   # N
            8: 2,   # O
            15: 3,  # P
            16: 2,  # S
        }
    
    cutoffs = [cutoff/2] * len(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    
    hanging_atoms = []
    for i, atom in enumerate(atoms):
        atomic_number = atom.number
        if atomic_number in expected_coordination:
            indices, offsets = nl.get_neighbors(i)
            num_neighbors = len(indices)
            expected = expected_coordination[atomic_number]
            if num_neighbors < expected:
                hanging_atoms.append(i)
    
    return hanging_atoms


def build_diamond_surface(size, highlight=False, output_file='diamond_surface.xyz', passivate_element=None):
    """Build diamond surface with boundaries and coordination highlighting"""
    try:
        size_int = int(size)
        if size_int <= 0:
            raise ValueError("Size must be positive")
    except ValueError:
        print(f"Error: Invalid size '{size}'. Please provide a positive integer.")
        return False
    
    # Create diamond structure with vacuum and no PBC for surface analysis
    diamond = diamond111('C', size=(size_int, size_int, size_int), vacuum=6.0)
    diamond.set_pbc([False, False, False])  # No PBC for proper surface analysis
    print(f"Created diamond surface with {len(diamond)} atoms")
    
    # Apply hydrogen passivation if requested
    if passivate_element:
        hanging_atoms = find_atoms_with_hanging_bonds(diamond)
        if hanging_atoms:
            print(f"Found {len(hanging_atoms)} atoms with hanging bonds - passivating with {passivate_element}")
            from .build_diamond import create_passivated_structure
            diamond = create_passivated_structure(diamond, hanging_atoms, passivate_element)
        else:
            print("No hanging bonds found - structure is fully coordinated")
    
    if highlight:
        # Find atoms with hanging bonds using no PBC
        hanging_atoms = find_atoms_with_hanging_bonds(diamond)
        print(f"Found {len(hanging_atoms)} atoms with hanging bonds")
        
        # Create a copy for highlighting
        highlighted_diamond = diamond.copy()
        
        # Add RGB color columns based on coordination
        cutoffs = [0.9] * len(highlighted_diamond)
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(highlighted_diamond)
        
        # Initialize RGB arrays
        rgb_colors = []
        
        for i, atom in enumerate(highlighted_diamond):
            indices, offsets = nl.get_neighbors(i)
            num_neighbors = len(indices)
            
            if num_neighbors == 4:
                # Fully coordinated - Green (0, 255, 0)
                rgb_colors.append([0, 255, 0])
            elif num_neighbors == 1:
                # 1 bond, 3 missing - Yellow (255, 255, 0)
                rgb_colors.append([255, 255, 0])
            elif num_neighbors == 0:
                # No bonds - Red (255, 0, 0)
                rgb_colors.append([255, 0, 0])
            else:
                # Other cases - Orange (255, 165, 0)
                rgb_colors.append([255, 165, 0])
        
        # Add RGB as arrays to the atoms object
        rgb_array = np.array(rgb_colors)
        highlighted_diamond.set_array('R', rgb_array[:, 0])
        highlighted_diamond.set_array('G', rgb_array[:, 1])
        highlighted_diamond.set_array('B', rgb_array[:, 2])
        
        # Save the highlighted structure
        from ase.io import write
        write(output_file, highlighted_diamond, format='extxyz')
        print(f"Diamond surface saved to '{output_file}' with RGB color highlighting")
        print("Colors: Green=4 bonds, Yellow=1 bond, Red=0 bonds, Orange=2-3 bonds")
        
        # Print coordination statistics
        coord_counts = {}
        for i, atom in enumerate(highlighted_diamond):
            indices, offsets = nl.get_neighbors(i)
            num_neighbors = len(indices)
            coord_counts[num_neighbors] = coord_counts.get(num_neighbors, 0) + 1
        
        print("Coordination statistics:")
        for coord, count in sorted(coord_counts.items()):
            print(f"  {coord} bonds: {count} atoms")
        
    else:
        # Save without highlighting
        from ase.io import write
        write(output_file, diamond)
        print(f"Diamond surface saved to '{output_file}'")
    
    return True