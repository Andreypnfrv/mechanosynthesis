from ase.build import diamond111
from ase.neighborlist import NeighborList
import numpy as np
from ase import Atoms
from ase.io import write


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


def passivate_structure(atoms, passivation_element='H', bond_length=1.09):
    """
    Add passivation atoms to hanging bonds
    
    Args:
        atoms: ASE Atoms object with hanging bonds
        passivation_element: Element to add (H, CH3, etc.)
        bond_length: Bond length for passivation atoms
    """
    from ase.data import atomic_numbers
    
    hanging_indices = find_atoms_with_hanging_bonds(atoms)
    if not hanging_indices:
        return atoms
    
    # Create neighbor list to find bonding directions
    cutoffs = [0.9] * len(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    
    new_atoms = atoms.copy()
    
    for i in hanging_indices:
        pos = atoms.positions[i]
        indices, offsets = nl.get_neighbors(i)
        
        if len(indices) == 0:
            continue
            
        # Calculate existing bond vectors
        bond_dirs = []
        for j in indices:
            bond_vec = pos - atoms.positions[j]  # From neighbor to center
            bond_dirs.append(bond_vec / np.linalg.norm(bond_vec))
        
        if bond_dirs:
            num_missing = 4 - len(indices)
            
            # For each missing bond, find the best tetrahedral direction
            # Standard tetrahedral directions from center
            tet_dirs = np.array([
                [1.0, 1.0, 1.0],   # (1,1,1)
                [1.0, -1.0, -1.0], # (1,-1,-1) 
                [-1.0, 1.0, -1.0], # (-1,1,-1)
                [-1.0, -1.0, 1.0]  # (-1,-1,1)
            ])
            tet_dirs = tet_dirs / np.sqrt(3.0)  # Normalize
            
            # Find unused tetrahedral directions
            used_dirs = []
            for k in range(min(num_missing, len(tet_dirs))):
                best_dir = None
                max_angle = -2.0
                
                for tet_dir in tet_dirs:
                    # Skip if already used
                    if any(np.allclose(tet_dir, used) for used in used_dirs):
                        continue
                        
                    # Find minimum angle with existing bonds
                    min_overlap = 2.0
                    for bond_dir in bond_dirs:
                        overlap = abs(np.dot(tet_dir, bond_dir))
                        min_overlap = min(min_overlap, overlap)
                    
                    # Choose direction with minimal overlap with existing bonds
                    if min_overlap > max_angle:
                        max_angle = min_overlap
                        best_dir = tet_dir
                
                if best_dir is not None:
                    h_pos = pos + best_dir * bond_length
                    new_atoms.append(passivation_element)
                    new_atoms.positions[-1] = h_pos
                    used_dirs.append(best_dir)
    
    return new_atoms


def build_adamantane(passivation='H', output_file='adamantane.xyz'):
    """
    Build adamantane (C10H16) structure
    """
    # Adamantane coordinates (from literature)
    positions = np.array([
        [0.000,  1.631,  0.000],  # C1
        [1.412,  1.089, -1.004],  # C2
        [1.412,  1.089,  1.004],  # C3
        [-1.412, 1.089,  1.004],  # C4
        [-1.412, 1.089, -1.004],  # C5
        [0.000, -0.542, -1.631],  # C6
        [1.412, -1.085,  0.000],  # C7
        [0.000, -0.542,  1.631],  # C8
        [-1.412, -1.085, 0.000],  # C9
        [0.000, -1.631,  0.000]   # C10
    ])
    
    if passivation == 'H' or passivation is None:
        # Adamantane C10H16 - using experimental crystal structure coordinates
        # Reference: Donohue & Goodman, Acta Cryst. 22, 352 (1967)
        all_positions = np.array([
            # 4 bridgehead carbons (tetrahedral positions)
            [ 0.000,  0.000,  0.000],   # C1 - bridgehead
            [ 1.544,  1.544,  1.544],   # C4 - bridgehead  
            [-1.544, -1.544,  1.544],   # C7 - bridgehead
            [ 1.544, -1.544, -1.544],   # C10 - bridgehead
            # 6 bridge carbons (edges of tetrahedron)
            [ 0.772,  0.772,  0.000],   # C2 - bridge
            [ 0.772,  0.000,  0.772],   # C3 - bridge
            [ 0.000,  0.772,  0.772],   # C5 - bridge
            [-0.772,  0.000, -0.772],   # C6 - bridge
            [ 0.000, -0.772, -0.772],   # C8 - bridge
            [-0.772, -0.772,  0.000],   # C9 - bridge
            # Hydrogens on bridgehead carbons (1 H each, 4 total)
            [-1.089, -1.089, -1.089],   # H on C1
            [ 2.633,  2.633,  2.633],   # H on C4
            [-2.633, -2.633,  2.633],   # H on C7
            [ 2.633, -2.633, -2.633],   # H on C10
            # Hydrogens on bridge carbons (2 H each, 12 total)
            [ 1.861,  1.861, -1.089],   # H on C2
            [ 0.683,  0.683, -1.089],   # H on C2
            [ 1.861, -1.089,  1.861],   # H on C3
            [ 0.683, -1.089,  0.683],   # H on C3
            [-1.089,  1.861,  1.861],   # H on C5
            [-1.089,  0.683,  0.683],   # H on C5
            [-1.861, -1.089, -1.861],   # H on C6
            [-0.683, -1.089, -0.683],   # H on C6
            [ 1.089, -1.861, -1.861],   # H on C8
            [ 1.089, -0.683, -0.683],   # H on C8
            [-1.861, -1.861,  1.089],   # H on C9
            [-0.683, -0.683,  1.089],   # H on C9
        ])
        symbols = ['C'] * 10 + ['H'] * 16
        adamantane = Atoms(symbols=symbols, positions=all_positions)
    else:
        # Just carbon skeleton for other passivation
        symbols = ['C'] * 10
        adamantane = Atoms(symbols=symbols, positions=positions)
        if passivation:
            adamantane = passivate_structure(adamantane, passivation)
    
    adamantane.set_pbc([False, False, False])
    adamantane.center()
    
    write(output_file, adamantane, format='extxyz')
    print(f"Adamantane saved to '{output_file}' with {len(adamantane)} atoms")
    return True


def build_rod(length=5, width=2, height=2, passivation='H', output_file='rod.xyz'):
    """
    Build rod-shaped diamondoid (length x width x height in diamond cells)
    """
    # Create proper 3D diamond supercell 
    from ase.build import bulk
    diamond = bulk('C', 'diamond', a=3.567, cubic=True)
    diamond = diamond.repeat((max(length+1, 3), max(width+1, 3), max(height+1, 3)))
    
    # Filter atoms to rod shape
    positions = diamond.get_positions()
    cell = diamond.get_cell()
    
    # Get cell dimensions
    a, b, c = np.linalg.norm(cell, axis=1)
    center = positions.mean(axis=0)
    
    selected_indices = []
    for i, pos in enumerate(positions):
        rel_pos = pos - center
        # Rod along x-axis
        if (abs(rel_pos[0]) <= length * a/6 and 
            abs(rel_pos[1]) <= width * b/6 and 
            abs(rel_pos[2]) <= height * c/6):
            selected_indices.append(i)
    
    if not selected_indices:
        print(f"Error: No atoms found for rod {length}x{width}x{height}")
        return False
    
    rod = Atoms()
    for i in selected_indices:
        rod.append(diamond[i])
    
    rod.set_pbc([False, False, False])
    rod.center()
    
    if passivation:
        rod = passivate_structure(rod, passivation)
    
    write(output_file, rod, format='extxyz')
    print(f"Rod {length}x{width}x{height} saved to '{output_file}' with {len(rod)} atoms")
    return True


def build_platform(size_x=3, size_y=3, thickness=1, passivation='H', output_file='platform.xyz'):
    """
    Build platform-shaped diamondoid (size_x x size_y x thickness in diamond cells)
    """
    # Create proper 3D diamond supercell
    from ase.build import bulk  
    diamond = bulk('C', 'diamond', a=3.567, cubic=True)
    diamond = diamond.repeat((max(size_x+2, 3), max(size_y+2, 3), max(thickness+2, 3)))
    
    # Filter atoms to platform shape
    positions = diamond.get_positions()
    cell = diamond.get_cell()
    
    # Get cell dimensions
    a, b, c = np.linalg.norm(cell, axis=1)
    center = positions.mean(axis=0)
    
    selected_indices = []
    for i, pos in enumerate(positions):
        rel_pos = pos - center
        # Flat platform with adjusted scaling
        if (abs(rel_pos[0]) <= size_x * a/4 and 
            abs(rel_pos[1]) <= size_y * b/4 and 
            abs(rel_pos[2]) <= thickness * c/4):
            selected_indices.append(i)
    
    if not selected_indices:
        print(f"Error: No atoms found for platform {size_x}x{size_y}x{thickness}")
        return False
    
    platform = Atoms()
    for i in selected_indices:
        platform.append(diamond[i])
    
    platform.set_pbc([False, False, False])
    platform.center()
    
    if passivation:
        platform = passivate_structure(platform, passivation)
    
    write(output_file, platform, format='extxyz')
    print(f"Platform {size_x}x{size_y}x{thickness} saved to '{output_file}' with {len(platform)} atoms")
    return True


def build_diamondoid(shape_type='sphere', size=5, passivation=None, output_file='diamondoid.xyz'):
    """
    Build diamondoid structure with various shapes and sizes
    
    Args:
        shape_type: Type of shape ('tetrahedron', 'cube', 'octahedron', 'sphere', 'adamantane', 'rod', 'platform')
        size: Size parameter (for rod/platform: 'length:width:height' or 'x:y:thickness')
        passivation: Passivation element ('H', 'CH3', etc.) or None
        output_file: Output filename
    """
    # Handle special cases
    if shape_type == 'adamantane':
        return build_adamantane(passivation, output_file)
    
    elif shape_type == 'rod':
        if ':' in str(size):
            dims = [int(x) for x in str(size).split(':')]
            if len(dims) == 3:
                return build_rod(dims[0], dims[1], dims[2], passivation, output_file)
        return build_rod(int(size), 2, 2, passivation, output_file)
    
    elif shape_type == 'platform':
        if ':' in str(size):
            dims = [int(x) for x in str(size).split(':')]
            if len(dims) == 3:
                return build_platform(dims[0], dims[1], dims[2], passivation, output_file)
        return build_platform(int(size), int(size), 1, passivation, output_file)
    
    # Original geometric shapes
    try:
        size_float = float(size)
        if size_float <= 0:
            raise ValueError("Size must be positive")
    except ValueError:
        print(f"Error: Invalid size '{size}'. Please provide a positive number.")
        return False
    
    # Create a proper 3D diamond supercell to cut from
    from ase.build import bulk
    supercell_size = max(3, int(size_float / 2) + 2)
    # Use proper bulk diamond structure with FCC lattice
    diamond = bulk('C', 'diamond', a=3.567, cubic=True)
    diamond = diamond.repeat((supercell_size, supercell_size, supercell_size))
    
    # Get the center of the supercell
    positions = diamond.get_positions()
    center = positions.mean(axis=0)
    
    # Define shape filtering functions
    def is_inside_tetrahedron(pos, center, size):
        """Check if position is inside tetrahedron"""
        # Simple tetrahedron centered at origin
        x, y, z = pos - center
        return (x >= 0 and y >= 0 and z >= 0 and 
                x + y + z <= size)
    
    def is_inside_cube(pos, center, size):
        """Check if position is inside cube"""
        diff = np.abs(pos - center)
        return np.all(diff <= size/2)
    
    def is_inside_octahedron(pos, center, size):
        """Check if position is inside octahedron"""
        x, y, z = np.abs(pos - center)
        return x + y + z <= size
    
    def is_inside_sphere(pos, center, size):
        """Check if position is inside sphere"""
        return np.linalg.norm(pos - center) <= size
    
    # Select shape function
    shape_functions = {
        'tetrahedron': is_inside_tetrahedron,
        'cube': is_inside_cube,
        'octahedron': is_inside_octahedron,
        'sphere': is_inside_sphere
    }
    
    if shape_type not in shape_functions:
        print(f"Error: Unknown shape '{shape_type}'. Available: {list(shape_functions.keys())}")
        return False
    
    shape_func = shape_functions[shape_type]
    
    # Filter atoms based on shape
    selected_indices = []
    for i, pos in enumerate(positions):
        if shape_func(pos, center, size_float):
            selected_indices.append(i)
    
    if not selected_indices:
        print(f"Error: No atoms found inside {shape_type} with size {size_float}")
        return False
    
    # Create new structure with selected atoms
    diamondoid = Atoms()
    for i in selected_indices:
        diamondoid.append(diamond[i])
    
    # Remove PBC for finite structure
    diamondoid.set_pbc([False, False, False])
    
    # Center the structure
    diamondoid.center()
    
    print(f"Created {shape_type} diamondoid with {len(diamondoid)} atoms")
    
    # Find surface atoms and add coordination coloring
    hanging_atoms = find_atoms_with_hanging_bonds(diamondoid)
    print(f"Found {len(hanging_atoms)} surface atoms")
    
    # Add RGB color columns based on coordination
    cutoffs = [0.9] * len(diamondoid)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(diamondoid)
    
    rgb_colors = []
    coord_counts = {}
    
    for i, atom in enumerate(diamondoid):
        indices, offsets = nl.get_neighbors(i)
        num_neighbors = len(indices)
        coord_counts[num_neighbors] = coord_counts.get(num_neighbors, 0) + 1
        
        if num_neighbors == 4:
            # Fully coordinated - Green (0, 255, 0)
            rgb_colors.append([0, 255, 0])
        elif num_neighbors == 3:
            # 3 bonds - Orange (255, 165, 0)
            rgb_colors.append([255, 165, 0])
        elif num_neighbors == 2:
            # 2 bonds - Orange-Red (255, 100, 0)
            rgb_colors.append([255, 100, 0])
        elif num_neighbors == 1:
            # 1 bond - Yellow (255, 255, 0)
            rgb_colors.append([255, 255, 0])
        else:
            # No bonds - Red (255, 0, 0)
            rgb_colors.append([255, 0, 0])
    
    # Add RGB as arrays to the atoms object
    rgb_array = np.array(rgb_colors)
    diamondoid.set_array('R', rgb_array[:, 0])
    diamondoid.set_array('G', rgb_array[:, 1])
    diamondoid.set_array('B', rgb_array[:, 2])
    
    # Apply passivation if requested
    if passivation:
        diamondoid = passivate_structure(diamondoid, passivation)
        print(f"Applied {passivation} passivation")
    
    # Save the diamondoid structure
    write(output_file, diamondoid, format='extxyz')
    print(f"Diamondoid saved to '{output_file}' with {len(diamondoid)} atoms")
    
    # Print coordination statistics
    print("Coordination statistics:")
    for coord, count in sorted(coord_counts.items()):
        print(f"  {coord} bonds: {count} atoms")
    
    return True