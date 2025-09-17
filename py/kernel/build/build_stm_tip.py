from ase.build import bulk
from ase.neighborlist import NeighborList
from ase.io import write
from ase import Atoms
import numpy as np


def create_bcc_tungsten_supercell(size=(5, 5, 8)):
    """
    Create BCC tungsten supercell
    
    Args:
        size: Tuple (nx, ny, nz) for supercell dimensions
        
    Returns:
        ASE Atoms object with BCC tungsten structure
    """
    # BCC tungsten lattice parameter: 3.165 Å
    w_bulk = bulk('W', 'bcc', a=3.165, cubic=True)
    w_supercell = w_bulk.repeat(size)
    
    # Remove PBC for finite structure
    w_supercell.set_pbc([False, False, False])
    w_supercell.center()
    
    return w_supercell


def filter_atoms_to_tip_shape(atoms, height_ratio=0.8, tip_sharpness=3):
    """
    Filter atoms to create tip shape (pyramid/cone)
    
    Args:
        atoms: ASE Atoms object
        height_ratio: Fraction of total height to use for tip
        tip_sharpness: Sharpening factor (1-5, higher = sharper)
        
    Returns:
        List of indices for atoms to keep
    """
    positions = atoms.get_positions()
    center = positions.mean(axis=0)
    
    # Find z-extent
    z_min, z_max = positions[:, 2].min(), positions[:, 2].max()
    z_total = z_max - z_min
    z_tip_height = z_total * height_ratio
    
    # Find xy-extent for base radius calculation
    xy_positions = positions[:, :2]
    xy_center = center[:2]
    max_radius = np.max(np.linalg.norm(xy_positions - xy_center, axis=1))
    
    selected_indices = []
    
    for i, pos in enumerate(positions):
        # Distance from center axis (z-axis)
        xy_dist = np.linalg.norm(pos[:2] - xy_center)
        
        # Relative height (0 = bottom, 1 = top)
        z_rel = (pos[2] - z_min) / z_total
        
        # Calculate allowed radius at this height
        # Linear taper: radius decreases linearly from base to tip
        allowed_radius = max_radius * (1 - z_rel) ** (tip_sharpness / 2.0)
        
        # Keep atom if within allowed radius
        if xy_dist <= allowed_radius:
            selected_indices.append(i)
    
    return selected_indices


def sharpen_tip(atoms, tip_atom_count=1):
    """
    Further sharpen the tip by keeping only the specified number of topmost atoms
    
    Args:
        atoms: ASE Atoms object
        tip_atom_count: Number of atoms to keep at the very tip (1-7)
        
    Returns:
        Modified atoms object with sharpened tip
    """
    positions = atoms.get_positions()
    
    # Find the topmost atoms
    z_max = positions[:, 2].max()
    z_threshold = z_max - 0.5  # Atoms within 0.5 Å of the top
    
    top_indices = []
    other_indices = []
    
    for i, pos in enumerate(positions):
        if pos[2] >= z_threshold:
            top_indices.append(i)
        else:
            other_indices.append(i)
    
    # Sort top atoms by z-coordinate (highest first)
    top_indices.sort(key=lambda i: positions[i, 2], reverse=True)
    
    # Keep only the specified number of tip atoms
    final_indices = other_indices + top_indices[:tip_atom_count]
    
    # Create new atoms object with selected atoms
    tip_atoms = Atoms()
    for i in final_indices:
        tip_atoms.append(atoms[i])
    
    tip_atoms.set_pbc([False, False, False])
    tip_atoms.center()
    
    return tip_atoms


def add_base_fixation_tags(atoms, fixed_layers=2):
    """
    Add tags to fix the bottom layers of the tip
    
    Args:
        atoms: ASE Atoms object
        fixed_layers: Number of bottom layers to fix
        
    Returns:
        Modified atoms object with fixation tags
    """
    positions = atoms.get_positions()
    z_min = positions[:, 2].min()
    z_max = positions[:, 2].max()
    z_range = z_max - z_min
    
    # Calculate layer thickness (approximate)
    layer_thickness = z_range / (len(set(positions[:, 2])) - 1) if len(set(positions[:, 2])) > 1 else z_range / 5
    
    # Assign tags based on z-position
    tags = []
    for pos in positions:
        layer_num = int((pos[2] - z_min) / layer_thickness)
        if layer_num < fixed_layers:
            tags.append(1)  # Fixed atoms
        else:
            tags.append(0)  # Mobile atoms
    
    atoms.set_tags(tags)
    return atoms


def build_tungsten_stm_tip(height=8, base_size=5, tip_sharpness=3, 
                          tip_atoms=1, freeze_base=True, output_file='stm_tip.xyz'):
    """
    Build tungsten STM tip structure
    
    Args:
        height: Height of the tip in BCC unit cells (default: 8)
        base_size: Size of the base in BCC unit cells (default: 5) 
        tip_sharpness: Sharpening factor 1-5 (higher = sharper, default: 3)
        tip_atoms: Number of atoms at the very tip (1-7, default: 1)
        freeze_base: Add fixation tags for base atoms (default: True)
        output_file: Output filename (default: 'stm_tip.xyz')
        
    Returns:
        bool: Success status
    """
    try:
        # Validate input parameters
        height = max(3, int(height))
        base_size = max(2, int(base_size))
        tip_sharpness = max(1, min(5, int(tip_sharpness)))
        tip_atoms = max(1, min(7, int(tip_atoms)))
        
        print(f"Building tungsten STM tip:")
        print(f"  Height: {height} BCC cells")
        print(f"  Base size: {base_size} BCC cells") 
        print(f"  Tip sharpness: {tip_sharpness}/5")
        print(f"  Tip atoms: {tip_atoms}")
        print(f"  Freeze base: {freeze_base}")
        
        # Step 1: Create BCC tungsten supercell
        print("\n1. Creating BCC tungsten supercell...")
        w_supercell = create_bcc_tungsten_supercell((base_size, base_size, height))
        print(f"   Created supercell with {len(w_supercell)} atoms")
        
        # Step 2: Filter atoms to tip shape
        print("\n2. Filtering atoms to tip shape...")
        tip_indices = filter_atoms_to_tip_shape(w_supercell, height_ratio=0.9, 
                                               tip_sharpness=tip_sharpness)
        
        # Create tip structure from filtered atoms
        tip = Atoms()
        for i in tip_indices:
            tip.append(w_supercell[i])
        
        tip.set_pbc([False, False, False])
        tip.center()
        print(f"   Filtered to {len(tip)} atoms")
        
        # Step 3: Sharpen the tip
        print(f"\n3. Sharpening tip to {tip_atoms} atom(s)...")
        tip = sharpen_tip(tip, tip_atoms)
        print(f"   Final tip has {len(tip)} atoms")
        
        # Step 4: Add base fixation if requested
        if freeze_base:
            print("\n4. Adding base fixation tags...")
            tip = add_base_fixation_tags(tip, fixed_layers=2)
            fixed_count = sum(1 for tag in tip.get_tags() if tag == 1)
            print(f"   Fixed {fixed_count} base atoms")
        
        # Step 5: Add coordination analysis for visualization
        print("\n5. Analyzing coordination...")
        cutoffs = [1.8] * len(tip)  # W-W bond cutoff ~3.5 Å, use 1.8 Å as half-cutoff
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(tip)
        
        coord_counts = {}
        rgb_colors = []
        
        for i in range(len(tip)):
            indices, offsets = nl.get_neighbors(i)
            num_neighbors = len(indices)
            coord_counts[num_neighbors] = coord_counts.get(num_neighbors, 0) + 1
            
            # Color coding by coordination
            if num_neighbors >= 8:
                # Bulk coordination - Blue
                rgb_colors.append([0, 0, 255])
            elif num_neighbors >= 6:
                # Surface - Green
                rgb_colors.append([0, 255, 0])
            elif num_neighbors >= 4:
                # Edge - Orange
                rgb_colors.append([255, 165, 0])
            elif num_neighbors >= 2:
                # Corner - Yellow
                rgb_colors.append([255, 255, 0])
            else:
                # Single atom - Red
                rgb_colors.append([255, 0, 0])
        
        # Add RGB arrays for visualization
        rgb_array = np.array(rgb_colors)
        tip.set_array('R', rgb_array[:, 0])
        tip.set_array('G', rgb_array[:, 1])
        tip.set_array('B', rgb_array[:, 2])
        
        # Step 6: Save the structure
        write(output_file, tip, format='extxyz')
        print(f"\n✅ STM tip saved to '{output_file}' with {len(tip)} atoms")
        
        # Print coordination statistics
        print("\nCoordination statistics:")
        for coord in sorted(coord_counts.keys()):
            count = coord_counts[coord]
            print(f"  {coord} neighbors: {count} atoms")
        
        # Print tip information
        positions = tip.get_positions()
        tip_height = positions[:, 2].max() - positions[:, 2].min()
        base_width = 2 * np.max(np.linalg.norm(positions[:, :2] - positions[:, :2].mean(axis=0), axis=1))
        
        print(f"\nTip dimensions:")
        print(f"  Height: {tip_height:.2f} Å")
        print(f"  Base width: {base_width:.2f} Å")
        print(f"  Aspect ratio: {tip_height/base_width:.2f}")
        
        return True
        
    except Exception as e:
        print(f"Error building STM tip: {e}")
        return False


def main():
    """Example usage of the STM tip builder"""
    print("Building example tungsten STM tips...")
    
    # Sharp single-atom tip
    build_tungsten_stm_tip(
        height=6, 
        base_size=4, 
        tip_sharpness=5, 
        tip_atoms=1,
        output_file='stm_tip_sharp.xyz'
    )
    
    # Blunter multi-atom tip  
    build_tungsten_stm_tip(
        height=8,
        base_size=5, 
        tip_sharpness=2,
        tip_atoms=3,
        output_file='stm_tip_blunt.xyz'
    )


if __name__ == "__main__":
    main()