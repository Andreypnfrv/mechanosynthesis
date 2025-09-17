"""
Silicon (100)-(2×1) reconstructed surface builder for mechanosynthesis calculations.

Creates clean Si(100) surfaces with (2×1) reconstruction, H-passivation, 
and proper dangling bond sites for TIMe-Ge molecular tool attachment.
"""

import numpy as np
from ase import Atoms
from ase.build import bulk, surface, add_adsorbate
from ase.visualize import view
import os

class Si100Surface:
    """Builder for Si(100)-(2×1) reconstructed surfaces."""
    
    # Silicon lattice parameter (Å)
    LATTICE_PARAM = 5.431
    
    def __init__(self, nx=4, ny=4, layers=3, passivate=True):
        """
        Initialize Si(100) surface builder.
        
        Parameters:
        -----------
        nx : int
            Number of dimer rows along x-direction
        ny : int  
            Number of dimers per row along y-direction
        layers : int
            Number of Si layers (minimum 3 recommended)
        passivate : bool
            Add H atoms to passivate edge dangling bonds
        """
        self.nx = nx
        self.ny = ny 
        self.layers = layers
        self.passivate = passivate
        self.atoms = None
        self.active_sites = []
        
    def build_bulk_structure(self):
        """Create proper Si(100) surface using ASE surface builder."""
        # Create bulk Si with diamond structure
        si_bulk = bulk('Si', 'diamond', a=self.LATTICE_PARAM)
        
        # Create (100) surface - this gives proper crystallographic structure
        # The surface function creates the correct Miller index cut
        si_surface = surface(si_bulk, (1, 0, 0), layers=self.layers, vacuum=10.0)
        
        # Repeat the surface to get desired size
        si_surface = si_surface.repeat((self.nx, self.ny, 1))
        
        return si_surface
        
    def create_surface_cut(self, si_surface):
        """Surface is already properly cut by ASE, just return it."""
        return si_surface
        
    def apply_reconstruction(self, si_surface):
        """Apply (2×1) reconstruction to surface atoms."""
        atoms = si_surface.copy()
        positions = atoms.positions
        
        # Find surface layer (highest z)
        z_coords = positions[:, 2]
        z_max = np.max(z_coords)
        z_threshold = z_max - 0.5  # More generous threshold
        
        surface_indices = np.where(z_coords > z_threshold)[0]
        
        if len(surface_indices) == 0:
            return atoms
        
        # Create (2×1) dimers by pairing surface atoms
        dimer_bond = 2.25  # Å, Si-Si dimer bond length
        dimer_tilt = 0.1   # Å, asymmetric dimer tilt
        
        # Group surface atoms by rows (y-coordinate)
        surface_positions = positions[surface_indices]
        y_coords = surface_positions[:, 1]
        
        # Round y-coordinates to find rows
        y_unique = np.unique(np.round(y_coords, 1))
        
        new_positions = positions.copy()
        
        for y_coord in y_unique:
            # Find atoms in this row
            row_mask = np.abs(y_coords - y_coord) < 0.2
            row_indices = surface_indices[row_mask]
            
            if len(row_indices) < 2:
                continue
                
            # Sort by x-coordinate
            row_positions = positions[row_indices]
            x_sorted = np.argsort(row_positions[:, 0])
            sorted_indices = row_indices[x_sorted]
            
            # Create dimers from pairs
            for i in range(0, len(sorted_indices) - 1, 2):
                idx1, idx2 = sorted_indices[i], sorted_indices[i + 1]
                pos1, pos2 = positions[idx1], positions[idx2]
                
                # Calculate dimer center
                center = (pos1 + pos2) / 2
                
                # Create asymmetric dimer
                new_positions[idx1] = center + np.array([-dimer_bond/2, 0, dimer_tilt])
                new_positions[idx2] = center + np.array([dimer_bond/2, 0, -dimer_tilt])
        
        atoms.positions = new_positions
        return atoms
        
    def add_hydrogen_passivation(self, positions):
        """Add H atoms to passivate edge dangling bonds."""
        symbols = ['Si'] * len(positions)
        
        if not self.passivate:
            return positions, symbols
            
        h_positions = []
        Si_H_bond = 1.48  # Å, Si-H bond length
        
        # Find edge atoms and add proper H passivation
        x_min, x_max = np.min(positions[:, 0]), np.max(positions[:, 0])
        y_min, y_max = np.min(positions[:, 1]), np.max(positions[:, 1])
        z_min = np.min(positions[:, 2])
        
        for i, pos in enumerate(positions):
            x, y, z = pos
            
            # Passivate edge atoms
            if abs(x - x_min) < 0.5:  # Left edge
                h_pos = pos + np.array([-Si_H_bond, 0, 0])
                h_positions.append(h_pos)
                symbols.append('H')
                
            if abs(x - x_max) < 0.5:  # Right edge
                h_pos = pos + np.array([Si_H_bond, 0, 0])
                h_positions.append(h_pos)
                symbols.append('H')
                
            if abs(y - y_min) < 0.5:  # Front edge
                h_pos = pos + np.array([0, -Si_H_bond, 0])
                h_positions.append(h_pos)
                symbols.append('H')
                
            if abs(y - y_max) < 0.5:  # Back edge
                h_pos = pos + np.array([0, Si_H_bond, 0])
                h_positions.append(h_pos)
                symbols.append('H')
                
            # Passivate top two layers (surface atoms)
            z_coords = positions[:, 2]
            unique_z = np.unique(np.round(z_coords, 1))
            top_2_layers = sorted(unique_z)[-2:] if len(unique_z) >= 2 else [np.max(z_coords)]
            
            for layer_z in top_2_layers:
                if abs(z - layer_z) < 0.1:
                    h_pos = pos + np.array([0, 0, Si_H_bond])
                    h_positions.append(h_pos)
                    symbols.append('H')
                    break
        
        if h_positions:
            all_positions = np.vstack([positions, np.array(h_positions)])
        else:
            all_positions = positions
            
        return all_positions, symbols
        
    def identify_active_sites(self, positions, symbols):
        """Identify dangling bond sites for molecular tool attachment."""
        # Find surface Si atoms without full coordination
        z_max = np.max(positions[np.array(symbols) == 'Si'][:, 2])
        z_threshold = z_max - 0.1
        
        active_sites = []
        for i, (pos, sym) in enumerate(zip(positions, symbols)):
            if sym == 'Si' and pos[2] > z_threshold:
                # This is a surface Si atom - potential active site
                active_sites.append(i)
                
        self.active_sites = active_sites
        return active_sites
        
    def build(self):
        """Build complete Si(100)-(2×1) surface."""
        # Create bulk structure using ASE
        si_surface = self.build_bulk_structure()
        
        # Apply reconstruction
        si_surface = self.apply_reconstruction(si_surface)
        
        # Add passivation
        positions, symbols = self.add_hydrogen_passivation(si_surface.positions)
        
        # Create final ASE atoms object
        self.atoms = Atoms(symbols=symbols, positions=positions)
        
        # Center the structure
        self.atoms.center()
        
        # Identify active sites
        self.identify_active_sites(positions, symbols)
        
        return self.atoms
        
    def save(self, project_path, filename='si_surface'):
        """Save surface structure in multiple formats."""
        if self.atoms is None:
            self.build()
            
        inputs_dir = os.path.join(project_path, 'inputs')
        os.makedirs(inputs_dir, exist_ok=True)
        
        # Save in multiple formats
        xyz_file = os.path.join(inputs_dir, f'{filename}.xyz')
        gen_file = os.path.join(inputs_dir, f'{filename}.gen')
        
        # XYZ format
        self.atoms.write(xyz_file)
        
        # GEN format for DFTB+
        with open(gen_file, 'w') as f:
            f.write(f"{len(self.atoms)} C\n")
            unique_symbols = list(set(self.atoms.get_chemical_symbols()))
            f.write(" ".join(unique_symbols) + "\n")
            
            symbol_to_index = {sym: i+1 for i, sym in enumerate(unique_symbols)}
            
            for i, (atom, pos) in enumerate(zip(self.atoms, self.atoms.positions)):
                symbol_idx = symbol_to_index[atom.symbol]
                f.write(f"{i+1:4d} {symbol_idx:2d} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
        
        # Save metadata
        symbols, counts = np.unique(self.atoms.get_chemical_symbols(), return_counts=True)
        composition = {sym: int(count) for sym, count in zip(symbols, counts)}
        
        metadata = {
            'nx': int(self.nx),
            'ny': int(self.ny), 
            'layers': int(self.layers),
            'passivated': bool(self.passivate),
            'n_atoms': len(self.atoms),
            'active_sites': [int(x) for x in self.active_sites],
            'composition': composition
        }
        
        import json
        with open(os.path.join(inputs_dir, 'parameters.json'), 'w') as f:
            json.dump(metadata, f, indent=2)
            
        print(f"Si(100)-(2×1) surface saved:")
        print(f"  Size: {self.nx}×{self.ny}×{self.layers}")
        print(f"  Atoms: {len(self.atoms)} ({metadata['composition']})")
        print(f"  Active sites: {len(self.active_sites)}")
        print(f"  Files: {xyz_file}, {gen_file}")
        
        return xyz_file, gen_file


def build_si_surface(nx=4, ny=4, layers=3, passivate=True, project_name=None):
    """
    Convenience function to build Si(100) surface.
    
    Parameters:
    -----------
    nx, ny, layers : int
        Surface dimensions
    passivate : bool
        H-passivation of edges
    project_name : str
        Project directory name
        
    Returns:
    --------
    tuple : (xyz_file, gen_file, metadata)
    """
    if project_name is None:
        project_name = f"si_surface_{nx}x{ny}x{layers}"
        
    # Create project directory
    from kernel.project_structure import create_project_structure
    project_path = create_project_structure(project_name)
    
    # Build surface
    surface = Si100Surface(nx=nx, ny=ny, layers=layers, passivate=passivate)
    surface.build()
    
    # Save files
    xyz_file, gen_file = surface.save(project_path)
    
    metadata = {
        'project_path': project_path,
        'surface_builder': surface,
        'xyz_file': xyz_file,
        'gen_file': gen_file
    }
    
    return metadata