from ase.build import bulk


def build_diamond_supercell(size, output_file='diamond_supercell.xyz', passivate=None):
    """Build bulk diamond supercell with full periodic boundary conditions"""
    try:
        size_int = int(size)
        if size_int <= 0:
            raise ValueError("Size must be positive")
    except ValueError:
        print(f"Error: Invalid size '{size}'. Please provide a positive integer.")
        return False
    
    # Create bulk diamond structure as isolated particle (no PBC)
    diamond = bulk('C', 'diamond', a=3.567, cubic=True)
    diamond = diamond.repeat((size_int, size_int, size_int))
    diamond.set_pbc([False, False, False])  # Remove periodic boundary conditions
    print(f"Created diamond supercell with {len(diamond)} atoms ({size_int}x{size_int}x{size_int})")
    
    if passivate:
        from .build_diamond import find_atoms_with_hanging_bonds, create_passivated_structure
        hanging_atoms = find_atoms_with_hanging_bonds(diamond)
        if hanging_atoms:
            print(f"Found {len(hanging_atoms)} atoms with hanging bonds - passivating with {passivate}")
            diamond = create_passivated_structure(diamond, hanging_atoms, passivate)
        else:
            print("No hanging bonds found - structure is fully coordinated")
    
    # Save the supercell structure
    from ase.io import write
    write(output_file, diamond)
    print(f"Diamond supercell saved to '{output_file}'")
    
    return True