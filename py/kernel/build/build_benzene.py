from ase.build import molecule
from kernel.project_structure import ase_to_xyz_input

benzene = molecule("C6H6", vacuum=6.0)

ase_to_xyz_input(benzene, "2025-08-22_benzene", "benzene")








