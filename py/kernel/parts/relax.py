"""
Geometry relaxation module for mechanosynthesis research.

Available methods:
- relax(): DFTB+ with ASE/BFGS optimizer
- relax_native_dftb(): Native DFTB+ optimization (fast)
- relax_passivated(): Multi-stage DFTB+ for passivated structures
- relax_orca_dft(): Orca DFT with ASE optimizer and video generation
- relax_orca_native(): Native Orca DFT optimization with video generation
- relax_xtb(): xTB semiempirical quantum chemistry optimization
- relax_gxtb(): g-xTB semiempirical quantum chemistry optimization (Docker)
"""

import os
import glob
import numpy as np
import time
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent.parent))
import settings

from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.optimize import BFGS
from ase.calculators.dftb import Dftb
from ase.calculators.orca import ORCA
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.patches import Circle
import subprocess
import shutil
import platform
from ..utils.timer import calculate_and_format_time, print_completion_message, print_failure_message, print_timeout_message
from ase.constraints import FixAtoms


def run_structural_analysis_if_enabled(initial_structure_path, final_structure_path, 
                                      project_folder, method_name, enable_analysis=None):
    """
    Run structural analysis if enabled
    
    Args:
        initial_structure_path: Path to initial structure
        final_structure_path: Path to final relaxed structure
        project_folder: Project folder name
        method_name: Name of relaxation method used
        enable_analysis: Whether to run analysis (None = check environment/config)
    """
    # Check if analysis is enabled
    if enable_analysis is None:
        # Use settings.py configuration
        try:
            import sys
            from pathlib import Path
            sys.path.append(str(Path(__file__).parent.parent.parent))
            import settings
            enable_analysis = getattr(settings, 'STRUCTURAL_ANALYSIS', True)
        except:
            # Fallback to environment variable or default True
            enable_analysis = os.environ.get('STRUCTURAL_ANALYSIS', 'true').lower() in ['true', '1', 'yes', 'on']
    
    if not enable_analysis:
        return
    
    try:
        from .structural_analysis import analyze_structural_changes, print_structural_summary
        
        # Set up output directory for analysis within method-specific folder
        # Try to determine the method-specific output directory from the final structure path
        if os.path.exists(final_structure_path):
            method_output_dir = os.path.dirname(final_structure_path)
            analysis_dir = os.path.join(method_output_dir, "structural_analysis")
        else:
            # Fallback to general outputs directory
            script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            project_root = os.path.dirname(script_dir)
            base_path = os.path.join(project_root, "data", project_folder)
            analysis_dir = os.path.join(base_path, "outputs", "structural_analysis")
        
        os.makedirs(analysis_dir, exist_ok=True)
        
        print(f"\n{'='*50}")
        print(f"RUNNING STRUCTURAL ANALYSIS ({method_name.upper()})")
        print(f"{'='*50}")
        
        # Run analysis
        analysis = analyze_structural_changes(
            initial_structure_path, 
            final_structure_path, 
            analysis_dir
        )
        
        # Print summary to console
        print_structural_summary(analysis)
        
    except ImportError as e:
        print(f"Warning: Could not import structural analysis module: {e}")
    except Exception as e:
        print(f"Warning: Structural analysis failed: {e}")


def apply_surface_constraints(atoms, fix_bottom_layers=True, n_layers_to_fix=1):
    """
    Apply constraints to surface structures by fixing bottom layers.
    
    Args:
        atoms: ASE atoms object
        fix_bottom_layers: bool, whether to fix bottom layers
        n_layers_to_fix: int, number of bottom layers to fix
        
    Returns:
        atoms: ASE atoms object with constraints applied
    """
    if not fix_bottom_layers:
        return atoms
        
    # Get only Si atoms for layer analysis (ignore H atoms)
    symbols = atoms.get_chemical_symbols()
    si_indices = [i for i, s in enumerate(symbols) if s == 'Si']
    
    if not si_indices:
        print("No Si atoms found for constraint analysis")
        return atoms
    
    # Analyze Si atoms z-coordinates
    si_positions = atoms.positions[si_indices]
    z_coords = si_positions[:, 2]
    
    # Sort Si atoms by z-coordinate
    z_sorted_indices = np.argsort(z_coords)
    z_sorted = z_coords[z_sorted_indices]
    
    # Use clustering to identify distinct layers
    # Silicon layers typically have ~1.36 Å spacing
    layer_tolerance = 0.8  # Angstrom - more generous for surface structures
    
    # Find layer boundaries using clustering
    layers = []
    current_layer = [z_sorted[0]]
    
    for i in range(1, len(z_sorted)):
        if abs(z_sorted[i] - z_sorted[i-1]) < layer_tolerance:
            current_layer.append(z_sorted[i])
        else:
            layers.append(current_layer)
            current_layer = [z_sorted[i]]
    layers.append(current_layer)  # Don't forget the last layer
    
    print(f"Found {len(layers)} Si layers:")
    layer_centers = []
    for i, layer in enumerate(layers):
        center_z = np.mean(layer)
        layer_centers.append(center_z)
        print(f"  Layer {i}: z={center_z:.2f} Å ({len(layer)} atoms)")
    
    # Fix atoms in bottom n_layers_to_fix layers
    atoms_to_fix = []
    layers_to_fix = min(n_layers_to_fix, len(layers))
    
    for layer_idx in range(layers_to_fix):
        layer_center = layer_centers[layer_idx]
        for i, si_idx in enumerate(si_indices):
            z = atoms.positions[si_idx, 2]
            if abs(z - layer_center) <= layer_tolerance:
                atoms_to_fix.append(si_idx)
    
    if atoms_to_fix:
        constraint = FixAtoms(indices=atoms_to_fix)
        atoms.set_constraint(constraint)
        print(f"Fixed {len(atoms_to_fix)} Si atoms in bottom {layers_to_fix} layer(s)")
        
        # Also fix any H atoms attached to fixed Si atoms (within 2 Å)
        h_indices = [i for i, s in enumerate(symbols) if s == 'H']
        additional_h_fixes = []
        
        for h_idx in h_indices:
            h_pos = atoms.positions[h_idx]
            for si_idx in atoms_to_fix:
                si_pos = atoms.positions[si_idx]
                distance = np.linalg.norm(h_pos - si_pos)
                if distance < 2.0:  # H bonded to fixed Si
                    additional_h_fixes.append(h_idx)
                    break
        
        if additional_h_fixes:
            all_fixed = atoms_to_fix + additional_h_fixes
            constraint = FixAtoms(indices=all_fixed)
            atoms.set_constraint(constraint)
            print(f"Also fixed {len(additional_h_fixes)} H atoms bonded to fixed Si")
    
    # Additional edge stabilization for surface structures
    edge_fixes = []
    all_positions = atoms.positions
    
    # Find structure boundaries
    x_coords = all_positions[:, 0]
    y_coords = all_positions[:, 1]
    z_coords = all_positions[:, 2]
    
    x_min, x_max = x_coords.min(), x_coords.max()
    y_min, y_max = y_coords.min(), y_coords.max()
    z_min = z_coords.min()
    
    # Fix edge atoms in lower half of structure for additional stability
    z_middle = z_min + (z_coords.max() - z_min) * 0.5
    
    for si_idx in si_indices:
        if si_idx in atoms_to_fix:  # Already fixed
            continue
            
        pos = all_positions[si_idx]
        x, y, z = pos
        
        # Check if atom is at edge and in lower half
        is_edge = (abs(x - x_min) < 3.0 or abs(x - x_max) < 3.0 or 
                   abs(y - y_min) < 3.0 or abs(y - y_max) < 3.0)
        is_lower = z < z_middle
        
        if is_edge and is_lower:
            edge_fixes.append(si_idx)
    
    if edge_fixes:
        all_fixed = atoms_to_fix + additional_h_fixes + edge_fixes
        constraint = FixAtoms(indices=all_fixed)
        atoms.set_constraint(constraint)
        print(f"Additionally fixed {len(edge_fixes)} edge Si atoms for stability")
    
    return atoms


def detect_unpaired_electrons(atoms):
    """
    Detect number of unpaired electrons in a molecule based on even/odd electron count.
    
    Args:
        atoms: ASE atoms object
        
    Returns:
        int: Number of unpaired electrons (0 for closed-shell, 1+ for radicals)
    """
    # Count total number of electrons
    total_electrons = sum(atom.number for atom in atoms)  # atomic numbers sum to electron count for neutral molecule
    
    # If odd number of electrons, we have 1 unpaired electron (doublet radical)
    # If even number, assume closed shell (singlet) unless specified otherwise
    unpaired = total_electrons % 2
    
    return unpaired


def estimate_time_range(atoms, calc_type, method='B3LYP', nprocs=1):
    """
    Estimate calculation time range based on real performance data
    
    Args:
        atoms: ASE atoms object
        calc_type: 'dftb_relax', 'dftb_native', 'dftb_passivated', 'orca_relax', 'orca_hessian'
        method: DFT method (affects scaling for Orca)
        nprocs: Number of processors (for parallel speedup estimation)
        
    Returns:
        str: Formatted time range estimate with step prediction
    """
    n = len(atoms)
    
    # Data-driven prediction formulas based on performance_logs.csv analysis
    def predict_dftb_time(atoms_count):
        # Based on real data: time ≈ 0.15 + 0.08 * atoms^1.75
        # Accounts for variation in DFTB calculations
        base_time = 0.15 + 0.08 * (atoms_count ** 1.75)
        return (base_time * 0.6, base_time * 2.5)  # min-max range from data variance
    
    def predict_xtb_time(atoms_count):
        # xTB is consistently very fast from our data (0.01s for all tested molecules)
        return (0.01, max(0.05, 0.01 * atoms_count / 10))
    
    def predict_orca_time(atoms_count, method_name='B3LYP'):
        # Based on w_hexamethyl (25 atoms) = 1651s for Hessian (much slower than relaxation)
        # For relaxation, estimate ~1/10 of Hessian time  
        method_factor = {'HF': 0.5, 'PBE': 0.8, 'B3LYP': 1.0, 'M06-2X': 1.5}.get(method_name.upper(), 1.0)
        base_time = 2.0 * (atoms_count ** 2.2) * method_factor  # More conservative scaling
        return (base_time * 0.3, base_time * 3.0)
    
    # Get time predictions based on method
    if calc_type.startswith('dftb'):
        min_time, max_time = predict_dftb_time(n)
    elif calc_type.startswith('xtb'):
        min_time, max_time = predict_xtb_time(n)
    elif calc_type.startswith('orca'):
        min_time, max_time = predict_orca_time(n, method)
        if 'hessian' in calc_type:
            min_time *= 5  # Hessian calculations are much more expensive
            max_time *= 10
    elif calc_type.startswith('gxtb'):
        # g-xTB similar to xTB but slower due to Docker overhead
        xtb_min, xtb_max = predict_xtb_time(n)
        min_time, max_time = xtb_min * 5, xtb_max * 25
        if 'hessian' in calc_type:
            min_time *= 120  # Numerical hessian is extremely slow
            max_time *= 300
    else:
        # Fallback for unknown methods
        min_time, max_time = 1.0, 10.0
    
    # Apply parallel speedup for Orca calculations
    parallel_speedup = 1.0
    if nprocs > 1 and calc_type.startswith('orca'):
        effective_cores = min(nprocs, 4)  # Diminishing returns beyond 4 cores
        parallel_speedup = 1 + (effective_cores - 1) * 0.7  # 70% efficiency per additional core
        min_time /= parallel_speedup
        max_time /= parallel_speedup
    
    # Predict optimization steps based on system complexity
    def predict_steps(atoms_count, calc_method):
        """Predict min-max optimization steps based on molecule size and method"""
        if calc_method.startswith('xtb'):
            # xTB converges quickly but sometimes with issues (quality=2)
            return (10, 50)
        elif calc_method.startswith('dftb'):
            # DFTB+ typically needs more steps for good convergence (quality=4)
            if atoms_count <= 5:
                return (20, 80)
            elif atoms_count <= 30:
                return (50, 150)
            else:
                return (100, 300)  # Large systems need more steps
        elif calc_method.startswith('orca'):
            # DFT methods are more stable but slower per step
            return (30, 120)
        else:
            return (20, 100)
    
    min_steps, max_steps = predict_steps(n, calc_type)
    
    def sec_to_str(sec):
        if sec < 60: 
            return f"{sec:.1f}с"
        elif sec < 3600: 
            return f"{int(sec//60)}м {int(sec%60)}с"
        else: 
            hours = int(sec // 3600)
            minutes = int((sec % 3600) // 60)
            return f"{hours}ч {minutes}м"
    
    speedup_note = f" (на {nprocs} ядрах)" if nprocs > 1 and calc_type.startswith('orca') else ""
    steps_info = f" (~{min_steps}-{max_steps} шагов)"
    
    return f"{sec_to_str(min_time)} - {sec_to_str(max_time)}{speedup_note}{steps_info}"


def get_angular_momentum_params(atoms):
    """Generate SlaterKoster angular momentum parameters from atoms"""
    unique_elements = set(atoms.get_chemical_symbols())
    params = {}
    
    # Default angular momentum for common elements
    angular_momentum = {
        'H': 's',
        'C': 'p',
        'N': 'p', 
        'O': 'p',
        'Si': 'd',
        'S': 'd',
        'P': 'd'
    }
    
    for element in unique_elements:
        if element in angular_momentum:
            params[f'Hamiltonian_MaxAngularMomentum_{element}'] = angular_momentum[element]
        else:
            # Default to 'p' for unknown elements
            params[f'Hamiltonian_MaxAngularMomentum_{element}'] = 'p'
    
    return params


def relax_native_dftb(project_folder):
    """
    Relax using native DFTB+ geometry optimization (much faster)
    
    Args:
        project_folder (str): Name of the project folder or molecule name
    """
    base_path = f"data/{project_folder}"
    
    # If exact folder doesn't exist, try to find a folder containing the molecule name
    if not os.path.exists(base_path):
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
        else:
            raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    if not os.path.exists(base_path):
        raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    # Find all XYZ files in base directory and inputs subdirectory
    xyz_files = []
    xyz_pattern = os.path.join(base_path, "*.xyz")
    xyz_files.extend(glob.glob(xyz_pattern))
    
    inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
    xyz_files.extend(glob.glob(inputs_pattern))
    
    # Convert to absolute paths to handle directory changes
    xyz_files = [os.path.abspath(f) for f in xyz_files]
    
    if not xyz_files:
        print(f"No XYZ files found in {base_path}")
        return False
    
    # Create logs, outputs, and temp directories
    logs_dir = os.path.join(base_path, "logs")
    outputs_dir = os.path.join(base_path, "outputs")
    native_outputs_dir = os.path.join(base_path, "outputs", "native_dftb")
    temp_dir = os.path.join(base_path, "temp")
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    os.makedirs(native_outputs_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Store original working directory
    original_cwd = os.getcwd()
    
    try:
        # Change to temp directory for DFTB+ execution
        os.chdir(temp_dir)
        
        for xyz_file in xyz_files:
            filename = os.path.basename(xyz_file)
            name_without_ext = os.path.splitext(filename)[0]
            
            print(f"Processing {filename} with native DFTB+ optimization...")
            
            # Read atoms and convert to DFTB+ gen format
            atoms = read(xyz_file)
            
            # Display time estimate
            time_estimate = estimate_time_range(atoms, 'dftb_native')
            print(f"⏱️  Оценка времени оптимизации: {time_estimate}")
            print(f"    Формула: {atoms.get_chemical_formula()} ({len(atoms)} атомов)")
            gen_file = f"{name_without_ext}.gen"
            write(gen_file, atoms, format='gen')
            
            # Generate angular momentum parameters
            angular_params = get_angular_momentum_params(atoms)
            
            # Create DFTB+ input file with geometry optimization
            hsd_content = f"""
Geometry = GenFormat {{
    <<< "{gen_file}"
}}

Driver = GeometryOptimization {{
    Optimizer = Rational {{}}
    MaxSteps = 500
    OutputPrefix = "{name_without_ext}_opt"
}}

Hamiltonian = DFTB {{
    SCC = Yes
    SCCTolerance = 1e-6
    MaxSCCIterations = 500
    Filling = Fermi {{Temperature[K] = 300.0}}
    HCorrection = Damping {{Exponent = 4.0}}
    Mixer = Anderson {{MixingParameter = 0.05}}
    SlaterKosterFiles = Type2FileNames {{
        Prefix = "{os.environ.get('SKDIR', '/Users/andreypanferov/opt/dftb+/slakos')}/mio-1-1/"
        Separator = "-"
        Suffix = ".skf"
    }}
    MaxAngularMomentum = {{
"""
            
            # Add angular momentum parameters
            for element in set(atoms.get_chemical_symbols()):
                angular_momentum = {'H': 's', 'C': 'p', 'N': 'p', 'O': 'p', 'Si': 'd', 'S': 'd', 'P': 'd'}
                am = angular_momentum.get(element, 'p')
                hsd_content += f"        {element} = {am}\n"
            
            hsd_content += """    }
}

Analysis {
    PrintForces = Yes
}

Options {
    WriteResultsTag = Yes
}
"""
            
            # Write input file
            with open("dftb_in.hsd", "w") as f:
                f.write(hsd_content)
            
            print(f"Starting native DFTB+ optimization for {filename} with {len(atoms)} atoms")
            print(f"Formula: {atoms.get_chemical_formula()}")
            
            # Run DFTB+ with real-time output
            import subprocess
            print("Running DFTB+ optimization...")
            process = subprocess.Popen(["dftb+"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, 
                                     universal_newlines=True, bufsize=1)
            
            with open("dftb.out", "w") as logfile:
                step_count = 0
                for line in process.stdout:
                    logfile.write(line)
                    if "Geometry step:" in line:
                        step_count += 1
                        print(f"  Step {step_count}: {line.strip()}")
                    elif "Total energy:" in line:
                        print(f"  {line.strip()}")
                    elif "Geometry converged" in line:
                        print(f"  {line.strip()}")
            
            result = process.wait()
            
            if result != 0:
                print(f"DFTB+ failed for {filename}")
                continue
            
            # Copy results to outputs
            xyz_traj = f"{name_without_ext}_opt.xyz"
            gen_final = f"{name_without_ext}_opt.gen"
            
            if os.path.exists(xyz_traj):
                final_output = os.path.join("..", "outputs", "native_dftb", f"{name_without_ext}_relaxed.xyz")
                extxyz_output = os.path.join("..", "outputs", "native_dftb", f"{name_without_ext}_relaxed.extxyz")
                os.system(f"cp {xyz_traj} {final_output}")
                os.system(f"cp {xyz_traj} {extxyz_output}")
                print(f"Native DFTB+ optimization completed for {filename} in {step_count} steps")
                print(f"Results saved to outputs/native_dftb/")
            else:
                print(f"Warning: No trajectory output found for {filename}")
            
            # Copy log to logs directory
            log_output = os.path.join("..", "logs", f"{name_without_ext}_dftb_native.log")
            os.system(f"cp dftb.out {log_output}")
            
    finally:
        # Always return to original working directory
        os.chdir(original_cwd)
    
    # All files processed successfully
    return True


def relax_passivated(project_folder, stages=3):
    """
    Multi-stage relaxation for passivated structures (with H or other surface atoms)
    
    Stage 1: Very gentle optimization (maxstep=0.02, fmax=0.1)
    Stage 2: Medium optimization (maxstep=0.05, fmax=0.05) 
    Stage 3: Final optimization (maxstep=0.1, fmax=0.02)
    
    Args:
        project_folder (str): Name of the project folder or molecule name
        stages (int): Number of relaxation stages (default: 3)
    """
    base_path = f"data/{project_folder}"
    
    # If exact folder doesn't exist, try to find a folder containing the molecule name
    if not os.path.exists(base_path):
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
        else:
            raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    if not os.path.exists(base_path):
        raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    # Find all XYZ files in base directory and inputs subdirectory
    xyz_files = []
    xyz_pattern = os.path.join(base_path, "*.xyz")
    xyz_files.extend(glob.glob(xyz_pattern))
    
    inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
    xyz_files.extend(glob.glob(inputs_pattern))
    
    # Convert to absolute paths to handle directory changes
    xyz_files = [os.path.abspath(f) for f in xyz_files]
    
    if not xyz_files:
        print(f"No XYZ files found in {base_path}")
        return False
    
    # Create logs, outputs, and temp directories
    logs_dir = os.path.join(base_path, "logs")
    outputs_dir = os.path.join(base_path, "outputs")
    staged_outputs_dir = os.path.join(base_path, "outputs", "staged_ase")
    temp_dir = os.path.join(base_path, "temp")
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    os.makedirs(staged_outputs_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Store original working directory
    original_cwd = os.getcwd()
    
    # Define relaxation stages
    stage_params = [
        {'maxstep': 0.02, 'fmax': 0.1, 'steps': 100, 'name': 'gentle'},
        {'maxstep': 0.05, 'fmax': 0.05, 'steps': 200, 'name': 'medium'},
        {'maxstep': 0.1, 'fmax': 0.02, 'steps': 300, 'name': 'final'}
    ][:stages]
    
    try:
        # Change to project directory to contain all DFTB+ files
        os.chdir(base_path)
    
        for xyz_file in xyz_files:
            filename = os.path.basename(xyz_file)
            name_without_ext = os.path.splitext(filename)[0]
            
            print(f"Processing {filename} with {stages}-stage relaxation...")
            
            # Read atoms using absolute path since we changed directory
            atoms = read(xyz_file)
            
            # Display time estimate
            time_estimate = estimate_time_range(atoms, 'dftb_passivated')
            print(f"⏱️  Оценка времени многоэтапной релаксации: {time_estimate}")
            print(f"    Формула: {atoms.get_chemical_formula()} ({len(atoms)} атомов)")
            
            # Check if structure is passivated (contains H atoms)
            symbols = atoms.get_chemical_symbols()
            has_hydrogen = 'H' in symbols
            if not has_hydrogen:
                print(f"Warning: {filename} does not contain hydrogen - consider using regular relax() instead")
            
            # Generate angular momentum parameters
            angular_params = get_angular_momentum_params(atoms)
            
            # Set up DFTB calculator with label pointing to temp/ folder
            calc_params = {
                'atoms': atoms,
                'label': f'temp/{name_without_ext}',
                'kpts': (1,1,1),
                'Hamiltonian_SCC': 'Yes',
                'Hamiltonian_SCCTolerance': 1e-6,
                'Hamiltonian_MaxSCCIterations': 500,
                'Hamiltonian_SlaterKosterFiles_Prefix': os.environ.get('SKDIR', '/Users/andreypanferov/opt/dftb+/slakos') + '/mio-1-1/',
                'Hamiltonian_SlaterKosterFiles_Separator': '-',
                'Hamiltonian_SlaterKosterFiles_Suffix': '.skf',
                'Hamiltonian_MaxAngularMomentum_': '',
                'Hamiltonian_HCorrection': 'Damping {Exponent = 4.0}',
                'Hamiltonian_Filling': 'Fermi {Temperature[K] = 300.0}',
                'Hamiltonian_Mixer': 'Anderson {MixingParameter = 0.05}'
            }
            
            # Add angular momentum parameters
            calc_params.update(angular_params)
            
            calc = Dftb(**calc_params)
            atoms.calc = calc
            
            # Set up trajectory and logs
            traj_file = os.path.join("outputs", "staged_ase", f"{name_without_ext}_relaxed.traj")
            overall_log = os.path.join("logs", f"{name_without_ext}_staged_relax.log")
            
            traj = Trajectory(traj_file, 'w', atoms)
            
            print(f"Starting staged relaxation of {filename} with {len(atoms)} atoms")
            print(f"Formula: {atoms.get_chemical_formula()}")
            
            # Run staged optimization
            for stage_idx, stage in enumerate(stage_params, 1):
                print(f"Stage {stage_idx}/{len(stage_params)}: {stage['name']} relaxation (maxstep={stage['maxstep']}, fmax={stage['fmax']})")
                
                stage_log = os.path.join("logs", f"{name_without_ext}_stage{stage_idx}_{stage['name']}.log")
                
                opt = BFGS(atoms, logfile=stage_log, maxstep=stage['maxstep'])
                opt.attach(traj.write, interval=1)
                
                try:
                    opt.run(fmax=stage['fmax'], steps=stage['steps'])
                    
                    # Get intermediate results
                    energy = atoms.get_potential_energy()
                    max_force = np.max(np.linalg.norm(atoms.get_forces(), axis=1))
                    print(f"  Stage {stage_idx} completed: E={energy:.6f} eV, max_force={max_force:.4f} eV/Å")
                    
                except Exception as e:
                    print(f"  Stage {stage_idx} failed: {e}")
                    break
            
            # Write final relaxed structure
            relaxed_xyz = os.path.join("outputs", "staged_ase", f"{name_without_ext}_relaxed.xyz")
            write(relaxed_xyz, atoms)
            
            # Export trajectory in multiple formats
            try:
                imgs = Trajectory(traj_file)
                trajectory_list = [atoms for atoms in imgs]
                
                # OVITO-compatible XYZ
                ovito_xyz = os.path.join("outputs", "staged_ase", f"{name_without_ext}_trajectory.xyz")
                create_ovito_trajectory(trajectory_list, ovito_xyz, name_without_ext)
                
                # EXTXYZ format
                extxyz_file = os.path.join("outputs", "staged_ase", f"{name_without_ext}_trajectory.extxyz")
                write(extxyz_file, trajectory_list, format='extxyz')
                
                print(f"Results saved to outputs/staged_ase/")
                print(f"Trajectory: {name_without_ext}_trajectory.xyz (OVITO compatible)")
            except Exception as e:
                print(f"Could not create trajectory files for {filename}: {e}")
            
            # Print final results
            try:
                final_energy = atoms.get_potential_energy()
                max_force = np.max(np.linalg.norm(atoms.get_forces(), axis=1))
                print(f"Final energy: {final_energy:.6f} eV")
                print(f"Max force: {max_force:.4f} eV/Å")
            except Exception as e:
                print(f"Could not get final energy/forces for {filename}: {e}")
            
            print(f"Staged relaxation completed for {filename}")
            print()
    
    finally:
        # Always return to original working directory
        os.chdir(original_cwd)


def relax(project_folder):
    """
    Relax all XYZ files in data/PROJECT_FOLDER/
    
    Args:
        project_folder (str): Name of the project folder or molecule name
    """
    base_path = f"data/{project_folder}"
    
    # If exact folder doesn't exist, try to find a folder containing the molecule name
    if not os.path.exists(base_path):
        # Look for folders that end with the molecule name
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
        else:
            raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    if not os.path.exists(base_path):
        raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    # Find all XYZ files in base directory and inputs subdirectory
    xyz_files = []
    xyz_pattern = os.path.join(base_path, "*.xyz")
    xyz_files.extend(glob.glob(xyz_pattern))
    
    inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
    xyz_files.extend(glob.glob(inputs_pattern))
    
    # Convert to absolute paths to handle directory changes
    xyz_files = [os.path.abspath(f) for f in xyz_files]
    
    if not xyz_files:
        print(f"No XYZ files found in {base_path} but there's inputs/benzene.xyz. always check inputs")
        return
    
    # Create logs, outputs, and temp directories
    logs_dir = os.path.join(base_path, "logs")
    outputs_dir = os.path.join(base_path, "outputs")
    ase_outputs_dir = os.path.join(base_path, "outputs", "ase_bfgs")
    temp_dir = os.path.join(base_path, "temp")
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    os.makedirs(ase_outputs_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Store original working directory
    original_cwd = os.getcwd()
    
    try:
        # Change to project directory to contain all DFTB+ files
        os.chdir(base_path)
    
        for xyz_file in xyz_files:
            filename = os.path.basename(xyz_file)
            name_without_ext = os.path.splitext(filename)[0]
            
            print(f"Processing {filename}...")
            
            # Read atoms using absolute path since we changed directory
            atoms = read(xyz_file)
            
            # Apply surface constraints for surface structures
            if 'surface' in project_folder.lower():
                atoms = apply_surface_constraints(atoms, fix_bottom_layers=True, n_layers_to_fix=2)
            
            # Display time estimate
            time_estimate = estimate_time_range(atoms, 'dftb_relax')
            print(f"⏱️  Оценка времени релаксации: {time_estimate}")
            print(f"    Формула: {atoms.get_chemical_formula()} ({len(atoms)} атомов)")
            
            # Generate angular momentum parameters
            angular_params = get_angular_momentum_params(atoms)
            
            # Set up DFTB calculator with label pointing to temp/ folder
            calc_params = {
                'atoms': atoms,
                'label': f'temp/{name_without_ext}',
                'kpts': (1,1,1),
                'Hamiltonian_SCC': 'Yes',
                'Hamiltonian_SCCTolerance': 1e-6,
                'Hamiltonian_MaxSCCIterations': 500,
                'Hamiltonian_SlaterKosterFiles_Prefix': os.environ.get('SKDIR', '/Users/andreypanferov/opt/dftb+/slakos') + '/mio-1-1/',
                'Hamiltonian_SlaterKosterFiles_Separator': '-',
                'Hamiltonian_SlaterKosterFiles_Suffix': '.skf',
                'Hamiltonian_MaxAngularMomentum_': '',
                'Hamiltonian_HCorrection': 'Damping {Exponent = 4.0}',
                'Hamiltonian_Filling': 'Fermi {Temperature[K] = 300.0}',
                'Hamiltonian_Mixer': 'Anderson {MixingParameter = 0.05}'
            }
            
            # Add angular momentum parameters
            calc_params.update(angular_params)
            
            calc = Dftb(**calc_params)
            atoms.calc = calc
            
            # Set up trajectory and optimizer
            traj_file = os.path.join("outputs", "ase_bfgs", f"{name_without_ext}_relaxed.traj")
            log_file = os.path.join("logs", f"{name_without_ext}_relax.log")
        
            traj = Trajectory(traj_file, 'w', atoms)
            opt = BFGS(atoms, logfile=log_file, maxstep=0.1)
            opt.attach(traj.write, interval=1)
            
            print(f"Starting relaxation of {filename} with {len(atoms)} atoms")
            print(f"Formula: {atoms.get_chemical_formula()}")
            
            # Run optimization
            opt.run(fmax=0.02, steps=300)
            
            # Write relaxed structure
            relaxed_xyz = os.path.join("outputs", "ase_bfgs", f"{name_without_ext}_relaxed.xyz")
            write(relaxed_xyz, atoms)
            
            # Export trajectory in multiple formats
            try:
                imgs = Trajectory(traj_file)
                trajectory_list = [atoms for atoms in imgs]
                
                # OVITO-compatible XYZ
                ovito_xyz = os.path.join("outputs", "ase_bfgs", f"{name_without_ext}_trajectory.xyz")
                create_ovito_trajectory(trajectory_list, ovito_xyz, name_without_ext)
                
                # EXTXYZ format
                extxyz_file = os.path.join("outputs", "ase_bfgs", f"{name_without_ext}_trajectory.extxyz")
                write(extxyz_file, trajectory_list, format='extxyz')
                
                print(f"Results saved to outputs/ase_bfgs/")
                print(f"Trajectory: {name_without_ext}_trajectory.xyz (OVITO compatible)")
            except Exception as e:
                print(f"Could not create trajectory files for {filename}: {e}")
            
            # Print final results
            try:
                final_energy = atoms.get_potential_energy()
                max_force = np.max(np.linalg.norm(atoms.get_forces(), axis=1))
                print(f"Final energy: {final_energy:.6f} eV")
                print(f"Max force: {max_force:.4f} eV/Å")
            except Exception as e:
                print(f"Could not get final energy/forces for {filename}: {e}")
            
            print(f"Relaxation completed for {filename}")
            print()
    
    finally:
        # Always return to original working directory
        os.chdir(original_cwd)
    
    # All files processed successfully
    return True


def create_ovito_trajectory(trajectory, output_file, molecule_name):
    """
    Create OVITO-compatible multi-frame XYZ file with energy information
    
    Args:
        trajectory (list): List of ASE atoms objects
        output_file (str): Path to output XYZ file  
        molecule_name (str): Name of the molecule
    """
    with open(output_file, 'w') as f:
        for i, atoms in enumerate(trajectory):
            f.write(f"{len(atoms)}\n")
            
            # Try to get energy from atoms object
            try:
                energy = atoms.get_potential_energy()
                f.write(f"Frame {i+1}, Energy = {energy:.6f} eV, {molecule_name} optimization\n")
            except:
                f.write(f"Frame {i+1}, {molecule_name} optimization step\n")
            
            # Write coordinates
            positions = atoms.get_positions()
            symbols = atoms.get_chemical_symbols()
            
            for symbol, pos in zip(symbols, positions):
                f.write(f"{symbol:2s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
    
    print(f"OVITO-compatible trajectory saved to {output_file}")


def create_pdb_trajectory(trajectory, output_file, molecule_name):
    """
    Create PDB trajectory file for visualization programs like VMD, PyMOL
    
    Args:
        trajectory (list): List of ASE atoms objects
        output_file (str): Path to output PDB file  
        molecule_name (str): Name of the molecule
    """
    with open(output_file, 'w') as f:
        for i, atoms in enumerate(trajectory):
            f.write(f"MODEL {i+1:8d}\n")
            f.write(f"REMARK Frame {i+1}, {molecule_name} optimization\n")
            
            # Try to get energy from atoms object
            try:
                energy = atoms.get_potential_energy()
                f.write(f"REMARK Energy = {energy:.6f} eV\n")
            except:
                f.write(f"REMARK Optimization step {i+1}\n")
            
            positions = atoms.get_positions()
            symbols = atoms.get_chemical_symbols()
            
            for j, (symbol, pos) in enumerate(zip(symbols, positions), 1):
                f.write(f"ATOM  {j:5d}  {symbol:<2s}   MOL     1    {pos[0]:8.3f}{pos[1]:8.3f}{pos[2]:8.3f}  1.00  0.00           {symbol:>2s}\n")
            
            f.write("ENDMDL\n")
    
    print(f"PDB trajectory saved to {output_file}")


def generate_optimization_video(trajectory, output_file, molecule_name):
    """
    Generate a video showing the geometry optimization trajectory
    
    Args:
        trajectory (list): List of ASE atoms objects from optimization
        output_file (str): Path to output MP4 file
        molecule_name (str): Name of the molecule for title
    """
    if len(trajectory) < 2:
        print("Not enough frames for video generation")
        return
    
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Color map for different elements
    element_colors = {
        'H': 'white', 'C': 'gray', 'N': 'blue', 'O': 'red',
        'S': 'yellow', 'P': 'orange', 'Si': 'tan'
    }
    
    def animate(frame_num):
        ax.clear()
        atoms = trajectory[frame_num]
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        # Plot atoms
        for i, (pos, symbol) in enumerate(zip(positions, symbols)):
            color = element_colors.get(symbol, 'purple')
            size = 100 if symbol == 'H' else 200
            ax.scatter(pos[0], pos[1], pos[2], c=color, s=size, alpha=0.8)
        
        # Set equal aspect ratio and labels
        ax.set_xlabel('X (Å)')
        ax.set_ylabel('Y (Å)')
        ax.set_zlabel('Z (Å)')
        ax.set_title(f'{molecule_name} Optimization - Step {frame_num+1}/{len(trajectory)}')
        
        # Set consistent limits
        all_pos = np.concatenate([frame.get_positions() for frame in trajectory])
        margin = 2.0
        ax.set_xlim(all_pos[:, 0].min() - margin, all_pos[:, 0].max() + margin)
        ax.set_ylim(all_pos[:, 1].min() - margin, all_pos[:, 1].max() + margin)
        ax.set_zlim(all_pos[:, 2].min() - margin, all_pos[:, 2].max() + margin)
    
    # Create animation
    anim = animation.FuncAnimation(fig, animate, frames=len(trajectory), 
                                 interval=500, blit=False, repeat=True)
    
    # Save video
    try:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=2, metadata=dict(artist='Mechanosynthesis'), bitrate=1800)
        anim.save(output_file, writer=writer)
        print(f"Optimization video saved to {output_file}")
    except Exception as e:
        print(f"Could not save video: {e}")
        # Try saving as GIF instead
        gif_file = output_file.replace('.mp4', '.gif')
        try:
            anim.save(gif_file, writer='pillow', fps=2)
            print(f"Optimization GIF saved to {gif_file}")
        except Exception as e2:
            print(f"Could not save GIF either: {e2}")
    
    plt.close(fig)


# def relax_orca_native_DISABLED(project_folder, method='B3LYP', basis='6-31G(d)', nprocs=4, calc_hessian=False, hessian_method='phonopy'):
#     """
#     Relax using native Orca DFT optimization with video generation
#     
#     Args:
#         project_folder (str): Name of the project folder or molecule name
#         method (str): DFT method (default: B3LYP)
#         basis (str): Basis set (default: 6-31G(d) - faster than def2-SVP with similar geometry accuracy)
#         nprocs (int): Number of processors (default: 4)
#         calc_hessian (bool): Calculate Hessian after relaxation completion (default: False)
#         hessian_method (str): Hessian calculation method - 'phonopy' or 'orca' (default: 'phonopy')
#     """
#     base_path = f"data/{project_folder}"
#     
#     # Find project folder
#     if not os.path.exists(base_path):
#         pattern = f"data/*_{project_folder}"
#         matching_folders = glob.glob(pattern)
#         
#         if len(matching_folders) == 1:
#             base_path = matching_folders[0]
#         elif len(matching_folders) > 1:
#             raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
#         else:
#             raise FileNotFoundError(f"Project folder {base_path} does not exist")
#     
#     # Find all XYZ files
#     xyz_files = []
#     xyz_pattern = os.path.join(base_path, "*.xyz")
#     xyz_files.extend(glob.glob(xyz_pattern))
#     
#     inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
#     xyz_files.extend(glob.glob(inputs_pattern))
#     
#     xyz_files = [os.path.abspath(f) for f in xyz_files]
#     
#     if not xyz_files:
#         print(f"No XYZ files found in {base_path}")
#         return False
#     
#     # Create directories
#     logs_dir = os.path.join(base_path, "logs")
#     outputs_dir = os.path.join(base_path, "outputs")
#     orca_native_dir = os.path.join(base_path, "outputs", "orca_native")
#     videos_dir = os.path.join(base_path, "outputs", "videos")
#     temp_dir = os.path.join(base_path, "temp")
#     
#     for dir_path in [logs_dir, outputs_dir, orca_native_dir, videos_dir, temp_dir]:
#         os.makedirs(dir_path, exist_ok=True)
#     
#     original_cwd = os.getcwd()
#     
#     try:
#         os.chdir(temp_dir)
#         
#         for xyz_file in xyz_files:
#             filename = os.path.basename(xyz_file)
#             name_without_ext = os.path.splitext(filename)[0]
#             
#             print(f"Processing {filename} with native Orca DFT ({method}/{basis})...")
#             
#             atoms = read(xyz_file)
#             
#             # Determine number of processors
#             if platform.system() == 'Darwin':
#                 actual_nprocs = min(nprocs, 4)  # Cap at 4 cores for reasonable performance
#             else:
#                 actual_nprocs = nprocs
#             
#             # Display time estimate
#             time_estimate = estimate_time_range(atoms, 'orca_relax', method, actual_nprocs)
#             print(f"⏱️  Оценка времени релаксации (~10-50 шагов): {time_estimate}")
#             print(f"    Формула: {atoms.get_chemical_formula()} ({len(atoms)} атомов)")
#             
#             # Ask for confirmation
#             response = input("Продолжить расчет? (y/n): ").lower().strip()
#             if response not in ['y', 'yes', 'да', '']:
#                 print(f"Пропускаем расчет для {filename}")
#                 continue
#             # Calculate number of electrons and determine multiplicity
#             n_electrons = sum(atoms.get_atomic_numbers()) - 0  # assuming neutral molecule
#             multiplicity = 2 if n_electrons % 2 == 1 else 1  # odd electrons = doublet, even = singlet
#             print("Actual processors actual_nprocs", actual_nprocs)
#             
#             inp_content = f"""! {method} {basis} OPT
# 
# %pal nprocs {actual_nprocs} end
# %maxcore 1000
# 
# %scf
#     MaxIter 500
#     ConvForced true
# end
# 
# * xyzfile 0 {multiplicity} {os.path.abspath(xyz_file)}
# 
# """
#             
#             inp_file = f"{name_without_ext}.inp"
#             with open(inp_file, 'w') as f:
#                 f.write(inp_content)
#             
#             print(f"Starting native Orca DFT optimization of {filename} with {len(atoms)} atoms")
#             print(f"Formula: {atoms.get_chemical_formula()}")
#             
#             # Start timing
#             start_time = time.time()
#             
#             # Run Orca with full path (required for parallel execution)
#             try:
#                 orca_path = shutil.which('orca') or '/Applications/orca/orca'
#                 result = subprocess.run([orca_path, inp_file], 
#                                       capture_output=True, text=True, timeout=3600)
#                 
#                 # Calculate elapsed time
#                 elapsed_time, _ = calculate_and_format_time(start_time)
#                 
#                 if result.returncode != 0:
#                     print_failure_message("Orca optimization", filename, elapsed_time, result.stderr)
#                     continue
#                 
#                 print_completion_message("Orca optimization", filename, elapsed_time)
#                 
#                 # Read results
#                 opt_xyz = f"{name_without_ext}.xyz"
#                 trj_xyz = f"{name_without_ext}_trj.xyz"
#                 
#                 if os.path.exists(opt_xyz):
#                     # Copy final structure
#                     final_output = os.path.join("..", "outputs", "orca_native", f"{name_without_ext}_relaxed.xyz")
#                     shutil.copy2(opt_xyz, final_output)
#                         
#                         # Process trajectory if available
#                         if os.path.exists(trj_xyz):
#                             trajectory = read(trj_xyz, ":")
#                             
#                             # Create OVITO-compatible trajectory file
#                             traj_xyz_output = os.path.join("..", "outputs", "orca_native", f"{name_without_ext}_trajectory.xyz")
#                             create_ovito_trajectory(trajectory, traj_xyz_output, name_without_ext)
#                             
#                             # Also save as EXTXYZ with energy information
#                             extxyz_output = os.path.join("..", "outputs", "orca_native", f"{name_without_ext}_trajectory.extxyz")
#                             write(extxyz_output, trajectory, format='extxyz')
#                             
#                             # Create PDB trajectory for VMD/PyMOL
#                             pdb_output = os.path.join("..", "outputs", "orca_native", f"{name_without_ext}_trajectory.pdb")
#                             create_pdb_trajectory(trajectory, pdb_output, name_without_ext)
#                             
#                             # Save ASE trajectory for compatibility
#                             traj_file = os.path.join("..", "outputs", "orca_native", f"{name_without_ext}_relaxed.traj")
#                             traj_writer = Trajectory(traj_file, 'w')
#                             for frame in trajectory:
#                                 traj_writer.write(frame)
#                             traj_writer.close()
#                             
#                             # Generate video
#                             video_file = os.path.join("..", "outputs", "videos", f"{name_without_ext}_orca_optimization.mp4")
#                             generate_optimization_video(trajectory, video_file, name_without_ext)
#                             
#                             print(f"Optimization completed in {len(trajectory)} steps")
#                             print(f"Results saved to outputs/orca_native/")
#                             print(f"Trajectory: {name_without_ext}_trajectory.xyz (OVITO compatible)")
#                             print(f"Video saved to outputs/videos/{name_without_ext}_orca_optimization.mp4")
#                         
#                         else:
#                             print(f"Warning: No trajectory file found for {filename}")
#                     
#                     # Copy log files
#                     log_output = os.path.join("..", "logs", f"{name_without_ext}_orca.out")
#                     out_file = f"{name_without_ext}.out"
#                     if os.path.exists(out_file):
#                         shutil.copy2(out_file, log_output)
#                         
#                     # Set quality rating for successful optimization
#                     if os.path.exists(opt_xyz):
#                         trajectory_steps = len(trajectory) if os.path.exists(trj_xyz) and 'trajectory' in locals() else None
#                         if trajectory_steps:
#                             tracker.set_results(convergence_steps=trajectory_steps)
#                             tracker.set_quality(4, f"Successful optimization in {trajectory_steps} steps")
#                         else:
#                             tracker.set_quality(3, "Optimization completed, no trajectory info")
#                 
#                 except subprocess.TimeoutExpired:
#                     elapsed_time, _ = calculate_and_format_time(start_time)
#                     print_timeout_message("Orca calculation", filename, elapsed_time)
#                     tracker.set_quality(1, "Calculation timed out")
#                 except Exception as e:
#                     elapsed_time, _ = calculate_and_format_time(start_time)
#                     print_failure_message("Orca calculation", filename, elapsed_time, str(e))
#                     tracker.set_quality(1, f"Calculation failed: {str(e)}")
#             
#             print()
#         
#         # Calculate Hessian after all relaxations are complete if requested
#         if calc_hessian:
#             print("="*60)
#             print("Starting Hessian calculation after Orca relaxation completion...")
#             print(f"Using method: {hessian_method}")
#             print("="*60)
#             
#             # Find all relaxed structures created during this run
#             relaxed_files = []
#             for xyz_file in xyz_files:
#                 filename = os.path.basename(xyz_file)
#                 name_without_ext = os.path.splitext(filename)[0]
#                 relaxed_path = os.path.join(orca_native_dir, f"{name_without_ext}_relaxed.xyz")
#                 if os.path.exists(relaxed_path):
#                     relaxed_files.append(relaxed_path)
#             
#             if hessian_method == 'orca':
#                 # Use Orca for Hessian calculation with same method/basis as relaxation
#                 for relaxed_file in relaxed_files:
#                     calculate_orca_hessian(relaxed_file, method, basis, actual_nprocs)
#             else:
#                 # Use phonopy/DFTB+ for Hessian calculation (default)
#                 project_name = os.path.basename(base_path)
#                 for relaxed_file in relaxed_files:
#                     filename = os.path.basename(relaxed_file)
#                     name_without_ext = os.path.splitext(filename)[0]
#                     # Use the existing calculate_hessian function with specific structure file
#                     try:
#                         relative_path = os.path.relpath(relaxed_file, os.path.join(base_path, "outputs"))
#                         calculate_hessian(project_name, relative_path)
#                     except Exception as e:
#                         print(f"Error calculating Hessian for {filename}: {e}")
#             
#             print("="*60)
#             print("Hessian calculations completed!")
#             print("="*60)
#     
#     finally:
#         os.chdir(original_cwd)
#     
#     # All files processed successfully
#     return True
# 
# 
def relax_orca_dft(project_folder, method='B3LYP', basis='6-31G(d)', nprocs=4):
    """
    Relax using Orca DFT with ASE optimizer and trajectory video generation
    
    Args:
        project_folder (str): Name of the project folder or molecule name
        method (str): DFT method (default: B3LYP)
        basis (str): Basis set (default: 6-31G(d) - faster than def2-SVP with similar geometry accuracy)
        nprocs (int): Number of processors for Orca (default: 4)
    """
    base_path = f"data/{project_folder}"
    
    # Find project folder
    if not os.path.exists(base_path):
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
        else:
            raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    # Find all XYZ files
    xyz_files = []
    xyz_pattern = os.path.join(base_path, "*.xyz")
    xyz_files.extend(glob.glob(xyz_pattern))
    
    inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
    xyz_files.extend(glob.glob(inputs_pattern))
    
    xyz_files = [os.path.abspath(f) for f in xyz_files]
    
    if not xyz_files:
        print(f"No XYZ files found in {base_path}")
        return False
    
    # Create directories
    logs_dir = os.path.join(base_path, "logs")
    outputs_dir = os.path.join(base_path, "outputs")
    orca_outputs_dir = os.path.join(base_path, "outputs", "orca_dft")
    videos_dir = os.path.join(base_path, "outputs", "videos")
    temp_dir = os.path.join(base_path, "temp")
    
    for dir_path in [logs_dir, outputs_dir, orca_outputs_dir, videos_dir, temp_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    original_cwd = os.getcwd()
    
    try:
        os.chdir(base_path)
    
        for xyz_file in xyz_files:
            filename = os.path.basename(xyz_file)
            name_without_ext = os.path.splitext(filename)[0]
            
            print(f"Processing {filename} with Orca DFT ({method}/{basis})...")
            
            atoms = read(xyz_file)
            
            # Determine number of processors
            if platform.system() == 'Darwin':
                actual_nprocs = min(nprocs, 4)  # Cap at 4 cores for reasonable performance
            else:
                actual_nprocs = nprocs
            
            # Display time estimate
            time_estimate = estimate_time_range(atoms, 'orca_relax', method, actual_nprocs)
            print(f"⏱️  Оценка времени релаксации с ASE: {time_estimate}")
            print(f"    Формула: {atoms.get_chemical_formula()} ({len(atoms)} атомов)")
            
            # Calculate number of electrons and determine multiplicity
            n_electrons = sum(atoms.get_atomic_numbers()) - 0  # assuming neutral molecule
            multiplicity = 2 if n_electrons % 2 == 1 else 1  # odd electrons = doublet, even = singlet
            
            orca_simple = f"{method} {basis}"
            orca_blocks = f"""%pal nprocs {actual_nprocs} end
%maxcore 1000
%scf
    MaxIter 500
    ConvForced true
end
"""
            
            calc = ORCA(label=f"temp/{name_without_ext}_orca",
                       orcasimpleinput=orca_simple,
                       orcablocks=orca_blocks,
                       charge=0,
                       mult=multiplicity)
            atoms.calc = calc
            
            # Set up trajectory and optimizer
            traj_file = os.path.join("outputs", "orca_dft", f"{name_without_ext}_relaxed.traj")
            log_file = os.path.join("logs", f"{name_without_ext}_orca_relax.log")
        
            traj = Trajectory(traj_file, 'w', atoms)
            opt = BFGS(atoms, logfile=log_file, maxstep=0.1)
            opt.attach(traj.write, interval=1)
            
            print(f"Starting Orca DFT optimization of {filename} with {len(atoms)} atoms")
            print(f"Formula: {atoms.get_chemical_formula()}")
            
            # Start timing
            start_time = time.time()
            
            # Run optimization
            try:
                opt.run(fmax=0.02, steps=300)
                
                # Write final structure
                final_xyz = os.path.join("outputs", "orca_dft", f"{name_without_ext}_relaxed.xyz")
                write(final_xyz, atoms)
                
                # Convert trajectory to OVITO-compatible formats
                trajectory = Trajectory(traj_file)
                trajectory_list = [atoms for atoms in trajectory]
                
                # Create OVITO-compatible trajectory file
                traj_xyz_output = os.path.join("outputs", "orca_dft", f"{name_without_ext}_trajectory.xyz")
                create_ovito_trajectory(trajectory_list, traj_xyz_output, name_without_ext)
                
                # Save as EXTXYZ with energy information
                extxyz_output = os.path.join("outputs", "orca_dft", f"{name_without_ext}_trajectory.extxyz")
                write(extxyz_output, trajectory_list, format='extxyz')
                
                # Generate video from trajectory
                video_file = os.path.join("outputs", "videos", f"{name_without_ext}_orca_ase_optimization.mp4")
                generate_optimization_video(trajectory_list, video_file, name_without_ext)
                
                print(f"Results saved to outputs/orca_dft/")
                print(f"Trajectory: {name_without_ext}_trajectory.xyz (OVITO compatible)")
                print(f"Video saved to outputs/videos/{name_without_ext}_orca_ase_optimization.mp4")
                
                # Calculate elapsed time
                elapsed_time, _ = calculate_and_format_time(start_time)
                
                # Print final results
                final_energy = atoms.get_potential_energy()
                max_force = np.max(np.linalg.norm(atoms.get_forces(), axis=1))
                print(f"Final energy: {final_energy:.6f} eV")
                print(f"Max force: {max_force:.4f} eV/Å")
                print_completion_message("Orca DFT optimization", filename, elapsed_time)
                
            except Exception as e:
                elapsed_time, _ = calculate_and_format_time(start_time)
                print_failure_message("Orca DFT optimization", filename, elapsed_time, str(e))
            
            print()
    
    finally:
        os.chdir(original_cwd)
    
    # All files processed successfully
    return True


def relax_xtb(project_folder, method='GFN2-xTB', nprocs=1, backend='native'):
    """
    Relax using xTB semiempirical quantum chemistry optimization
    
    Args:
        project_folder (str): Name of the project folder or molecule name
        method (str): xTB method (GFN2-xTB, GFN1-xTB, GFN0-xTB) (default: GFN2-xTB)
        nprocs (int): Number of processors (default: 1, xtb has limited parallel scaling)
        backend (str): xTB backend - 'native' (system xtb, has v6.7.1 bug), 
                      'gxtb' (Docker), 'fixed' (our patched version, recommended)
    
    Note: Native xTB v6.7.1 has a Fortran format bug causing "Missing comma between descriptors".
          For reliable calculations, use --xtb-backend gxtb (Docker container with g-xTB).
    """
    base_path = f"data/{project_folder}"
    
    # Find project folder
    if not os.path.exists(base_path):
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
        else:
            raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    # Find all XYZ files
    xyz_files = []
    xyz_pattern = os.path.join(base_path, "*.xyz")
    xyz_files.extend(glob.glob(xyz_pattern))
    
    inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
    xyz_files.extend(glob.glob(inputs_pattern))
    
    xyz_files = [os.path.abspath(f) for f in xyz_files]
    
    if not xyz_files:
        print(f"No XYZ files found in {base_path}")
        return False
    
    # Create directories
    logs_dir = os.path.join(base_path, "logs")
    outputs_dir = os.path.join(base_path, "outputs")
    xtb_outputs_dir = os.path.join(base_path, "outputs", "xtb")
    videos_dir = os.path.join(base_path, "outputs", "videos")
    temp_dir = os.path.join(base_path, "temp")
    
    for dir_path in [logs_dir, outputs_dir, xtb_outputs_dir, videos_dir, temp_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    original_cwd = os.getcwd()
    
    # Method mapping
    method_flags = {
        'GFN2-xTB': '--gfn 2',
        'GFN1-xTB': '--gfn 1', 
        'GFN0-xTB': '--gfn 0'
    }
    
    method_flag = method_flags.get(method, '--gfn 2')
    
    try:
        os.chdir(temp_dir)
        
        # Warn user about native backend issues
        if backend == 'native':
            print("⚠️  Warning: Native xTB v6.7.1 has a known Fortran format bug")
            print("   Error: 'Missing comma between descriptors' in optimizer.f90:852")
            print("   💡 Recommendation: Use --xtb-backend fixed (our patched version) or --xtb-backend gxtb")
            print()
        
        for xyz_file in xyz_files:
            filename = os.path.basename(xyz_file)
            name_without_ext = os.path.splitext(filename)[0]
            
            print(f"Processing {filename} with xTB ({method}, backend: {backend})...")
            
            atoms = read(xyz_file)
            
            # Display time estimate
            time_estimate = estimate_time_range(atoms, 'xtb_relax')
            print(f"⏱️  Оценка времени xTB оптимизации: {time_estimate}")
            print(f"    Формула: {atoms.get_chemical_formula()} ({len(atoms)} атомов)")
            
            # Copy input file to temp directory
            local_xyz = f"{name_without_ext}.xyz"
            shutil.copy2(xyz_file, local_xyz)
            
            print(f"Starting xTB optimization of {filename} with {len(atoms)} atoms")
            print(f"Formula: {atoms.get_chemical_formula()}")
            
            # Detect unpaired electrons for radical systems
            unpaired_electrons = detect_unpaired_electrons(atoms)
            if unpaired_electrons > 0:
                print(f"🔬 Detected radical system with {unpaired_electrons} unpaired electron(s)")
                print(f"   Using unrestricted Hartree-Fock (UHF) calculation")

            # Start timing
            start_time = time.time()
            
            # Choose backend for xTB calculation
            if backend == 'native':
                success = _run_native_xtb(xyz_file, filename, name_without_ext, method_flag, nprocs, start_time, project_folder, unpaired_electrons)
            elif backend == 'fixed':
                success = _run_fixed_xtb(xyz_file, filename, name_without_ext, method_flag, nprocs, start_time, project_folder, unpaired_electrons)
            elif backend == 'gxtb':
                success = _run_gxtb_backend(xyz_file, filename, name_without_ext, atoms, start_time)
            else:
                print(f"Unknown backend '{backend}'. Available backends: native, fixed, gxtb")
                continue
                    
            # If any file fails, return False immediately
            if not success:
                print(f"❌ Relaxation failed for {filename}. Stopping pipeline.")
                return False
            
            print()
    
    finally:
        os.chdir(original_cwd)
    
    # If we get here, all files were successful
    return True


def _run_native_xtb(xyz_file, filename, name_without_ext, method_flag, nprocs, start_time, project_folder, unpaired_electrons=0):
    """
    Run native system xTB optimization
    
    Args:
        xyz_file (str): Path to initial structure file
        project_folder (str): Project folder name
        unpaired_electrons (int): Number of unpaired electrons (0 for closed-shell)
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Set environment variable for OpenMP threads if nprocs > 1
        env = os.environ.copy()
        if nprocs > 1:
            env['OMP_NUM_THREADS'] = str(nprocs)
        
        # Ensure temp directory is clean and temp files stay contained
        temp_files_before = set(os.listdir('.'))
        
        # Use our wrapper to handle v6.7.1 bug automatically
        wrapper_path = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'xtb_wrapper.py'))
        if os.path.exists(wrapper_path):
            xtb_cmd = ['python3', wrapper_path, filename, '--opt', method_flag, '--cycles', '500']
        else:
            # Fallback to direct xtb call
            xtb_cmd = ['xtb', filename, '--opt', method_flag, '--cycles', '500']
        
        # Add UHF flag for radical systems
        if unpaired_electrons > 0:
            xtb_cmd.extend(['--uhf', str(unpaired_electrons)])
        
        result = subprocess.run(xtb_cmd, capture_output=True, text=True, timeout=1800, env=env)
        
        # Clean up any temporary files created by xTB in current directory
        temp_files_after = set(os.listdir('.'))
        temp_files_created = temp_files_after - temp_files_before
        for temp_file in temp_files_created:
            if temp_file not in [filename, 'xtbopt.xyz']:  # Keep input and main output
                try:
                    os.remove(temp_file)
                except OSError:
                    pass
        
        elapsed_time, _ = calculate_and_format_time(start_time)
        
        if result.returncode != 0:
            # Check for the known v6.7.1 bug
            if "Missing comma between descriptors" in result.stderr or "(1x," in result.stderr:
                print(f"❌ Native xTB failed due to v6.7.1 bug: Fortran format string error")
                print(f"   Error: 'Missing comma between descriptors' in optimizer.f90:852")
                print(f"   This is a known bug in xTB v6.7.1")
                print(f"   💡 Suggestion: Use --xtb-backend gxtb or --xtb-backend fixed")
                return False
            else:
                print_failure_message("Native xTB optimization", filename, elapsed_time, result.stderr)
                return False
        
        print_completion_message("Native xTB optimization", filename, elapsed_time)
        
        # Process native xTB output
        opt_xyz = 'xtbopt.xyz'  # Standard xTB output file
        
        if os.path.exists(opt_xyz):
            final_output = os.path.join("..", "outputs", "xtb", f"{name_without_ext}_relaxed.xyz")
            shutil.copy2(opt_xyz, final_output)
            
            print(f"Results saved to outputs/xtb/")
            print(f"Final structure: {name_without_ext}_relaxed.xyz")
            
            # Run structural analysis
            run_structural_analysis_if_enabled(
                xyz_file, final_output, project_folder, "xTB (native)"
            )
            
            # Copy optimization log
            log_dest = os.path.join("..", "logs", f"{name_without_ext}_xtb_native.log")
            with open(log_dest, 'w') as f:
                f.write("Native xTB Optimization Log\n")
                f.write("="*50 + "\n")
                f.write(f"Command: {' '.join(xtb_cmd)}\n")
                f.write(f"Return code: {result.returncode}\n")
                f.write(f"Elapsed time: {elapsed_time}\n")
                f.write("\nSTDOUT:\n")
                f.write(result.stdout)
                f.write("\n\nSTDERR:\n")
                f.write(result.stderr)
            
            # Try to extract energy from stdout
            if 'TOTAL ENERGY' in result.stdout:
                import re
                energy_match = re.search(r'TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh', result.stdout)
                if energy_match:
                    energy_hartree = float(energy_match.group(1))
                    print(f"Final energy: {energy_hartree:.6f} Hartree")
            
            return True
        else:
            print(f"Warning: No optimized structure found from native xTB")
            return False
            
    except subprocess.TimeoutExpired:
        elapsed_time, _ = calculate_and_format_time(start_time)
        print_timeout_message("Native xTB calculation", filename, elapsed_time)
        return False
    except Exception as e:
        elapsed_time, _ = calculate_and_format_time(start_time)
        print_failure_message("Native xTB calculation", filename, elapsed_time, str(e))
        return False


def _run_fixed_xtb(xyz_file, filename, name_without_ext, method_flag, nprocs, start_time, project_folder, unpaired_electrons=0):
    """
    Run patched xTB optimization from our fixed xTB fork at 3rdparty/xtb/fork/build_cmake/xtb
    
    Args:
        xyz_file (str): Path to initial structure file
        project_folder (str): Project folder name
        unpaired_electrons (int): Number of unpaired electrons (0 for closed-shell)
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Check if fixed xtb binary exists in our fork location
        # Navigate back to the root of the project from temp directory 
        script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # Gets mechanosynthesis root
        fixed_xtb_path = os.path.join(script_dir, "3rdparty", "xtb", "fork", "build_cmake", "xtb")
        
        if not os.path.exists(fixed_xtb_path):
            print(f"🔧 Fixed xTB binary not found. Building from patched source...")
            
            # Try to build the fixed version using CMake
            original_cwd = os.getcwd()
            try:
                xtb_fork_dir = os.path.join(script_dir, "3rdparty", "xtb", "fork")
                os.chdir(xtb_fork_dir)
                
                # Set up environment for macOS build  
                env = os.environ.copy()
                env['MACOSX_DEPLOYMENT_TARGET'] = '14.0'
                env['SDKROOT'] = subprocess.check_output(['xcrun', '--show-sdk-path'], text=True).strip()
                
                # Configure with CMake using homebrew OpenBLAS
                cmake_cmd = [
                    'cmake', '-B', 'build_cmake',
                    '-DCMAKE_Fortran_COMPILER=gfortran',
                    '-DBLAS_LIBRARIES=/opt/homebrew/opt/openblas/lib/libopenblas.dylib',
                    '-DLAPACK_LIBRARIES=/opt/homebrew/opt/openblas/lib/libopenblas.dylib'
                ]
                
                build_result = subprocess.run(cmake_cmd, capture_output=True, text=True, timeout=300, env=env)
                
                if build_result.returncode == 0:
                    # Build the xTB executable
                    compile_result = subprocess.run(['cmake', '--build', 'build_cmake', '--parallel', '4', '--target', 'xtb-exe'], 
                                                  capture_output=True, text=True, timeout=600, env=env)
                    
                    if compile_result.returncode == 0:
                        # Check if binary was created
                        fixed_xtb_path = os.path.join(xtb_fork_dir, "build_cmake", "xtb")
                        if os.path.exists(fixed_xtb_path):
                            print(f"✅ Successfully built fixed xTB at {fixed_xtb_path}")
                        else:
                            print(f"❌ Build succeeded but binary not found at expected location")
                            return False
                    else:
                        print(f"❌ Failed to compile fixed xTB: {compile_result.stderr}")
                        return False
                else:
                    print(f"❌ Failed to configure build for fixed xTB: {build_result.stderr}")
                    print(f"💡 Suggestion: Use --xtb-backend gxtb as alternative")
                    return False
                    
            finally:
                os.chdir(original_cwd)
        
        if not os.path.exists(fixed_xtb_path):
            print(f"Error: Fixed xTB binary still not found at {fixed_xtb_path}")
            return False
        
        # Set environment variable for OpenMP threads if nprocs > 1
        env = os.environ.copy()
        if nprocs > 1:
            env['OMP_NUM_THREADS'] = str(nprocs)
        
        # Track temporary files for cleanup
        temp_files_before = set(os.listdir('.'))
        
        # Run fixed xtb optimization
        xtb_cmd = [fixed_xtb_path, filename, '--opt', method_flag, '--cycles', '500']
        
        # Add UHF flag for radical systems
        if unpaired_electrons > 0:
            xtb_cmd.extend(['--uhf', str(unpaired_electrons)])
        
        result = subprocess.run(xtb_cmd, capture_output=True, text=True, timeout=1800, env=env)
        
        # Clean up temporary files created by xTB
        temp_files_after = set(os.listdir('.'))
        temp_files_created = temp_files_after - temp_files_before
        for temp_file in temp_files_created:
            if temp_file not in [filename, 'xtbopt.xyz']:  # Keep input and main output
                try:
                    os.remove(temp_file)
                except OSError:
                    pass
        
        elapsed_time, _ = calculate_and_format_time(start_time)
        
        if result.returncode != 0:
            print_failure_message("Fixed xTB optimization", filename, elapsed_time, result.stderr)
            return False
        
        print_completion_message("Fixed xTB optimization", filename, elapsed_time)
        
        # Process fixed xTB output
        opt_xyz = 'xtbopt.xyz'
        
        if os.path.exists(opt_xyz):
            final_output = os.path.join("..", "outputs", "xtb", f"{name_without_ext}_relaxed.xyz")
            shutil.copy2(opt_xyz, final_output)
            
            print(f"Results saved to outputs/xtb/")
            print(f"Final structure: {name_without_ext}_relaxed.xyz")
            
            # Run structural analysis
            run_structural_analysis_if_enabled(
                xyz_file, final_output, project_folder, "xTB (fixed)"
            )
            
            # Copy optimization log
            log_dest = os.path.join("..", "logs", f"{name_without_ext}_xtb_fixed.log")
            with open(log_dest, 'w') as f:
                f.write("Fixed xTB Optimization Log\n")
                f.write("="*50 + "\n")
                f.write(f"Command: {' '.join(xtb_cmd)}\n")
                f.write(f"Return code: {result.returncode}\n")
                f.write(f"Elapsed time: {elapsed_time}\n")
                f.write("\nSTDOUT:\n")
                f.write(result.stdout)
                f.write("\n\nSTDERR:\n")
                f.write(result.stderr)
            
            # Try to extract energy from stdout
            if 'TOTAL ENERGY' in result.stdout:
                import re
                energy_match = re.search(r'TOTAL ENERGY\s+(-?\d+\.\d+)\s+Eh', result.stdout)
                if energy_match:
                    energy_hartree = float(energy_match.group(1))
                    print(f"Final energy: {energy_hartree:.6f} Hartree")
            
            return True
        else:
            print(f"Warning: No optimized structure found from fixed xTB")
            return False
            
    except subprocess.TimeoutExpired:
        elapsed_time, _ = calculate_and_format_time(start_time)
        print_timeout_message("Fixed xTB calculation", filename, elapsed_time)
        return False
    except Exception as e:
        elapsed_time, _ = calculate_and_format_time(start_time)
        print_failure_message("Fixed xTB calculation", filename, elapsed_time, str(e))
        return False


def _run_gxtb_backend(xyz_file, filename, name_without_ext, atoms, start_time):
    """
    Run g-xTB via Docker backend
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Create output directory for g-xTB backend 
        docker_output_dir = os.path.join("..", "outputs", "xtb", name_without_ext)
        os.makedirs(docker_output_dir, exist_ok=True)
        abs_docker_output_dir = os.path.abspath(docker_output_dir)
        
        # Read atoms for Docker input
        n_atoms = len(atoms)
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        
        # Create XYZ content for Docker
        xyz_content = f"{n_atoms}\n{name_without_ext} molecule\n"
        for sym, pos in zip(symbols, positions):
            xyz_content += f"{sym:2s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n"
        
        # Run g-xTB via Docker container
        docker_cmd = [
            'docker', 'run', '--rm',
            '-v', f'{abs_docker_output_dir}:/data',
            'gxtb:latest',
            'bash', '-c',
            f'''cat > /tmp/{filename} << 'EOF'
{xyz_content}EOF
python3 /usr/local/bin/run_gxtb.py /tmp/{filename} --relax --output-dir /data'''
        ]
        
        result = subprocess.run(docker_cmd, capture_output=True, text=True, timeout=1800)
        
        elapsed_time, _ = calculate_and_format_time(start_time)
        
        if result.returncode != 0:
            print_failure_message("xTB optimization (g-xTB backend)", filename, elapsed_time, f"Docker error: {result.stderr}")
            return False
        
        print_completion_message("xTB optimization (g-xTB backend)", filename, elapsed_time)
        
        # Process Docker output files (g-xTB backend)
        docker_files = os.listdir(docker_output_dir)
        
        # Look for optimized geometry from g-xTB backend (in coord format)
        final_output = os.path.join("..", "outputs", "xtb", f"{name_without_ext}_relaxed.xyz")
        optimized_found = False
        
        # g-xTB creates a "coord" file in Turbomole format - convert it to XYZ
        if 'coord' in docker_files:
            coord_file = os.path.join(docker_output_dir, 'coord')
            
            # Convert coord (Turbomole) to XYZ format
            try:
                optimized_atoms = []
                with open(coord_file, 'r') as f:
                    lines = f.readlines()
                
                for line in lines:
                    if line.strip().startswith('$') or not line.strip():
                        continue
                    parts = line.strip().split()
                    if len(parts) >= 4:
                        # Convert from Bohr to Angstrom (1 Bohr = 0.529177 Angstrom)
                        x_ang = float(parts[0]) * 0.529177249
                        y_ang = float(parts[1]) * 0.529177249
                        z_ang = float(parts[2]) * 0.529177249
                        element = parts[3].capitalize()
                        optimized_atoms.append((element, x_ang, y_ang, z_ang))
                
                # Write XYZ file
                with open(final_output, 'w') as f:
                    f.write(f"{len(optimized_atoms)}\n")
                    f.write(f"Optimized {name_without_ext} via g-xTB\n")
                    for element, x, y, z in optimized_atoms:
                        f.write(f"{element:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
                
                optimized_found = True
                
            except Exception as e:
                print(f"Error converting coord to XYZ: {e}")
                pass
        
        if optimized_found:
            print(f"Results saved to outputs/xtb/")
            print(f"Final structure: {name_without_ext}_relaxed.xyz")
            
            # Copy optimization log from g-xTB backend output if available
            if 'optimization.log' in docker_files:
                log_source = os.path.join(docker_output_dir, 'optimization.log')
                log_dest = os.path.join("..", "logs", f"{name_without_ext}_xtb_gxtb.log")
                shutil.copy2(log_source, log_dest)
                
                # Parse log for energy information
                try:
                    with open(log_source, 'r') as f:
                        log_content = f.read()
                        
                    # Extract final energy from log
                    if 'total' in log_content:
                        import re
                        energy_match = re.search(r'total\s+(-?\d+\.\d+)', log_content)
                        if energy_match:
                            energy_hartree = float(energy_match.group(1))
                            print(f"Final energy: {energy_hartree:.6f} Hartree")
                except Exception:
                    pass
        else:
            print(f"Warning: No optimized structure found in Docker output")
        
        # Save Docker execution log
        docker_log = os.path.join("..", "logs", f"{name_without_ext}_xtb_gxtb_docker.log")
        with open(docker_log, 'w') as f:
            f.write("xTB Optimization via g-xTB Docker Backend\n")
            f.write("="*50 + "\n")
            f.write(f"Command: {' '.join(docker_cmd[:5])} ... [truncated]\n")
            f.write(f"Return code: {result.returncode}\n")
            f.write(f"Elapsed time: {elapsed_time}\n")
            f.write("\nSTDOUT:\n")
            f.write(result.stdout)
            f.write("\n\nSTDERR:\n")
            f.write(result.stderr)
        
        return optimized_found
        
    except subprocess.TimeoutExpired:
        elapsed_time, _ = calculate_and_format_time(start_time)
        print_timeout_message("xTB calculation (g-xTB backend)", filename, elapsed_time)
        return False
    except Exception as e:
        elapsed_time, _ = calculate_and_format_time(start_time)
        print_failure_message("xTB calculation (g-xTB backend)", filename, elapsed_time, str(e))
        return False


def relax_gxtb(project_folder):
    """
    Relax using g-xTB semiempirical quantum chemistry optimization via Docker
    
    Args:
        project_folder (str): Name of the project folder or molecule name
    """
    base_path = f"data/{project_folder}"
    
    # Find project folder
    if not os.path.exists(base_path):
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
        else:
            raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    # Find all XYZ files
    xyz_files = []
    xyz_pattern = os.path.join(base_path, "*.xyz")
    xyz_files.extend(glob.glob(xyz_pattern))
    
    inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
    xyz_files.extend(glob.glob(inputs_pattern))
    
    xyz_files = [os.path.abspath(f) for f in xyz_files]
    
    if not xyz_files:
        print(f"No XYZ files found in {base_path}")
        return False
    
    # Create directories
    logs_dir = os.path.join(base_path, "logs")
    outputs_dir = os.path.join(base_path, "outputs")
    gxtb_outputs_dir = os.path.join(base_path, "outputs", "gxtb")
    videos_dir = os.path.join(base_path, "outputs", "videos")
    
    for dir_path in [logs_dir, outputs_dir, gxtb_outputs_dir, videos_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    original_cwd = os.getcwd()
    
    try:
        for xyz_file in xyz_files:
            filename = os.path.basename(xyz_file)
            name_without_ext = os.path.splitext(filename)[0]
            
            print(f"Processing {filename} with g-xTB (Docker)...")
            
            atoms = read(xyz_file)
            
            # Display time estimate (similar to xTB)
            time_estimate = estimate_time_range(atoms, 'xtb_relax')
            print(f"⏱️  Estimated g-xTB optimization time: {time_estimate}")
            print(f"    Formula: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")
            
            print(f"Running g-xTB geometry optimization via Docker...")
            
            start_time = time.time()
            
            # Use g-xTB Docker container for actual g-xTB calculation
            docker_output_dir = os.path.join(gxtb_outputs_dir, name_without_ext)
            os.makedirs(docker_output_dir, exist_ok=True)
            
            # Convert to absolute path for Docker
            abs_docker_output_dir = os.path.abspath(docker_output_dir)
            
            # Read the XYZ file content to embed in Docker command
            atoms = read(xyz_file)
            n_atoms = len(atoms)
            symbols = atoms.get_chemical_symbols()
            positions = atoms.get_positions()
            
            # Create XYZ content
            xyz_content = f"{n_atoms}\n{name_without_ext} molecule\n"
            for sym, pos in zip(symbols, positions):
                xyz_content += f"{sym:2s} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n"
            
            # Run g-xTB via Docker container
            docker_cmd = [
                'docker', 'run', '--rm',
                '-v', f'{abs_docker_output_dir}:/data',
                'gxtb:latest',
                'bash', '-c',
                f'''cat > /tmp/{filename} << 'EOF'
{xyz_content}EOF
python3 /usr/local/bin/run_gxtb.py /tmp/{filename} --relax --output-dir /data'''
            ]
            
            try:
                result = subprocess.run(docker_cmd, capture_output=True, text=True, timeout=1800)
                
                elapsed_time, _ = calculate_and_format_time(start_time)
                
                if result.returncode == 0:
                    print_completion_message("g-xTB optimization", filename, elapsed_time)
                    
                    # Process Docker output files
                    docker_files = os.listdir(docker_output_dir)
                    
                    # Copy optimization results
                    final_output = os.path.join(gxtb_outputs_dir, f"{name_without_ext}_gxtb_relaxed.xyz")
                    
                    # Look for optimized geometry from g-xTB - it might be under different names
                    optimized_found = False
                    for possible_file in [f'{filename}', 'xtbopt.xyz', 'opt.xyz']:
                        if possible_file in docker_files:
                            source_file = os.path.join(docker_output_dir, possible_file)
                            shutil.copy2(source_file, final_output)
                            optimized_found = True
                            
                            # Read final structure to get energy info
                            try:
                                final_atoms = read(source_file)
                                if hasattr(final_atoms, 'info') and 'energy' in final_atoms.info:
                                    print(f"Final energy: {final_atoms.info['energy']:.6f} Hartree")
                            except Exception:
                                pass
                            break
                    
                    if optimized_found:
                        print(f"Results saved to outputs/gxtb/")
                        print(f"Final structure: {name_without_ext}_gxtb_relaxed.xyz")
                    else:
                        print(f"Warning: No optimized structure found in Docker output")
                    
                    # Copy optimization log
                    if 'optimization.log' in docker_files:
                        log_source = os.path.join(docker_output_dir, 'optimization.log')
                        log_dest = os.path.join(logs_dir, f"{name_without_ext}_gxtb.log")
                        shutil.copy2(log_source, log_dest)
                        
                        # Parse log for additional info
                        try:
                            with open(log_source, 'r') as f:
                                log_content = f.read()
                                
                            # Extract final energy from log
                            if 'total' in log_content:
                                import re
                                energy_match = re.search(r'total\s+(-?\d+\.\d+)', log_content)
                                if energy_match:
                                    energy_hartree = float(energy_match.group(1))
                                    print(f"Final energy: {energy_hartree:.6f} Hartree")
                        except Exception:
                            pass
                        
                else:
                    print_failure_message("g-xTB optimization", filename, elapsed_time)
                    print(f"Docker error: {result.stderr}")
                
                # Save Docker execution log
                docker_log = os.path.join(logs_dir, f"{name_without_ext}_gxtb_docker.log")
                with open(docker_log, 'w') as f:
                    f.write("g-xTB Docker Container Output\n")
                    f.write("="*50 + "\n")
                    f.write(f"Command: {' '.join(docker_cmd[:5])} ... [truncated]\n")
                    f.write(f"Return code: {result.returncode}\n")
                    f.write(f"Elapsed time: {elapsed_time}\n")
                    f.write("\nSTDOUT:\n")
                    f.write(result.stdout)
                    f.write("\n\nSTDERR:\n")
                    f.write(result.stderr)
                
            except subprocess.TimeoutExpired:
                elapsed_time, _ = calculate_and_format_time(start_time)
                print_timeout_message("g-xTB Docker optimization", filename, elapsed_time)
            except Exception as e:
                elapsed_time, _ = calculate_and_format_time(start_time)
                print_failure_message("g-xTB Docker optimization", filename, elapsed_time, str(e))
                
    finally:
        os.chdir(original_cwd)
    
    # All files processed successfully
    return True


def parse_xtb_trajectory(log_file):
    """
    Parse xTB optimization log to extract trajectory frames
    
    Args:
        log_file (str): Path to xtb optimization log file
        
    Returns:
        list: List of ASE atoms objects representing trajectory frames
    """
    trajectory_frames = []
    
    try:
        with open(log_file, 'r') as f:
            content = f.read()
        
        # Try to find coordinate blocks in the log
        # This is a simplified parser - xTB log format can vary
        lines = content.split('\n')
        
        current_coords = []
        current_symbols = []
        in_coord_block = False
        
        for line in lines:
            line = line.strip()
            
            # Look for coordinate sections (this may need adjustment based on xTB version)
            if 'molecular geometry' in line.lower() or 'coordinates' in line.lower():
                in_coord_block = True
                current_coords = []
                current_symbols = []
                continue
            
            # End of coordinate block
            if in_coord_block and (line == '' or line.startswith('---') or 'energy' in line.lower()):
                if current_coords and current_symbols:
                    atoms = create_atoms_from_coords(current_symbols, current_coords)
                    if atoms is not None:
                        trajectory_frames.append(atoms)
                in_coord_block = False
                continue
            
            # Parse coordinate line
            if in_coord_block:
                parts = line.split()
                if len(parts) >= 4:
                    try:
                        symbol = parts[0]
                        x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                        current_symbols.append(symbol)
                        current_coords.append([x, y, z])
                    except (ValueError, IndexError):
                        continue
        
        return trajectory_frames
        
    except Exception as e:
        print(f"Warning: Could not parse xTB trajectory: {e}")
        return []


def create_atoms_from_coords(symbols, coordinates):
    """Create ASE atoms object from symbols and coordinates"""
    try:
        from ase import Atoms
        atoms = Atoms(symbols=symbols, positions=coordinates)
        return atoms
    except Exception:
        return None


def relax_native_dftb_ptbp(project_folder):
    """
    Relax using native DFTB+ geometry optimization with PTBP parameters (supports tungsten)
    
    Args:
        project_folder (str): Name of the project folder or molecule name
        
    Returns:
        bool: True if successful, False if failed
    """
    base_path = f"data/{project_folder}"
    
    # If exact folder doesn't exist, try to find a folder containing the molecule name
    if not os.path.exists(base_path):
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
        else:
            raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    if not os.path.exists(base_path):
        raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    # Convert to absolute path 
    base_path = os.path.abspath(base_path)
    
    # Check if PTBP parameters are available
    ptbp_dir = os.path.join(os.environ.get('SKDIR', '/Users/andreypanferov/opt/dftb+/slakos'), 'ptbp-2024')
    if not os.path.exists(ptbp_dir):
        print(f"❌ PTBP parameters not found at {ptbp_dir}")
        print("Please download and extract PTBP parameters first.")
        return False
    
    # Find all XYZ files in base directory and inputs subdirectory
    xyz_files = []
    xyz_pattern = os.path.join(base_path, "*.xyz")
    xyz_files.extend(glob.glob(xyz_pattern))
    
    inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
    xyz_files.extend(glob.glob(inputs_pattern))
    
    # Convert to absolute paths to handle directory changes
    xyz_files = [os.path.abspath(f) for f in xyz_files]
    
    if not xyz_files:
        print(f"No XYZ files found in {base_path}")
        return False
    
    # Create logs, outputs, and temp directories
    logs_dir = os.path.join(base_path, "logs")
    outputs_dir = os.path.join(base_path, "outputs")
    ptbp_outputs_dir = os.path.join(base_path, "outputs", "ptbp_dftb")
    temp_dir = os.path.join(base_path, "temp")
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    os.makedirs(ptbp_outputs_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Store original working directory
    original_cwd = os.getcwd()
    
    try:
        # Change to temp directory for DFTB+ execution
        os.chdir(temp_dir)
        
        for xyz_file in xyz_files:
            filename = os.path.basename(xyz_file)
            name_without_ext = os.path.splitext(filename)[0]
            
            print(f"Processing {filename} with native DFTB+ optimization (PTBP parameters)...")
            
            # Read atoms and convert to DFTB+ gen format
            atoms = read(xyz_file)
            
            # Display time estimate
            time_estimate = estimate_time_range(atoms, 'dftb_native')
            print(f"⏱️  Оценка времени оптимизации: {time_estimate}")
            print(f"    Формула: {atoms.get_chemical_formula()} ({len(atoms)} атомов)")
            gen_file = f"{name_without_ext}.gen"
            write(gen_file, atoms, format='gen')
            
            # Generate angular momentum parameters for PTBP
            angular_params = get_angular_momentum_params_ptbp(atoms)
            
            # Create DFTB+ input file with geometry optimization using PTBP
            hsd_content = f"""
Geometry = GenFormat {{
    <<< "{gen_file}"
}}

Driver = GeometryOptimization {{
    Optimizer = Rational {{}}
    MaxSteps = 500
    OutputPrefix = "{name_without_ext}_opt"
}}

Hamiltonian = DFTB {{
    SCC = Yes
    SCCTolerance = 1e-6
    MaxSCCIterations = 500
    Filling = Fermi {{Temperature[K] = 300.0}}
    HCorrection = Damping {{Exponent = 4.0}}
    Mixer = Anderson {{MixingParameter = 0.05}}
    SlaterKosterFiles = Type2FileNames {{
        Prefix = "{ptbp_dir}/"
        Separator = "-"
        Suffix = ".skf"
    }}
    MaxAngularMomentum = {{
{angular_params}
    }}
}}

Options {{
    WriteResultsTag = Yes
}}
"""
            
            # Write DFTB+ input file
            hsd_file = f"dftb_in.hsd"
            with open(hsd_file, 'w') as f:
                f.write(hsd_content)
            
            print(f"Starting native DFTB+ optimization for {filename} with {len(atoms)} atoms")
            print(f"Formula: {atoms.get_chemical_formula()}")
            
            start_time = time.time()
            print("Running DFTB+ optimization...")
            
            # Run DFTB+
            try:
                dftb_executable = '/Users/andreypanferov/opt/dftb+/bin/dftb+'
                result = subprocess.run([dftb_executable], 
                                     capture_output=True, 
                                     text=True, 
                                     timeout=3600)  # 1 hour timeout
                
                elapsed_time = time.time() - start_time
                
                if result.returncode == 0:
                    print_completion_message("Native DFTB+ (PTBP)", filename, elapsed_time)
                    
                    # Process outputs
                    opt_gen_file = f"{name_without_ext}_opt.gen"
                    if os.path.exists(opt_gen_file):
                        # Convert optimized structure to XYZ
                        opt_atoms = read(opt_gen_file, format='gen')
                        output_xyz = os.path.join(ptbp_outputs_dir, f"{name_without_ext}_relaxed.xyz")
                        write(output_xyz, opt_atoms)
                        print(f"Results saved to {ptbp_outputs_dir}/")
                        print(f"Final structure: {name_without_ext}_relaxed.xyz")
                        
                        # Get final energy from results.tag
                        if os.path.exists('results.tag'):
                            try:
                                with open('results.tag', 'r') as f:
                                    content = f.read()
                                    if 'total_energy' in content:
                                        for line in content.split('\n'):
                                            if 'total_energy' in line and ':real:0:' in line:
                                                energy = float(line.split()[1])
                                                print(f"Final energy: {energy:.6f} Hartree")
                                                break
                            except:
                                pass
                    else:
                        print("⚠️ Warning: Optimized geometry file not found")
                        return False
                        
                else:
                    print_failure_message("Native DFTB+ (PTBP)", filename, elapsed_time)
                    print("STDOUT:", result.stdout[-1000:])  # Last 1000 chars
                    print("STDERR:", result.stderr[-1000:])  # Last 1000 chars
                    return False
                    
            except subprocess.TimeoutExpired:
                print_timeout_message("Native DFTB+ (PTBP)", filename, 3600)
                return False
            except FileNotFoundError:
                print(f"❌ DFTB+ executable not found. Make sure DFTB+ is installed and in PATH")
                return False
            
    finally:
        # Always return to original directory
        os.chdir(original_cwd)
    
    return True


def get_angular_momentum_params_ptbp(atoms):
    """
    Generate MaxAngularMomentum parameters for PTBP parameter set
    
    PTBP supports elements H-Rn with appropriate angular momentum:
    - s-block (H, Li-Cs): s orbital
    - p-block (B-F, Al-Cl, etc): p orbital  
    - d-block (Sc-Zn, Y-Cd, La, Hf-Hg): d orbital
    - f-block (Ce-Lu, Th-Lr): f orbital
    """
    # Define angular momentum mapping for PTBP
    angular_momentum_map = {
        # s-block elements
        'H': 's', 'Li': 's', 'Be': 's', 'Na': 's', 'Mg': 's',
        'K': 's', 'Ca': 's', 'Rb': 's', 'Sr': 's', 'Cs': 's', 'Ba': 's',
        'Fr': 's', 'Ra': 's',
        
        # p-block elements  
        'B': 'p', 'C': 'p', 'N': 'p', 'O': 'p', 'F': 'p',
        'Al': 'p', 'Si': 'p', 'P': 'p', 'S': 'p', 'Cl': 'p',
        'Ga': 'p', 'Ge': 'p', 'As': 'p', 'Se': 'p', 'Br': 'p',
        'In': 'p', 'Sn': 'p', 'Sb': 'p', 'Te': 'p', 'I': 'p',
        'Tl': 'p', 'Pb': 'p', 'Bi': 'p', 'Po': 'p', 'At': 'p', 'Rn': 'p',
        
        # d-block elements (transition metals)
        'Sc': 'd', 'Ti': 'd', 'V': 'd', 'Cr': 'd', 'Mn': 'd', 'Fe': 'd', 'Co': 'd', 'Ni': 'd', 'Cu': 'd', 'Zn': 'd',
        'Y': 'd', 'Zr': 'd', 'Nb': 'd', 'Mo': 'd', 'Tc': 'd', 'Ru': 'd', 'Rh': 'd', 'Pd': 'd', 'Ag': 'd', 'Cd': 'd',
        'La': 'd', 'Hf': 'd', 'Ta': 'd', 'W': 'd', 'Re': 'd', 'Os': 'd', 'Ir': 'd', 'Pt': 'd', 'Au': 'd', 'Hg': 'd',
        'Ac': 'd', 'Rf': 'd', 'Db': 'd', 'Sg': 'd', 'Bh': 'd', 'Hs': 'd', 'Mt': 'd', 'Ds': 'd', 'Rg': 'd', 'Cn': 'd',
        
        # f-block elements (lanthanides and actinides)
        'Ce': 'f', 'Pr': 'f', 'Nd': 'f', 'Pm': 'f', 'Sm': 'f', 'Eu': 'f', 'Gd': 'f', 'Tb': 'f', 
        'Dy': 'f', 'Ho': 'f', 'Er': 'f', 'Tm': 'f', 'Yb': 'f', 'Lu': 'f',
        'Th': 'f', 'Pa': 'f', 'U': 'f', 'Np': 'f', 'Pu': 'f', 'Am': 'f', 'Cm': 'f', 'Bk': 'f',
        'Cf': 'f', 'Es': 'f', 'Fm': 'f', 'Md': 'f', 'No': 'f', 'Lr': 'f',
    }
    
    unique_symbols = set(atoms.get_chemical_symbols())
    angular_params = []
    
    for symbol in sorted(unique_symbols):
        angular_momentum = angular_momentum_map.get(symbol, 'p')  # Default to p if unknown
        angular_params.append(f'        {symbol} = "{angular_momentum}"')
    
    return '\n'.join(angular_params)


def relax_dftb_unified(project_folder, params='auto'):
    """
    Unified DFTB+ relaxation with automatic parameter selection
    
    Args:
        project_folder (str): Name of the project folder or molecule name
        params (str): Parameter set - 'auto', 'ptbp', '3ob'
        
    Returns:
        bool: True if successful, False if failed
    """
    from ..backends.dftb_backend import DFTBBackend
    
    base_path = f"data/{project_folder}"
    
    # If exact folder doesn't exist, try to find a folder containing the molecule name
    if not os.path.exists(base_path):
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
        else:
            raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    if not os.path.exists(base_path):
        raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    # Convert to absolute path 
    base_path = os.path.abspath(base_path)
    
    # Find all XYZ files in base directory and inputs subdirectory
    xyz_files = []
    xyz_pattern = os.path.join(base_path, "*.xyz")
    xyz_files.extend(glob.glob(xyz_pattern))
    
    inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
    xyz_files.extend(glob.glob(inputs_pattern))
    
    # Convert to absolute paths to handle directory changes
    xyz_files = [os.path.abspath(f) for f in xyz_files]
    
    if not xyz_files:
        print(f"No XYZ files found in {base_path}")
        return False
    
    # Create backend instance
    backend = DFTBBackend(params=params)
    
    if not backend.is_available():
        print("❌ DFTB+ backend not available")
        return False
    
    # Create logs, outputs, and temp directories
    logs_dir = os.path.join(base_path, "logs")
    outputs_dir = os.path.join(base_path, "outputs")
    temp_dir = os.path.join(base_path, "temp")
    os.makedirs(logs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    os.makedirs(temp_dir, exist_ok=True)
    
    # Store original working directory
    original_cwd = os.getcwd()
    
    try:
        # Change to temp directory for DFTB+ execution
        os.chdir(temp_dir)
        
        for xyz_file in xyz_files:
            filename = os.path.basename(xyz_file)
            name_without_ext = os.path.splitext(filename)[0]
            
            # Read atoms to determine parameter set
            atoms = read(xyz_file)
            
            # Auto-detect parameter set if needed
            if params == 'auto':
                param_set = backend._detect_parameter_set(atoms)
            else:
                param_set = params
            
            param_path = backend._get_parameter_path(param_set)
            if not os.path.exists(param_path):
                print(f"❌ Parameter set '{param_set}' not found at {param_path}")
                continue
            
            # Create parameter-specific output directory
            param_outputs_dir = os.path.join(base_path, "outputs", f"{param_set}_dftb")
            os.makedirs(param_outputs_dir, exist_ok=True)
            
            print(f"Processing {filename} with DFTB+ ({param_set.upper()} parameters)...")
            
            # Display time estimate
            time_estimate = estimate_time_range(atoms, 'dftb_native')
            print(f"⏱️  Оценка времени оптимизации: {time_estimate}")
            print(f"    Формула: {atoms.get_chemical_formula()} ({len(atoms)} атомов)")
            
            # Write gen file
            gen_file = f"{name_without_ext}.gen"
            write(gen_file, atoms, format='gen')
            
            # Generate angular momentum parameters
            angular_params = backend._get_angular_momentum_params(atoms, param_set)
            
            # Create DFTB+ input file
            hsd_content = f"""
Geometry = GenFormat {{
    <<< "{gen_file}"
}}

Driver = GeometryOptimization {{
    Optimizer = Rational {{}}
    MaxSteps = 500
    OutputPrefix = "{name_without_ext}_opt"
}}

Hamiltonian = DFTB {{
    SCC = Yes
    SCCTolerance = 1e-6
    MaxSCCIterations = 500
    Filling = Fermi {{Temperature[K] = 300.0}}
    HCorrection = Damping {{Exponent = 4.0}}
    Mixer = Anderson {{MixingParameter = 0.05}}
    SlaterKosterFiles = Type2FileNames {{
        Prefix = "{param_path}/"
        Separator = "-"
        Suffix = ".skf"
    }}
    MaxAngularMomentum = {{
{angular_params}
    }}
}}

Options {{
    WriteResultsTag = Yes
}}
"""
            
            # Write DFTB+ input file
            hsd_file = f"dftb_in.hsd"
            with open(hsd_file, 'w') as f:
                f.write(hsd_content)
            
            print(f"Starting DFTB+ optimization for {filename} with {len(atoms)} atoms")
            print(f"Formula: {atoms.get_chemical_formula()}")
            print(f"Parameters: {param_set.upper()}")
            
            start_time = time.time()
            print("Running DFTB+ optimization...")
            
            # Run DFTB+
            try:
                dftb_executable = '/Users/andreypanferov/opt/dftb+/bin/dftb+'
                result = subprocess.run([dftb_executable], 
                                     capture_output=True, 
                                     text=True, 
                                     timeout=3600)  # 1 hour timeout
                
                elapsed_time = time.time() - start_time
                
                if result.returncode == 0:
                    print_completion_message(f"DFTB+ ({param_set.upper()})", filename, elapsed_time)
                    
                    # Process outputs
                    opt_gen_file = f"{name_without_ext}_opt.gen"
                    if os.path.exists(opt_gen_file):
                        # Convert optimized structure to XYZ
                        optimized_atoms = read(opt_gen_file, format='gen')
                        output_xyz = os.path.join(param_outputs_dir, f"{name_without_ext}_relaxed.xyz")
                        write(output_xyz, optimized_atoms)
                        
                        print(f"Results saved to {param_outputs_dir}/")
                        
                        # Run structural analysis
                        run_structural_analysis_if_enabled(
                            xyz_file, output_xyz, project_folder, f"DFTB+ ({param_set.upper()})"
                        )
                    else:
                        print(f"Warning: Optimized structure file not found")
                else:
                    print_failure_message(f"DFTB+ ({param_set.upper()})", filename, elapsed_time)
                    if result.stderr:
                        print(f"DFTB+ error: {result.stderr}")
                
            except subprocess.TimeoutExpired:
                print_timeout_message(f"DFTB+ ({param_set.upper()})", filename, 3600)
            except Exception as e:
                elapsed_time = time.time() - start_time
                print_failure_message(f"DFTB+ ({param_set.upper()})", filename, elapsed_time, str(e))
            
            # Clean up temporary files
            for temp_file in glob.glob("*"):
                if temp_file.endswith(('.gen', '.hsd', '.out', '.xyz')):
                    try:
                        os.remove(temp_file)
                    except:
                        pass
            
            print()
    
    finally:
        os.chdir(original_cwd)
    
    return True


