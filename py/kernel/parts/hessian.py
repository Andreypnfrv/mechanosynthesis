"""
Hessian calculation module for mechanosynthesis research with automatic method matching.

Main functions:
- calculate_hessian_with_method(project_folder, relax_method, ...): Automatically matches Hessian method to relaxation method
- calculate_hessian(project_folder, structure_file=None): DFTB+ with phonopy (legacy)
- calculate_hessian_orca_standalone(...): Orca DFT Hessian
- calculate_hessian_gxtb_standalone(...): g-xTB Hessian via Docker

Method matching:
- dftb/passivated → DFTB+/phonopy
- xtb/gxtb → g-xTB via Docker  
- orca/orca-native → Orca DFT
"""

import os
import glob
import numpy as np
import subprocess
import shutil
import time
from ase.io import read, write
from ase.calculators.dftb import Dftb
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
import yaml
from ..utils.timer import calculate_and_format_time, print_completion_message, print_failure_message, print_timeout_message
from ..utils.performance_tracker import PerformanceTracker


def extract_frequencies_from_dftb_hessian(hessian_file, atoms):
    """
    Extract vibrational frequencies from DFTB+ hessian.out file.
    
    Args:
        hessian_file: Path to hessian.out file
        atoms: ASE Atoms object for mass information
        
    Returns:
        List of frequencies in cm⁻¹
    """
    import numpy as np
    from ase.units import Hartree, Bohr
    
    try:
        # Read hessian matrix from file
        with open(hessian_file, 'r') as f:
            lines = f.readlines()
        
        # Parse the hessian matrix
        # DFTB+ hessian.out format: matrix elements in atomic units
        n_atoms = len(atoms)
        n_dof = 3 * n_atoms
        
        hessian_matrix = np.zeros((n_dof, n_dof))
        
        # Parse the hessian matrix - DFTB+ writes it as continuous data
        # All numeric values are matrix elements in row-major order
        all_values = []
        
        for line in lines:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Split line and try to parse all numbers
            elements = line.split()
            for element in elements:
                try:
                    value = float(element)
                    all_values.append(value)
                except ValueError:
                    continue
        
        # Fill the hessian matrix from the flat array
        if len(all_values) >= n_dof * n_dof:
            hessian_matrix = np.array(all_values[:n_dof*n_dof]).reshape(n_dof, n_dof)
        else:
            print(f"Warning: Expected {n_dof*n_dof} values, got {len(all_values)}")
            return []
        
        # Convert Hessian from Hartree/Bohr² to proper units and mass-weight
        # DFTB+ hessian is in Hartree/Bohr²
        masses = atoms.get_masses()  # in amu
        
        # Mass-weight the Hessian matrix first
        for i in range(n_dof):
            for j in range(n_dof):
                atom_i = i // 3
                atom_j = j // 3
                hessian_matrix[i, j] /= np.sqrt(masses[atom_i] * masses[atom_j])
        
        # Diagonalize mass-weighted Hessian
        eigenvalues, eigenvectors = np.linalg.eigh(hessian_matrix)
        
        # Convert eigenvalues to frequencies in cm⁻¹
        # Formula: ν = (1/2π) * √(eigenvalue * Hartree / (amu * Bohr²)) / c
        # Constants in atomic units to cm⁻¹
        frequencies = []
        
        # Conversion factor from Hartree/(amu*Bohr²) to cm⁻¹
        hartree_to_joule = 4.359744e-18  # J
        amu_to_kg = 1.66054e-27  # kg  
        bohr_to_m = 5.29177e-11  # m
        c_light = 2.998e10  # cm/s
        
        conversion = np.sqrt(hartree_to_joule / (amu_to_kg * (bohr_to_m * 100)**2)) / (2 * np.pi * c_light)
        
        for eigenval in eigenvalues:
            if eigenval > 1e-10:  # Positive eigenvalue
                freq_cm = np.sqrt(eigenval) * conversion
                frequencies.append(freq_cm)
            elif eigenval < -1e-10:  # Negative eigenvalue (imaginary frequency)
                freq_cm = np.sqrt(abs(eigenval)) * conversion
                frequencies.append(-freq_cm)
            else:
                frequencies.append(0.0)  # Nearly zero eigenvalue
        
        return sorted(frequencies)
        
    except Exception as e:
        print(f"Error parsing DFTB+ hessian.out: {e}")
        return []


def save_frequencies_for_thermodynamics(frequencies, project_dir, backend_name='dftb'):
    """
    Save frequencies in standard format for thermodynamic analysis.
    
    Args:
        frequencies: List of frequencies in cm⁻¹
        project_dir: Project directory path
        backend_name: Name of backend used
    """
    import numpy as np
    
    # Create outputs directory structure
    outputs_dir = os.path.join(project_dir, "outputs")
    backend_dir = os.path.join(outputs_dir, backend_name)
    os.makedirs(backend_dir, exist_ok=True)
    
    # Save frequencies in multiple formats for compatibility
    freq_file = os.path.join(backend_dir, "frequencies.dat")
    freq_file_legacy = os.path.join(outputs_dir, "frequencies.dat")
    
    # Save as simple text file
    np.savetxt(freq_file, frequencies, fmt='%.6f', header='Vibrational frequencies (cm-1)')
    np.savetxt(freq_file_legacy, frequencies, fmt='%.6f', header='Vibrational frequencies (cm-1)')
    
    print(f"Frequencies saved for thermodynamic analysis: {freq_file}")
    
    return freq_file


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


def ase_to_phonopy_atoms(ase_atoms):
    """Convert ASE atoms to PhonopyAtoms"""
    return PhonopyAtoms(
        symbols=ase_atoms.get_chemical_symbols(),
        positions=ase_atoms.get_positions(),
        cell=ase_atoms.get_cell()
    )


def analyze_frequencies(frequencies):
    """Analyze frequencies to determine structure stability"""
    if len(frequencies) == 0:
        return "No frequencies found"
    
    negative_freqs = [f for f in frequencies if f < 0]
    small_freqs = [f for f in frequencies if 0 <= f < 50]
    
    analysis = []
    analysis.append(f"Total modes: {len(frequencies)}")
    analysis.append(f"Negative frequencies: {len(negative_freqs)}")
    analysis.append(f"Small frequencies (0-50 cm⁻¹): {len(small_freqs)}")
    
    if len(negative_freqs) == 0:
        if len(small_freqs) <= 6:  # 6 translational/rotational modes expected
            analysis.append("RESULT: True minimum - structure is stable")
        else:
            analysis.append("RESULT: Nearly stable - consider tighter optimization")
    elif len(negative_freqs) == 1:
        analysis.append("RESULT: Transition state (1st order saddle point)")
    else:
        analysis.append(f"RESULT: Higher order saddle point ({len(negative_freqs)} imaginary modes)")
        analysis.append("RECOMMENDATION: Re-optimize with different starting geometry")
    
    # Add frequency ranges
    if len(frequencies) > 0:
        analysis.append(f"Frequency range: {min(frequencies):.1f} to {max(frequencies):.1f} cm⁻¹")
    
    return "\n".join(analysis)


def calculate_orca_hessian(structure_path, method='B3LYP', basis='def2-SVP', nprocs=1):
    """
    Calculate Hessian and vibrational frequencies using Orca
    
    Args:
        structure_path (str): Path to the optimized structure file
        method (str): DFT method (default: B3LYP)
        basis (str): Basis set (default: def2-SVP)  
        nprocs (int): Number of processors (default: 1)
    
    Returns:
        tuple: (frequencies_cm, success)
    """
    filename = os.path.basename(structure_path)
    name_without_ext = os.path.splitext(filename)[0]
    
    # Create hessian directory
    base_dir = os.path.dirname(structure_path)
    hessian_dir = os.path.join(base_dir, f"{name_without_ext}_hessian_orca")
    os.makedirs(hessian_dir, exist_ok=True)
    
    original_cwd = os.getcwd()
    
    try:
        os.chdir(hessian_dir)
        
        # Clean working directory to avoid conflicts with previous calculations
        for ext in ['*.gbw', '*.ges', '*.tmp', '*.int.tmp', '*.loc', '*.qro', '*.uno', '*.unso', '*.uco', '*.hess', '*.cis']:
            for f in glob.glob(ext):
                try:
                    os.remove(f)
                except:
                    pass
        
        # Read atoms
        atoms = read(structure_path)
        
        # Determine number of processors
        import platform
        if platform.system() == 'Darwin':
            actual_nprocs = min(nprocs, 4)  # Cap at 4 cores for reasonable performance
        else:
            actual_nprocs = nprocs
        
        print(f"DEBUG: method={method}, basis={basis}, nprocs={nprocs}, actual_nprocs={actual_nprocs}")
        
        # Display time estimate
        from .relax import estimate_time_range
        time_estimate = estimate_time_range(atoms, 'orca_hessian', method, actual_nprocs)
        print(f"⏱️  Оценка времени расчета гессиана: {time_estimate}")
        print(f"    Формула: {atoms.get_chemical_formula()} ({len(atoms)} атомов)")
        
        # Ask for confirmation (auto-proceed if not interactive)
        try:
            response = input("Продолжить расчет гессиана? (y/n): ").lower().strip()
            if response not in ['y', 'yes', 'да', '']:
                print(f"Пропускаем расчет гессиана для {filename}")
                return [], False
        except (EOFError, KeyboardInterrupt):
            # Auto-proceed if not in interactive mode
            print("Автоматически продолжаем расчет гессиана...")
            pass
        
        # Calculate number of electrons and determine multiplicity
        n_electrons = sum(atoms.get_atomic_numbers()) - 0  # assuming neutral molecule
        multiplicity = 2 if n_electrons % 2 == 1 else 1  # odd electrons = doublet, even = singlet
        
        # Use numerical frequencies for problematic functionals
        freq_method = "NUMFREQ" if method in ['B3LYP', 'PBE0'] else "FREQ"
        
        # Configure parallelization - use shared memory on macOS to avoid MPI issues
        if platform.system() == 'Darwin' and actual_nprocs > 1:
            pal_settings = f"%pal nprocs {actual_nprocs} end"
        else:
            pal_settings = f"%pal nprocs {actual_nprocs} end"
            
        inp_content = f"""! {method} {basis} {freq_method}

{pal_settings}
%maxcore 2000

%scf
    MaxIter 500
    ConvForced true
end"""

        # Add numerical frequency settings if needed
        if freq_method == "NUMFREQ":
            inp_content += """

%freq
    CentralDiff true
    increment 0.005
end"""
        
        inp_content += f"""

* xyzfile 0 {multiplicity} {os.path.abspath(structure_path)}

"""
        
        inp_file = f"{name_without_ext}_freq.inp"
        with open(inp_file, 'w') as f:
            f.write(inp_content)
        
        print(f"Starting Orca frequency calculation for {filename}")
        print(f"Method: {method}/{basis}")
        print(f"Formula: {atoms.get_chemical_formula()}")
        
        # Start timing
        start_time = time.time()
        
        # Run Orca frequency calculation
        orca_path = shutil.which('orca') or '/Applications/orca/orca'
        print(f"Using Orca at: {orca_path}")
        print(f"Using frequency method: {freq_method}")
        
        # Create output file name
        out_file = f"{name_without_ext}_freq.out"
        
        # Longer timeout for numerical frequencies
        timeout = 14400 if freq_method == "NUMFREQ" else 7200  # 4 hours vs 2 hours
        
        # Set environment variables for MPI on macOS
        env = os.environ.copy()
        hostfile_path = None
        if platform.system() == 'Darwin':
            # Create a simple hostfile for localhost only
            hostfile_path = os.path.join(os.getcwd(), 'hostfile')
            with open(hostfile_path, 'w') as hf:
                hf.write('localhost slots=4\n')
            
            # Configure Open MPI for local execution
            env['OMPI_MCA_btl_vader_single_copy_mechanism'] = 'none'
            env['OMPI_MCA_btl'] = 'vader,self'
            env['OMPI_MCA_plm'] = 'rsh'  
            env['OMPI_MCA_plm_rsh_agent'] = 'sh'
            env['OMPI_MCA_orte_default_hostfile'] = hostfile_path
            env['OMPI_MCA_rmaps_base_oversubscribe'] = '1'  # Allow oversubscription
            env['RSH_COMMAND'] = 'sh'  # Set remote shell command
            
            # Add MPI library path - create symlink if needed
            lib_path = '/usr/local/lib/libmpi.40.dylib'
            if not os.path.exists(lib_path):
                try:
                    os.makedirs('/usr/local/lib', exist_ok=True)
                    if not os.path.exists(lib_path):
                        os.symlink('/opt/homebrew/lib/libmpi.40.dylib', lib_path)
                except:
                    pass
            
            # Also try DYLD_LIBRARY_PATH as fallback
            env['DYLD_LIBRARY_PATH'] = '/opt/homebrew/lib:/usr/local/lib'
            env['DYLD_FALLBACK_LIBRARY_PATH'] = '/opt/homebrew/lib:/usr/local/lib'
        
        # Run Orca with output redirection to ensure .out file is created
        with open(out_file, 'w') as outf:
            result = subprocess.run([orca_path, inp_file], 
                                  stdout=outf, stderr=subprocess.PIPE, 
                                  text=True, timeout=timeout, env=env)
        
        # Also capture stderr for debugging
        if result.stderr:
            with open(f"{name_without_ext}_freq.err", 'w') as errf:
                errf.write(result.stderr)
        
        # Calculate elapsed time
        elapsed_time, _ = calculate_and_format_time(start_time)
        
        if result.returncode != 0:
            print_failure_message("Orca frequency calculation", filename, elapsed_time, f"return code: {result.returncode}")
            if result.stderr:
                print(f"STDERR: {result.stderr}")
            # Don't return immediately - try to parse partial results
            print("Attempting to parse partial results...")
        else:
            print_completion_message("Orca frequency calculation", filename, elapsed_time)
        
        # Parse frequencies from output file
        frequencies_cm = []
        analysis_created = False
        
        if os.path.exists(out_file):
            with open(out_file, 'r') as f:
                content = f.read()
            
            print(f"Parsing output file: {out_file}")
            
            # Extract frequencies from output - multiple possible formats
            freq_section = False
            for line in content.split('\n'):
                # Try different frequency section markers
                if any(marker in line for marker in ['VIBRATIONAL FREQUENCIES', 'NORMAL MODES', 'FREQUENCY']):
                    freq_section = True
                    continue
                elif freq_section and line.strip() == '':
                    continue
                elif freq_section and ':' in line:
                    try:
                        parts = line.split()
                        if len(parts) >= 2 and parts[1].replace('-', '').replace('.', '').isdigit():
                            freq = float(parts[1])
                            frequencies_cm.append(freq)
                    except (ValueError, IndexError):
                        continue
                # Alternative frequency format: look for "cm-1" or "cm^-1"
                elif 'cm-1' in line or 'cm^-1' in line or 'wavenumber' in line.lower():
                    parts = line.split()
                    for i, part in enumerate(parts):
                        if part.replace('-', '').replace('.', '').isdigit():
                            try:
                                freq = float(part)
                                if abs(freq) < 10000:  # reasonable frequency range
                                    frequencies_cm.append(freq)
                            except ValueError:
                                continue
            
            # Remove duplicates and sort
            frequencies_cm = sorted(list(set(frequencies_cm)))
            
            # Write analysis file only if frequencies were found and calculation was successful
            if frequencies_cm and result.returncode == 0:
                analysis = analyze_frequencies(frequencies_cm)
                
                with open(f"{name_without_ext}_orca_analysis.txt", "w") as f:
                    f.write(f"Orca {method}/{basis} Frequency Analysis for {filename}\n")
                    f.write("="*60 + "\n\n")
                    f.write(analysis + "\n\n")
                    f.write("All frequencies (cm⁻¹):\n")
                    for i, freq in enumerate(frequencies_cm):
                        f.write(f"{i+1:3d}: {freq:8.2f}\n")
            elif result.returncode != 0:
                # Only write error analysis if calculation failed
                with open(f"{name_without_ext}_orca_error.txt", "w") as f:
                    f.write(f"Orca {method}/{basis} Error Report for {filename}\n")
                    f.write("="*60 + "\n\n")
                    f.write(f"Calculation FAILED with return code: {result.returncode}\n")
                    if result.stderr:
                        f.write(f"Error message: {result.stderr}\n")
                    f.write("\nCheck .out and .err files for detailed error information.\n")
            
            analysis_created = True
            
            if frequencies_cm and result.returncode == 0:
                print(f"Found {len(frequencies_cm)} frequencies")
                print("Analysis:")
                print(analyze_frequencies(frequencies_cm))
                return frequencies_cm, True
            else:
                print(f"Warning: No frequencies found or calculation failed for {filename}")
                return [], False
        else:
            print(f"Error: Orca output file not found for {filename}")
            
            # Create error file if .out is missing
            with open(f"{name_without_ext}_orca_error.txt", "w") as f:
                f.write(f"Orca {method}/{basis} Error Report for {filename}\n")
                f.write("="*60 + "\n\n")
                f.write("ERROR: Output file not created\n")
                f.write(f"Return code: {result.returncode}\n")
                if result.stderr:
                    f.write(f"Error message: {result.stderr}\n")
                f.write("\nPossible issues:\n")
                f.write("- Orca installation problem\n")
                f.write("- Input file format error\n") 
                f.write("- Insufficient memory or disk space\n")
                f.write("- Structure convergence issues\n")
            
            return [], False
            
    except subprocess.TimeoutExpired:
        elapsed_time, _ = calculate_and_format_time(start_time)
        print_timeout_message("Orca frequency calculation", filename, elapsed_time)
        return [], False
    except Exception as e:
        try:
            elapsed_time, _ = calculate_and_format_time(start_time)
            print_failure_message("Orca frequency calculation", filename, elapsed_time, str(e))
        except:
            print(f"Error in Orca frequency calculation for {filename}: {e}")
        return [], False
    finally:
        # Clean up hostfile if created
        if hostfile_path and os.path.exists(hostfile_path):
            try:
                os.remove(hostfile_path)
            except:
                pass
        os.chdir(original_cwd)


def calculate_gxtb_hessian(structure_path):
    """
    Calculate Hessian and vibrational frequencies using g-xTB Docker container
    
    Args:
        structure_path (str): Path to the optimized structure file
    
    Returns:
        tuple: (frequencies_cm, success)
    """
    filename = os.path.basename(structure_path)
    name_without_ext = os.path.splitext(filename)[0]
    
    # Create hessian directory
    base_dir = os.path.dirname(structure_path)
    hessian_dir = os.path.join(base_dir, f"{name_without_ext}_hessian_gxtb")
    os.makedirs(hessian_dir, exist_ok=True)
    
    original_cwd = os.getcwd()
    
    try:
        os.chdir(hessian_dir)
        
        # Read atoms for info
        atoms = read(structure_path)
        
        print(f"Starting g-xTB Hessian calculation for {filename}")
        print(f"Formula: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")
        
        # Display time estimate and warning
        from .relax import estimate_time_range
        time_estimate = estimate_time_range(atoms, 'gxtb_hessian')
        print(f"⏱️  Оценка времени расчета гессиана: {time_estimate}")
        print(f"⚠️  WARNING: g-xTB uses NUMERICAL Hessian (3N+1 = {3*len(atoms)+1} evaluations)")
        print(f"    For large molecules, consider using --hessian-orca for analytical Hessian")
        
        # Ask for confirmation (auto-proceed if not interactive)
        try:
            response = input("Продолжить расчет гессиана g-xTB? (y/n): ").lower().strip()
            if response not in ['y', 'yes', 'да', '']:
                print(f"Пропускаем расчет гессиана для {filename}")
                return [], False
        except (EOFError, KeyboardInterrupt):
            # Auto-proceed if not in interactive mode
            print("Автоматически продолжаем расчет гессиана g-xTB...")
            pass
        
        # Start timing
        start_time = time.time()
        
        # Run g-xTB Hessian calculation using Docker
        script_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'run_gxtb.sh')
        
        # Run the command and capture output (longer timeout for Hessian)
        result = subprocess.run([script_path, structure_path, '--hessian', '--output-dir', hessian_dir], 
                              capture_output=True, text=True, timeout=7200)  # 2 hour timeout
        
        # Calculate elapsed time
        elapsed_time, _ = calculate_and_format_time(start_time)
        
        if result.returncode != 0:
            print_failure_message("g-xTB Hessian calculation", filename, elapsed_time, f"return code: {result.returncode}")
            if result.stderr:
                print(f"STDERR: {result.stderr}")
            return [], False
        else:
            print_completion_message("g-xTB Hessian calculation", filename, elapsed_time)
        
        # Parse frequencies from g-xTB output files
        frequencies_cm = []
        
        # Look for frequency output files
        freq_files = glob.glob(os.path.join(hessian_dir, "*.frequencies"))
        if not freq_files:
            freq_files = glob.glob(os.path.join(hessian_dir, "vibspectrum"))
            
        for freq_file in freq_files:
            if os.path.exists(freq_file):
                with open(freq_file, 'r') as f:
                    content = f.read()
                
                # Parse g-xTB frequency format
                for line in content.split('\n'):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                # First column is usually mode number, second is frequency
                                freq = float(parts[1])
                                if abs(freq) < 10000:  # reasonable frequency range
                                    frequencies_cm.append(freq)
                            except (ValueError, IndexError):
                                continue
        
        # Also check standard output for frequencies
        if not frequencies_cm and result.stdout:
            lines = result.stdout.split('\n')
            freq_section = False
            for line in lines:
                if 'frequencies' in line.lower() or 'vibration' in line.lower():
                    freq_section = True
                    continue
                elif freq_section and line.strip():
                    parts = line.split()
                    for part in parts:
                        try:
                            freq = float(part)
                            if abs(freq) < 10000:  # reasonable frequency range
                                frequencies_cm.append(freq)
                        except ValueError:
                            continue
        
        # Remove duplicates and sort
        frequencies_cm = sorted(list(set(frequencies_cm)))
        
        # Write analysis file
        if frequencies_cm:
            analysis = analyze_frequencies(frequencies_cm)
            
            with open(f"{name_without_ext}_gxtb_analysis.txt", "w") as f:
                f.write(f"g-xTB Hessian Analysis for {filename}\n")
                f.write("="*60 + "\n\n")
                f.write(analysis + "\n\n")
                f.write("All frequencies (cm⁻¹):\n")
                for i, freq in enumerate(frequencies_cm):
                    f.write(f"{i+1:3d}: {freq:8.2f}\n")
            
            print(f"Found {len(frequencies_cm)} frequencies")
            print("Analysis:")
            print(analyze_frequencies(frequencies_cm))
            return frequencies_cm, True
        else:
            print(f"Warning: No frequencies found for {filename}")
            
            # Create error file
            with open(f"{name_without_ext}_gxtb_error.txt", "w") as f:
                f.write(f"g-xTB Hessian Error Report for {filename}\n")
                f.write("="*60 + "\n\n")
                f.write("ERROR: No frequencies found in output\n")
                if result.stdout:
                    f.write(f"STDOUT:\n{result.stdout}\n\n")
                if result.stderr:
                    f.write(f"STDERR:\n{result.stderr}\n\n")
            
            return [], False
            
    except subprocess.TimeoutExpired:
        elapsed_time, _ = calculate_and_format_time(start_time)
        print_timeout_message("g-xTB Hessian calculation", filename, elapsed_time)
        return [], False
    except Exception as e:
        try:
            elapsed_time, _ = calculate_and_format_time(start_time)
            print_failure_message("g-xTB Hessian calculation", filename, elapsed_time, str(e))
        except:
            print(f"Error in g-xTB Hessian calculation for {filename}: {e}")
        return [], False
    finally:
        os.chdir(original_cwd)


def calculate_hessian_gxtb_standalone(project_folder):
    """
    Calculate Hessian and vibrational frequencies using g-xTB for all relaxed structures in project
    
    Args:
        project_folder (str): Name of the project folder
    """
    base_path = f"data/{project_folder}"
    
    # Find project folder (same logic as other functions)
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
    
    outputs_dir = os.path.join(base_path, "outputs")
    if not os.path.exists(outputs_dir):
        raise FileNotFoundError(f"Outputs directory {outputs_dir} does not exist")
    
    # Find all relaxed structure files to analyze
    structure_files = []
    
    # Look in all output subdirectories for relaxed structures
    for subdir in ['ase_bfgs', 'staged_ase', 'native_dftb', 'orca_native', 'orca_dft', 'gxtb']:
        subdir_path = os.path.join(outputs_dir, subdir)
        if os.path.exists(subdir_path):
            pattern = os.path.join(subdir_path, "*_relaxed.xyz")
            structure_files.extend(glob.glob(pattern))
    
    # Also check main outputs directory
    pattern = os.path.join(outputs_dir, "*_relaxed.xyz")
    structure_files.extend(glob.glob(pattern))
    
    if not structure_files:
        print(f"No relaxed structure files found in {outputs_dir}")
        print("Looking for files matching pattern: *_relaxed.xyz")
        return
    
    # Convert to absolute paths and remove duplicates
    structure_files = list(set([os.path.abspath(f) for f in structure_files]))
    
    print(f"Found {len(structure_files)} relaxed structures for g-xTB Hessian calculation")
    print("Using g-xTB semiempirical method via Docker")
    print("="*60)
    
    # Calculate Hessian for each structure
    for structure_path in structure_files:
        frequencies, success = calculate_gxtb_hessian(structure_path)
        
        if success:
            filename = os.path.basename(structure_path)
            print(f"g-xTB Hessian calculation successful for {filename}")
        else:
            filename = os.path.basename(structure_path) 
            print(f"g-xTB Hessian calculation failed for {filename}")
        print()
    
    print("="*60)
    print("All g-xTB Hessian calculations completed!")
    print("="*60)


def calculate_hessian_xtb_standalone(project_folder, xtb_method='GFN2-xTB', backend='fixed'):
    """
    Calculate Hessian and vibrational frequencies using xTB for all relaxed structures in project
    
    Args:
        project_folder (str): Name of the project folder
        xtb_method (str): xTB method (GFN2-xTB, GFN1-xTB, GFN0-xTB) (default: GFN2-xTB)
        backend (str): xTB backend - 'native' (system xtb), 'fixed' (patched version), 'gxtb' (Docker) (default: fixed)
    """
    base_path = f"data/{project_folder}"
    
    # Find project folder (same logic as other functions)
    if not os.path.exists(base_path):
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            raise FileNotFoundError(f"Multiple project folders found for '{project_folder}': {[os.path.basename(f) for f in matching_folders]}")
        else:
            raise FileNotFoundError(f"Project folder {base_path} does not exist")
    
    outputs_dir = os.path.join(base_path, "outputs")
    if not os.path.exists(outputs_dir):
        raise FileNotFoundError(f"Outputs directory {outputs_dir} does not exist")
    
    # Find all relaxed structure files to analyze
    structure_files = []
    
    # Look for xTB relaxed structures first
    xtb_subdir = os.path.join(outputs_dir, "xtb")
    if os.path.exists(xtb_subdir):
        pattern = os.path.join(xtb_subdir, "*_relaxed.xyz")
        structure_files.extend(glob.glob(pattern))
    
    # Also look in other output subdirectories for relaxed structures
    for subdir in ['ase_bfgs', 'staged_ase', 'native_dftb', 'orca_native', 'orca_dft', 'gxtb']:
        subdir_path = os.path.join(outputs_dir, subdir)
        if os.path.exists(subdir_path):
            pattern = os.path.join(subdir_path, "*_relaxed.xyz")
            structure_files.extend(glob.glob(pattern))
    
    # Also check main outputs directory
    pattern = os.path.join(outputs_dir, "*_relaxed.xyz")
    structure_files.extend(glob.glob(pattern))
    
    if not structure_files:
        print(f"No relaxed structure files found in {outputs_dir}")
        print("Looking for files matching pattern: *_relaxed.xyz")
        return
    
    # Convert to absolute paths and remove duplicates
    structure_files = list(set([os.path.abspath(f) for f in structure_files]))
    
    print(f"Found {len(structure_files)} relaxed structures for xTB Hessian calculation")
    print(f"Using xTB method: {xtb_method} (backend: {backend})")
    print("="*60)
    
    # Create hessian output directory
    hessian_dir = os.path.join(base_path, "outputs", "hessian_xtb")
    os.makedirs(hessian_dir, exist_ok=True)
    
    # Calculate Hessian for each structure
    for structure_path in structure_files:
        filename = os.path.basename(structure_path)
        name_without_ext = os.path.splitext(filename)[0]
        
        # Read atoms for tracking
        atoms = read(structure_path)
        
        # Track performance of this calculation
        method_params = {
            'method': xtb_method,
            'backend': backend,
            'calculation_type': 'hessian'
        }
        
        with PerformanceTracker(
            method=f"xtb-hessian-{backend}",
            molecule_name=name_without_ext,
            method_params=method_params,
            n_cores=1  # xTB hessian typically single-threaded
        ) as tracker:
            tracker.set_atoms(atoms)
            
            frequencies, success = calculate_xtb_hessian_single(structure_path, xtb_method, backend, hessian_dir)
            
            if success:
                tracker.set_quality(4, "Successful Hessian calculation")
                print(f"xTB Hessian calculation successful for {filename}")
            else:
                tracker.set_quality(1, "Hessian calculation failed")
                print(f"xTB Hessian calculation failed for {filename}")
        print()
    
    print("="*60)
    print("All xTB Hessian calculations completed!")
    print("="*60)


def calculate_hessian_orca_standalone(project_folder, method='B3LYP', basis='def2-SVP', nprocs=4):
    """
    Calculate Hessian and vibrational frequencies using Orca for all relaxed structures in project
    
    Args:
        project_folder (str): Name of the project folder
        method (str): DFT method (default: B3LYP)
        basis (str): Basis set (default: def2-SVP)  
        nprocs (int): Number of processors (default: 4)
    """
    base_path = f"data/{project_folder}"
    
    # Find project folder (same logic as other functions)
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
    
    outputs_dir = os.path.join(base_path, "outputs")
    if not os.path.exists(outputs_dir):
        raise FileNotFoundError(f"Outputs directory {outputs_dir} does not exist")
    
    # Find all relaxed structure files to analyze
    structure_files = []
    
    # Look in all output subdirectories for relaxed structures
    for subdir in ['ase_bfgs', 'staged_ase', 'native_dftb', 'orca_native', 'orca_dft']:
        subdir_path = os.path.join(outputs_dir, subdir)
        if os.path.exists(subdir_path):
            pattern = os.path.join(subdir_path, "*_relaxed.xyz")
            structure_files.extend(glob.glob(pattern))
    
    # Also check main outputs directory
    pattern = os.path.join(outputs_dir, "*_relaxed.xyz")
    structure_files.extend(glob.glob(pattern))
    
    if not structure_files:
        print(f"No relaxed structure files found in {outputs_dir}")
        print("Looking for files matching pattern: *_relaxed.xyz")
        return
    
    # Convert to absolute paths and remove duplicates
    structure_files = list(set([os.path.abspath(f) for f in structure_files]))
    
    # Use multiple processors on macOS
    import platform
    if platform.system() == 'Darwin':
        actual_nprocs = min(nprocs, 4)  # Cap at 4 cores for reasonable performance
    else:
        actual_nprocs = nprocs
    
    print(f"DEBUG: method={method}, basis={basis}, nprocs={nprocs}, actual_nprocs={actual_nprocs}")
    
    print(f"Found {len(structure_files)} relaxed structures for Orca Hessian calculation")
    print(f"Using Orca {method}/{basis} with {actual_nprocs} processors")
    print("="*60)
    
    # Calculate Hessian for each structure
    for structure_path in structure_files:
        frequencies, success = calculate_orca_hessian(structure_path, method, basis, actual_nprocs)
        
        if success:
            filename = os.path.basename(structure_path)
            print(f"Orca Hessian calculation successful for {filename}")
        else:
            filename = os.path.basename(structure_path) 
            print(f"Orca Hessian calculation failed for {filename}")
        print()
    
    print("="*60)
    print("All Orca Hessian calculations completed!")
    print("="*60)


def calculate_hessian(project_folder, structure_file=None, displacement=0.01):
    """
    Calculate Hessian and vibrational frequencies using phonopy+DFTB+
    
    Args:
        project_folder (str): Name of the project folder
        structure_file (str): Specific structure file to analyze (optional)
        displacement (float): Displacement for finite difference (Angstrom)
    """
    base_path = f"data/{project_folder}"
    
    # Find project folder (same logic as relax.py)
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
    
    outputs_dir = os.path.join(base_path, "outputs")
    if not os.path.exists(outputs_dir):
        raise FileNotFoundError(f"Outputs directory {outputs_dir} does not exist")
    
    # Find structure files to analyze
    structure_files = []
    if structure_file:
        # Use specific file
        full_path = os.path.join(outputs_dir, structure_file)
        if os.path.exists(full_path):
            structure_files.append(full_path)
        else:
            raise FileNotFoundError(f"Structure file {full_path} does not exist")
    else:
        # Find all relaxed structures
        pattern = os.path.join(outputs_dir, "*_relaxed.xyz")
        structure_files = glob.glob(pattern)
    
    if not structure_files:
        print(f"No structure files found in {outputs_dir}")
        return
    
    # Convert to absolute paths
    structure_files = [os.path.abspath(f) for f in structure_files]
    
    # Store original working directory
    original_cwd = os.getcwd()
    
    try:
        for structure_path in structure_files:
            filename = os.path.basename(structure_path)
            name_without_ext = os.path.splitext(filename)[0]
            
            # Create hessian directory for this structure
            hessian_dir = os.path.join(outputs_dir, f"{name_without_ext}_hessian")
            os.makedirs(hessian_dir, exist_ok=True)
            
            print(f"Calculating Hessian for {filename} using phonopy+DFTB+...")
            
            # Start timing
            start_time = time.time()
            
            # Change to hessian directory for phonopy calculation
            os.chdir(hessian_dir)
            
            # Read atoms
            atoms = read(structure_path)
            print(f"Starting phonopy calculation for {filename} with {len(atoms)} atoms")
            print(f"Formula: {atoms.get_chemical_formula()}")
            
            # Ensure periodic cell for phonopy (needed even for molecules)
            if np.allclose(atoms.get_cell(), 0):
                # Add large box for isolated molecule
                center = atoms.get_center_of_mass()
                atoms.center(vacuum=15.0)  # 15 Angstrom vacuum
                atoms.set_pbc([True, True, True])  # Enable periodicity
                print(f"Added periodic box with 15 Å vacuum")
            
            # Convert to phonopy atoms
            phonopy_atoms = ase_to_phonopy_atoms(atoms)
            
            # Create phonopy object for 1x1x1 supercell (gamma point only)
            phonopy = Phonopy(phonopy_atoms, [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            
            # Generate displacements
            phonopy.generate_displacements(distance=displacement)
            
            print(f"Generated {len(phonopy.supercells_with_displacements)} displaced structures")
            
            # Set up DFTB+ calculator
            angular_params = get_angular_momentum_params(atoms)
            
            calc_params = {
                'label': 'phonopy_calc',
                'kpts': (1, 1, 1),
                'Hamiltonian_SCC': 'Yes',
                'Hamiltonian_SCCTolerance': 1e-8,
                'Hamiltonian_MaxSCCIterations': 500,
                'Hamiltonian_SlaterKosterFiles_Prefix': os.environ.get('SKDIR', '/Users/andreypanferov/opt/dftb+/slakos') + '/mio-1-1/',
                'Hamiltonian_SlaterKosterFiles_Separator': '-',
                'Hamiltonian_SlaterKosterFiles_Suffix': '.skf',
                'Hamiltonian_MaxAngularMomentum_': '',
                'Hamiltonian_HCorrection': 'Damping {Exponent = 4.0}',
                'Hamiltonian_Filling': 'Fermi {Temperature[K] = 300.0}',
                'Hamiltonian_Mixer': 'Anderson {MixingParameter = 0.05}',
                'Options_WriteResultsTag': 'Yes',
                'Analysis_CalculateForces': 'Yes'
            }
            
            # Add angular momentum parameters
            calc_params.update(angular_params)
            
            # Calculate forces for displaced structures
            forces = []
            
            for i, supercell in enumerate(phonopy.supercells_with_displacements):
                print(f"Calculating forces for displacement {i+1}/{len(phonopy.supercells_with_displacements)}")
                
                # Convert phonopy supercell to ASE atoms
                displaced_atoms = read('displaced_structure.xyz') if os.path.exists('displaced_structure.xyz') else atoms.copy()
                displaced_atoms.set_positions(supercell.positions)
                displaced_atoms.set_cell(supercell.cell)
                displaced_atoms.set_chemical_symbols(supercell.symbols)
                
                # Set calculator and calculate forces
                calc = Dftb(**calc_params)
                displaced_atoms.calc = calc
                
                try:
                    force = displaced_atoms.get_forces()
                    forces.append(force)
                    print(f"  Max force: {np.max(np.linalg.norm(force, axis=1)):.4f} eV/Å")
                except Exception as e:
                    print(f"  Error calculating forces for displacement {i+1}: {e}")
                    forces.append(np.zeros((len(displaced_atoms), 3)))
            
            # Set calculated forces to phonopy
            phonopy.forces = forces
            
            # Calculate dynamical matrix and frequencies
            phonopy.produce_force_constants()
            
            # Get frequencies at Gamma point
            phonopy.run_mesh([1, 1, 1])
            frequencies = phonopy.get_frequencies([0, 0, 0])  # Gamma point
            
            # Convert from THz to cm^-1
            frequencies_cm = frequencies * 33.356  # THz to cm^-1 conversion
            
            # Write results
            if len(frequencies_cm) > 0:
                analysis = analyze_frequencies(frequencies_cm)
                
                with open(f"{name_without_ext}_analysis.txt", "w") as f:
                    f.write(f"Phonopy+DFTB+ Hessian Analysis for {filename}\n")
                    f.write("="*50 + "\n\n")
                    f.write(f"Displacement: {displacement} Å\n")
                    f.write(f"Number of displaced structures: {len(phonopy.supercells_with_displacements)}\n\n")
                    f.write(analysis + "\n\n")
                    f.write("All frequencies (cm⁻¹):\n")
                    for i, freq in enumerate(frequencies_cm):
                        f.write(f"{i+1:3d}: {freq:8.2f}\n")
                
                # Save phonopy results
                phonopy.save("phonopy_disp.yaml")
                
                # Calculate elapsed time
                elapsed_time, _ = calculate_and_format_time(start_time)
                
                print_completion_message("Phonopy Hessian calculation", filename, elapsed_time)
                print("Analysis:")
                print(analysis)
                print()
            else:
                # Calculate elapsed time for failed calculation too
                elapsed_time, _ = calculate_and_format_time(start_time)
                print(f"Warning: No frequencies calculated for {filename}")
                print_failure_message("Phonopy Hessian calculation", filename, elapsed_time, "No frequencies calculated")
                
    finally:
        # Always return to original working directory
        os.chdir(original_cwd)


def calculate_hessian_with_method(project_folder, relax_method, orca_method='B3LYP', orca_basis='def2-SVP', nprocs=4):
    """
    Calculate Hessian using method that matches the relaxation method used.
    
    Args:
        project_folder (str): Name of the project folder
        relax_method (str): Relaxation method used: dftb|xtb|gxtb|orca|orca-native|passivated
        orca_method (str): DFT method for Orca calculations (default: B3LYP)
        orca_basis (str): Basis set for Orca calculations (default: def2-SVP)
        nprocs (int): Number of processors for parallel calculations
    """
    
    # Method mapping from relaxation method to Hessian method
    method_mapping = {
        'dftb': 'dftb',
        'passivated': 'dftb',  # DFTB+ relaxed structures
        'xtb': 'xtb',          # xTB → xTB Hessian (native)
        'gxtb': 'gxtb',        # g-xTB → g-xTB 
        'orca': 'orca',        # Orca DFT → Orca DFT
        'orca-native': 'orca'  # Native Orca → Orca DFT
    }
    
    if relax_method not in method_mapping:
        raise ValueError(f"Unknown relaxation method '{relax_method}'. "
                        f"Supported methods: {list(method_mapping.keys())}")
    
    hessian_method = method_mapping[relax_method]
    
    print(f"Using {hessian_method.upper()} Hessian calculation for {relax_method} relaxed structures")
    print(f"Project: {project_folder}")
    print()
    
    # Call appropriate Hessian calculation function
    if hessian_method == 'dftb':
        print("Calculating Hessian with DFTB+/phonopy...")
        calculate_hessian(project_folder)
        
    elif hessian_method == 'xtb':
        print("Calculating Hessian with xTB semiempirical method...")
        calculate_hessian_xtb_standalone(project_folder, 'GFN2-xTB', 'fixed')
        
    elif hessian_method == 'gxtb':
        print("Calculating Hessian with g-xTB via Docker...")
        print("⚠️  Warning: g-xTB Hessian is VERY SLOW (numerical derivatives with 3N+1 evaluations)")
        calculate_hessian_gxtb_standalone(project_folder)
        
    elif hessian_method == 'orca':
        print(f"Calculating Hessian with Orca DFT ({orca_method}/{orca_basis})...")
        calculate_hessian_orca_standalone(project_folder, orca_method, orca_basis, nprocs)
        
    else:
        raise RuntimeError(f"Internal error: unknown Hessian method '{hessian_method}'")


def calculate_xtb_hessian_single(structure_path, xtb_method='GFN2-xTB', backend='fixed', output_dir=None):
    """
    Calculate Hessian and vibrational frequencies using xTB for a single structure
    
    Args:
        structure_path (str): Path to the optimized structure file
        xtb_method (str): xTB method (GFN2-xTB, GFN1-xTB, GFN0-xTB)
        backend (str): xTB backend ('native', 'fixed', 'gxtb')
        output_dir (str): Directory to save results (optional)
    
    Returns:
        tuple: (frequencies_cm, success)
    """
    filename = os.path.basename(structure_path)
    name_without_ext = os.path.splitext(filename)[0].replace('_relaxed', '')
    
    # Create hessian directory
    if output_dir is None:
        base_dir = os.path.dirname(structure_path)
        hessian_dir = os.path.join(base_dir, f"{name_without_ext}_hessian_xtb")
    else:
        hessian_dir = output_dir
        
    os.makedirs(hessian_dir, exist_ok=True)
    
    original_cwd = os.getcwd()
    
    try:
        os.chdir(hessian_dir)
        
        # Read atoms for info
        atoms = read(structure_path)
        
        print(f"Starting xTB Hessian calculation for {filename}")
        print(f"Formula: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")
        print(f"Method: {xtb_method} (backend: {backend})")
        
        # Copy structure to working directory
        local_xyz = f"{name_without_ext}.xyz"
        shutil.copy2(structure_path, local_xyz)
        
        # Determine xTB executable based on backend
        if backend == 'fixed':
            xtb_exe = os.path.join(os.getcwd(), "../../../3rdparty/xtb/fork/build_cmake/xtb")
            if not os.path.exists(xtb_exe):
                xtb_exe = "xtb"  # fallback to system xtb
        elif backend == 'gxtb':
            # Use g-xTB via Docker - redirect to g-xTB function
            os.chdir(original_cwd)
            return calculate_gxtb_hessian(structure_path)
        else:  # native
            xtb_exe = "xtb"
            
        # Build xTB command for Hessian calculation
        xtb_cmd = [xtb_exe, local_xyz, "--hess"]
        
        # Add method parameter  
        if xtb_method == 'GFN1-xTB':
            xtb_cmd.append("--gfn1")
        elif xtb_method == 'GFN0-xTB':
            xtb_cmd.append("--gfn0")
        # GFN2-xTB is default, no flag needed
        
        print(f"Running: {' '.join(xtb_cmd)}")
        print("This may take several minutes for larger molecules...")
        
        # Run xTB Hessian calculation
        start_time = time.time()
        
        try:
            result = subprocess.run(xtb_cmd, capture_output=True, text=True, timeout=3600)  # 1 hour timeout
            
            if result.returncode != 0:
                print(f"xTB Hessian failed with return code {result.returncode}")
                if result.stderr:
                    print(f"Error: {result.stderr}")
                return None, False
                
            # Parse output for frequencies
            output_lines = result.stdout.split('\n')
            frequencies = []
            
            # Look for vibrational frequencies in output
            freq_section = False
            for line in output_lines:
                if 'vibrational frequencies' in line.lower():
                    freq_section = True
                    continue
                    
                if freq_section and line.strip():
                    # Try to parse frequency values
                    parts = line.strip().split()
                    for part in parts:
                        try:
                            freq = float(part)
                            if freq > 0:  # Only positive frequencies
                                frequencies.append(freq)
                        except ValueError:
                            continue
                            
                if freq_section and not line.strip():
                    break
            
            elapsed = time.time() - start_time
            
            print(f"xTB Hessian completed in {elapsed:.1f} seconds")
            
            if frequencies:
                print(f"Found {len(frequencies)} vibrational frequencies")
                print(f"Lowest frequency: {min(frequencies):.1f} cm⁻¹")
                print(f"Highest frequency: {max(frequencies):.1f} cm⁻¹")
                
                # Save frequencies to file
                freq_file = f"{name_without_ext}_frequencies.txt"
                with open(freq_file, 'w') as f:
                    f.write(f"xTB {xtb_method} Vibrational Frequencies\n")
                    f.write(f"Structure: {filename}\n")
                    f.write(f"Formula: {atoms.get_chemical_formula()} ({len(atoms)} atoms)\n")
                    f.write(f"Calculation time: {elapsed:.1f} seconds\n")
                    f.write("\nFrequencies (cm⁻¹):\n")
                    for i, freq in enumerate(frequencies, 1):
                        f.write(f"{i:3d}  {freq:8.1f}\n")
                        
                print(f"Frequencies saved to: {freq_file}")
                
                # Save frequencies for thermodynamic analysis
                try:
                    project_path = os.path.dirname(os.path.dirname(hessian_dir))  # Get project directory
                    # Save directly to outputs directory in standard format
                    import numpy as np
                    outputs_dir = os.path.join(project_path, "outputs")
                    os.makedirs(outputs_dir, exist_ok=True)
                    freq_file_thermo = os.path.join(outputs_dir, "frequencies.dat")
                    np.savetxt(freq_file_thermo, frequencies, fmt='%.6f', header='Vibrational frequencies (cm-1)')
                    print(f"📊 Frequencies saved for thermodynamic analysis: {freq_file_thermo}")
                except Exception as e:
                    print(f"⚠️  Warning: Could not save frequencies for thermodynamics: {e}")
                
            # Save full output
            with open(f"{name_without_ext}_xtb_hess.out", 'w') as f:
                f.write(result.stdout)
            
            return frequencies, True
            
        except subprocess.TimeoutExpired:
            print(f"xTB Hessian calculation timed out after 1 hour")
            return None, False
            
        except Exception as e:
            print(f"Error during xTB Hessian calculation: {e}")
            return None, False
            
    finally:
        os.chdir(original_cwd)
        
    return None, False


def calculate_hessian_dftb_native(project_folder):
    """
    Calculate Hessian and vibrational frequencies using native DFTB+ analytical Hessian
    
    This function uses DFTB+ built-in Hessian calculation instead of phonopy finite differences.
    Much faster and more accurate than phonopy approach.
    
    Args:
        project_folder (str): Name of the project folder
        
    Returns:
        bool: True if successful, False if failed
    """
    import subprocess
    import os
    import glob
    import time
    from ase.io import read, write
    
    # Import angular momentum function from relax.py
    from .relax import get_angular_momentum_params_ptbp
    
    base_path = f"data/{project_folder}"
    
    # Find project folder (same logic as relax.py)
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
    
    outputs_dir = os.path.join(base_path, "outputs")
    if not os.path.exists(outputs_dir):
        raise FileNotFoundError(f"Outputs directory {outputs_dir} does not exist")
    
    # Find relaxed structures to analyze
    structure_files = []
    
    # Look for DFTB+ PTBP relaxed structures first (preferred)
    ptbp_pattern = os.path.join(outputs_dir, "ptbp_dftb", "*_relaxed.xyz")
    ptbp_files = glob.glob(ptbp_pattern)
    if ptbp_files:
        structure_files.extend(ptbp_files)
    
    # Also look in root outputs directory
    root_pattern = os.path.join(outputs_dir, "*_relaxed.xyz") 
    root_files = glob.glob(root_pattern)
    structure_files.extend(root_files)
    
    if not structure_files:
        print(f"No relaxed structures found in {outputs_dir}")
        print("Please run DFTB+ relaxation first")
        return False
    
    # Remove duplicates
    structure_files = list(set(structure_files))
    
    # Check if PTBP parameters are available
    skdir = os.environ.get('SKDIR', '/Users/andreypanferov/opt/dftb+/slakos/')
    ptbp_dir = os.path.join(skdir, 'ptbp-2024')
    
    if not os.path.exists(ptbp_dir):
        raise FileNotFoundError(f"PTBP parameter directory not found: {ptbp_dir}")
    
    print(f"Found {len(structure_files)} structures for DFTB+ native Hessian calculation")
    print(f"Using PTBP parameters from: {ptbp_dir}")
    
    success_count = 0
    
    for structure_path in structure_files:
        filename = os.path.basename(structure_path)
        name_without_ext = os.path.splitext(filename)[0]
        
        print(f"\nCalculating native DFTB+ Hessian for {filename}...")
        
        # Create hessian calculation directory
        hessian_dir = os.path.join(outputs_dir, f"{name_without_ext}_hessian_dftb_native")
        os.makedirs(hessian_dir, exist_ok=True)
        
        # Store original working directory and absolute paths
        original_cwd = os.getcwd()
        abs_structure_path = os.path.abspath(structure_path)
        abs_hessian_dir = os.path.abspath(hessian_dir)
        
        try:
            # Change to hessian calculation directory
            os.chdir(hessian_dir)
            
            # Read structure using absolute path
            atoms = read(abs_structure_path)
            
            print(f"Starting DFTB+ native Hessian for {filename} with {len(atoms)} atoms")
            print(f"Formula: {atoms.get_chemical_formula()}")
            
            # Convert to simple XYZ format for DFTB+ compatibility
            temp_xyz = "structure.xyz"
            with open(temp_xyz, 'w') as f:
                f.write(f"{len(atoms)}\n")
                f.write(f"{atoms.get_chemical_formula()}\n")
                for atom in atoms:
                    f.write(f"{atom.symbol:2s} {atom.position[0]:12.8f} {atom.position[1]:12.8f} {atom.position[2]:12.8f}\n")
            
            print(f"📝 Structure converted to simple XYZ format: {temp_xyz}")
            
            # Write structure in gen format
            gen_file = f"{name_without_ext}.gen"
            write(gen_file, atoms, format='gen')
            
            # Generate angular momentum parameters for PTBP
            angular_params = get_angular_momentum_params_ptbp(atoms)
            
            # Create DFTB+ input file for Hessian calculation
            hsd_content = f"""
Geometry = GenFormat {{
    <<< "{gen_file}"
}}

Driver = SecondDerivatives {{}}

Hamiltonian = DFTB {{
    SCC = Yes
    SCCTolerance = 1e-8
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


Options = {{
    WriteResultsTag = Yes
}}
"""
            
            # Write DFTB+ input file
            hsd_file = "dftb_in.hsd"
            with open(hsd_file, 'w') as f:
                f.write(hsd_content)
            
            print("Running DFTB+ native Hessian calculation...")
            start_time = time.time()
            
            # Run DFTB+
            dftb_executable = '/Users/andreypanferov/opt/dftb+/bin/dftb+'
            result = subprocess.run([dftb_executable], 
                                 capture_output=True, 
                                 text=True, 
                                 timeout=3600)  # 1 hour timeout
            
            elapsed_time = time.time() - start_time
            
            if result.returncode == 0:
                print(f"✅ Native DFTB+ Hessian completed for {filename} за {elapsed_time:.0f}с")
                
                # Check for output files
                vibrations_file = "vibrations.tag"
                hessian_file = "hessian.out"
                
                if os.path.exists(vibrations_file):
                    print(f"Vibrational frequencies saved to: {hessian_dir}/vibrations.tag")
                    
                    # Parse and display frequencies
                    try:
                        with open(vibrations_file, 'r') as f:
                            content = f.read()
                            
                        # Extract frequencies
                        frequencies = []
                        in_freq_section = False
                        
                        for line in content.split('\n'):
                            line = line.strip()
                            if 'vibration_frequencies' in line:
                                in_freq_section = True
                                continue
                            elif in_freq_section and line:
                                if line.startswith('vibration_frequencies'):
                                    continue
                                elif ':real:' in line:
                                    # Parse frequency value
                                    parts = line.split(':')
                                    if len(parts) >= 3:
                                        try:
                                            freq_cm1 = float(parts[2]) * 219474.63  # Convert from Hartree to cm^-1
                                            frequencies.append(freq_cm1)
                                        except (ValueError, IndexError):
                                            continue
                                elif line == '}':
                                    break
                        
                        if frequencies:
                            print(f"Found {len(frequencies)} vibrational frequencies:")
                            
                            # Separate imaginary and real frequencies
                            imaginary_freqs = [f for f in frequencies if f < 0]
                            real_freqs = [f for f in frequencies if f >= 0]
                            
                            if imaginary_freqs:
                                print(f"⚠️  Imaginary frequencies ({len(imaginary_freqs)}): {imaginary_freqs[:5]}...")
                            
                            if real_freqs:
                                print(f"✅ Real frequencies ({len(real_freqs)}): {real_freqs[:10]}... cm⁻¹")
                                
                        else:
                            print("⚠️  No frequencies found in vibrations.tag")
                            
                    except Exception as e:
                        print(f"⚠️  Error parsing frequencies: {e}")
                
                if os.path.exists(hessian_file):
                    print(f"Hessian matrix saved to: {hessian_dir}/hessian.out")
                    
                    # Post-process Hessian to get frequencies
                    print("\n🔄 Post-processing Hessian matrix to extract frequencies...")
                    frequencies, modes, freq_success = process_dftb_hessian_to_frequencies(abs_hessian_dir, abs_structure_path)
                    
                    if freq_success and frequencies:
                        print("✅ Frequency analysis completed")
                        
                        # Save frequencies for thermodynamic analysis
                        try:
                            save_frequencies_for_thermodynamics(frequencies, abs_project_path, 'dftb')
                            print("📊 Frequencies saved for thermodynamic analysis")
                        except Exception as e:
                            print(f"⚠️  Warning: Could not save frequencies for thermodynamics: {e}")
                    else:
                        print("⚠️  Frequency analysis failed, but Hessian matrix is available")
                
                print(f"All results saved to: {hessian_dir}/")
                success_count += 1
                
            else:
                print(f"❌ DFTB+ Hessian failed for {filename}")
                print("STDOUT:", result.stdout)
                print("STDERR:", result.stderr)
                
        except subprocess.TimeoutExpired:
            print(f"❌ DFTB+ Hessian calculation timed out after 1 hour for {filename}")
            
        except Exception as e:
            print(f"❌ Error during DFTB+ Hessian calculation for {filename}: {e}")
            
        finally:
            # Return to original directory
            os.chdir(original_cwd)
    
    if success_count > 0:
        print(f"\n✅ Successfully calculated native DFTB+ Hessian for {success_count}/{len(structure_files)} structures")
        return True
    else:
        print(f"\n❌ Failed to calculate DFTB+ Hessian for any structures")
        return False


def process_dftb_hessian_to_frequencies(hessian_dir, structure_path):
    """
    Post-process DFTB+ Hessian matrix to extract vibrational frequencies
    
    Args:
        hessian_dir (str): Directory containing DFTB+ Hessian results
        structure_path (str): Path to the molecular structure file
        
    Returns:
        tuple: (frequencies_cm1, normal_modes, success)
    """
    import numpy as np
    import os
    from ase.io import read
    from ase.units import Hartree, Bohr
    
    try:
        # Physical constants
        hartree_to_cm1 = 219474.63  # Hartree to cm^-1 conversion
        amu_to_au = 1822.888486209  # atomic mass units to atomic units
        
        # Read molecular structure
        atoms = read(structure_path)
        n_atoms = len(atoms)
        masses = atoms.get_masses()
        
        print(f"Processing Hessian for {n_atoms} atoms...")
        
        # Read Hessian matrix from hessian.out
        hessian_file = os.path.join(hessian_dir, "hessian.out")
        if not os.path.exists(hessian_file):
            print(f"❌ Hessian file not found: {hessian_file}")
            return None, None, False
        
        # Parse Hessian matrix (3N x 3N)
        hessian_data = []
        with open(hessian_file, 'r') as f:
            for line in f:
                # Each line contains multiple Hessian matrix elements
                elements = line.strip().split()
                for element in elements:
                    try:
                        hessian_data.append(float(element))
                    except ValueError:
                        continue
        
        # Reshape to 3N x 3N matrix
        n_dof = 3 * n_atoms  # degrees of freedom
        expected_size = n_dof * n_dof
        
        if len(hessian_data) != expected_size:
            print(f"❌ Hessian matrix size mismatch: expected {expected_size}, got {len(hessian_data)}")
            return None, None, False
        
        hessian = np.array(hessian_data).reshape(n_dof, n_dof)
        print(f"✅ Read {n_dof}x{n_dof} Hessian matrix")
        
        # Mass-weight the Hessian matrix
        # H_mw[i,j] = H[i,j] / sqrt(m_i * m_j)
        mass_weights = np.zeros(n_dof)
        for i in range(n_atoms):
            mass_weights[3*i:3*i+3] = 1.0 / np.sqrt(masses[i] * amu_to_au)
        
        # Create mass-weighted Hessian
        mass_weighted_hessian = np.zeros_like(hessian)
        for i in range(n_dof):
            for j in range(n_dof):
                mass_weighted_hessian[i,j] = hessian[i,j] * mass_weights[i] * mass_weights[j]
        
        print("✅ Mass-weighted Hessian matrix")
        
        # Diagonalize mass-weighted Hessian
        eigenvalues, eigenvectors = np.linalg.eigh(mass_weighted_hessian)
        
        # Convert eigenvalues to frequencies
        frequencies = np.zeros(n_dof)
        for i, eigenval in enumerate(eigenvalues):
            if eigenval >= 0:
                frequencies[i] = np.sqrt(eigenval) * hartree_to_cm1
            else:
                # Imaginary frequency
                frequencies[i] = -np.sqrt(-eigenval) * hartree_to_cm1
        
        print(f"✅ Calculated {n_dof} vibrational frequencies")
        
        # Sort frequencies (lowest to highest)
        sorted_indices = np.argsort(frequencies)
        frequencies = frequencies[sorted_indices]
        eigenvectors = eigenvectors[:, sorted_indices]
        
        # Classify frequencies
        # First 6 modes should be translations (3) and rotations (3) ≈ 0 cm^-1
        # Remaining modes are vibrational
        
        # Count near-zero frequencies (within 50 cm^-1 of zero)
        near_zero_threshold = 50.0
        translation_rotation_modes = np.sum(np.abs(frequencies) < near_zero_threshold)
        vibrational_modes = n_dof - translation_rotation_modes
        
        print(f"📊 Frequency analysis:")
        print(f"   Translation/rotation modes: {translation_rotation_modes} (≈0 cm⁻¹)")
        print(f"   Vibrational modes: {vibrational_modes}")
        
        # Separate imaginary and real frequencies
        imaginary_freqs = frequencies[frequencies < -near_zero_threshold]
        real_freqs = frequencies[frequencies > near_zero_threshold]
        
        if len(imaginary_freqs) > 0:
            freq_str = ', '.join([f'{f:.1f}' for f in imaginary_freqs[:5]])
            print(f"⚠️  Imaginary frequencies ({len(imaginary_freqs)}): {freq_str}... cm⁻¹")
        
        if len(real_freqs) > 0:
            print(f"✅ Real vibrational frequencies ({len(real_freqs)}):")
            low_str = ', '.join([f'{f:.1f}' for f in real_freqs[:5]])
            print(f"   Lowest: {low_str} cm⁻¹")
            if len(real_freqs) > 5:
                high_str = ', '.join([f'{f:.1f}' for f in real_freqs[-5:]])
                print(f"   Highest: {high_str} cm⁻¹")
        
        # Save frequencies to file
        freq_file = os.path.join(hessian_dir, "frequencies.txt")
        with open(freq_file, 'w') as f:
            f.write(f"# Vibrational frequencies for {n_atoms} atom system\n")
            f.write(f"# {len(imaginary_freqs)} imaginary, {len(real_freqs)} real frequencies\n")
            f.write("# Mode  Frequency(cm^-1)  Type\n")
            
            for i, freq in enumerate(frequencies):
                if abs(freq) < near_zero_threshold:
                    freq_type = "trans/rot"
                elif freq < -near_zero_threshold:
                    freq_type = "imaginary"
                else:
                    freq_type = "vibrational"
                f.write(f"{i+1:4d}  {freq:12.2f}      {freq_type}\n")
        
        print(f"💾 Frequencies saved to: {freq_file}")
        
        return frequencies, eigenvectors, True
        
    except Exception as e:
        print(f"❌ Error processing Hessian: {e}")
        return None, None, False


def analyze_geometry_accuracy(structure1_path, structure2_path, method1_name="Method 1", method2_name="Method 2"):
    """
    Compare geometric parameters between two optimized structures
    
    Args:
        structure1_path (str): Path to first structure (e.g., DFTB+)
        structure2_path (str): Path to second structure (e.g., Orca)
        method1_name (str): Name of first method for display
        method2_name (str): Name of second method for display
        
    Returns:
        dict: Analysis results with bond lengths, angles, and errors
    """
    import numpy as np
    from ase.io import read
    from ase import Atoms
    import itertools
    
    try:
        # Read structures
        atoms1 = read(structure1_path)
        atoms2 = read(structure2_path)
        
        print(f"\n🔍 Geometry Comparison: {method1_name} vs {method2_name}")
        print(f"Structure 1: {len(atoms1)} atoms from {structure1_path}")
        print(f"Structure 2: {len(atoms2)} atoms from {structure2_path}")
        
        if len(atoms1) != len(atoms2):
            print("❌ Different number of atoms in structures")
            return None
        
        # Get positions and chemical symbols
        pos1 = atoms1.get_positions()
        pos2 = atoms2.get_positions()
        symbols = atoms1.get_chemical_symbols()
        
        results = {
            'bond_lengths': {'method1': [], 'method2': [], 'errors': []},
            'bond_angles': {'method1': [], 'method2': [], 'errors': []},
            'dihedral_angles': {'method1': [], 'method2': [], 'errors': []}
        }
        
        # Analyze W-C bond lengths
        w_indices = [i for i, s in enumerate(symbols) if s == 'W']
        c_indices = [i for i, s in enumerate(symbols) if s == 'C']
        
        print(f"\n📏 W-C Bond Length Analysis:")
        print(f"Found {len(w_indices)} W atoms and {len(c_indices)} C atoms")
        
        wc_bonds1 = []
        wc_bonds2 = []
        
        for w_idx in w_indices:
            for c_idx in c_indices:
                dist1 = np.linalg.norm(pos1[w_idx] - pos1[c_idx])
                dist2 = np.linalg.norm(pos2[w_idx] - pos2[c_idx])
                
                # Consider as bonded if distance < 3.0 Å
                if dist1 < 3.0 and dist2 < 3.0:
                    wc_bonds1.append(dist1)
                    wc_bonds2.append(dist2)
                    error = abs(dist1 - dist2)
                    results['bond_lengths']['method1'].append(dist1)
                    results['bond_lengths']['method2'].append(dist2)
                    results['bond_lengths']['errors'].append(error)
                    
                    print(f"   W{w_idx}-C{c_idx}: {dist1:.4f} Å vs {dist2:.4f} Å (Δ = {error:.4f} Å)")
        
        # Analyze C-H bond lengths  
        h_indices = [i for i, s in enumerate(symbols) if s == 'H']
        print(f"\n📏 C-H Bond Length Analysis:")
        
        ch_bonds1 = []
        ch_bonds2 = []
        
        for c_idx in c_indices:
            for h_idx in h_indices:
                dist1 = np.linalg.norm(pos1[c_idx] - pos1[h_idx])
                dist2 = np.linalg.norm(pos2[c_idx] - pos2[h_idx])
                
                # Consider as bonded if distance < 1.3 Å
                if dist1 < 1.3 and dist2 < 1.3:
                    ch_bonds1.append(dist1)
                    ch_bonds2.append(dist2)
                    error = abs(dist1 - dist2)
                    results['bond_lengths']['method1'].append(dist1)
                    results['bond_lengths']['method2'].append(dist2)
                    results['bond_lengths']['errors'].append(error)
        
        # Print C-H statistics instead of all individual bonds (too many)
        if ch_bonds1:
            ch_mean_error = np.mean([abs(b1-b2) for b1,b2 in zip(ch_bonds1, ch_bonds2)])
            ch_max_error = max([abs(b1-b2) for b1,b2 in zip(ch_bonds1, ch_bonds2)])
            print(f"   {len(ch_bonds1)} C-H bonds: mean Δ = {ch_mean_error:.4f} Å, max Δ = {ch_max_error:.4f} Å")
        
        # Analyze C-W-C bond angles
        print(f"\n📐 C-W-C Bond Angle Analysis:")
        
        for w_idx in w_indices:
            # Find carbons bonded to this tungsten
            bonded_carbons = []
            for c_idx in c_indices:
                dist = np.linalg.norm(pos1[w_idx] - pos1[c_idx])
                if dist < 3.0:  # bonded
                    bonded_carbons.append(c_idx)
            
            # Calculate angles between carbon pairs
            for i, c1_idx in enumerate(bonded_carbons):
                for j, c2_idx in enumerate(bonded_carbons):
                    if i < j:  # avoid duplicates
                        # Calculate angle C1-W-C2
                        vec1_1 = pos1[c1_idx] - pos1[w_idx]
                        vec2_1 = pos1[c2_idx] - pos1[w_idx]
                        vec1_2 = pos2[c1_idx] - pos2[w_idx]
                        vec2_2 = pos2[c2_idx] - pos2[w_idx]
                        
                        # Calculate angles
                        cos_angle1 = np.dot(vec1_1, vec2_1) / (np.linalg.norm(vec1_1) * np.linalg.norm(vec2_1))
                        cos_angle2 = np.dot(vec1_2, vec2_2) / (np.linalg.norm(vec1_2) * np.linalg.norm(vec2_2))
                        
                        # Clamp to [-1, 1] to avoid numerical errors
                        cos_angle1 = np.clip(cos_angle1, -1, 1)
                        cos_angle2 = np.clip(cos_angle2, -1, 1)
                        
                        angle1 = np.degrees(np.arccos(cos_angle1))
                        angle2 = np.degrees(np.arccos(cos_angle2))
                        error = abs(angle1 - angle2)
                        
                        results['bond_angles']['method1'].append(angle1)
                        results['bond_angles']['method2'].append(angle2)
                        results['bond_angles']['errors'].append(error)
                        
                        print(f"   C{c1_idx}-W{w_idx}-C{c2_idx}: {angle1:.2f}° vs {angle2:.2f}° (Δ = {error:.2f}°)")
        
        # Summary statistics
        print(f"\n📊 Summary Statistics:")
        
        if results['bond_lengths']['errors']:
            bond_errors = np.array(results['bond_lengths']['errors'])
            print(f"Bond Length Errors:")
            print(f"   Mean: {np.mean(bond_errors):.4f} Å")
            print(f"   Max:  {np.max(bond_errors):.4f} Å")
            print(f"   RMS:  {np.sqrt(np.mean(bond_errors**2)):.4f} Å")
        
        if results['bond_angles']['errors']:
            angle_errors = np.array(results['bond_angles']['errors'])
            print(f"Bond Angle Errors:")
            print(f"   Mean: {np.mean(angle_errors):.2f}°")
            print(f"   Max:  {np.max(angle_errors):.2f}°") 
            print(f"   RMS:  {np.sqrt(np.mean(angle_errors**2)):.2f}°")
        
        return results
        
    except Exception as e:
        print(f"❌ Error analyzing geometry: {e}")
        return None


def escape_saddle_point(project_folder, displacement_scale=0.1, backend='dftb'):
    """
    Escape from saddle point by following imaginary modes and re-optimizing
    
    Args:
        project_folder (str): Project folder name
        displacement_scale (float): Scale factor for displacement along imaginary modes
        backend (str): Backend to use for re-optimization ('dftb', 'xtb', 'orca')
    
    Returns:
        bool: True if successful, False otherwise
    """
    import numpy as np
    from ase.io import read, write
    import os
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
    import settings
    
    try:
        print(f"\n🔄 Escaping saddle point for {project_folder}")
        
        # Setup paths
        project_data_dir = os.path.join(settings.DATA_FOLDER, project_folder)
        
        # Find the relaxed structure (try different backends)
        possible_paths = [
            os.path.join(project_data_dir, "outputs", "ptbp_dftb", f"{project_folder}_relaxed.xyz"),
            os.path.join(project_data_dir, "outputs", "xtb", f"{project_folder}_relaxed.xyz"),
            os.path.join(project_data_dir, "outputs", "orca", "PBE_def2-SVP", f"{project_folder}_relaxed.xyz")
        ]
        
        structure_path = None
        for path in possible_paths:
            if os.path.exists(path):
                structure_path = path
                break
                
        if not structure_path:
            print("❌ Could not find relaxed structure")
            return False
            
        # Find the Hessian output directory
        hessian_dirs = [
            os.path.join(project_data_dir, "outputs", f"{project_folder}_relaxed_hessian_dftb_native"),
            os.path.join(project_data_dir, "outputs", f"{project_folder}_relaxed_hessian_orca"),
            os.path.join(project_data_dir, "outputs", f"{project_folder}_relaxed_hessian_xtb")
        ]
        
        hessian_dir = None
        for hdir in hessian_dirs:
            if os.path.exists(os.path.join(hdir, "frequencies.txt")):
                hessian_dir = hdir
                break
                
        if not hessian_dir:
            print("❌ Could not find Hessian calculation results")
            return False
            
        # Read frequencies and modes
        freq_file = os.path.join(hessian_dir, "frequencies.txt")
        with open(freq_file, 'r') as f:
            lines = f.readlines()
            
        # Find imaginary frequencies
        imaginary_modes = []
        for line in lines:
            if 'imaginary' in line and not 'trans/rot' in line and not line.strip().startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        mode_idx = int(parts[0]) - 1  # Convert to 0-based
                        frequency = float(parts[1])
                        imaginary_modes.append((mode_idx, frequency))
                    except (ValueError, IndexError):
                        continue
                
        if not imaginary_modes:
            print("✅ No imaginary frequencies found - structure is at minimum")
            return True
            
        print(f"Found {len(imaginary_modes)} imaginary modes:")
        for mode_idx, freq in imaginary_modes:
            print(f"   Mode {mode_idx+1}: {freq:.2f} cm⁻¹")
            
        # Read the structure
        atoms = read(structure_path)
        original_positions = atoms.get_positions().copy()
        
        # Try to read eigenvectors from hessian calculation
        # For now, use random displacement along imaginary modes
        # In a full implementation, we'd read the actual normal mode vectors
        
        # Apply small displacement to escape saddle point
        # Use the most negative imaginary mode
        worst_mode = min(imaginary_modes, key=lambda x: x[1])
        print(f"Using mode {worst_mode[0]+1} (freq: {worst_mode[1]:.2f} cm⁻¹) for displacement")
        
        # Apply random displacement scaled by the magnitude of imaginary frequency
        np.random.seed(42)  # For reproducibility
        n_atoms = len(atoms)
        displacement = np.random.randn(n_atoms, 3) * displacement_scale
        
        # Scale displacement by frequency magnitude (larger displacement for more negative frequencies)
        scale_factor = abs(worst_mode[1]) / 100.0  # Scale by frequency magnitude
        displacement *= scale_factor
        
        new_positions = original_positions + displacement
        atoms.set_positions(new_positions)
        
        # Create output directory for continued relaxation
        output_dir = os.path.join(project_data_dir, "outputs", f"continued_{backend}")
        os.makedirs(output_dir, exist_ok=True)
        
        # Save perturbed structure
        perturbed_path = os.path.join(output_dir, f"{project_folder}_perturbed.xyz")
        write(perturbed_path, atoms)
        
        print(f"💾 Saved perturbed structure to {perturbed_path}")
        print(f"   RMS displacement: {np.sqrt(np.mean(displacement**2)):.4f} Å")
        
        # Now re-optimize with the chosen backend
        print(f"\n🔄 Re-optimizing with {backend} backend...")
        
        if backend == 'dftb':
            from .relax import relax_native_dftb_ptbp
            # Copy perturbed structure to project folder as initial structure
            import shutil
            project_input_path = os.path.join(project_data_dir, f"{project_folder}_initial.xyz")
            shutil.copy(perturbed_path, project_input_path)
            success = relax_native_dftb_ptbp(project_folder)
        elif backend == 'xtb':
            from .relax import relax_xtb
            # Copy perturbed structure to project folder as initial structure  
            import shutil
            project_input_path = os.path.join(project_data_dir, f"{project_folder}_initial.xyz")
            shutil.copy(perturbed_path, project_input_path)
            success = relax_xtb(project_folder, method='GFN2-xTB', backend='fixed')
        elif backend == 'orca':
            from .relax import relax_orca_dft
            # Copy perturbed structure to project folder as initial structure
            import shutil
            project_input_path = os.path.join(project_data_dir, f"{project_folder}_initial.xyz")
            shutil.copy(perturbed_path, project_input_path)
            success = relax_orca_dft(project_folder, method='PBE', basis='def2-SVP', nprocs=4)
        else:
            print(f"❌ Unknown backend: {backend}")
            return False
            
        if success:
            print("✅ Successfully continued optimization from saddle point")
            
            # Check if we need to calculate new Hessian to verify we found a minimum
            print("\n🎯 Recommendation: Calculate Hessian for the new structure to verify it's a minimum")
            return True
        else:
            print("❌ Failed to continue optimization")
            return False
            
    except Exception as e:
        print(f"❌ Error escaping saddle point: {str(e)}")
        return False


def get_relaxed_filename(project_name, iteration=None):
    """Get filename for relaxed structure with iteration support"""
    if iteration is None or iteration == 1:
        return f"{project_name}_relaxed.xyz"
    else:
        return f"{project_name}_relax_iter{iteration}.xyz"


def get_perturbed_filename(project_name, iteration):
    """Get filename for perturbed structure"""
    return f"{project_name}_perturbed_iter{iteration}.xyz"


def get_frequencies_filename(iteration=None):
    """Get filename for frequencies with iteration support"""
    if iteration is None or iteration == 1:
        return "frequencies.txt"
    else:
        return f"frequencies_iter{iteration}.txt"


def get_backend_output_dir(project_folder, backend):
    """Get output directory for specific backend"""
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
    import settings
    
    project_data_dir = os.path.join(settings.DATA_FOLDER, project_folder)
    
    backend_map = {
        'dftb': 'ptbp_dftb',
        'dftb-ptbp': 'ptbp_dftb',
        'dftb-native': 'dftb_native',
        'dftb-passivated': 'dftb_passivated', 
        'xtb': 'xtb',
        'gxtb': 'gxtb',
        'orca': 'orca'
    }
    
    backend_dir = backend_map.get(backend, backend)
    return os.path.join(project_data_dir, "outputs", backend_dir)


def count_imaginary_frequencies(project_folder, backend, iteration=None):
    """Count imaginary frequencies from frequencies file"""
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
    import settings
    
    try:
        project_data_dir = os.path.join(settings.DATA_FOLDER, project_folder)
        
        # Find Hessian directory - try different naming patterns
        hessian_suffix = "dftb_native" if "dftb" in backend else backend
        outputs_dir = os.path.join(project_data_dir, "outputs")
        
        # Try different possible hessian directory names
        possible_hessian_dirs = [
            f"{project_folder}_relaxed_hessian_{hessian_suffix}",  # project-based name
            f"*_relaxed_hessian_{hessian_suffix}"  # any structure-based name
        ]
        
        hessian_dir = None
        for pattern in possible_hessian_dirs:
            if '*' in pattern:
                # Use glob to find matching directories
                import glob
                matches = glob.glob(os.path.join(outputs_dir, pattern))
                if matches:
                    hessian_dir = matches[0]  # Use first match
                    break
            else:
                candidate = os.path.join(outputs_dir, pattern)
                if os.path.exists(candidate):
                    hessian_dir = candidate
                    break
        
        if hessian_dir is None:
            print(f"❌ Hessian directory not found in {outputs_dir}")
            return -1
        
        freq_file = os.path.join(hessian_dir, get_frequencies_filename(iteration))
        
        if not os.path.exists(freq_file):
            print(f"❌ Frequencies file not found: {freq_file}")
            return -1
            
        imaginary_count = 0
        with open(freq_file, 'r') as f:
            for line in f:
                if 'imaginary' in line and not 'trans/rot' in line and not line.strip().startswith('#'):
                    imaginary_count += 1
                    
        return imaginary_count
        
    except Exception as e:
        print(f"❌ Error counting imaginary frequencies: {str(e)}")
        return -1


def rename_relaxed_result(project_folder, backend, iteration):
    """Rename relaxed result to iteration-specific filename"""
    try:
        backend_dir = get_backend_output_dir(project_folder, backend)
        
        # Source: standard relaxed filename
        src_file = os.path.join(backend_dir, f"{project_folder}_relaxed.xyz")
        
        # Destination: iteration-specific filename  
        dst_file = os.path.join(backend_dir, get_relaxed_filename(project_folder, iteration))
        
        if os.path.exists(src_file) and src_file != dst_file:
            import shutil
            shutil.copy2(src_file, dst_file)
            print(f"📁 Saved iteration {iteration} result: {os.path.basename(dst_file)}")
            return True
            
    except Exception as e:
        print(f"❌ Error renaming relaxed result: {str(e)}")
        
    return False


def rename_hessian_result(project_folder, backend, iteration):
    """Rename hessian frequencies to iteration-specific filename"""
    try:
        import sys
        sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
        import settings
        
        project_data_dir = os.path.join(settings.DATA_FOLDER, project_folder)
        hessian_suffix = "dftb_native" if "dftb" in backend else backend
        hessian_dir = os.path.join(project_data_dir, "outputs", 
                                  f"{project_folder}_relaxed_hessian_{hessian_suffix}")
        
        # Source: standard frequencies filename
        src_file = os.path.join(hessian_dir, "frequencies.txt")
        
        # Destination: iteration-specific filename
        dst_file = os.path.join(hessian_dir, get_frequencies_filename(iteration))
        
        if os.path.exists(src_file) and src_file != dst_file:
            import shutil
            shutil.copy2(src_file, dst_file)
            print(f"📁 Saved iteration {iteration} frequencies: {os.path.basename(dst_file)}")
            return True
            
    except Exception as e:
        print(f"❌ Error renaming hessian result: {str(e)}")
        
    return False


def finalize_results(project_folder, backend, final_iteration):
    """Finalize results by ensuring final files have standard names"""
    try:
        print(f"\n🎯 Finalizing results from iteration {final_iteration}...")
        
        backend_dir = get_backend_output_dir(project_folder, backend)
        
        # Copy final relaxed structure to standard name
        final_relaxed = os.path.join(backend_dir, get_relaxed_filename(project_folder, final_iteration))
        standard_relaxed = os.path.join(backend_dir, f"{project_folder}_relaxed.xyz")
        
        if os.path.exists(final_relaxed) and final_relaxed != standard_relaxed:
            import shutil
            shutil.copy2(final_relaxed, standard_relaxed)
            print(f"✅ Final structure: {os.path.basename(standard_relaxed)}")
            
        # Copy final frequencies to standard name
        import sys
        sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
        import settings
        
        project_data_dir = os.path.join(settings.DATA_FOLDER, project_folder)
        hessian_suffix = "dftb_native" if "dftb" in backend else backend
        hessian_dir = os.path.join(project_data_dir, "outputs", 
                                  f"{project_folder}_relaxed_hessian_{hessian_suffix}")
        
        final_frequencies = os.path.join(hessian_dir, get_frequencies_filename(final_iteration))
        standard_frequencies = os.path.join(hessian_dir, "frequencies.txt")
        
        if os.path.exists(final_frequencies) and final_frequencies != standard_frequencies:
            import shutil
            shutil.copy2(final_frequencies, standard_frequencies)
            print(f"✅ Final frequencies: frequencies.txt")
            
        return True
        
    except Exception as e:
        print(f"❌ Error finalizing results: {str(e)}")
        return False


def escape_saddle_point_with_iteration(project_folder, backend, iteration):
    """Escape from saddle point with iteration support and proper file naming"""
    import numpy as np
    from ase.io import read, write
    import os
    import sys
    sys.path.append(os.path.join(os.path.dirname(__file__), '..', '..'))
    import settings
    
    try:
        print(f"\n🔄 Escaping saddle point for {project_folder} (iteration {iteration})")
        
        # Setup paths
        project_data_dir = os.path.join(settings.DATA_FOLDER, project_folder)
        
        # Find the relaxed structure from current iteration
        backend_dir = get_backend_output_dir(project_folder, backend)
        current_relaxed = os.path.join(backend_dir, get_relaxed_filename(project_folder, iteration))
        
        if not os.path.exists(current_relaxed):
            # Fallback to standard name if iteration-specific doesn't exist
            current_relaxed = os.path.join(backend_dir, f"{project_folder}_relaxed.xyz")
            
        if not os.path.exists(current_relaxed):
            print(f"❌ Could not find relaxed structure: {current_relaxed}")
            return False
            
        # Find the Hessian results from current iteration
        hessian_suffix = "dftb_native" if "dftb" in backend else backend
        hessian_dir = os.path.join(project_data_dir, "outputs", 
                                  f"{project_folder}_relaxed_hessian_{hessian_suffix}")
        
        freq_file = os.path.join(hessian_dir, get_frequencies_filename(iteration))
        
        if not os.path.exists(freq_file):
            # Fallback to standard name
            freq_file = os.path.join(hessian_dir, "frequencies.txt")
            
        if not os.path.exists(freq_file):
            print(f"❌ Could not find frequencies file")
            return False
            
        # Read frequencies and find imaginary modes
        with open(freq_file, 'r') as f:
            lines = f.readlines()
            
        imaginary_modes = []
        for line in lines:
            if 'imaginary' in line and not 'trans/rot' in line and not line.strip().startswith('#'):
                parts = line.strip().split()
                if len(parts) >= 3:
                    try:
                        mode_idx = int(parts[0]) - 1  # Convert to 0-based
                        frequency = float(parts[1])
                        imaginary_modes.append((mode_idx, frequency))
                    except (ValueError, IndexError):
                        continue
                        
        if not imaginary_modes:
            print("✅ No imaginary frequencies found - structure is at minimum")
            return True
            
        print(f"Found {len(imaginary_modes)} imaginary modes:")
        for mode_idx, freq in imaginary_modes:
            print(f"   Mode {mode_idx+1}: {freq:.2f} cm⁻¹")
            
        # Read the structure
        atoms = read(current_relaxed)
        original_positions = atoms.get_positions().copy()
        
        # Apply displacement based on most negative imaginary mode
        worst_mode = min(imaginary_modes, key=lambda x: x[1])
        print(f"Using mode {worst_mode[0]+1} (freq: {worst_mode[1]:.2f} cm⁻¹) for displacement")
        
        # Apply random displacement scaled by frequency magnitude  
        np.random.seed(42 + iteration)  # Different seed for each iteration
        n_atoms = len(atoms)
        displacement_scale = 0.1  # Default scale
        displacement = np.random.randn(n_atoms, 3) * displacement_scale
        
        # Scale by frequency magnitude
        scale_factor = abs(worst_mode[1]) / 100.0
        displacement *= scale_factor
        
        new_positions = original_positions + displacement
        atoms.set_positions(new_positions)
        
        # Save perturbed structure with iteration-specific name
        perturbed_filename = get_perturbed_filename(project_folder, iteration)
        perturbed_path = os.path.join(backend_dir, perturbed_filename)
        
        # Ensure backend directory exists
        os.makedirs(backend_dir, exist_ok=True)
        
        write(perturbed_path, atoms)
        
        print(f"💾 Saved perturbed structure: {os.path.basename(perturbed_path)}")
        print(f"   RMS displacement: {np.sqrt(np.mean(displacement**2)):.4f} Å")
        
        # Copy perturbed structure as initial structure for next iteration
        initial_path = os.path.join(project_data_dir, f"{project_folder}_initial.xyz")
        import shutil
        shutil.copy2(perturbed_path, initial_path)
        
        return True
        
    except Exception as e:
        print(f"❌ Error escaping saddle point: {str(e)}")
        return False


def calculate_hessian_dftb_unified(project_folder, params='auto'):
    """
    Unified DFTB+ Hessian calculation with automatic parameter selection
    
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
    project_data_dir = base_path
    
    # Create backend instance
    backend = DFTBBackend(params=params)
    
    if not backend.is_available():
        print("❌ DFTB+ backend not available")
        return False
    
    # Look for relaxed structures in various backend output directories
    structure_files = []
    outputs_dir = os.path.join(project_data_dir, "outputs")
    
    # Search order: PTBP -> native -> passivated -> any
    search_patterns = [
        os.path.join(outputs_dir, "ptbp_dftb", "*_relaxed.xyz"),
        os.path.join(outputs_dir, "native", "*_relaxed.xyz"), 
        os.path.join(outputs_dir, "passivated", "*_relaxed.xyz"),
        os.path.join(outputs_dir, "*", "*_relaxed.xyz"),
        os.path.join(project_data_dir, "*_relaxed.xyz")
    ]
    
    for pattern in search_patterns:
        files = glob.glob(pattern)
        if files:
            structure_files.extend(files)
            break
    
    if not structure_files:
        print(f"❌ No relaxed structures found in {project_data_dir}")
        print("Run relaxation first with: python main.py relax --backend dftb <project_name>")
        return False
    
    # Use the first structure found
    structure_file = structure_files[0]
    print(f"📁 Using relaxed structure: {os.path.basename(structure_file)}")
    
    # Read atoms to determine parameter set
    atoms = read(structure_file)
    
    # Auto-detect parameter set if needed
    if params == 'auto':
        param_set = backend._detect_parameter_set(atoms)
    else:
        param_set = params
    
    param_path = backend._get_parameter_path(param_set)
    if not os.path.exists(param_path):
        print(f"❌ Parameter set '{param_set}' not found at {param_path}")
        return False
    
    # Convert to absolute path
    param_path = os.path.abspath(param_path)
    
    print(f"🔧 Using {param_set.upper()} parameters for Hessian calculation")
    print(f"📊 Structure: {atoms.get_chemical_formula()} ({len(atoms)} atoms)")
    
    # Create temp directory for calculation
    temp_dir = os.path.join(project_data_dir, "temp")
    os.makedirs(temp_dir, exist_ok=True)
    
    # Store original working directory
    original_cwd = os.getcwd()
    
    try:
        os.chdir(temp_dir)
        
        # Copy structure to temp directory
        name_without_ext = os.path.splitext(os.path.basename(structure_file))[0]
        temp_structure = f"{name_without_ext}.xyz"
        import shutil
        shutil.copy2(structure_file, temp_structure)
        
        # Write gen file
        gen_file = f"{name_without_ext}.gen"
        write(gen_file, atoms, format='gen')
        
        # Generate angular momentum parameters
        angular_params = backend._get_angular_momentum_params(atoms, param_set)
        
        # Create DFTB+ input file for Hessian calculation
        hsd_content = f"""
Geometry = GenFormat {{
    <<< "{gen_file}"
}}

Driver = SecondDerivatives {{
    Atoms = 1:-1
    Delta = 1.0e-4
}}

Hamiltonian = DFTB {{
    SCC = Yes
    SCCTolerance = 1e-7
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
    WriteDetailedOut = Yes
}}

ParserOptions {{
    ParserVersion = 7
}}
"""
        
        # Write DFTB+ input file
        hsd_file = f"dftb_in.hsd"
        with open(hsd_file, 'w') as f:
            f.write(hsd_content)
        
        print(f"🚀 Starting DFTB+ Hessian calculation...")
        print(f"   Parameters: {param_set.upper()}")
        print(f"   Formula: {atoms.get_chemical_formula()}")
        
        start_time = time.time()
        
        # Run DFTB+
        try:
            dftb_executable = '/Users/andreypanferov/opt/dftb+/bin/dftb+'
            
            # Set environment variable for Slater-Koster files
            env = os.environ.copy()
            env['SKDIR'] = param_path
            
            print(f"💻 Running DFTB+ with environment SKDIR={param_path}")
            
            result = subprocess.run([dftb_executable], 
                                 capture_output=True, 
                                 text=True, 
                                 timeout=7200,  # 2 hour timeout
                                 env=env)
            
            elapsed_time = time.time() - start_time
            
            # Always show DFTB+ output for debugging
            if result.stdout:
                print("📋 DFTB+ STDOUT:")
                print(result.stdout[-1000:])  # Last 1000 chars
            if result.stderr:
                print("❌ DFTB+ STDERR:")
                print(result.stderr[-1000:])  # Last 1000 chars
            
            if result.returncode == 0:
                print_completion_message(f"DFTB+ Hessian ({param_set.upper()})", name_without_ext, elapsed_time)
                
                # Process results - look for hessian output
                hessian_file = "hessian.out"
                if os.path.exists(hessian_file):
                    # Copy results back to project directory
                    output_dir = os.path.join(project_data_dir, "outputs", f"{param_set}_dftb")
                    os.makedirs(output_dir, exist_ok=True)
                    
                    shutil.copy2(hessian_file, os.path.join(output_dir, f"{name_without_ext}_hessian.out"))
                    
                    # Copy other output files
                    for outfile in ["results.tag", "detailed.out", "band.out"]:
                        if os.path.exists(outfile):
                            shutil.copy2(outfile, os.path.join(output_dir, f"{name_without_ext}_{outfile}"))
                    
                    # Extract frequencies from hessian.out
                    print("🔄 Extracting vibrational frequencies from Hessian matrix...")
                    try:
                        frequencies = extract_frequencies_from_dftb_hessian(hessian_file, atoms)
                        if frequencies:
                            print(f"✅ Found {len(frequencies)} vibrational frequencies")
                            print(f"   Range: {min(frequencies):.1f} - {max(frequencies):.1f} cm⁻¹")
                            
                            # Save frequencies for thermodynamic analysis
                            save_frequencies_for_thermodynamics(frequencies, project_data_dir, 'dftb')
                            print("📊 Frequencies saved for thermodynamic analysis")
                        else:
                            print("⚠️  No frequencies extracted from Hessian matrix")
                    except Exception as e:
                        print(f"⚠️  Error extracting frequencies: {e}")
                    
                    print(f"💾 Results saved to {output_dir}/")
                    return True
                else:
                    print("❌ Hessian output file not found")
                    return False
            else:
                print_failure_message(f"DFTB+ Hessian ({param_set.upper()})", name_without_ext, elapsed_time)
                if result.stderr:
                    print(f"DFTB+ error: {result.stderr}")
                return False
                
        except subprocess.TimeoutExpired:
            print_timeout_message(f"DFTB+ Hessian ({param_set.upper()})", name_without_ext, 7200)
            return False
        except Exception as e:
            elapsed_time = time.time() - start_time
            print_failure_message(f"DFTB+ Hessian ({param_set.upper()})", name_without_ext, elapsed_time, str(e))
            return False
    
    finally:
        os.chdir(original_cwd)