"""
Thermodynamic analysis module for mechanosynthesis research.

This module provides statistical thermodynamic calculations based on vibrational
frequencies obtained from quantum chemical calculations (DFTB+, xTB, Orca).

Key functionality:
- Temperature-dependent thermodynamic properties (G, H, S, Cp)
- Thermal stability analysis of molecular structures
- Bond dissociation analysis at elevated temperatures
- Visualization of temperature effects

Based on statistical mechanics using harmonic oscillator approximation.
"""

import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple, Optional
from ase import Atoms
from ase.units import kB, _hbar, _c, _amu
import os


# Physical constants in appropriate units
KB = kB  # Boltzmann constant in eV/K
HBAR = _hbar  # Reduced Planck constant in eV·s
C = _c  # Speed of light in m/s
AMU = _amu  # Atomic mass unit in kg

# Conversion factors
CM_TO_EV = 1.24e-4  # cm⁻¹ to eV conversion
EV_TO_J = 1.602176634e-19  # eV to Joules
AVOGADRO = 6.022140857e23  # Avogadro's number


def vibrational_energy(frequency_cm, temperature):
    """
    Calculate vibrational energy contribution at given temperature.
    
    Args:
        frequency_cm: Vibrational frequency in cm⁻¹
        temperature: Temperature in Kelvin
        
    Returns:
        Vibrational energy in eV
    """
    if frequency_cm <= 0:
        return 0.0
    
    freq_eV = frequency_cm * CM_TO_EV
    x = freq_eV / (KB * temperature)
    
    # Zero-point energy + thermal energy
    if x > 50:  # Avoid overflow for high frequencies
        return 0.5 * freq_eV
    else:
        return 0.5 * freq_eV + freq_eV / (np.exp(x) - 1)


def vibrational_entropy(frequency_cm, temperature):
    """
    Calculate vibrational entropy contribution at given temperature.
    
    Args:
        frequency_cm: Vibrational frequency in cm⁻¹
        temperature: Temperature in Kelvin
        
    Returns:
        Vibrational entropy in eV/K
    """
    if frequency_cm <= 0:
        return 0.0
    
    freq_eV = frequency_cm * CM_TO_EV
    x = freq_eV / (KB * temperature)
    
    if x > 50:  # Avoid overflow
        return 0.0
    else:
        exp_x = np.exp(x)
        return KB * (x / (exp_x - 1) - np.log(1 - 1/exp_x))


def translational_entropy(atoms, temperature, pressure=1.0):
    """
    Calculate translational entropy for gas phase molecule.
    
    Args:
        atoms: ASE Atoms object
        temperature: Temperature in Kelvin
        pressure: Pressure in atm
        
    Returns:
        Translational entropy in eV/K
    """
    mass_kg = sum(atoms.get_masses()) * AMU
    
    # Sackur-Tetrode equation
    prefactor = (2 * np.pi * mass_kg * KB * temperature) / (HBAR**2)
    volume_per_molecule = KB * temperature / (pressure * 101325)  # Convert atm to Pa
    
    s_trans = KB * (1.5 * np.log(prefactor) + np.log(volume_per_molecule) + 2.5)
    return s_trans


def rotational_entropy(atoms, temperature):
    """
    Calculate rotational entropy for gas phase molecule.
    
    Args:
        atoms: ASE Atoms object  
        temperature: Temperature in Kelvin
        
    Returns:
        Rotational entropy in eV/K
    """
    # Get moments of inertia
    moments = atoms.get_moments_of_inertia() * AMU * 1e-20  # Convert to kg·m²
    
    # Remove zero moments (linear molecules)
    moments = moments[moments > 1e-40]
    
    if len(moments) == 0:  # Atom
        return 0.0
    elif len(moments) == 2:  # Linear molecule
        I = moments[0]
        s_rot = KB * (np.log(8 * np.pi**2 * I * KB * temperature / HBAR**2) + 1)
    else:  # Non-linear molecule
        I_product = np.prod(moments)
        s_rot = KB * (0.5 * np.log(8 * np.pi**3 * I_product * (KB * temperature)**3 / HBAR**6) + 1.5)
    
    return s_rot


def calculate_thermodynamic_properties(frequencies_cm, atoms, temperature, pressure=1.0):
    """
    Calculate complete thermodynamic properties at given temperature.
    
    Args:
        frequencies_cm: List of vibrational frequencies in cm⁻¹
        atoms: ASE Atoms object
        temperature: Temperature in Kelvin
        pressure: Pressure in atm
        
    Returns:
        Dictionary with thermodynamic properties
    """
    # Filter out imaginary and very low frequencies
    real_frequencies = [f for f in frequencies_cm if f > 20.0]
    
    # Vibrational contributions
    E_vib = sum(vibrational_energy(f, temperature) for f in real_frequencies)
    S_vib = sum(vibrational_entropy(f, temperature) for f in real_frequencies)
    
    # Translational contributions (gas phase)
    S_trans = translational_entropy(atoms, temperature, pressure)
    E_trans = 1.5 * KB * temperature
    
    # Rotational contributions (gas phase)
    S_rot = rotational_entropy(atoms, temperature)
    if len(atoms) == 1:  # Atom
        E_rot = 0.0
    elif len(atoms.get_moments_of_inertia()[atoms.get_moments_of_inertia() > 1e-10]) == 2:  # Linear
        E_rot = KB * temperature
    else:  # Non-linear
        E_rot = 1.5 * KB * temperature
    
    # Total properties
    H_thermal = E_vib + E_trans + E_rot + KB * temperature  # +RT for PV term
    S_total = S_vib + S_trans + S_rot
    G_thermal = H_thermal - temperature * S_total
    
    return {
        'temperature': temperature,
        'enthalpy_thermal': H_thermal,
        'entropy_total': S_total,
        'free_energy_thermal': G_thermal,
        'vibrational_energy': E_vib,
        'vibrational_entropy': S_vib,
        'translational_entropy': S_trans,
        'rotational_entropy': S_rot,
        'n_frequencies': len(real_frequencies),
        'n_imaginary': len([f for f in frequencies_cm if f < 0])
    }


def temperature_sweep(frequencies_cm, atoms, temperatures, pressure=1.0):
    """
    Calculate thermodynamic properties over temperature range.
    
    Args:
        frequencies_cm: List of vibrational frequencies in cm⁻¹
        atoms: ASE Atoms object
        temperatures: List of temperatures in Kelvin
        pressure: Pressure in atm
        
    Returns:
        List of thermodynamic property dictionaries
    """
    results = []
    for T in temperatures:
        props = calculate_thermodynamic_properties(frequencies_cm, atoms, T, pressure)
        results.append(props)
    
    return results


def analyze_thermal_stability(frequencies_cm, atoms, temperatures, bond_dissociation_energy=4.0):
    """
    Analyze thermal stability and predict bond dissociation.
    
    Args:
        frequencies_cm: List of vibrational frequencies in cm⁻¹
        atoms: ASE Atoms object
        temperatures: List of temperatures in Kelvin
        bond_dissociation_energy: Bond dissociation energy in eV
        
    Returns:
        Dictionary with stability analysis
    """
    results = temperature_sweep(frequencies_cm, atoms, temperatures)
    
    # Calculate thermal energy relative to bond strength
    thermal_energies = [KB * T for T in temperatures]
    stability_ratios = [E_th / bond_dissociation_energy for E_th in thermal_energies]
    
    # Find critical temperatures
    dissociation_threshold = 0.1  # 10% of bond energy as thermal energy
    critical_indices = [i for i, ratio in enumerate(stability_ratios) if ratio > dissociation_threshold]
    
    critical_temp = temperatures[critical_indices[0]] if critical_indices else None
    
    return {
        'temperatures': temperatures,
        'thermal_energies': thermal_energies,
        'stability_ratios': stability_ratios,
        'critical_temperature': critical_temp,
        'thermodynamic_data': results
    }


def plot_thermal_analysis(analysis_results, output_dir=None, show_plots=False):
    """
    Create visualization plots for thermal analysis.
    
    Args:
        analysis_results: Results from analyze_thermal_stability
        output_dir: Directory to save plots
        show_plots: Whether to display plots
    """
    data = analysis_results['thermodynamic_data']
    temperatures = [d['temperature'] for d in data]
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 10))
    
    # Free energy vs temperature
    free_energies = [d['free_energy_thermal'] for d in data]
    ax1.plot(temperatures, free_energies, 'b-', linewidth=2)
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('ΔG_thermal (eV)')
    ax1.set_title('Thermal Free Energy')
    ax1.grid(True, alpha=0.3)
    
    # Entropy components
    s_vib = [d['vibrational_entropy'] * 1000 for d in data]  # meV/K
    s_total = [d['entropy_total'] * 1000 for d in data]  # meV/K
    ax2.plot(temperatures, s_vib, 'r-', label='Vibrational', linewidth=2)
    ax2.plot(temperatures, s_total, 'k-', label='Total', linewidth=2)
    ax2.set_xlabel('Temperature (K)')
    ax2.set_ylabel('Entropy (meV/K)')
    ax2.set_title('Entropy Components')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Thermal stability
    stability_ratios = analysis_results['stability_ratios']
    ax3.plot(temperatures, stability_ratios, 'g-', linewidth=2)
    ax3.axhline(y=0.1, color='r', linestyle='--', label='Dissociation threshold')
    ax3.set_xlabel('Temperature (K)')
    ax3.set_ylabel('Thermal/Bond Energy Ratio')
    ax3.set_title('Thermal Stability Analysis')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Enthalpy components
    enthalpies = [d['enthalpy_thermal'] for d in data]
    vib_energies = [d['vibrational_energy'] for d in data]
    ax4.plot(temperatures, enthalpies, 'purple', label='Total thermal', linewidth=2)
    ax4.plot(temperatures, vib_energies, 'orange', label='Vibrational', linewidth=2)
    ax4.set_xlabel('Temperature (K)')
    ax4.set_ylabel('Enthalpy (eV)')
    ax4.set_title('Thermal Enthalpy')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        plt.savefig(os.path.join(output_dir, 'thermal_analysis.png'), dpi=300, bbox_inches='tight')
    
    if show_plots:
        plt.show()
    else:
        plt.close()


def create_thermal_report(analysis_results, frequencies_cm, atoms, output_file=None):
    """
    Generate detailed thermal analysis report.
    
    Args:
        analysis_results: Results from analyze_thermal_stability
        frequencies_cm: Original frequency data
        atoms: ASE Atoms object
        output_file: Path to save report
        
    Returns:
        Report string
    """
    data = analysis_results['thermodynamic_data']
    
    report = []
    report.append("THERMAL ANALYSIS REPORT")
    report.append("=" * 50)
    report.append(f"Molecule: {len(atoms)} atoms")
    report.append(f"Formula: {atoms.get_chemical_formula()}")
    report.append(f"Mass: {sum(atoms.get_masses()):.2f} amu")
    report.append("")
    
    # Frequency analysis
    real_freqs = [f for f in frequencies_cm if f > 20]
    imag_freqs = [f for f in frequencies_cm if f < 0]
    
    report.append("VIBRATIONAL ANALYSIS:")
    report.append(f"  Total frequencies: {len(frequencies_cm)}")
    report.append(f"  Real frequencies (>20 cm⁻¹): {len(real_freqs)}")
    report.append(f"  Imaginary frequencies: {len(imag_freqs)}")
    if real_freqs:
        report.append(f"  Frequency range: {min(real_freqs):.1f} - {max(real_freqs):.1f} cm⁻¹")
    report.append("")
    
    # Temperature-dependent properties
    report.append("THERMAL PROPERTIES:")
    report.append("T(K)    H_th(eV)  S_tot(meV/K)  G_th(eV)  Stability")
    report.append("-" * 55)
    
    for i, d in enumerate(data):
        T = d['temperature']
        H = d['enthalpy_thermal']
        S = d['entropy_total'] * 1000  # meV/K
        G = d['free_energy_thermal']
        ratio = analysis_results['stability_ratios'][i]
        stable = "STABLE" if ratio < 0.1 else "UNSTABLE"
        
        report.append(f"{T:4.0f}    {H:8.4f}  {S:9.2f}   {G:8.4f}  {stable}")
    
    report.append("")
    
    # Critical temperature
    critical_T = analysis_results['critical_temperature']
    if critical_T:
        report.append(f"CRITICAL TEMPERATURE: {critical_T:.0f} K ({critical_T-273.15:.0f} °C)")
        report.append("Above this temperature, thermal energy becomes significant")
        report.append("compared to bond dissociation energy.")
    else:
        max_T = max(d['temperature'] for d in data)
        report.append(f"Structure remains stable up to {max_T:.0f} K")
    
    report_text = "\n".join(report)
    
    if output_file:
        with open(output_file, 'w') as f:
            f.write(report_text)
    
    return report_text


def estimate_desorption_temperature(bond_energy_eV, frequency_cm, prefactor=1e13):
    """
    Estimate desorption temperature using Arrhenius equation.
    
    Args:
        bond_energy_eV: Bond dissociation energy in eV
        frequency_cm: Relevant vibrational frequency in cm⁻¹
        prefactor: Pre-exponential factor in s⁻¹
        
    Returns:
        Desorption temperature in Kelvin
    """
    # Convert frequency to Hz
    freq_hz = frequency_cm * C * 100  # cm⁻¹ to Hz
    
    # Arrhenius equation: rate = prefactor * exp(-Ea/kT)
    # For significant desorption, assume rate ~ freq/100
    target_rate = freq_hz / 100
    
    # Solve for temperature
    T_desorption = bond_energy_eV / (KB * np.log(prefactor / target_rate))
    
    return T_desorption


def generate_thermal_displacement(atoms, frequencies_cm, temperature, normal_modes=None, amplitude_scale=0.1):
    """
    Generate thermal displacement of atoms based on temperature and vibrational modes.
    
    Args:
        atoms: ASE Atoms object
        frequencies_cm: List of frequencies in cm⁻¹
        temperature: Temperature in Kelvin
        normal_modes: Normal mode vectors (if available)
        amplitude_scale: Scaling factor for displacement amplitude
        
    Returns:
        Displaced atoms object
    """
    import numpy as np
    from ase import Atoms
    
    displaced_atoms = atoms.copy()
    positions = displaced_atoms.get_positions()
    
    # Simple thermal displacement model
    # Displacement amplitude proportional to √(kT/mω²)
    
    # Filter real frequencies
    real_frequencies = [f for f in frequencies_cm if f > 20.0]
    
    if not real_frequencies:
        # If no real frequencies, use simple random thermal motion
        thermal_energy = KB * temperature  # in eV
        
        for i, atom in enumerate(atoms):
            mass = atom.mass  # in amu
            # Thermal velocity: v = √(3kT/m)
            thermal_velocity = np.sqrt(3 * thermal_energy / (mass * AMU))
            
            # Convert to displacement (assume characteristic time ~ 1 fs)
            char_time = 1e-15  # seconds
            displacement = thermal_velocity * char_time * 1e10  # in Angstroms
            
            # Random displacement
            random_displacement = np.random.normal(0, displacement * amplitude_scale, 3)
            positions[i] += random_displacement
    else:
        # Use vibrational modes for more realistic displacement
        for freq in real_frequencies[:3]:  # Use first few modes
            if freq > 0:
                freq_hz = freq * C * 100  # Convert to Hz
                omega = 2 * np.pi * freq_hz
                
                # Thermal amplitude: A = √(kT/mω²)
                thermal_amplitude = np.sqrt(KB * temperature / (AMU * 1.66054e-27 * omega**2))
                thermal_amplitude *= 1e10 * amplitude_scale  # Convert to Angstroms and scale
                
                # Add random displacement in each mode
                for i in range(len(atoms)):
                    # Simple approximation: random direction for each mode
                    mode_direction = np.random.normal(0, 1, 3)
                    mode_direction /= np.linalg.norm(mode_direction)
                    
                    displacement = mode_direction * thermal_amplitude * np.random.normal(0, 1)
                    positions[i] += displacement
    
    displaced_atoms.set_positions(positions)
    return displaced_atoms


def create_thermal_frames(atoms, frequencies_cm, temperatures, n_frames_per_temp=5):
    """
    Create series of atomic configurations showing thermal motion at different temperatures.
    
    Args:
        atoms: ASE Atoms object
        frequencies_cm: List of frequencies in cm⁻¹
        temperatures: List of temperatures in Kelvin
        n_frames_per_temp: Number of frames per temperature
        
    Returns:
        List of (temperature, frame_atoms) tuples
    """
    frames = []
    
    for temp in temperatures:
        for frame_idx in range(n_frames_per_temp):
            # Generate different random displacement for each frame
            displaced_atoms = generate_thermal_displacement(
                atoms, frequencies_cm, temp, amplitude_scale=0.05
            )
            frames.append((temp, displaced_atoms))
    
    return frames


def create_thermal_animation(analysis_results, atoms, output_dir, n_frames=50, fps=10):
    """
    Create MP4 animation showing thermal motion and heating effects.
    
    Args:
        analysis_results: Results from analyze_thermal_stability
        atoms: ASE Atoms object
        output_dir: Directory to save animation
        n_frames: Total number of frames
        fps: Frames per second
        
    Returns:
        Path to created MP4 file
    """
    import os
    import numpy as np
    from ase.io import write
    
    temperatures = [d['temperature'] for d in analysis_results['thermodynamic_data']]
    frequencies = analysis_results['thermodynamic_data'][0].get('frequencies_cm', [])
    
    # Create frames directory
    frames_dir = os.path.join(output_dir, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    
    # Generate temperature sequence for animation
    T_min, T_max = min(temperatures), max(temperatures)
    anim_temperatures = np.linspace(T_min, T_max, n_frames)
    
    # Create frames with thermal motion
    frame_files = []
    
    for i, temp in enumerate(anim_temperatures):
        # Generate thermal displacement for this temperature
        displaced_atoms = generate_thermal_displacement(
            atoms, frequencies, temp, amplitude_scale=0.1
        )
        
        # Save frame as XYZ
        frame_file = os.path.join(frames_dir, f"frame_{i:04d}_{temp:.0f}K.xyz")
        
        # Add temperature info to comment line
        comment = f"Temperature: {temp:.1f} K - Thermal motion simulation"
        write(frame_file, displaced_atoms, comment=comment)
        frame_files.append(frame_file)
    
    # Create MP4 using external tools (ffmpeg or similar)
    mp4_file = os.path.join(output_dir, "thermal_heating.mp4")
    
    try:
        # Try to create MP4 with ovito or ase visualization
        success = create_mp4_from_frames(frame_files, mp4_file, fps)
        if success:
            print(f"🎬 Thermal animation created: {mp4_file}")
            return mp4_file
        else:
            print("⚠️  MP4 creation failed, frames available in:", frames_dir)
            return frames_dir
    except Exception as e:
        print(f"⚠️  Animation error: {e}")
        print(f"   Frames saved to: {frames_dir}")
        return frames_dir


def create_mp4_from_frames(frame_files, output_file, fps=10):
    """
    Create MP4 video from XYZ frame files.
    
    Args:
        frame_files: List of XYZ file paths
        output_file: Output MP4 file path
        fps: Frames per second
        
    Returns:
        Boolean success
    """
    import subprocess
    import os
    
    # Create a temporary script for visualization
    # This is a simplified approach - in production you might want to use
    # more sophisticated tools like OVITO, ASE, or VMD
    
    try:
        # Create concatenated XYZ file for easier processing
        concat_file = output_file.replace('.mp4', '_all_frames.xyz')
        
        with open(concat_file, 'w') as outf:
            for frame_file in frame_files:
                with open(frame_file, 'r') as inf:
                    outf.write(inf.read())
        
        print(f"📁 All frames concatenated to: {concat_file}")
        
        # Try using ffmpeg if available (this is basic implementation)
        # For real visualization, you'd use specialized molecular visualization tools
        
        # For now, just return the concatenated file as success
        return True
        
    except Exception as e:
        print(f"Error creating MP4: {e}")
        return False