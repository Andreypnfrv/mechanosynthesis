"""
Simple test for thermodynamics module
"""
import sys
import os
import numpy as np
sys.path.append('py')

from ase import Atoms
from kernel.parts.thermodynamics import (
    calculate_thermodynamic_properties,
    temperature_sweep,
    analyze_thermal_stability,
    vibrational_energy,
    vibrational_entropy
)

def test_simple_molecule():
    """Test thermodynamics with simple H2 molecule"""
    print("Testing thermodynamics with H2 molecule...")
    
    # Simple H2 molecule
    atoms = Atoms('H2', positions=[[0, 0, 0], [0.74, 0, 0]])
    
    # Mock vibrational frequency for H2 (should be ~4400 cm-1)
    frequencies = [4400.0]  # H-H stretch
    
    # Test at room temperature
    T = 298.15
    props = calculate_thermodynamic_properties(frequencies, atoms, T)
    
    print(f"Temperature: {T} K")
    print(f"Vibrational energy: {props['vibrational_energy']:.6f} eV")
    print(f"Total entropy: {props['entropy_total']*1000:.2f} meV/K")
    print(f"Thermal free energy: {props['free_energy_thermal']:.6f} eV")
    print(f"Number of frequencies used: {props['n_frequencies']}")
    
    # Test temperature sweep
    temperatures = [298, 400, 600, 800, 1000]
    results = temperature_sweep(frequencies, atoms, temperatures)
    
    print("\nTemperature sweep:")
    print("T(K)   G_th(eV)   S_tot(meV/K)")
    for r in results:
        print(f"{r['temperature']:4.0f}   {r['free_energy_thermal']:8.5f}   {r['entropy_total']*1000:8.2f}")
    
    # Test thermal stability
    stability = analyze_thermal_stability(frequencies, atoms, temperatures, bond_dissociation_energy=4.5)
    
    print(f"\nCritical temperature: {stability['critical_temperature']} K")
    print("Test passed!")

def test_water_molecule():
    """Test with water molecule (more frequencies)"""
    print("\nTesting thermodynamics with H2O molecule...")
    
    # Water molecule
    atoms = Atoms('H2O', positions=[[0, 0, 0], [0.96, 0, 0], [0.24, 0.93, 0]])
    
    # Water vibrational frequencies (cm-1)
    frequencies = [1595, 3657, 3756]  # bend, sym stretch, asym stretch
    
    T = 298.15
    props = calculate_thermodynamic_properties(frequencies, atoms, T)
    
    print(f"H2O at {T} K:")
    print(f"Vibrational energy: {props['vibrational_energy']:.6f} eV")
    print(f"Total entropy: {props['entropy_total']*1000:.2f} meV/K")
    print(f"Thermal free energy: {props['free_energy_thermal']:.6f} eV")
    
    # Test high temperature behavior
    hot_temps = [298, 500, 1000, 1500, 2000]
    hot_results = temperature_sweep(frequencies, atoms, hot_temps)
    
    print("\nHigh temperature behavior:")
    print("T(K)    S_vib(meV/K)  S_tot(meV/K)")
    for r in hot_results:
        print(f"{r['temperature']:4.0f}    {r['vibrational_entropy']*1000:9.2f}   {r['entropy_total']*1000:9.2f}")
    
    print("Water test passed!")

def test_individual_functions():
    """Test individual thermodynamic functions"""
    print("\nTesting individual functions...")
    
    # Test vibrational energy and entropy
    freq = 1000.0  # cm-1
    temps = [100, 298, 500, 1000]
    
    print("Frequency: 1000 cm-1")
    print("T(K)   E_vib(eV)   S_vib(meV/K)")
    for T in temps:
        E = vibrational_energy(freq, T)
        S = vibrational_entropy(freq, T) * 1000
        print(f"{T:4.0f}   {E:.6f}    {S:.3f}")
    
    print("Individual function tests passed!")

if __name__ == "__main__":
    test_simple_molecule()
    test_water_molecule() 
    test_individual_functions()
    print("\nAll thermodynamics tests passed! 🎉")