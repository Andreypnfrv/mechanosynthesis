"""
Structural analysis module for mechanosynthesis research.

Provides detailed analysis of structural changes during relaxation:
- Bond type analysis (C-C, C-H, Si-C, W-C, etc.)
- Bond formation/breaking statistics
- Average bond length changes
- Valence angle analysis
- Coordination number analysis
"""

import os
import numpy as np
import json
from typing import Dict, List, Tuple, Optional, Any
from collections import defaultdict, Counter
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from ase import Atoms
from ase.io import read, write
from ase.neighborlist import NeighborList
from ase.data import covalent_radii, atomic_numbers
# import pandas as pd  # Optional dependency - only needed for advanced CSV export


class StructuralAnalyzer:
    """Main class for structural analysis"""
    
    def __init__(self, cutoff_factor: float = 1.2):
        """
        Initialize structural analyzer
        
        Args:
            cutoff_factor: Multiplier for covalent radii to determine bonds
        """
        self.cutoff_factor = cutoff_factor
        
    def get_bond_cutoffs(self, atoms: Atoms) -> List[float]:
        """Get bond cutoff distances for each atom"""
        cutoffs = []
        for symbol in atoms.get_chemical_symbols():
            atomic_num = atomic_numbers[symbol]
            cutoffs.append(covalent_radii[atomic_num] * self.cutoff_factor)
        return cutoffs
    
    def analyze_bonds(self, atoms: Atoms) -> Dict[str, Any]:
        """
        Analyze all bonds in the structure
        
        Returns:
            Dictionary with bond information including types, lengths, and connectivity
        """
        cutoffs = self.get_bond_cutoffs(atoms)
        nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
        nl.update(atoms)
        
        bonds = []
        bond_types = defaultdict(list)
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        for i in range(len(atoms)):
            indices, offsets = nl.get_neighbors(i)
            for j in indices:
                if i < j:  # Avoid double counting
                    # Calculate bond length
                    bond_vector = positions[j] - positions[i]
                    bond_length = np.linalg.norm(bond_vector)
                    
                    # Determine bond type
                    atom1, atom2 = symbols[i], symbols[j]
                    bond_type = f"{min(atom1, atom2)}-{max(atom1, atom2)}"
                    
                    bond_info = {
                        'atoms': (i, j),
                        'symbols': (atom1, atom2),
                        'type': bond_type,
                        'length': bond_length,
                        'vector': bond_vector
                    }
                    
                    bonds.append(bond_info)
                    bond_types[bond_type].append(bond_length)
        
        return {
            'bonds': bonds,
            'bond_types': dict(bond_types),
            'total_bonds': len(bonds)
        }
    
    def analyze_angles(self, atoms: Atoms, bonds: List[Dict]) -> Dict[str, Any]:
        """
        Analyze valence angles in the structure
        
        Args:
            atoms: ASE atoms object
            bonds: List of bond dictionaries from analyze_bonds
            
        Returns:
            Dictionary with angle information
        """
        positions = atoms.get_positions()
        symbols = atoms.get_chemical_symbols()
        
        # Build adjacency list from bonds
        adjacency = defaultdict(list)
        for bond in bonds:
            i, j = bond['atoms']
            adjacency[i].append(j)
            adjacency[j].append(i)
        
        angles = []
        angle_types = defaultdict(list)
        
        for center_atom in range(len(atoms)):
            neighbors = adjacency[center_atom]
            if len(neighbors) < 2:
                continue
                
            # Calculate all angles with center_atom as vertex
            for i, neighbor1 in enumerate(neighbors):
                for neighbor2 in neighbors[i+1:]:
                    # Calculate angle
                    vec1 = positions[neighbor1] - positions[center_atom]
                    vec2 = positions[neighbor2] - positions[center_atom]
                    
                    cos_angle = np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2))
                    # Clamp to valid range to avoid numerical errors
                    cos_angle = np.clip(cos_angle, -1.0, 1.0)
                    angle = np.degrees(np.arccos(cos_angle))
                    
                    # Determine angle type
                    center_symbol = symbols[center_atom]
                    neighbor1_symbol = symbols[neighbor1]
                    neighbor2_symbol = symbols[neighbor2]
                    
                    # Sort neighbor symbols for consistent naming
                    sorted_neighbors = sorted([neighbor1_symbol, neighbor2_symbol])
                    angle_type = f"{sorted_neighbors[0]}-{center_symbol}-{sorted_neighbors[1]}"
                    
                    angle_info = {
                        'atoms': (neighbor1, center_atom, neighbor2),
                        'symbols': (neighbor1_symbol, center_symbol, neighbor2_symbol),
                        'type': angle_type,
                        'angle': angle
                    }
                    
                    angles.append(angle_info)
                    angle_types[angle_type].append(angle)
        
        return {
            'angles': angles,
            'angle_types': dict(angle_types),
            'total_angles': len(angles)
        }
    
    def analyze_coordination(self, atoms: Atoms, bonds: List[Dict]) -> Dict[str, Any]:
        """
        Analyze coordination numbers for each atom type
        
        Args:
            atoms: ASE atoms object  
            bonds: List of bond dictionaries from analyze_bonds
            
        Returns:
            Dictionary with coordination analysis
        """
        symbols = atoms.get_chemical_symbols()
        coordination = defaultdict(int)
        
        # Count bonds for each atom
        for bond in bonds:
            i, j = bond['atoms']
            coordination[i] += 1
            coordination[j] += 1
        
        # Group by element type
        element_coordination = defaultdict(list)
        for atom_idx, coord_num in coordination.items():
            element = symbols[atom_idx]
            element_coordination[element].append(coord_num)
        
        # Calculate statistics
        coord_stats = {}
        for element, coord_nums in element_coordination.items():
            coord_stats[element] = {
                'mean': np.mean(coord_nums),
                'std': np.std(coord_nums),
                'min': min(coord_nums),
                'max': max(coord_nums),
                'distribution': Counter(coord_nums)
            }
        
        return {
            'atom_coordination': dict(coordination),
            'element_coordination': dict(element_coordination),
            'coordination_stats': coord_stats
        }
    
    def compare_structures(self, initial_atoms: Atoms, final_atoms: Atoms) -> Dict[str, Any]:
        """
        Compare two structures and analyze changes
        
        Args:
            initial_atoms: Initial structure
            final_atoms: Final structure after relaxation
            
        Returns:
            Comprehensive comparison analysis
        """
        # Analyze both structures
        initial_analysis = {
            'bonds': self.analyze_bonds(initial_atoms),
            'angles': self.analyze_angles(initial_atoms, self.analyze_bonds(initial_atoms)['bonds']),
            'coordination': self.analyze_coordination(initial_atoms, self.analyze_bonds(initial_atoms)['bonds'])
        }
        
        final_analysis = {
            'bonds': self.analyze_bonds(final_atoms),
            'angles': self.analyze_angles(final_atoms, self.analyze_bonds(final_atoms)['bonds']),
            'coordination': self.analyze_coordination(final_atoms, self.analyze_bonds(final_atoms)['bonds'])
        }
        
        # Compare bond types and lengths
        bond_changes = self._compare_bonds(
            initial_analysis['bonds']['bond_types'], 
            final_analysis['bonds']['bond_types']
        )
        
        # Compare angles
        angle_changes = self._compare_angles(
            initial_analysis['angles']['angle_types'],
            final_analysis['angles']['angle_types'] 
        )
        
        # Compare coordination
        coord_changes = self._compare_coordination(
            initial_analysis['coordination']['coordination_stats'],
            final_analysis['coordination']['coordination_stats']
        )
        
        # Calculate RMSD
        rmsd = self._calculate_rmsd(initial_atoms, final_atoms)
        
        return {
            'initial': initial_analysis,
            'final': final_analysis,
            'changes': {
                'bonds': bond_changes,
                'angles': angle_changes,
                'coordination': coord_changes,
                'rmsd': rmsd
            }
        }
    
    def _compare_bonds(self, initial_bonds: Dict, final_bonds: Dict) -> Dict[str, Any]:
        """Compare bond types between structures"""
        all_bond_types = set(initial_bonds.keys()) | set(final_bonds.keys())
        
        bond_changes = {}
        for bond_type in all_bond_types:
            initial_lengths = initial_bonds.get(bond_type, [])
            final_lengths = final_bonds.get(bond_type, [])
            
            initial_count = len(initial_lengths)
            final_count = len(final_lengths)
            count_change = final_count - initial_count
            
            initial_avg = np.mean(initial_lengths) if initial_lengths else 0
            final_avg = np.mean(final_lengths) if final_lengths else 0
            length_change = final_avg - initial_avg if initial_lengths and final_lengths else None
            
            bond_changes[bond_type] = {
                'initial_count': initial_count,
                'final_count': final_count, 
                'count_change': count_change,
                'initial_avg_length': initial_avg,
                'final_avg_length': final_avg,
                'avg_length_change': length_change,
                'initial_lengths': initial_lengths,
                'final_lengths': final_lengths
            }
        
        return bond_changes
    
    def _compare_angles(self, initial_angles: Dict, final_angles: Dict) -> Dict[str, Any]:
        """Compare angle types between structures"""
        all_angle_types = set(initial_angles.keys()) | set(final_angles.keys())
        
        angle_changes = {}
        for angle_type in all_angle_types:
            initial_vals = initial_angles.get(angle_type, [])
            final_vals = final_angles.get(angle_type, [])
            
            initial_count = len(initial_vals)
            final_count = len(final_vals)
            count_change = final_count - initial_count
            
            initial_avg = np.mean(initial_vals) if initial_vals else 0
            final_avg = np.mean(final_vals) if final_vals else 0
            angle_change = final_avg - initial_avg if initial_vals and final_vals else None
            
            angle_changes[angle_type] = {
                'initial_count': initial_count,
                'final_count': final_count,
                'count_change': count_change, 
                'initial_avg_angle': initial_avg,
                'final_avg_angle': final_avg,
                'avg_angle_change': angle_change,
                'initial_angles': initial_vals,
                'final_angles': final_vals
            }
        
        return angle_changes
    
    def _compare_coordination(self, initial_coord: Dict, final_coord: Dict) -> Dict[str, Any]:
        """Compare coordination numbers between structures"""
        all_elements = set(initial_coord.keys()) | set(final_coord.keys())
        
        coord_changes = {}
        for element in all_elements:
            initial_stats = initial_coord.get(element, {'mean': 0})
            final_stats = final_coord.get(element, {'mean': 0})
            
            coord_changes[element] = {
                'initial_avg_coordination': initial_stats['mean'],
                'final_avg_coordination': final_stats['mean'],
                'coordination_change': final_stats['mean'] - initial_stats['mean'],
                'initial_distribution': initial_stats.get('distribution', {}),
                'final_distribution': final_stats.get('distribution', {})
            }
        
        return coord_changes
    
    def _calculate_rmsd(self, atoms1: Atoms, atoms2: Atoms) -> float:
        """Calculate RMSD between two structures"""
        if len(atoms1) != len(atoms2):
            return float('inf')
        
        pos1 = atoms1.get_positions()
        pos2 = atoms2.get_positions()
        
        # Center both structures
        pos1_centered = pos1 - np.mean(pos1, axis=0)
        pos2_centered = pos2 - np.mean(pos2, axis=0)
        
        diff = pos2_centered - pos1_centered
        rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
        
        return rmsd


def analyze_structural_changes(initial_structure_path: str, 
                             final_structure_path: str,
                             output_dir: Optional[str] = None) -> Dict[str, Any]:
    """
    Main function to analyze structural changes between two structures
    
    Args:
        initial_structure_path: Path to initial structure file
        final_structure_path: Path to final structure file  
        output_dir: Directory to save analysis results
        
    Returns:
        Complete structural analysis
    """
    # Read structures
    initial_atoms = read(initial_structure_path)
    final_atoms = read(final_structure_path)
    
    # Initialize analyzer
    analyzer = StructuralAnalyzer()
    
    # Perform analysis
    analysis = analyzer.compare_structures(initial_atoms, final_atoms)
    
    # Save results if output directory provided
    if output_dir:
        save_analysis_results(analysis, output_dir, 
                            Path(initial_structure_path).stem,
                            Path(final_structure_path).stem)
    
    return analysis


def save_analysis_results(analysis: Dict[str, Any], 
                         output_dir: str,
                         initial_name: str,
                         final_name: str):
    """
    Save analysis results to files
    
    Args:
        analysis: Analysis results dictionary
        output_dir: Output directory path
        initial_name: Name of initial structure
        final_name: Name of final structure
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Save JSON results
    json_path = output_path / f"structural_analysis_{initial_name}_to_{final_name}.json"
    with open(json_path, 'w') as f:
        # Convert numpy arrays to lists for JSON serialization
        json_analysis = _prepare_for_json(analysis)
        json.dump(json_analysis, f, indent=2)
    
    # Generate text report
    report_path = output_path / f"structural_report_{initial_name}_to_{final_name}.txt"
    with open(report_path, 'w') as f:
        f.write(generate_text_report(analysis, initial_name, final_name))
    
    # Generate visualizations
    create_analysis_plots(analysis, output_path, initial_name, final_name)
    
    print(f"Analysis results saved to {output_dir}")


def _prepare_for_json(obj):
    """Recursively prepare object for JSON serialization"""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, (np.integer, np.int64, np.int32)):
        return int(obj)
    elif isinstance(obj, (np.floating, np.float64, np.float32)):
        return float(obj)
    elif isinstance(obj, dict):
        return {str(key): _prepare_for_json(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [_prepare_for_json(item) for item in obj]
    elif isinstance(obj, tuple):
        return [_prepare_for_json(item) for item in obj]
    elif hasattr(obj, 'dtype'):  # Handle other numpy types
        return obj.item() if hasattr(obj, 'item') else str(obj)
    else:
        return obj


def generate_text_report(analysis: Dict[str, Any], 
                        initial_name: str,
                        final_name: str) -> str:
    """Generate a human-readable text report"""
    report = []
    report.append(f"STRUCTURAL ANALYSIS REPORT")
    report.append(f"=" * 50)
    report.append(f"Initial structure: {initial_name}")
    report.append(f"Final structure: {final_name}")
    report.append(f"RMSD: {analysis['changes']['rmsd']:.3f} Å")
    report.append("")
    
    # Bond analysis
    report.append("BOND ANALYSIS")
    report.append("-" * 20)
    bond_changes = analysis['changes']['bonds']
    
    for bond_type, changes in bond_changes.items():
        if changes['count_change'] != 0 or changes['avg_length_change'] is not None:
            report.append(f"\n{bond_type} bonds:")
            report.append(f"  Count: {changes['initial_count']} → {changes['final_count']} ({changes['count_change']:+d})")
            
            if changes['avg_length_change'] is not None:
                report.append(f"  Avg length: {changes['initial_avg_length']:.3f} → {changes['final_avg_length']:.3f} Å ({changes['avg_length_change']:+.3f} Å)")
    
    # Angle analysis  
    report.append("\n\nANGLE ANALYSIS")
    report.append("-" * 20)
    angle_changes = analysis['changes']['angles']
    
    for angle_type, changes in angle_changes.items():
        if changes['count_change'] != 0 or changes['avg_angle_change'] is not None:
            report.append(f"\n{angle_type} angles:")
            report.append(f"  Count: {changes['initial_count']} → {changes['final_count']} ({changes['count_change']:+d})")
            
            if changes['avg_angle_change'] is not None:
                report.append(f"  Avg angle: {changes['initial_avg_angle']:.1f}° → {changes['final_avg_angle']:.1f}° ({changes['avg_angle_change']:+.1f}°)")
    
    # Coordination analysis
    report.append("\n\nCOORDINATION ANALYSIS") 
    report.append("-" * 20)
    coord_changes = analysis['changes']['coordination']
    
    for element, changes in coord_changes.items():
        if abs(changes['coordination_change']) > 0.01:
            report.append(f"\n{element} atoms:")
            report.append(f"  Avg coordination: {changes['initial_avg_coordination']:.2f} → {changes['final_avg_coordination']:.2f} ({changes['coordination_change']:+.2f})")
    
    return "\n".join(report)


def create_analysis_plots(analysis: Dict[str, Any], 
                         output_path: Path,
                         initial_name: str, 
                         final_name: str):
    """Create visualization plots for the analysis"""
    
    # Plot 1: Bond length changes
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    bond_changes = analysis['changes']['bonds']
    bond_types = []
    length_changes = []
    
    for bond_type, changes in bond_changes.items():
        if changes['avg_length_change'] is not None:
            bond_types.append(bond_type)
            length_changes.append(changes['avg_length_change'])
    
    if bond_types:
        colors = ['red' if x < 0 else 'green' for x in length_changes]
        ax1.barh(bond_types, length_changes, color=colors, alpha=0.7)
        ax1.set_xlabel('Average Bond Length Change (Å)')
        ax1.set_title('Bond Length Changes')
        ax1.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    
    # Plot 2: Angle changes
    angle_changes = analysis['changes']['angles'] 
    angle_types = []
    angle_deltas = []
    
    for angle_type, changes in angle_changes.items():
        if changes['avg_angle_change'] is not None:
            angle_types.append(angle_type)
            angle_deltas.append(changes['avg_angle_change'])
    
    if angle_types:
        colors = ['red' if x < 0 else 'green' for x in angle_deltas]
        ax2.barh(angle_types, angle_deltas, color=colors, alpha=0.7)
        ax2.set_xlabel('Average Angle Change (°)')
        ax2.set_title('Valence Angle Changes')
        ax2.axvline(x=0, color='black', linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path / f"structural_changes_{initial_name}_to_{final_name}.png", 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Bond count changes
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    
    bond_types_count = []
    count_changes = []
    
    for bond_type, changes in bond_changes.items():
        if changes['count_change'] != 0:
            bond_types_count.append(bond_type)
            count_changes.append(changes['count_change'])
    
    if bond_types_count:
        colors = ['red' if x < 0 else 'green' for x in count_changes]
        ax.barh(bond_types_count, count_changes, color=colors, alpha=0.7)
        ax.set_xlabel('Bond Count Change')
        ax.set_title('Bond Formation/Breaking')
        ax.axvline(x=0, color='black', linestyle='--', alpha=0.5)
        
        # Add count labels
        for i, (bond_type, count) in enumerate(zip(bond_types_count, count_changes)):
            initial = bond_changes[bond_type]['initial_count']
            final = bond_changes[bond_type]['final_count'] 
            ax.text(count + 0.1 if count > 0 else count - 0.1, i, 
                   f"{initial}→{final}", ha='left' if count > 0 else 'right', va='center')
    
    plt.tight_layout()
    plt.savefig(output_path / f"bond_count_changes_{initial_name}_to_{final_name}.png",
                dpi=300, bbox_inches='tight')
    plt.close()


def print_structural_summary(analysis: Dict[str, Any]):
    """Print a brief summary of structural changes to console"""
    changes = analysis['changes']
    
    print(f"\n{'='*50}")
    print("STRUCTURAL ANALYSIS SUMMARY")
    print(f"{'='*50}")
    print(f"RMSD: {changes['rmsd']:.3f} Å")
    
    # Count significant changes
    significant_bond_changes = sum(1 for changes in changes['bonds'].values() 
                                 if changes['count_change'] != 0)
    significant_angle_changes = sum(1 for changes in changes['angles'].values()
                                  if changes['avg_angle_change'] is not None and abs(changes['avg_angle_change']) > 1.0)
    
    print(f"Bond types with count changes: {significant_bond_changes}")
    print(f"Angle types with significant changes: {significant_angle_changes}")
    
    # Show most significant changes
    if changes['bonds']:
        max_bond_change = max(changes['bonds'].items(), 
                            key=lambda x: abs(x[1]['avg_length_change']) if x[1]['avg_length_change'] else 0)
        if max_bond_change[1]['avg_length_change']:
            print(f"Largest bond length change: {max_bond_change[0]} ({max_bond_change[1]['avg_length_change']:+.3f} Å)")
    
    print(f"{'='*50}")