#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import shutil
import re
from datetime import datetime
from pathlib import Path
import settings
from kernel.project_structure import create_project_structure
from kernel.parts.relax import relax, relax_native_dftb, relax_passivated, relax_xtb, relax_gxtb, relax_native_dftb_ptbp, relax_dftb_unified
# Temporarily commented out due to indentation issues: relax_orca_native, relax_orca_dft
from kernel.parts.hessian import calculate_hessian_orca_standalone, calculate_hessian_gxtb_standalone, calculate_hessian_xtb_standalone, calculate_hessian_with_method, calculate_hessian_dftb_native, calculate_hessian_dftb_unified, escape_saddle_point, count_imaginary_frequencies, rename_relaxed_result, rename_hessian_result, finalize_results, escape_saddle_point_with_iteration
from kernel.build.build_molecule import build_molecule
from kernel.build.build_diamond import build_diamond_with_highlighting
from kernel.build.build_diamond_supercell import build_diamond_supercell
from kernel.build.build_diamond_surface import build_diamond_surface
from kernel.build.build_diamondoid import build_diamondoid
from kernel.build.build_sdf import build_from_sdf
from kernel.build.build_stm_tip import build_tungsten_stm_tip
from kernel.build.si.surface import build_si_surface
from kernel.utils.performance_tracker import PerformanceTracker
from kernel.sequences.scene import SceneGenerator
from ase.io import read
import glob

def get_project_info_for_tracking(project_folder):
    """Get project information for performance tracking"""
    base_path = f"data/{project_folder}"
    
    # Find project folder
    if not os.path.exists(base_path):
        pattern = f"data/*_{project_folder}"
        matching_folders = glob.glob(pattern)
        
        if len(matching_folders) == 1:
            base_path = matching_folders[0]
        elif len(matching_folders) > 1:
            # Take the first match for tracking purposes
            base_path = matching_folders[0]
        else:
            # Return basic info if folder not found
            return project_folder, None, 0, 0
    
    # Find XYZ files
    xyz_files = []
    xyz_pattern = os.path.join(base_path, "*.xyz")
    xyz_files.extend(glob.glob(xyz_pattern))
    
    inputs_pattern = os.path.join(base_path, "inputs", "*.xyz")
    xyz_files.extend(glob.glob(inputs_pattern))
    
    if not xyz_files:
        return project_folder, None, 0, 0
    
    # Read the first molecule for basic info
    try:
        atoms = read(xyz_files[0])
        n_atoms = len(atoms)
        n_electrons = sum(atoms.get_atomic_numbers())
        return project_folder, atoms, n_atoms, n_electrons
    except:
        return project_folder, None, 0, 0

def get_project_name(base_name, structure_type=None, size=None):
    """Generate project name without date prefix, including size for relevant structures"""
    if structure_type in ['supercell', 'surface'] and size:
        return f"{base_name}_{structure_type}_{size}"
    elif structure_type == 'crystal' and size:
        return f"{base_name}_{size}"
    return base_name

def add_common_args(parser):
    """Add common arguments to subparsers"""
    parser.add_argument('--method', default='B3LYP', 
                       help='Method for calculations (default: B3LYP for Orca)')
    parser.add_argument('--basis', default='def2-SVP',
                       help='Basis set for Orca (default: def2-SVP)')
    parser.add_argument('--nprocs', type=int, default=4,
                       help='Number of processors (default: 4)')
    parser.add_argument('--mem', default='2GB',
                       help='Memory allocation (default: 2GB)')

def add_xtb_args(parser):
    """Add xTB-specific arguments"""
    parser.add_argument('--xtb-method', default='GFN2-xTB',
                       choices=['GFN2-xTB', 'GFN1-xTB', 'GFN0-xTB'],
                       help='xTB method (default: GFN2-xTB)')
    parser.add_argument('--xtb-backend', default='fixed',
                       choices=['native', 'gxtb', 'fixed'],
                       help='xTB backend (default: fixed - our patched version)')

def add_build_args(parser):
    """Add build-specific arguments"""
    parser.add_argument('--highlight', action='store_true',
                       help='Highlight atoms with hanging bonds')
    parser.add_argument('--passivate', metavar='ELEMENT',
                       help='Passivate hanging bonds with element (e.g., H)')
    parser.add_argument('--output', metavar='FILENAME',
                       help='Output filename (auto-generated if not specified)')
    parser.add_argument('--relax', action='store_true',
                       help='Run relaxation after building')
    parser.add_argument('--hessian', action='store_true',
                       help='Run Hessian calculation after relaxation')
    parser.add_argument('--continue-from-saddle', action='store_true',
                       help='Continue optimization from saddle point (requires prior Hessian calculation)')
    parser.add_argument('--relax-hessian', action='store_true',
                       help='Automated relax-hessian workflow until true minimum is found')
    parser.add_argument('--max-iterations', type=int, default=5,
                       help='Maximum iterations for saddle point escape (default: 5)')
    parser.add_argument('--backend', default='xtb',
                       choices=['dftb', 'dftb-native', 'dftb-passivated', 'dftb-ptbp', 'xtb', 'gxtb', 'orca'],
                       help='Method for both relaxation and Hessian (default: xtb). Use "dftb" for unified auto-selection.')
    parser.add_argument('--relax-backend',
                       choices=['dftb', 'dftb-native', 'dftb-passivated', 'dftb-ptbp', 'xtb', 'gxtb'],
                       help='Override relaxation method (overrides --backend). Use "dftb" for unified auto-selection.')
    parser.add_argument('--dftb-params', default='auto',
                       choices=['auto', 'ptbp', '3ob'],
                       help='DFTB parameter set selection (default: auto). Only used with --backend dftb.')
    parser.add_argument('--hessian-backend',
                       choices=['dftb', 'xtb', 'gxtb', 'orca', 'orca-simple'],
                       help='Override Hessian method (overrides --backend)')

def setup_parser():
    """Setup the main argument parser with subcommands"""
    parser = argparse.ArgumentParser(description='Mechanosynthesis pipeline management')
    subparsers = parser.add_subparsers(dest='command', help='Available commands')

    # CREATE subcommand
    create_parser = subparsers.add_parser('create', help='Create new project')
    create_parser.add_argument('name', help='Project name')

    # RELAX subcommand
    relax_parser = subparsers.add_parser('relax', help='Run geometry optimization')
    relax_parser.add_argument('project', help='Project folder name')
    relax_parser.add_argument('--backend', default='xtb',
                             choices=['dftb', 'dftb-native', 'dftb-passivated', 'dftb-ptbp', 'xtb', 'gxtb', 'orca', 'orca-ase'],
                             help='Relaxation backend (default: xtb)')
    relax_parser.add_argument('--dftb-params', default='auto',
                             choices=['auto', 'ptbp', '3ob'],
                             help='DFTB parameter set selection (default: auto). Only used with --backend dftb.')
    relax_parser.add_argument('--calc-hessian', action='store_true',
                             help='Calculate Hessian after relaxation')
    relax_parser.add_argument('--hessian-method', choices=['phonopy', 'orca'], default='phonopy',
                             help='Hessian calculation method')
    add_common_args(relax_parser)
    add_xtb_args(relax_parser)

    # HESSIAN subcommand  
    hessian_parser = subparsers.add_parser('hessian', help='Calculate Hessian/frequencies')
    hessian_parser.add_argument('project', help='Project folder name')
    hessian_parser.add_argument('--backend', default='xtb',
                               choices=['dftb', 'xtb', 'gxtb', 'orca', 'orca-simple'],
                               help='Hessian backend (default: xtb)')
    hessian_parser.add_argument('--dftb-params', default='auto',
                               choices=['auto', 'ptbp', '3ob'],
                               help='DFTB parameter set selection (default: auto). Only used with --backend dftb.')
    add_common_args(hessian_parser)
    add_xtb_args(hessian_parser)
    
    # THERMO subcommand
    thermo_parser = subparsers.add_parser('thermo', help='Thermodynamic analysis')
    thermo_parser.add_argument('project', help='Project folder name')
    thermo_parser.add_argument('--backend', default='dftb',
                              choices=['dftb', 'xtb', 'orca'],
                              help='Backend for thermodynamic analysis (default: dftb)')
    thermo_parser.add_argument('--temperature-range', default='298,1200',
                              help='Temperature range in K as min,max (default: 298,1200)')
    thermo_parser.add_argument('--steps', type=int, default=20,
                              help='Number of temperature points (default: 20)')
    thermo_parser.add_argument('--plot', action='store_true',
                              help='Create thermal analysis plots')
    thermo_parser.add_argument('--report', action='store_true',
                              help='Generate detailed thermal report')
    thermo_parser.add_argument('--bond-energy', type=float, default=4.0,
                              help='Bond dissociation energy in eV (default: 4.0)')
    thermo_parser.add_argument('--animate', action='store_true',
                              help='Generate MP4 animation of thermal heating')
    thermo_parser.add_argument('--frames', type=int, default=50,
                              help='Number of frames for animation (default: 50)')
    thermo_parser.add_argument('--fps', type=int, default=10,
                              help='Frames per second for animation (default: 10)')

    # BUILD subcommand
    build_parser = subparsers.add_parser('build', help='Build molecular structures')
    build_subparsers = build_parser.add_subparsers(dest='build_type', help='Structure types')
    
    # Build molecule
    mol_parser = build_subparsers.add_parser('molecule', help='Build simple molecules')
    mol_parser.add_argument('formula', help='Molecule formula (e.g., C6H6, H2O)')
    add_build_args(mol_parser)
    
    # Build diamond
    diamond_parser = build_subparsers.add_parser('diamond', help='Build diamond structures')
    diamond_parser.add_argument('size', help='Structure size')
    diamond_parser.add_argument('--type', default='crystal', choices=['crystal', 'supercell', 'surface'],
                               help='Diamond structure type (default: crystal)')
    add_build_args(diamond_parser)
    
    # Build diamondoid
    diamondoid_parser = build_subparsers.add_parser('diamondoid', help='Build diamondoid structures')
    diamondoid_parser.add_argument('shape', help='Shape:size (e.g., sphere:5, rod:5:2:2)')
    add_build_args(diamondoid_parser)
    
    # Build STM tip
    tip_parser = build_subparsers.add_parser('stm-tip', help='Build STM tips')
    tip_parser.add_argument('height', type=int, help='Tip height in BCC cells')
    tip_parser.add_argument('--base-size', type=int, default=5, help='Base size (default: 5)')
    tip_parser.add_argument('--sharpness', type=int, default=3, choices=[1,2,3,4,5],
                           help='Tip sharpness 1-5 (default: 3)')
    tip_parser.add_argument('--tip-atoms', type=int, default=1, choices=[1,2,3,4,5,6,7],
                           help='Atoms at tip apex (default: 1)')
    tip_parser.add_argument('--freeze-base', action='store_true', default=True,
                           help='Add fixation tags for base')
    add_build_args(tip_parser)
    add_common_args(tip_parser)
    add_xtb_args(tip_parser)
    
    # Build from SDF
    sdf_parser = build_subparsers.add_parser('sdf', help='Build from SDF file')
    sdf_parser.add_argument('file', help='SDF file path')
    sdf_parser.add_argument('--project', help='Project name (auto-generated if not specified)')
    
    # Build Si surface
    si_parser = build_subparsers.add_parser('si-surface', help='Build Si(100)-(2×1) surfaces')
    si_parser.add_argument('--size', default='4x4x3', help='Surface size nx×ny×layers (default: 4x4x3)')
    si_parser.add_argument('--no-passivate', action='store_true', help='Skip H-passivation of edges')
    add_build_args(si_parser)
    add_common_args(si_parser)
    add_xtb_args(si_parser)

    # TEST subcommand
    test_parser = subparsers.add_parser('test', help='Run pytest tests')
    test_parser.add_argument('pattern', nargs='?', help='Test pattern to match (e.g., "unified_dftb" for test_unified_dftb.py)')
    test_parser.add_argument('--verbose', '-v', action='store_true', help='Verbose test output')
    test_parser.add_argument('--markers', '-m', help='Run tests with specific markers (slow, integration, requires_xtb, requires_orca)')
    test_parser.add_argument('--direct', action='store_true', help='Run test directly without pytest')

    # SEQUENCE subcommand
    sequence_parser = subparsers.add_parser('sequence', help='Run molecular sequence simulations')
    sequence_parser.add_argument('scene_name', help='Scene name (YAML file in data/scenes/)')
    sequence_parser.add_argument('--generate', action='store_true', help='Generate input files without running calculations')
    sequence_parser.add_argument('--frame', type=int, help='Specific frame to process')
    sequence_parser.add_argument('--backend', default='dftb', 
                                choices=['dftb', 'xtb', 'orca'],
                                help='Calculation backend (default: dftb)')

    return parser

def execute_relax(args):
    """Execute relaxation with selected backend"""
    # Get project information for tracking
    project_name, atoms, n_atoms, n_electrons = get_project_info_for_tracking(args.project)
    
    # Prepare method parameters based on backend and args
    method_params = {
        'backend': args.backend,
    }
    
    # Add backend-specific parameters
    if args.backend in ['xtb', 'gxtb']:
        method_params.update({
            'xtb_method': getattr(args, 'xtb_method', 'GFN2-xTB'),
            'xtb_backend': getattr(args, 'xtb_backend', 'native')
        })
    elif args.backend in ['orca', 'orca-ase']:
        method_params.update({
            'method': getattr(args, 'method', 'B3LYP'),
            'basis': getattr(args, 'basis', 'def2-SVP'),
            'calc_hessian': getattr(args, 'calc_hessian', False)
        })
    
    # Track performance of this calculation
    with PerformanceTracker(
        method=f"relax-{args.backend}",
        molecule_name=project_name,
        method_params=method_params,
        n_cores=getattr(args, 'nprocs', 1)
    ) as tracker:
        if atoms is not None:
            tracker.set_atoms(atoms)
        
        backend_map = {
            'dftb': lambda: relax_dftb_unified(args.project, getattr(args, 'dftb_params', 'auto')),
            'dftb-native': lambda: relax_native_dftb(args.project),
            'dftb-passivated': lambda: relax_passivated(args.project),
            'dftb-ptbp': lambda: relax_native_dftb_ptbp(args.project),
            'xtb': lambda: relax_xtb(args.project, args.xtb_method, args.nprocs, args.xtb_backend),
            'gxtb': lambda: relax_gxtb(args.project),
            # Temporarily commented out due to indentation issues:
            # 'orca': lambda: relax_orca_native(args.project, args.method, args.basis, args.nprocs, args.calc_hessian, args.hessian_method),
            # 'orca-ase': lambda: relax_orca_dft(args.project, args.method, args.basis, args.nprocs)
        }
        
        try:
            success = backend_map[args.backend]()
            if success:
                tracker.set_quality(4, f"Successful {args.backend} relaxation")
            else:
                tracker.set_quality(2, f"{args.backend} relaxation completed with issues")
        except Exception as e:
            tracker.set_quality(1, f"{args.backend} relaxation failed: {str(e)}")
            raise

def execute_hessian(args):
    """Execute Hessian calculation with selected backend"""
    # Get project information for tracking
    project_name, atoms, n_atoms, n_electrons = get_project_info_for_tracking(args.project)
    
    # Prepare method parameters based on backend and args
    method_params = {
        'backend': args.backend,
        'calculation_type': 'hessian'
    }
    
    # Add backend-specific parameters
    if args.backend in ['xtb']:
        method_params.update({
            'xtb_method': getattr(args, 'xtb_method', 'GFN2-xTB'),
            'xtb_backend': getattr(args, 'xtb_backend', 'native')
        })
    elif args.backend in ['orca', 'orca-simple']:
        method_params.update({
            'method': getattr(args, 'method', 'B3LYP') if args.backend == 'orca' else 'PBE',
            'basis': getattr(args, 'basis', 'def2-SVP') if args.backend == 'orca' else 'def2-SVP'
        })
    
    # Track performance of this calculation
    with PerformanceTracker(
        method=f"hessian-{args.backend}",
        molecule_name=project_name,
        method_params=method_params,
        n_cores=getattr(args, 'nprocs', 1)
    ) as tracker:
        if atoms is not None:
            tracker.set_atoms(atoms)
        
        backend_map = {
            'dftb': lambda: calculate_hessian_dftb_unified(args.project, getattr(args, 'dftb_params', 'auto')),
            'xtb': lambda: calculate_hessian_xtb_standalone(args.project, args.xtb_method, args.xtb_backend),
            'gxtb': lambda: calculate_hessian_gxtb_standalone(args.project),
            'orca': lambda: calculate_hessian_orca_standalone(args.project, args.method, args.basis, args.nprocs),
            'orca-simple': lambda: calculate_hessian_orca_standalone(args.project, 'PBE', 'def2-SVP', args.nprocs)
        }
        
        try:
            backend_map[args.backend]()
            tracker.set_quality(4, f"Successful {args.backend} hessian calculation")
        except Exception as e:
            tracker.set_quality(1, f"{args.backend} hessian calculation failed: {str(e)}")
            raise

def run_relaxation_workflow(project_name, args):
    """Run relaxation with specified backend
    
    Returns:
        bool: True if relaxation succeeded, False otherwise
    """
    backend_map = {
        'dftb': lambda: relax_dftb_unified(project_name, getattr(args, 'dftb_params', 'auto')),
        'dftb-native': lambda: relax_native_dftb(project_name),
        'dftb-passivated': lambda: relax_passivated(project_name),
        'dftb-ptbp': lambda: relax_native_dftb_ptbp(project_name),
        'xtb': lambda: relax_xtb(project_name, args.xtb_method, args.nprocs, args.xtb_backend),
        'gxtb': lambda: relax_gxtb(project_name),
        'orca': lambda: relax_orca_native(project_name, args.method, args.basis, args.nprocs, False, 'phonopy')
    }
    
    # Determine backend: specific override or unified backend
    backend = args.relax_backend if hasattr(args, 'relax_backend') and args.relax_backend else args.backend
    
    # Smart default selection based on passivation
    if hasattr(args, 'passivate') and args.passivate and backend == 'xtb':
        print("Note: Using dftb-passivated for passivated structure")
        backend = 'dftb-passivated'
    
    # Execute the backend function and return its result
    return backend_map[backend]()

def run_hessian_workflow(project_name, args):
    """Run Hessian calculation with specified backend"""
    backend_map = {
        'dftb': lambda: calculate_hessian_dftb_unified(project_name, getattr(args, 'dftb_params', 'auto')),
        'dftb-native': lambda: calculate_hessian_dftb_native(project_name),
        'dftb-passivated': lambda: calculate_hessian_dftb_native(project_name),
        'dftb-ptbp': lambda: calculate_hessian_dftb_native(project_name),
        'xtb': lambda: calculate_hessian_xtb_standalone(project_name, args.xtb_method, args.xtb_backend),
        'gxtb': lambda: calculate_hessian_gxtb_standalone(project_name),
        'orca': lambda: calculate_hessian_orca_standalone(project_name, args.method, args.basis, args.nprocs),
        'orca-simple': lambda: calculate_hessian_orca_standalone(project_name, 'PBE', 'def2-SVP', args.nprocs)
    }
    
    # Determine backend: specific override or unified backend
    backend = args.hessian_backend if hasattr(args, 'hessian_backend') and args.hessian_backend else args.backend
    
    # Map some relaxation backends to appropriate Hessian backends
    # Map legacy DFTB backends to unified backend for hessian calculation compatibility
    if backend in ['dftb-native', 'dftb-passivated', 'dftb-ptbp']:
        backend = 'dftb'
    
    try:
        backend_map[backend]()
        return True
    except Exception as e:
        print(f"❌ Hessian calculation failed: {str(e)}")
        return False


def execute_thermo(args):
    """Execute thermodynamic analysis"""
    from kernel.backends.dftb_backend import DFTBBackend
    
    # Parse temperature range
    try:
        t_min, t_max = map(float, args.temperature_range.split(','))
        if t_min <= 0 or t_max <= t_min:
            raise ValueError("Invalid temperature range")
    except ValueError:
        print("❌ Invalid temperature range. Use format: min,max (e.g., 298,1200)")
        return
    
    # Create backend instance
    if args.backend == 'dftb':
        backend = DFTBBackend()
    else:
        print(f"❌ Backend {args.backend} not yet supported for thermodynamic analysis")
        print("   Currently supported: dftb")
        return
    
    print(f"🌡️  Running thermodynamic analysis for {args.project}")
    print(f"   Backend: {args.backend}")
    print(f"   Temperature range: {t_min:.0f}K - {t_max:.0f}K ({t_min-273.15:.0f}°C - {t_max-273.15:.0f}°C)")
    print(f"   Steps: {args.steps}")
    
    try:
        # Run thermodynamic analysis
        results = backend.calculate_thermodynamics(
            project_name=args.project,
            temperature_range=(t_min, t_max),
            n_points=args.steps,
            plot_results=args.plot,
            save_report=args.report,
            animate=args.animate,
            n_frames=args.frames,
            fps=args.fps
        )
        
        # Display key results
        if results['critical_temperature']:
            print(f"\n🔥 Critical temperature: {results['critical_temperature']:.0f}K ({results['critical_temperature']-273.15:.0f}°C)")
            print("   Above this temperature, thermal effects become significant")
        else:
            print(f"\n✅ Structure remains stable up to {t_max:.0f}K")
        
        print(f"\n📊 Results saved to: {results['output_dir']}")
        if args.plot:
            print("   - thermal_analysis.png (visualization)")
        if args.report:
            print("   - thermal_report.txt (detailed report)")
        if args.animate and results.get('animation_result'):
            if results['animation_result'].endswith('.mp4'):
                print("   - thermal_heating.mp4 (animation)")
            else:
                print(f"   - frames/ (animation frames)")
        print("   - thermal_data.json (raw data)")
        
    except FileNotFoundError as e:
        print(f"❌ {str(e)}")
        print("   Make sure to run Hessian calculation first:")
        print(f"   PYTHONPATH=py python py/main.py hessian {args.project} --backend {args.backend}")
    except Exception as e:
        print(f"❌ Thermodynamic analysis failed: {str(e)}")


def run_automated_relax_hessian_workflow(project_name, args):
    """Automated relax-hessian workflow until true minimum is found"""
    backend = args.backend
    iteration = 1
    max_iterations = args.max_iterations
    
    print(f"\n🔄 Starting automated relax-hessian workflow for {project_name}")
    print(f"   Backend: {backend}, Max iterations: {max_iterations}")
    
    while iteration <= max_iterations:
        print(f"\n--- Iteration {iteration} ---")
        
        # 1. Relaxation
        print(f"🔧 Running relaxation (iteration {iteration})...")
        relax_success = run_relaxation_workflow(project_name, args)
        
        if not relax_success:
            print(f"❌ Relaxation failed at iteration {iteration}")
            return False
            
        # Save relaxation result with iteration-specific name (if not first iteration)
        if iteration > 1:
            rename_relaxed_result(project_name, backend, iteration)
            
        # 2. Hessian calculation  
        print(f"🔍 Running Hessian calculation (iteration {iteration})...")
        hessian_success = run_hessian_workflow(project_name, args)
        
        if not hessian_success:
            print(f"❌ Hessian calculation failed at iteration {iteration}")
            return False
            
        # Save hessian result with iteration-specific name (if not first iteration)
        if iteration > 1:
            rename_hessian_result(project_name, backend, iteration)
            
        # 3. Check for imaginary frequencies
        imaginary_count = count_imaginary_frequencies(project_name, backend, iteration)
        
        if imaginary_count < 0:
            print(f"❌ Error checking frequencies at iteration {iteration}")
            return False
            
        if imaginary_count == 0:
            print(f"✅ True minimum found at iteration {iteration}!")
            print(f"   No imaginary frequencies detected")
            
            # Finalize results (copy to standard names)
            finalize_results(project_name, backend, iteration)
            return True
            
        print(f"⚠️  Found {imaginary_count} imaginary frequencies at iteration {iteration}")
        
        # Check if we've reached maximum iterations
        if iteration >= max_iterations:
            print(f"❌ Maximum iterations ({max_iterations}) reached")
            print(f"   Structure still contains {imaginary_count} imaginary frequencies")
            print(f"   Consider increasing --max-iterations or using different parameters")
            return False
            
        # 4. Escape from saddle point for next iteration
        print(f"🚀 Escaping saddle point for next iteration...")
        escape_success = escape_saddle_point_with_iteration(project_name, backend, iteration)
        
        if not escape_success:
            print(f"❌ Failed to escape saddle point at iteration {iteration}")
            return False
            
        iteration += 1
    
    # Should not reach here
    return False


def execute_build(args):
    """Execute build commands"""
    project_name = None  # Initialize project_name
    
    if args.build_type == 'molecule':
        project_name = f"2025-08-23_{args.formula.lower()}"  # Generate consistent project name
        build_molecule(args.formula)
        
    elif args.build_type == 'diamond':
        project_name = get_project_name('diamond', args.type, args.size)
        create_project_structure(project_name)
        
        if args.type == 'crystal':
            output_file = args.output or f'diamond_{args.size}.xyz'
            build_diamond_with_highlighting(args.size, args.highlight, output_file, args.passivate)
        elif args.type == 'supercell':
            output_file = args.output or f'diamond_supercell_{args.size}.xyz'
            build_diamond_supercell(args.size, output_file, args.passivate)
        elif args.type == 'surface':
            output_file = args.output or f'diamond_surface_{args.size}.xyz'
            build_diamond_surface(args.size, args.highlight, output_file, args.passivate)
            
        # Move to project inputs
        project_path = settings.get_project_data_path(project_name)
        inputs_path = project_path / 'inputs' / output_file
        os.rename(output_file, str(inputs_path))
        print(f"Moved {output_file} to {inputs_path}")
        
        if args.relax:
            backend = args.relax_backend if hasattr(args, 'relax_backend') and args.relax_backend else args.backend
            print(f"Running relaxation with {backend}...")
            relaxation_success = run_relaxation_workflow(project_name, args)
            if not relaxation_success:
                print("⚠️  Relaxation failed. Skipping Hessian calculation.")
                return
            
        if args.hessian:
            backend = args.hessian_backend if hasattr(args, 'hessian_backend') and args.hessian_backend else args.backend
            print(f"Running Hessian calculation with {backend}...")
            run_hessian_workflow(project_name, args)
            
        if args.continue_from_saddle:
            backend = args.backend
            print(f"Continuing optimization from saddle point with {backend}...")
            success = escape_saddle_point(project_name, displacement_scale=0.1, backend=backend)
            if not success:
                print("❌ Failed to continue from saddle point")
                return
                
                
    elif args.build_type == 'diamondoid':
        try:
            if ':' in args.shape:
                parts = args.shape.split(':', 1)
                shape = parts[0]
                size = parts[1]
            else:
                shape = args.shape
                size = '5' if shape in ['sphere', 'cube', 'octahedron', 'tetrahedron'] else None
            
            if shape in ['rod', 'platform']:
                project_name = f"diamondoid_{shape}_{size.replace(':', 'x')}"
            else:
                project_name = f"diamondoid_{shape}"
                
            create_project_structure(project_name)
            output_file = args.output or f'diamondoid_{shape}.xyz'
            build_diamondoid(shape, size, args.passivate, output_file)
            
            project_path = settings.get_project_data_path(project_name)
            inputs_path = project_path / 'inputs' / output_file
            os.rename(output_file, str(inputs_path))
            print(f"Moved {output_file} to {inputs_path}")
            
            if args.relax:
                backend = args.relax_backend if hasattr(args, 'relax_backend') and args.relax_backend else args.backend
                print(f"Running relaxation with {backend}...")
                relaxation_success = run_relaxation_workflow(project_name, args)
                if not relaxation_success:
                    print("⚠️  Relaxation failed. Skipping Hessian calculation.")
                    return
                
            if args.hessian:
                backend = args.hessian_backend if hasattr(args, 'hessian_backend') and args.hessian_backend else args.backend
                print(f"Running Hessian calculation with {backend}...")
                run_hessian_workflow(project_name, args)
                    
        except ValueError as e:
            print(f"Error parsing diamondoid parameters: {e}")
            sys.exit(1)
            
    elif args.build_type == 'stm-tip':
        project_name = f'w_stm_tip_{args.height}_{args.base_size}'
        create_project_structure(project_name)
        output_file = args.output or f'w_stm_tip_{args.height}_{args.base_size}.xyz'
        build_tungsten_stm_tip(
            height=args.height,
            base_size=args.base_size,
            tip_sharpness=args.sharpness,
            tip_atoms=args.tip_atoms,
            freeze_base=args.freeze_base,
            output_file=output_file
        )
        
        os.rename(output_file, f'data/{project_name}/inputs/{output_file}')
        print(f"Moved {output_file} to data/{project_name}/inputs/")
        
        if args.relax:
            backend = args.relax_backend if hasattr(args, 'relax_backend') and args.relax_backend else args.backend
            print(f"Running relaxation with {backend}...")
            relaxation_success = run_relaxation_workflow(project_name, args)
            if not relaxation_success:
                print("⚠️  Relaxation failed. Skipping Hessian calculation.")
                return
            
        if args.hessian:
            backend = args.hessian_backend if hasattr(args, 'hessian_backend') and args.hessian_backend else args.backend
            print(f"Running Hessian calculation with {backend}...")
            run_hessian_workflow(project_name, args)
            
            
    elif args.build_type == 'sdf':
        project_name = args.project if args.project else None
        build_from_sdf(args.file, project_name)
        
    elif args.build_type == 'si-surface':
        # Parse size parameter
        try:
            if 'x' in args.size:
                parts = args.size.split('x')
                if len(parts) == 3:
                    nx, ny, layers = map(int, parts)
                else:
                    raise ValueError("Size must be nx×ny×layers format")
            else:
                raise ValueError("Size must be nx×ny×layers format")
        except ValueError as e:
            print(f"Error parsing size parameter: {e}")
            print("Example: --size 4x4x3")
            sys.exit(1)
            
        passivate = not args.no_passivate
        project_name = f'si_surface_{args.size}'
        
        print(f"Building Si(100)-(2×1) surface: {nx}×{ny}×{layers}, passivate={passivate}")
        
        metadata = build_si_surface(nx=nx, ny=ny, layers=layers, 
                                   passivate=passivate, project_name=project_name)
        
        if args.relax:
            backend = args.relax_backend if hasattr(args, 'relax_backend') and args.relax_backend else args.backend
            print(f"Running relaxation with {backend}...")
            relaxation_success = run_relaxation_workflow(project_name, args)
            if not relaxation_success:
                print("⚠️  Relaxation failed. Skipping Hessian calculation.")
                return
                
        if args.hessian:
            backend = args.hessian_backend if hasattr(args, 'hessian_backend') and args.hessian_backend else args.backend
            print(f"Running Hessian calculation with {backend}...")
            run_hessian_workflow(project_name, args)
    
    # Universal post-processing workflow handlers
    # These work regardless of build type
    if hasattr(args, 'relax_hessian') and args.relax_hessian:
        if project_name is not None:
            print(f"\n🔄 Running automated relax-hessian workflow...")
            success = run_automated_relax_hessian_workflow(project_name, args)
            if not success:
                print("❌ Automated workflow failed to find true minimum")
                return
        else:
            print("❌ Cannot run automated workflow: no project created")


def execute_test(args):
    """Execute tests with optional pattern and marker filtering"""
    
    # Handle direct execution mode for specific tests
    if args.direct and args.pattern:
        # Convert pattern to test file path
        if not args.pattern.startswith('test_'):
            pattern = f'test_{args.pattern}'
        else:
            pattern = args.pattern
        
        if not pattern.endswith('.py'):
            pattern = f'{pattern}.py'
            
        test_path = f'py/tests/{pattern}'
        
        # Check if test file exists
        if not os.path.exists(test_path):
            print(f"❌ Test file not found: {test_path}")
            return False
        
        print(f"🚀 Running test directly: {test_path}")
        
        # Run test directly as Python script
        try:
            result = subprocess.run(['python', test_path], cwd=os.getcwd())
            return result.returncode == 0
        except Exception as e:
            print(f"❌ Error running test: {e}")
            return False
    
    # Standard pytest execution
    cmd = ['python', '-m', 'pytest']
    
    # Add verbose flag if requested
    if args.verbose:
        cmd.append('-v')
    
    # Add marker filtering if specified
    if args.markers:
        cmd.extend(['-m', args.markers])
    
    # Add test pattern if specified
    if args.pattern:
        # Convert pattern to test file path
        if not args.pattern.startswith('test_'):
            pattern = f'test_{args.pattern}'
        else:
            pattern = args.pattern
        
        if not pattern.endswith('.py'):
            pattern = f'{pattern}.py'
            
        test_path = f'py/tests/{pattern}'
        cmd.append(test_path)
    else:
        # Run all tests in tests directory
        cmd.append('py/tests/')
    
    print(f"Running: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, cwd=os.getcwd())
        return result.returncode == 0
    except Exception as e:
        print(f"Error running tests: {e}")
        print("💡 Try using --direct flag to run test directly")
        return False

def execute_sequence(args):
    """Execute sequence generation or calculation"""
    scene_name = args.scene_name
    scene_path = f"data/scenes/{scene_name}/{scene_name}.yaml"
    
    if not os.path.exists(scene_path):
        print(f"❌ Scene file not found: {scene_path}")
        return False
    
    try:
        generator = SceneGenerator(scene_path)
        
        if args.generate:
            if args.frame is not None:
                print(f"🔧 Generating frame {args.frame} for scene {scene_name}")
                success = generator.generate_frame(args.frame)
                if success:
                    print(f"✅ Frame {args.frame} generated successfully")
                else:
                    print(f"❌ Failed to generate frame {args.frame}")
                return success
            else:
                print(f"❌ --frame required with --generate")
                return False
        else:
            print(f"❌ Calculation mode not implemented yet")
            return False
            
    except Exception as e:
        print(f"❌ Error processing sequence: {e}")
        return False

def main():
    parser = setup_parser()
    
    if len(sys.argv) == 1:
        parser.print_help()
        return
        
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        return
    
    # Execute commands
    if args.command == 'create':
        project_name = get_project_name(args.name)
        create_project_structure(project_name)
    elif args.command == 'relax':
        execute_relax(args)
    elif args.command == 'hessian':
        execute_hessian(args)
    elif args.command == 'thermo':
        execute_thermo(args)
    elif args.command == 'build':
        if not args.build_type:
            print("Error: build command requires a structure type")
            parser.parse_args(['build', '--help'])
            return
        execute_build(args)
    elif args.command == 'test':
        success = execute_test(args)
        if not success:
            sys.exit(1)
    elif args.command == 'sequence':
        success = execute_sequence(args)
        if not success:
            sys.exit(1)

if __name__ == '__main__':
    main()