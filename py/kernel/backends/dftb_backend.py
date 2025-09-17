"""
Unified DFTB+ Backend Implementation

Handles DFTB+ density functional tight-binding calculations with automatic
parameter selection (PTBP/3OB) based on molecular composition
"""
import os
import glob
from typing import Dict, Any, Optional, Set, List
from ase import Atoms
from ase.io import read
from .base import RelaxationBackend, HessianBackend, BackendResult
from ..parts.thermodynamics import analyze_thermal_stability, create_thermal_report, plot_thermal_analysis, create_thermal_animation


class DFTBBackend(RelaxationBackend, HessianBackend):
    """Unified DFTB+ backend with automatic parameter selection"""
    
    def __init__(self, params: str = 'auto', **kwargs):
        super().__init__('dftb', **kwargs)
        self.params = params  # 'auto', 'ptbp', '3ob'
        self._heavy_metals = {
            'W', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
            'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Tc'
        }
        
    def _get_slakos_path(self) -> str:
        """Get path to slakos parameters"""
        # Check for local 3rdparty first
        local_slakos = "3rdparty/slakos"
        if os.path.exists(local_slakos):
            return local_slakos
        
        # Fallback to system installation
        return os.environ.get('SKDIR', '/Users/andreypanferov/opt/dftb+/slakos')
    
    def _detect_parameter_set(self, atoms: Atoms) -> str:
        """Automatically detect best parameter set for given atoms"""
        if self.params != 'auto':
            return self.params
            
        elements = set(atoms.get_chemical_symbols())
        
        # If contains heavy metals, use PTBP
        if elements.intersection(self._heavy_metals):
            return 'ptbp'
        
        # For organic molecules, could use 3OB but PTBP is more universal
        return 'ptbp'
    
    def _get_parameter_path(self, param_set: str) -> str:
        """Get full path to parameter set"""
        slakos_base = self._get_slakos_path()
        
        if param_set == 'ptbp':
            return os.path.join(slakos_base, 'ptbp-2024')
        elif param_set == '3ob':
            return os.path.join(slakos_base, '3ob-3-1')
        else:
            # Default to PTBP
            return os.path.join(slakos_base, 'ptbp-2024')
    
    def _get_angular_momentum_params(self, atoms: Atoms, param_set: str) -> str:
        """Generate angular momentum parameters for given atoms and parameter set"""
        if param_set == 'ptbp':
            return self._get_angular_momentum_params_ptbp(atoms)
        else:
            return self._get_angular_momentum_params_3ob(atoms)
    
    def _get_angular_momentum_params_ptbp(self, atoms: Atoms) -> str:
        """Generate MaxAngularMomentum parameters for PTBP parameter set"""
        # PTBP supports elements H-Rn with appropriate angular momentum
        angular_momentum_map = {
            'H': 's', 'He': 's',
            'Li': 's', 'Be': 's', 'B': 'p', 'C': 'p', 'N': 'p', 'O': 'p', 'F': 'p', 'Ne': 'p',
            'Na': 's', 'Mg': 's', 'Al': 'p', 'Si': 'p', 'P': 'p', 'S': 'p', 'Cl': 'p', 'Ar': 'p',
            'K': 's', 'Ca': 's', 'Sc': 'd', 'Ti': 'd', 'V': 'd', 'Cr': 'd', 'Mn': 'd', 'Fe': 'd',
            'Co': 'd', 'Ni': 'd', 'Cu': 'd', 'Zn': 'd', 'Ga': 'p', 'Ge': 'p', 'As': 'p', 'Se': 'p',
            'Br': 'p', 'Kr': 'p', 'Rb': 's', 'Sr': 's', 'Y': 'd', 'Zr': 'd', 'Nb': 'd', 'Mo': 'd',
            'Tc': 'd', 'Ru': 'd', 'Rh': 'd', 'Pd': 'd', 'Ag': 'd', 'Cd': 'd', 'In': 'p', 'Sn': 'p',
            'Sb': 'p', 'Te': 'p', 'I': 'p', 'Xe': 'p', 'Cs': 's', 'Ba': 's', 'La': 'd', 'Ce': 'f',
            'Hf': 'd', 'Ta': 'd', 'W': 'd', 'Re': 'd', 'Os': 'd', 'Ir': 'd', 'Pt': 'd', 'Au': 'd',
            'Hg': 'd', 'Tl': 'p', 'Pb': 'p', 'Bi': 'p', 'Po': 'p', 'At': 'p', 'Rn': 'p'
        }
        
        unique_elements = sorted(set(atoms.get_chemical_symbols()))
        angular_params = []
        
        for element in unique_elements:
            if element in angular_momentum_map:
                angular_params.append(f"{element} = \"{angular_momentum_map[element]}\"")
            else:
                # Default to 'p' for unknown elements
                angular_params.append(f"{element} = \"p\"")
        
        return "    " + "\n    ".join(angular_params)
    
    def _get_angular_momentum_params_3ob(self, atoms: Atoms) -> str:
        """Generate MaxAngularMomentum parameters for 3OB parameter set"""
        angular_momentum_map = {
            'H': 's', 'C': 'p', 'N': 'p', 'O': 'p', 'F': 'p', 'S': 'p', 'P': 'p',
            'Cl': 'p', 'Br': 'p', 'I': 'p', 'Zn': 'd', 'Na': 's', 'K': 's', 
            'Ca': 's', 'Mg': 's'
        }
        
        unique_elements = sorted(set(atoms.get_chemical_symbols()))
        angular_params = []
        
        for element in unique_elements:
            if element in angular_momentum_map:
                angular_params.append(f"{element} = \"{angular_momentum_map[element]}\"")
            else:
                raise ValueError(f"Element {element} not supported by 3OB parameter set")
        
        return "    " + "\n    ".join(angular_params)
    
    def is_available(self) -> bool:
        """Check if DFTB+ is available"""
        dftb_path = '/Users/andreypanferov/opt/dftb+/bin/dftb+'
        if not os.path.exists(dftb_path):
            return False
        
        # Check if at least one parameter set is available
        slakos_base = self._get_slakos_path()
        ptbp_path = os.path.join(slakos_base, 'ptbp-2024')
        ob3_path = os.path.join(slakos_base, '3ob-3-1')
        
        return os.path.exists(ptbp_path) or os.path.exists(ob3_path)
    
    def single_point(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """Single point DFTB+ calculation"""
        return {'energy': 0.0, 'forces': [], 'success': False}
    
    def optimize(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """DFTB+ geometry optimization"""
        return {'atoms': atoms, 'converged': False, 'success': False}
    
    def frequencies(self, atoms: Atoms, **kwargs) -> Dict[str, Any]:
        """DFTB+ frequency calculation"""
        return {'frequencies': [], 'success': False}
    
    def relax_structure(self, project_name: str, atoms: Optional[Atoms] = None, **kwargs) -> bool:
        """Relax structure using unified DFTB+ implementation"""
        try:
            from ..parts.relax import relax_dftb_unified
            return relax_dftb_unified(project_name, self.params)
        except ImportError:
            # Fallback to existing implementations during transition
            from ..parts.relax import relax_native_dftb_ptbp
            return relax_native_dftb_ptbp(project_name)
    
    def calculate_hessian(self, project_name: str, atoms: Optional[Atoms] = None, **kwargs) -> bool:
        """Calculate Hessian using unified DFTB+ implementation"""
        try:
            from ..parts.hessian import calculate_hessian_dftb_unified
            return calculate_hessian_dftb_unified(project_name, self.params)
        except ImportError:
            # Fallback to existing implementation during transition
            from ..parts.hessian import calculate_hessian_dftb_native
            return calculate_hessian_dftb_native(project_name)
    
    def calculate_thermodynamics(self, project_name: str, temperature_range=(298, 1200), 
                               n_points=20, plot_results=True, save_report=True, 
                               animate=False, n_frames=50, fps=10) -> Dict[str, Any]:
        """
        Calculate thermodynamic properties over temperature range after Hessian calculation.
        
        Args:
            project_name: Project folder name
            temperature_range: (T_min, T_max) in Kelvin
            n_points: Number of temperature points
            plot_results: Create visualization plots
            save_report: Save detailed report
            animate: Generate MP4 animation of thermal heating
            n_frames: Number of frames for animation
            fps: Frames per second for animation
            
        Returns:
            Dictionary with thermal analysis results
        """
        import numpy as np
        from ase.io import read
        
        # Find project directory
        project_dirs = [
            f"data/{project_name}",
            f"data/*_{project_name}",
        ]
        
        project_dir = None
        for pattern in project_dirs:
            matches = glob.glob(pattern)
            if matches:
                project_dir = matches[0]
                break
        
        if not project_dir or not os.path.exists(project_dir):
            raise FileNotFoundError(f"Project directory not found for {project_name}")
        
        # Look for frequency data from previous Hessian calculation
        freq_files = [
            os.path.join(project_dir, "outputs", "dftb", "frequencies.dat"),
            os.path.join(project_dir, "outputs", "frequencies.dat"),
            os.path.join(project_dir, "dftb_frequencies.dat"),
            os.path.join(project_dir, "frequencies.dat")
        ]
        
        frequencies = None
        for freq_file in freq_files:
            if os.path.exists(freq_file):
                try:
                    frequencies = np.loadtxt(freq_file)
                    if len(frequencies.shape) == 0:  # Single frequency
                        frequencies = [float(frequencies)]
                    else:
                        frequencies = frequencies.tolist()
                    break
                except:
                    continue
        
        if frequencies is None:
            raise FileNotFoundError(f"No frequency data found for {project_name}. Run Hessian calculation first.")
        
        # Load molecular structure
        structure_files = [
            os.path.join(project_dir, "outputs", "dftb", f"{project_name}_relaxed.xyz"),
            os.path.join(project_dir, "outputs", f"{project_name}_relaxed.xyz"),
            os.path.join(project_dir, f"{project_name}_relaxed.xyz"),
            os.path.join(project_dir, f"{project_name}.xyz")
        ]
        
        atoms = None
        for struct_file in structure_files:
            if os.path.exists(struct_file):
                try:
                    atoms = read(struct_file)
                    break
                except:
                    continue
        
        if atoms is None:
            raise FileNotFoundError(f"No structure file found for {project_name}")
        
        # Generate temperature points
        T_min, T_max = temperature_range
        temperatures = np.linspace(T_min, T_max, n_points)
        
        # Perform thermal analysis
        print(f"Calculating thermodynamic properties from {T_min}K to {T_max}K...")
        analysis_results = analyze_thermal_stability(frequencies, atoms, temperatures)
        
        # Create output directory
        output_dir = os.path.join(project_dir, "outputs", "thermodynamics")
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate plots
        if plot_results:
            print("Creating thermal analysis plots...")
            plot_thermal_analysis(analysis_results, output_dir=output_dir, show_plots=False)
        
        # Generate report
        if save_report:
            print("Generating thermal analysis report...")
            report_file = os.path.join(output_dir, "thermal_report.txt")
            report = create_thermal_report(analysis_results, frequencies, atoms, report_file)
            print(f"Report saved to: {report_file}")
        
        # Generate animation
        animation_result = None
        if animate:
            print("Creating thermal animation...")
            try:
                animation_result = create_thermal_animation(
                    analysis_results, atoms, output_dir, n_frames=n_frames, fps=fps
                )
            except Exception as e:
                print(f"⚠️  Animation creation failed: {e}")
                animation_result = None
        
        # Save thermal data
        import json
        data_file = os.path.join(output_dir, "thermal_data.json")
        with open(data_file, 'w') as f:
            # Convert numpy arrays to lists for JSON serialization
            json_data = {
                'temperatures': temperatures.tolist(),
                'frequencies_cm': frequencies,
                'thermal_properties': analysis_results['thermodynamic_data'],
                'stability_analysis': {
                    'thermal_energies': analysis_results['thermal_energies'],
                    'stability_ratios': analysis_results['stability_ratios'],
                    'critical_temperature': analysis_results['critical_temperature']
                }
            }
            json.dump(json_data, f, indent=2)
        
        print(f"Thermal analysis completed. Results saved to: {output_dir}")
        
        return {
            'project_dir': project_dir,
            'output_dir': output_dir,
            'temperatures': temperatures,
            'frequencies': frequencies,
            'analysis_results': analysis_results,
            'critical_temperature': analysis_results['critical_temperature'],
            'animation_result': animation_result
        }


# Legacy aliases for backward compatibility
DFTBRelaxationBackend = DFTBBackend
DFTBHessianBackend = DFTBBackend