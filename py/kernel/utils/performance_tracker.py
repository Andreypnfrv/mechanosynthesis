"""
Performance tracking system for mechanosynthesis calculations.

Provides context manager for automatic logging of relaxation performance metrics
to CSV file for analysis and benchmarking.
"""

import csv
import time
import platform
import psutil
import os
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, Any
import sys
sys.path.append(str(Path(__file__).parent.parent.parent))
import settings


class PerformanceTracker:
    """Context manager for tracking calculation performance metrics."""
    
    def __init__(self, method: str, molecule_name: str, method_params: Dict[str, Any] = None,
                 n_cores: int = 1, csv_path: Optional[Path] = None):
        """
        Initialize performance tracker.
        
        Args:
            method: Calculation method (xtb, dftb, orca, gxtb)
            molecule_name: Name/identifier of the molecule
            method_params: Dictionary of method parameters (basis, functional, etc.)
            n_cores: Number of CPU cores used
            csv_path: Path to CSV file (defaults to project root)
        """
        self.method = method
        self.molecule_name = molecule_name
        self.method_params = method_params or {}
        self.n_cores = n_cores
        
        # CSV file location
        if csv_path is None:
            # Check for environment variable override (useful for testing)
            env_path = os.environ.get('PERFORMANCE_LOG_PATH')
            if env_path:
                self.csv_path = Path(env_path)
            else:
                self.csv_path = settings.PROJECT_ROOT / "performance_logs.csv"
        else:
            self.csv_path = csv_path
            
        # Runtime tracking
        self.start_time = None
        self.end_time = None
        self.start_memory = None
        self.peak_memory = None
        
        # Results
        self.atoms = None
        self.success = False
        self.convergence_steps = None
        self.final_energy = None
        self.quality_rating = None
        self.notes = ""
        
        # Hardware info (cached)
        self.hardware_info = self._get_hardware_info()
        
    def _get_hardware_info(self) -> str:
        """Get hardware description."""
        system = platform.system()
        if system == "Darwin":
            # macOS
            try:
                import subprocess
                result = subprocess.run(['sysctl', '-n', 'machdep.cpu.brand_string'], 
                                      capture_output=True, text=True)
                cpu_model = result.stdout.strip() if result.returncode == 0 else "Unknown Mac"
                return f"macOS {cpu_model}"
            except:
                return f"macOS {platform.processor()}"
        elif system == "Linux":
            try:
                with open('/proc/cpuinfo', 'r') as f:
                    for line in f:
                        if line.startswith('model name'):
                            cpu_model = line.split(':', 1)[1].strip()
                            return f"Linux {cpu_model}"
            except:
                pass
            return f"Linux {platform.processor()}"
        else:
            return f"{system} {platform.processor()}"
    
    def __enter__(self):
        """Start tracking."""
        self.start_time = time.time()
        try:
            process = psutil.Process()
            self.start_memory = process.memory_info().rss / 1024 / 1024  # MB
        except:
            self.start_memory = None
        return self
        
    def __exit__(self, exc_type, exc_val, exc_tb):
        """End tracking and write to CSV."""
        self.end_time = time.time()
        
        # Determine success based on whether exception occurred
        self.success = exc_type is None
        
        # Get peak memory usage
        try:
            process = psutil.Process()
            current_memory = process.memory_info().rss / 1024 / 1024  # MB
            self.peak_memory = max(current_memory, self.start_memory or 0)
        except:
            self.peak_memory = None
            
        # Write to CSV
        self._write_to_csv()
        
    def set_atoms(self, atoms):
        """Set molecule atoms object for automatic property extraction."""
        self.atoms = atoms
        
    def set_results(self, convergence_steps: Optional[int] = None, 
                   final_energy: Optional[float] = None):
        """Set calculation results."""
        self.convergence_steps = convergence_steps
        self.final_energy = final_energy
        
    def set_quality(self, rating: int, notes: str = ""):
        """Set subjective quality assessment.
        
        Args:
            rating: Quality rating 1-5 (1=poor, 5=excellent)
            notes: Additional notes about quality/issues
        """
        if not 1 <= rating <= 5:
            raise ValueError("Quality rating must be between 1 and 5")
        self.quality_rating = rating
        self.notes = notes
        
    def _write_to_csv(self):
        """Write performance data to CSV file."""
        # Prepare data row
        wall_time = (self.end_time - self.start_time) if (self.start_time and self.end_time) else None
        
        # Extract molecule properties
        n_atoms = len(self.atoms) if self.atoms is not None else None
        n_electrons = sum(atom.number for atom in self.atoms) if self.atoms is not None else None
        
        # Format method parameters
        method_params_str = str(self.method_params) if self.method_params else ""
        
        row_data = {
            'timestamp': datetime.now().isoformat(),
            'molecule_name': self.molecule_name,
            'n_atoms': n_atoms,
            'n_electrons': n_electrons,
            'method': self.method,
            'method_params': method_params_str,
            'hardware': self.hardware_info,
            'n_cores': self.n_cores,
            'wall_time_sec': f"{wall_time:.2f}" if wall_time else None,
            'memory_peak_mb': f"{self.peak_memory:.1f}" if self.peak_memory else None,
            'success': self.success,
            'convergence_steps': self.convergence_steps,
            'final_energy': f"{self.final_energy:.6f}" if self.final_energy is not None else None,
            'quality_rating': self.quality_rating,
            'notes': self.notes
        }
        
        # Check if CSV exists and create header if needed
        file_exists = self.csv_path.exists()
        
        # Write to CSV
        with open(self.csv_path, 'a', newline='') as csvfile:
            fieldnames = list(row_data.keys())
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            # Write header if new file
            if not file_exists:
                writer.writeheader()
                
            # Write data row
            writer.writerow(row_data)
            
    @staticmethod
    def load_performance_data(csv_path: Optional[Path] = None) -> list:
        """Load performance data from CSV file."""
        if csv_path is None:
            # Check for environment variable override (useful for testing)
            env_path = os.environ.get('PERFORMANCE_LOG_PATH')
            if env_path:
                csv_path = Path(env_path)
            else:
                csv_path = settings.PROJECT_ROOT / "performance_logs.csv"
            
        if not csv_path.exists():
            return []
            
        with open(csv_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            return list(reader)


# Convenience function for quick tracking
def track_performance(method: str, molecule_name: str, **kwargs):
    """Convenience function to create PerformanceTracker with common defaults."""
    return PerformanceTracker(method=method, molecule_name=molecule_name, **kwargs)


if __name__ == "__main__":
    # Example usage
    print("Performance Tracker Test")
    
    # Simulate a calculation
    with track_performance("test", "benzene", method_params={"basis": "def2-SVP"}) as tracker:
        time.sleep(0.1)  # Simulate calculation
        tracker.set_quality(4, "Test run - good convergence")
        
    print(f"Logged to: {tracker.csv_path}")
    
    # Load and display data
    data = PerformanceTracker.load_performance_data()
    if data:
        print("\nRecent entries:")
        for entry in data[-3:]:  # Show last 3 entries
            print(f"  {entry['timestamp']}: {entry['method']} on {entry['molecule_name']} "
                  f"({entry['wall_time_sec']}s)")