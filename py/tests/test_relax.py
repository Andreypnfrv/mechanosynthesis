import pytest
import subprocess
import os
import tempfile
import shutil
import csv
from pathlib import Path


def create_h2o_xyz():
    """Create H2O molecule XYZ file for testing."""
    h2o_content = """3
Water molecule
O    0.000000    0.000000    0.000000
H    0.757000    0.587000    0.000000  
H   -0.757000    0.587000    0.000000
"""
    return h2o_content


class TestRelax:
    def setup_method(self):
        """Set up test environment before each test."""
        self.original_dir = os.getcwd()
        
        # Create test project in data/ directory
        self.project_name = "test_h2o_relax"
        self.project_path = Path(self.original_dir) / "data" / self.project_name
        self.project_path.mkdir(parents=True, exist_ok=True)
        
        # Create inputs directory
        inputs_dir = self.project_path / "inputs"
        inputs_dir.mkdir(exist_ok=True)
        
        # Create test input file in inputs directory
        self.h2o_file = inputs_dir / "h2o.xyz"
        with open(self.h2o_file, 'w') as f:
            f.write(create_h2o_xyz())
        
        # Create temporary analytics file for tests (don't pollute main one)
        self.temp_analytics_file = Path(self.original_dir) / "test_performance_logs.csv"
        
        # Backup main analytics file if it exists (copy, don't move!)
        self.main_analytics_file = Path(self.original_dir) / "performance_logs.csv"
        self.main_analytics_backup = None
        if self.main_analytics_file.exists():
            self.main_analytics_backup = Path(self.original_dir) / "performance_logs.csv.test_backup"
            shutil.copy2(self.main_analytics_file, self.main_analytics_backup)
    
    def teardown_method(self):
        """Clean up after each test."""
        os.chdir(self.original_dir)
        
        # Clean up test project
        if self.project_path.exists():
            shutil.rmtree(self.project_path)
        
        # Clean up temporary analytics file
        if self.temp_analytics_file.exists():
            self.temp_analytics_file.unlink()
        
        # Clean up backup file (main file should remain untouched)
        if self.main_analytics_backup and self.main_analytics_backup.exists():
            self.main_analytics_backup.unlink()
        
        # If main analytics file was modified during test, restore from backup
        # But since we use PERFORMANCE_LOG_PATH, this shouldn't happen
    
    def test_relax_dftb(self):
        """Test H2O relaxation using DFTB+ backend."""
        # Use our temporary analytics file for testing
        performance_log = self.temp_analytics_file
        
        cmd = [
            "python", "py/main.py", "relax", self.project_name,
            "--backend", "dftb"
        ]
        
        # Set environment variable to use temporary analytics file
        env = os.environ.copy()
        env['PERFORMANCE_LOG_PATH'] = str(performance_log)
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.original_dir, env=env)
        
        # Check that command succeeded
        assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"
        
        # Check that output files were created (DFTB uses "ptbp_dftb" directory)
        output_dir = self.project_path / "outputs" / "ptbp_dftb"
        assert output_dir.exists(), "Output directory was not created"
        
        # Check for expected output files
        expected_files = ["h2o_relaxed.xyz"]
        for filename in expected_files:
            filepath = output_dir / filename
            assert filepath.exists(), f"Expected output file {filename} was not created"
            assert filepath.stat().st_size > 0, f"Output file {filename} is empty"
        
        # Verify relaxed structure has correct number of atoms
        relaxed_file = output_dir / "h2o_relaxed.xyz"
        with open(relaxed_file, 'r') as f:
            lines = f.readlines()
            assert len(lines) >= 5, "Relaxed XYZ file has too few lines"
            assert lines[0].strip() == "3", "Relaxed structure should have 3 atoms"
            
        # Check that performance analytics were collected
        assert performance_log.exists(), "Performance analytics file was not created"
        assert performance_log.stat().st_size > 0, "Performance analytics file is empty"
        
        # Verify analytics content
        with open(performance_log, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) > 0, "No performance data was logged"
            
            latest_row = rows[-1]  # Get the last (most recent) entry
            assert latest_row['method'] == 'relax-dftb', f"Expected method 'relax-dftb', got '{latest_row['method']}'"
            assert latest_row['molecule_name'] == self.project_name, f"Expected molecule name '{self.project_name}', got '{latest_row['molecule_name']}'"
            assert float(latest_row['wall_time_sec']) > 0, "Duration should be positive"
            assert latest_row['success'] == 'True', "Relaxation should be marked as successful"
    
    def _test_relax_xtb(self):
        """Test H2O relaxation using xTB backend."""
        performance_log = Path(self.original_dir) / "performance_logs.csv"
        
        # Remove existing performance log to get clean test
        if performance_log.exists():
            performance_log.unlink()
            
        cmd = [
            "python", "py/main.py", "relax", self.project_name,
            "--backend", "xtb", "--xtb-backend", "native"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=self.original_dir)
        
        # Debug: print command output
        print(f"\nCommand: {' '.join(cmd)}")
        print(f"Return code: {result.returncode}")
        print(f"STDOUT:\n{result.stdout}")
        print(f"STDERR:\n{result.stderr}")
        
        # Check that command succeeded
        assert result.returncode == 0, f"Command failed with stderr: {result.stderr}"
        
        # Debug: check what files/dirs were actually created
        print(f"\nProject path contents:")
        if self.project_path.exists():
            for item in self.project_path.rglob("*"):
                print(f"  {item}")
        
        # Check that output files were created
        output_dir = self.project_path / "outputs" / "xtb"
        assert output_dir.exists(), "Output directory was not created"
        
        # Check for expected output files
        expected_files = ["h2o_relaxed.xyz", "h2o_trajectory.xyz"]
        for filename in expected_files:
            filepath = output_dir / filename
            assert filepath.exists(), f"Expected output file {filename} was not created"
            assert filepath.stat().st_size > 0, f"Output file {filename} is empty"
        
        # Verify relaxed structure has correct number of atoms
        relaxed_file = output_dir / "h2o_relaxed.xyz"
        with open(relaxed_file, 'r') as f:
            lines = f.readlines()
            assert len(lines) >= 5, "Relaxed XYZ file has too few lines"
            assert lines[0].strip() == "3", "Relaxed structure should have 3 atoms"
            
        # Check that performance analytics were collected
        assert performance_log.exists(), "Performance analytics file was not created"
        assert performance_log.stat().st_size > 0, "Performance analytics file is empty"
        
        # Verify analytics content
        with open(performance_log, 'r') as f:
            reader = csv.DictReader(f)
            rows = list(reader)
            assert len(rows) > 0, "No performance data was logged"
            
            latest_row = rows[-1]  # Get the last (most recent) entry
            assert latest_row['method'] == 'relax-xtb', f"Expected method 'relax-xtb', got '{latest_row['method']}'"
            assert latest_row['molecule_name'] == self.project_name, f"Expected molecule name '{self.project_name}', got '{latest_row['molecule_name']}'"
            assert float(latest_row['wall_time_sec']) > 0, "Duration should be positive"
            assert latest_row['success'] == 'True', "Relaxation should be marked as successful"