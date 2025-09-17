#!/usr/bin/env python3
"""
End-to-End Test for Unified DFTB Backend

Tests the complete workflow through command line interface:
1. Water molecule (H2O) - relaxation and Hessian
2. Tungsten hexamethyl (W(CH3)6) - relaxation with heavy metals

Uses pytest fixtures for test data management:
- Test fixtures stored in py/tests/fixtures/
- Temporary data copied to data/ during tests
- Automatic cleanup after tests

USAGE:
- Pytest mode: `python -m pytest py/tests/test_unified_dftb.py` (recommended)
- Direct mode: `python py/tests/test_unified_dftb.py` (legacy, requires existing data)
"""
import os
import sys
import subprocess
import time
import glob
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def run_command_with_validation(cmd, timeout=120):
    """
    Run command with timeout and capture output
    
    Args:
        cmd (list): Command and arguments
        timeout (int): Timeout in seconds
        
    Returns:
        tuple: (success, stdout, stderr, elapsed_time)
    """
    print(f"🚀 Running: {' '.join(cmd)}")
    start_time = time.time()
    
    try:
        # Set PYTHONPATH to parent directory
        parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        env = {**os.environ, 'PYTHONPATH': parent_dir}
        
        result = subprocess.run(
            cmd,
            cwd=parent_dir,  # Run from parent directory (mechanosynthesis/py)
            capture_output=True,
            text=True,
            timeout=timeout,
            env=env
        )
        elapsed_time = time.time() - start_time
        success = result.returncode == 0
        
        if success:
            print(f"✅ Command completed successfully in {elapsed_time:.2f}s")
        else:
            print(f"❌ Command failed with return code {result.returncode}")
            if result.stderr:
                print(f"   Error: {result.stderr[:200]}...")
        
        return success, result.stdout, result.stderr, elapsed_time
        
    except subprocess.TimeoutExpired:
        elapsed_time = time.time() - start_time
        print(f"⏰ Command timed out after {elapsed_time:.2f}s")
        return False, "", f"Command timed out after {timeout}s", elapsed_time
    except Exception as e:
        elapsed_time = time.time() - start_time
        print(f"💥 Command failed with exception: {e}")
        return False, "", str(e), elapsed_time


def validate_calculation_success(project_name, calc_type, expected_backend="ptbp_dftb"):
    """
    Validate that calculation completed successfully
    
    Args:
        project_name (str): Name of project folder
        calc_type (str): 'relax' or 'hessian'
        expected_backend (str): Expected output directory name
        
    Returns:
        bool: True if validation passed
    """
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.join(parent_dir, "data", project_name)
    output_dir = os.path.join(base_path, "outputs", expected_backend)
    
    print(f"🔍 Validating {calc_type} results for {project_name}...")
    
    # Check if output directory exists
    if not os.path.exists(output_dir):
        print(f"❌ Output directory not found: {output_dir}")
        return False
    
    if calc_type == 'relax':
        # Look for relaxed structure
        relaxed_files = glob.glob(os.path.join(output_dir, "*_relaxed.xyz"))
        if not relaxed_files:
            print(f"❌ No relaxed structure found in {output_dir}")
            return False
        
        # Check file size
        relaxed_file = relaxed_files[0]
        if os.path.getsize(relaxed_file) == 0:
            print(f"❌ Relaxed file is empty: {relaxed_file}")
            return False
        
        print(f"✅ Found relaxed structure: {os.path.basename(relaxed_file)} ({os.path.getsize(relaxed_file)} bytes)")
        
        # Check for structural analysis files
        analysis_dir = os.path.join(output_dir, "structural_analysis")
        if os.path.exists(analysis_dir):
            # Check for required analysis files
            analysis_files = {
                'json': glob.glob(os.path.join(analysis_dir, "*_structural_analysis_*.json")),
                'txt': glob.glob(os.path.join(analysis_dir, "*_structural_report_*.txt")),
                'png': glob.glob(os.path.join(analysis_dir, "*_structural_changes_*.png"))
            }
            
            missing_files = []
            for file_type, files in analysis_files.items():
                if not files:
                    missing_files.append(file_type)
                else:
                    # Check that files are not empty
                    for file_path in files:
                        if os.path.getsize(file_path) == 0:
                            missing_files.append(f"{file_type} (empty)")
            
            if missing_files:
                print(f"⚠️ Structural analysis incomplete. Missing: {', '.join(missing_files)}")
            else:
                print(f"✅ Found complete structural analysis files in {os.path.basename(analysis_dir)}/")
        else:
            print(f"⚠️ Structural analysis directory not found: {analysis_dir}")
        
    elif calc_type == 'hessian':
        # Look for Hessian output
        hessian_files = glob.glob(os.path.join(output_dir, "*_hessian.out"))
        if not hessian_files:
            print(f"❌ No Hessian file found in {output_dir}")
            return False
        
        # Check file size
        hessian_file = hessian_files[0]
        if os.path.getsize(hessian_file) == 0:
            print(f"❌ Hessian file is empty: {hessian_file}")
            return False
        
        print(f"✅ Found Hessian output: {os.path.basename(hessian_file)} ({os.path.getsize(hessian_file)} bytes)")
    
    return True


def setup_water_test():
    """Legacy setup function - now handled by pytest fixtures"""
    # This function is kept for backward compatibility with direct test execution
    # When run via pytest, the water_test_fixture handles setup/teardown
    project_name = "water_test"
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.join(parent_dir, "data", project_name)
    
    if not os.path.exists(base_path):
        print(f"⚠️  Water test project not found at {base_path}")
        print(f"   Run via pytest to use automatic fixture management")
        return None
    
    print(f"✅ Found water test project: {project_name}")
    return project_name


def test_water_e2e(project_name=None):
    """End-to-end test of water molecule"""
    print("\n" + "="*60)
    print("🌊 WATER MOLECULE E2E TEST")
    print("="*60)
    
    # Handle both pytest fixture and direct execution modes
    if project_name is None:
        project_name = setup_water_test()
        if project_name is None:
            return False
    
    # Test 1: Water relaxation
    print(f"\n1️⃣ Testing water relaxation...")
    cmd = ['python', 'main.py', 'relax', project_name, '--backend', 'dftb']
    success, stdout, stderr, elapsed = run_command_with_validation(cmd, timeout=60)
    
    if not success:
        print(f"❌ Water relaxation failed!")
        print(f"STDOUT: {stdout[:500]}...")
        print(f"STDERR: {stderr[:500]}...")
        return False
    
    # Check for time prediction in output
    if "Оценка времени оптимизации:" not in stdout:
        print(f"❌ Time prediction not found in relaxation output!")
        print(f"STDOUT: {stdout[:1000]}...")
        return False
    print(f"✅ Time prediction was performed during relaxation")
    
    # Check for structural analysis in output
    if "RUNNING STRUCTURAL ANALYSIS" not in stdout:
        print(f"❌ Structural analysis not found in relaxation output!")
        print(f"STDOUT: {stdout[:1000]}...")
        return False
    print(f"✅ Structural analysis was performed during relaxation")
    
    # Validate relaxation results
    if not validate_calculation_success(project_name, 'relax'):
        return False
    
    # Test 2: Water Hessian
    print(f"\n2️⃣ Testing water Hessian...")
    cmd = ['python', 'main.py', 'hessian', project_name, '--backend', 'dftb']
    success, stdout, stderr, elapsed = run_command_with_validation(cmd, timeout=120)
    
    if not success:
        print(f"❌ Water Hessian failed!")
        print(f"STDOUT: {stdout[:500]}...")
        print(f"STDERR: {stderr[:500]}...")
        return False
    
    # Validate Hessian results
    if not validate_calculation_success(project_name, 'hessian'):
        return False
    
    print(f"🎉 Water E2E test completed successfully!")
    return True


def test_tungsten_e2e(project_name=None):
    """End-to-end test of tungsten hexamethyl"""
    print("\n" + "="*60)
    print("⚡ TUNGSTEN HEXAMETHYL E2E TEST")
    print("="*60)
    
    # Handle both pytest fixture and direct execution modes
    if project_name is None:
        project_name = "w_hexamethyl"
    
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    base_path = os.path.join(parent_dir, "data", project_name)
    
    # Check if project exists
    if not os.path.exists(os.path.join(base_path, "inputs")):
        print(f"❌ Project {project_name} not found at {base_path}")
        print("   When run via pytest, this should be automatically set up from fixtures")
        return False
    
    # Check for input files
    input_files = glob.glob(os.path.join(base_path, "inputs", "*.xyz"))
    if not input_files:
        print(f"❌ No input XYZ files found in {base_path}/inputs/")
        return False
    
    print(f"✅ Found project with input: {os.path.basename(input_files[0])}")
    
    # Test: Tungsten relaxation (should auto-select PTBP for heavy metals)
    print(f"\n1️⃣ Testing tungsten compound relaxation...")
    cmd = ['python', 'main.py', 'relax', project_name, '--backend', 'dftb']
    success, stdout, stderr, elapsed = run_command_with_validation(cmd, timeout=120)
    
    if not success:
        print(f"❌ Tungsten relaxation failed!")
        print(f"STDOUT: {stdout[:500]}...")
        print(f"STDERR: {stderr[:500]}...")
        return False
    
    # Check for time prediction in output
    if "Оценка времени оптимизации:" not in stdout:
        print(f"❌ Time prediction not found in relaxation output!")
        print(f"STDOUT: {stdout[:1000]}...")
        return False
    print(f"✅ Time prediction was performed during relaxation")
    
    # Check for structural analysis in output
    if "RUNNING STRUCTURAL ANALYSIS" not in stdout:
        print(f"❌ Structural analysis not found in relaxation output!")
        print(f"STDOUT: {stdout[:1000]}...")
        return False
    print(f"✅ Structural analysis was performed during relaxation")
    
    # Validate relaxation results
    if not validate_calculation_success(project_name, 'relax'):
        return False
    
    # Check that PTBP parameters were used (should be automatic for tungsten)
    ptbp_output_dir = os.path.join(base_path, "outputs", "ptbp_dftb")
    if not os.path.exists(ptbp_output_dir):
        print(f"❌ Expected PTBP output directory not found: {ptbp_output_dir}")
        print("   Auto-detection of heavy metals may have failed")
        return False
    
    print(f"✅ Confirmed PTBP parameters were auto-selected for tungsten compound")
    print(f"🎉 Tungsten E2E test completed successfully!")
    return True


def cleanup_test_data():
    """Legacy cleanup function - now handled by pytest fixtures"""
    # This function is kept for backward compatibility with direct test execution
    # When run via pytest, fixtures handle cleanup automatically
    print(f"\n🧹 Cleaning up test data...")
    
    # Remove water test project
    parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    water_path = os.path.join(parent_dir, "data", "water_test")
    if os.path.exists(water_path):
        import shutil
        shutil.rmtree(water_path)
        print(f"✅ Removed {water_path}")
    
    # Remove tungsten test project (if it was created from fixture)
    tungsten_path = os.path.join(parent_dir, "data", "w_hexamethyl")
    # Only clean up if it looks like test data (check if it has only inputs/)
    if os.path.exists(tungsten_path):
        contents = os.listdir(tungsten_path)
        if contents == ['inputs'] or contents == []:
            import shutil
            shutil.rmtree(tungsten_path)
            print(f"✅ Removed test fixture: {tungsten_path}")


# Pytest functions for test framework integration
def test_unified_dftb_water(water_test_fixture):
    """Pytest wrapper for water E2E test"""
    project_name = water_test_fixture  # Fixture provides project name and handles setup/teardown
    success = test_water_e2e(project_name)
    assert success, "Water E2E test failed"


def test_unified_dftb_tungsten(w_hexamethyl_test_fixture):
    """Pytest wrapper for tungsten E2E test"""
    project_name = w_hexamethyl_test_fixture  # Fixture provides project name and handles setup/teardown
    success = test_tungsten_e2e(project_name)
    assert success, "Tungsten E2E test failed"


def run_all_tests():
    """Run all unified DFTB tests"""
    print("🧪 UNIFIED DFTB BACKEND - END-TO-END TESTS")
    print("=" * 70)
    print("Testing complete DFTB workflow through command line interface")
    print("=" * 70)
    
    all_tests_passed = True
    start_time = time.time()
    
    try:
        # Test 1: Water molecule (organic compound)
        if not test_water_e2e():
            all_tests_passed = False
        
        # Test 2: Tungsten compound (heavy metal) 
        if not test_tungsten_e2e():
            all_tests_passed = False
        
    except KeyboardInterrupt:
        print(f"\n⏹️  Tests interrupted by user")
        all_tests_passed = False
    except Exception as e:
        print(f"\n💥 Unexpected error during testing: {e}")
        import traceback
        traceback.print_exc()
        all_tests_passed = False
    finally:
        # Cleanup (for direct execution mode)
        cleanup_test_data()
    
    # Final results
    total_time = time.time() - start_time
    print(f"\n" + "="*70)
    if all_tests_passed:
        print(f"🎉 ALL TESTS PASSED! Total time: {total_time:.2f}s")
        print(f"✅ Unified DFTB backend is working correctly")
        print(f"✅ Auto-parameter selection is functional")
        print(f"✅ Both organic and heavy metal compounds are supported")
    else:
        print(f"❌ SOME TESTS FAILED! Total time: {total_time:.2f}s")
        print(f"🔧 Check the error messages above for details")
        sys.exit(1)
    
    print("="*70)
    return all_tests_passed


if __name__ == "__main__":
    success = run_all_tests()
    if not success:
        sys.exit(1)