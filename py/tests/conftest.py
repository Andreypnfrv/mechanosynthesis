#!/usr/bin/env python3
"""
Pytest configuration and fixtures for mechanosynthesis tests

Provides test fixtures and setup/teardown functionality to:
1. Copy test fixtures from py/tests/fixtures/ to data/ before tests
2. Clean up test data from data/ after tests
3. Ensure isolated test environment
"""
import os
import sys
import shutil
import pytest
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Test fixtures directory
FIXTURES_DIR = Path(__file__).parent / "fixtures"
# Data directory where tests run
DATA_DIR = Path(__file__).parent.parent.parent / "data"


def copy_fixture_to_data(fixture_name):
    """
    Copy a test fixture from fixtures/ to data/
    
    Args:
        fixture_name (str): Name of the fixture folder
    """
    fixture_path = FIXTURES_DIR / fixture_name
    data_path = DATA_DIR / fixture_name
    
    if not fixture_path.exists():
        raise FileNotFoundError(f"Fixture not found: {fixture_path}")
    
    # Remove existing data if it exists
    if data_path.exists():
        shutil.rmtree(data_path)
    
    # Copy fixture to data directory
    shutil.copytree(fixture_path, data_path)
    print(f"✅ Copied fixture {fixture_name} to data/")


def cleanup_test_data(project_names):
    """
    Clean up test data from data/ directory
    
    Args:
        project_names (list): List of project names to clean up
    """
    for project_name in project_names:
        project_path = DATA_DIR / project_name
        if project_path.exists():
            shutil.rmtree(project_path)
            print(f"🧹 Cleaned up test data: {project_name}")


@pytest.fixture(scope="function")
def water_test_fixture():
    """
    Fixture for water molecule test
    
    Sets up and tears down water_test project
    """
    project_name = "water_test"
    
    # Setup: Copy fixture to data directory
    copy_fixture_to_data(project_name)
    
    yield project_name
    
    # Teardown: Clean up test data
    cleanup_test_data([project_name])


@pytest.fixture(scope="function") 
def w_hexamethyl_test_fixture():
    """
    Fixture for tungsten hexamethyl test
    
    Sets up and tears down w_hexamethyl project for testing
    """
    project_name = "w_hexamethyl"
    
    # Setup: Copy fixture to data directory
    copy_fixture_to_data(project_name)
    
    yield project_name
    
    # Teardown: Clean up test data
    cleanup_test_data([project_name])


@pytest.fixture(scope="function")
def both_test_fixtures():
    """
    Fixture for tests that need both water and tungsten compounds
    
    Sets up and tears down both test projects
    """
    project_names = ["water_test", "w_hexamethyl"]
    
    # Setup: Copy both fixtures
    for project_name in project_names:
        copy_fixture_to_data(project_name)
    
    yield project_names
    
    # Teardown: Clean up all test data
    cleanup_test_data(project_names)


# Pytest hooks for session-level setup/teardown
def pytest_sessionstart(session):
    """Called after the Session object has been created"""
    print("🧪 Starting test session with fixture-based test data management")


def pytest_sessionfinish(session, exitstatus):
    """Called after whole test run finished"""
    # Final cleanup - remove any leftover test data
    test_projects = ["water_test", "w_hexamethyl"]
    cleanup_test_data(test_projects)
    print("🧹 Test session finished - all test data cleaned up")