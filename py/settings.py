"""
Settings module for mechanosynthesis pipeline

Loads configuration from .env file and provides centralized access to paths and settings.
"""
import os
from pathlib import Path

# Load environment variables from .env file in project root
project_root = Path(__file__).parent.parent
env_path = project_root / '.env'

# Simple .env loading without dotenv dependency
def load_env_file(path):
    env_vars = {}
    if path.exists():
        with open(path, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#') and '=' in line:
                    key, value = line.split('=', 1)
                    env_vars[key.strip()] = value.strip()
    return env_vars

env_vars = load_env_file(env_path)
# Update os.environ with loaded variables
os.environ.update(env_vars)

# Project structure
PROJECT_ROOT = project_root
PY_ROOT = Path(__file__).parent

# Data and documentation paths (outside py/ folder)
DATA_FOLDER = PROJECT_ROOT / os.getenv('DATA_FOLDER', 'data')
DOCS_FOLDER = PROJECT_ROOT / os.getenv('DOCS_FOLDER', 'docs')

# Python-specific paths
KERNEL_PATH = PY_ROOT / 'kernel'
TESTS_PATH = PY_ROOT / 'tests'
SCENES_PATH = DATA_FOLDER / 'scenes'

# External tool paths
DFTB_PATH = Path(os.getenv('DFTB_PATH', '/Users/andreypanferov/opt/dftb+/bin/dftb+'))
SLAKOS_PATH = Path(os.getenv('SKDIR', '/Users/andreypanferov/opt/dftb+/slakos/'))
XTB_FIXED_PATH = PROJECT_ROOT / '3rdparty' / 'xtb' / 'fork' / 'build_cmake' / 'xtb'
ORCA_PATH = Path(os.getenv('ORCA_PATH', '/Applications/orca-6.1.0/orca'))

# Cloud settings (if used)
CLOUD_PROJECT_ID = os.getenv('PROJECT_ID')
CLOUD_REGION = os.getenv('REGION', 'us-central1')
CLOUD_BUCKET = os.getenv('BUCKET_NAME')

# Analysis settings
STRUCTURAL_ANALYSIS = os.getenv('STRUCTURAL_ANALYSIS', 'true').lower() in ['true', '1', 'yes', 'on']

def get_project_data_path(project_name: str) -> Path:
    """Get data path for a specific project"""
    return DATA_FOLDER / project_name

def ensure_data_structure(project_name: str) -> dict:
    """Ensure project directory structure exists and return paths"""
    project_path = get_project_data_path(project_name)
    
    paths = {
        'project': project_path,
        'inputs': project_path / 'inputs',
        'outputs': project_path / 'outputs',
        'logs': project_path / 'logs',
        'temp': project_path / 'temp'
    }
    
    # Create directories
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    
    return paths

def validate_external_tools() -> dict:
    """Check availability of external calculation tools"""
    tools = {
        'dftb': DFTB_PATH.exists(),
        'xtb_fixed': XTB_FIXED_PATH.exists(),
        'orca': ORCA_PATH.exists(),
        'slakos': SLAKOS_PATH.exists() and SLAKOS_PATH.is_dir()
    }
    
    return tools

if __name__ == "__main__":
    print(f"Project Root: {PROJECT_ROOT}")
    print(f"Data Folder: {DATA_FOLDER}")
    print(f"Docs Folder: {DOCS_FOLDER}")
    print(f"Kernel Path: {KERNEL_PATH}")
    print("\nExternal Tools:")
    for tool, available in validate_external_tools().items():
        status = "✓" if available else "✗"
        print(f"  {status} {tool}")