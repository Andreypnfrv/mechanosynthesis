"""
Building functions for molecular parts

Consolidates all molecular building functionality from lib/build/ and individual build_*.py files
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Import functions that exist and handle missing ones gracefully
try:
    from ..build.build_diamond import *
except ImportError:
    pass

try:
    from ..build.build_molecule import *
except ImportError:
    pass

try:
    from ..build.build_stm_tip import *
except ImportError:
    pass

try:
    from ..build.build_adamantane import *
except ImportError:
    pass

try:
    from ..build.build_diamondoid import *
except ImportError:
    pass

# Import standalone build functions
try:
    from ..build_ge_ch2i4 import *
except ImportError:
    pass

try:
    from ..build_w_co6 import *
except ImportError:
    pass

# Only include functions that actually exist
__all__ = []