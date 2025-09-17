"""
Sequences Module - Scene composition, chaining, and animation

This module handles:
- Combining relaxed parts into scenes
- Positioning and orienting molecular pieces
- Scene calculations and optimization
- Artifact detection and extraction
- Animation and keyframe interpolation
- Quality assurance for multi-component systems
"""

from .scene import *
from .artifacts import *
from .animation import *

__all__ = [
    # Scene building
    'Scene', 'load_scene_yaml', 'position_pieces',
    # Artifact detection
    'detect_artifacts', 'extract_fragments', 'analyze_connectivity',
    # Animation
    'animate_sequence', 'interpolate_keyframes', 'render_animation'
]