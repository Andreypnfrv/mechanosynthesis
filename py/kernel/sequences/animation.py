"""
Animation and keyframe interpolation module

Handles creating animations from scene sequences and trajectories
"""
import numpy as np
from typing import List, Dict, Any, Optional
from ase import Atoms
from ase.io import write
import yaml


class KeyFrame:
    """Single keyframe in an animation sequence"""
    
    def __init__(self, time: float, atoms: Atoms, metadata: Optional[Dict] = None):
        self.time = time
        self.atoms = atoms
        self.metadata = metadata or {}


class SceneAnimator:
    """Handles animation of molecular scenes"""
    
    def __init__(self, name: str):
        self.name = name
        self.keyframes: List[KeyFrame] = []
        
    def add_keyframe(self, time: float, atoms: Atoms, metadata: Optional[Dict] = None):
        """Add a keyframe to the animation"""
        keyframe = KeyFrame(time, atoms, metadata)
        self.keyframes.append(keyframe)
        
        # Keep keyframes sorted by time
        self.keyframes.sort(key=lambda x: x.time)
    
    def interpolate_linear(self, start_frame: KeyFrame, end_frame: KeyFrame, 
                          num_steps: int) -> List[Atoms]:
        """Linear interpolation between two keyframes"""
        if len(start_frame.atoms) != len(end_frame.atoms):
            raise ValueError("Cannot interpolate between frames with different atom counts")
        
        start_pos = start_frame.atoms.get_positions()
        end_pos = end_frame.atoms.get_positions()
        symbols = start_frame.atoms.get_chemical_symbols()
        cell = start_frame.atoms.get_cell()
        
        interpolated_frames = []
        
        for i in range(num_steps):
            t = i / (num_steps - 1) if num_steps > 1 else 0
            
            # Linear interpolation of positions
            interp_pos = start_pos + t * (end_pos - start_pos)
            
            frame_atoms = Atoms(symbols=symbols, positions=interp_pos, cell=cell)
            interpolated_frames.append(frame_atoms)
        
        return interpolated_frames
    
    def generate_trajectory(self, fps: int = 24, total_duration: Optional[float] = None) -> List[Atoms]:
        """Generate complete animation trajectory"""
        if len(self.keyframes) < 2:
            raise ValueError("Need at least 2 keyframes for animation")
        
        trajectory = []
        
        for i in range(len(self.keyframes) - 1):
            start_frame = self.keyframes[i]
            end_frame = self.keyframes[i + 1]
            
            time_diff = end_frame.time - start_frame.time
            num_steps = max(1, int(time_diff * fps))
            
            interpolated = self.interpolate_linear(start_frame, end_frame, num_steps)
            
            # Don't duplicate the end frame (except for the last segment)
            if i < len(self.keyframes) - 2:
                interpolated = interpolated[:-1]
            
            trajectory.extend(interpolated)
        
        return trajectory
    
    def save_animation(self, output_path: str, format: str = 'xyz', fps: int = 24):
        """Save animation to file"""
        trajectory = self.generate_trajectory(fps=fps)
        
        if format.lower() == 'xyz':
            write(output_path, trajectory)
        else:
            raise ValueError(f"Unsupported animation format: {format}")
        
        print(f"Animation '{self.name}' saved to {output_path} ({len(trajectory)} frames)")


def load_animation_yaml(yaml_path: str) -> SceneAnimator:
    """Load animation definition from YAML file"""
    with open(yaml_path, 'r') as f:
        anim_data = yaml.safe_load(f)
    
    animator = SceneAnimator(anim_data['name'])
    
    # This is a placeholder - full implementation would handle:
    # - Loading keyframe definitions
    # - Interpolation settings
    # - Scene calculations at specific keyframes
    
    return animator


def create_approach_animation(tip_atoms: Atoms, surface_atoms: Atoms, 
                            approach_distance: float = 5.0,
                            contact_distance: float = 2.0,
                            num_frames: int = 50) -> SceneAnimator:
    """Create a simple approach animation between tip and surface"""
    
    animator = SceneAnimator("tip_approach")
    
    # Initial position - tip far from surface
    tip_high = tip_atoms.copy()
    tip_high.translate([0, 0, approach_distance])
    
    # Combine tip and surface
    initial_scene = tip_high + surface_atoms
    final_scene = tip_atoms.copy()
    final_scene.translate([0, 0, contact_distance])
    final_scene = final_scene + surface_atoms
    
    # Add keyframes
    animator.add_keyframe(0.0, initial_scene, {'description': 'initial_approach'})
    animator.add_keyframe(1.0, final_scene, {'description': 'contact'})
    
    return animator


# Utility functions for animation analysis
def analyze_trajectory_motion(trajectory: List[Atoms]) -> Dict[str, Any]:
    """Analyze motion patterns in a trajectory"""
    if len(trajectory) < 2:
        return {}
    
    # Calculate per-atom displacements
    initial_pos = trajectory[0].get_positions()
    final_pos = trajectory[-1].get_positions()
    
    displacements = np.linalg.norm(final_pos - initial_pos, axis=1)
    
    analysis = {
        'total_frames': len(trajectory),
        'max_displacement': displacements.max(),
        'mean_displacement': displacements.mean(),
        'highly_mobile_atoms': np.where(displacements > displacements.mean() + 2*displacements.std())[0].tolist()
    }
    
    return analysis