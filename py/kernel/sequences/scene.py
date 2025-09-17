"""
Scene composition and calculation module

Handles combining relaxed molecular parts into scenes for calculation
"""
import yaml
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from pathlib import Path
from ase import Atoms
from ase.io import read, write
import settings
from kernel.project_structure import create_project_structure


class Scene:
    """Scene builder for combining molecular parts"""
    
    def __init__(self, name: str):
        self.name = name
        self.pieces: Dict[str, Dict[str, Any]] = {}
        self.assembled_atoms: Optional[Atoms] = None
        
    def add_piece(self, name: str, source_path: str, position: Optional[List[float]] = None):
        """Add a molecular piece to the scene"""
        self.pieces[name] = {
            'source': source_path,
            'position': position or [0, 0, 0],
            'rotation': [0, 0, 0],
            'atoms': None  # Will be loaded when needed
        }
        
    def position_piece(self, piece_name: str, position: List[float]):
        """Set position for a piece"""
        if piece_name in self.pieces:
            self.pieces[piece_name]['position'] = position
        else:
            raise ValueError(f"Piece '{piece_name}' not found in scene")
    
    def _apply_rotation(self, positions: np.ndarray, rotation: List[float]) -> np.ndarray:
        """Apply rotation to positions (Euler angles in degrees)"""
        rx, ry, rz = np.radians(rotation)
        
        # Rotation matrices
        Rx = np.array([[1, 0, 0],
                       [0, np.cos(rx), -np.sin(rx)],
                       [0, np.sin(rx), np.cos(rx)]])
        
        Ry = np.array([[np.cos(ry), 0, np.sin(ry)],
                       [0, 1, 0],
                       [-np.sin(ry), 0, np.cos(ry)]])
        
        Rz = np.array([[np.cos(rz), -np.sin(rz), 0],
                       [np.sin(rz), np.cos(rz), 0],
                       [0, 0, 1]])
        
        # Combined rotation matrix (ZYX order)
        R = Rz @ Ry @ Rx
        
        # Center positions for rotation
        center = positions.mean(axis=0)
        centered_pos = positions - center
        
        # Apply rotation and recenter
        rotated_pos = (R @ centered_pos.T).T + center
        
        return rotated_pos
    
    def load_pieces(self):
        """Load atomic structures from source files"""
        for name, piece_info in self.pieces.items():
            if piece_info['atoms'] is None:
                try:
                    piece_info['atoms'] = read(piece_info['source'])
                except Exception as e:
                    raise FileNotFoundError(f"Could not load piece '{name}' from {piece_info['source']}: {e}")
    
    def assemble(self) -> Atoms:
        """Assemble all pieces into a single Atoms object"""
        self.load_pieces()
        
        if not self.pieces:
            raise ValueError("No pieces added to scene")
        
        all_positions = []
        all_symbols = []
        
        for name, piece_info in self.pieces.items():
            atoms = piece_info['atoms']
            positions = atoms.get_positions()
            symbols = atoms.get_chemical_symbols()
            
            # Apply rotation first (around origin)
            rotation = piece_info.get('rotation', [0, 0, 0])
            if any(r != 0 for r in rotation):
                positions = self._apply_rotation(positions, rotation)
            
            # Apply translation
            translation = np.array(piece_info['position'])
            positions += translation
            
            all_positions.extend(positions)
            all_symbols.extend(symbols)
        
        self.assembled_atoms = Atoms(symbols=all_symbols, positions=all_positions)
        
        # Set reasonable cell size
        pos_array = np.array(all_positions)
        margin = 5.0  # Å
        min_coords = pos_array.min(axis=0) - margin
        max_coords = pos_array.max(axis=0) + margin
        cell_size = max_coords - min_coords
        
        self.assembled_atoms.set_cell(cell_size)
        self.assembled_atoms.center()
        
        return self.assembled_atoms
    
    def save(self, output_path: str):
        """Save assembled scene to file"""
        if self.assembled_atoms is None:
            self.assemble()
        
        write(output_path, self.assembled_atoms)
        print(f"Scene '{self.name}' saved to {output_path}")


def load_scene_yaml(yaml_path: str) -> Scene:
    """Load scene definition from YAML file"""
    with open(yaml_path, 'r') as f:
        scene_data = yaml.safe_load(f)
    
    scene = Scene(scene_data['name'])
    
    for piece_data in scene_data.get('pieces', []):
        scene.add_piece(
            piece_data['name'],
            piece_data['source'],
            piece_data.get('position', [0, 0, 0])
        )
    
    return scene


def validate_scene_geometry(atoms: Atoms, min_distance: float = 0.8) -> List[str]:
    """Validate that scene geometry is reasonable"""
    warnings = []
    positions = atoms.get_positions()
    
    # Check for atom overlaps
    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            dist = np.linalg.norm(positions[i] - positions[j])
            if dist < min_distance:
                warnings.append(f"Atoms {i} and {j} too close: {dist:.3f} Å")
    
    return warnings


class SceneGenerator:
    """Scene generator from YAML animation definitions"""
    
    def __init__(self, yaml_path: str):
        self.yaml_path = yaml_path
        with open(yaml_path, 'r') as f:
            self.scene_data = yaml.safe_load(f)
        self.name = self.scene_data['name']
        
    def generate_frame(self, frame_number: int) -> bool:
        """Generate a specific frame from the animation"""
        try:
            # Find the keyframe
            keyframes = self.scene_data.get('animation', {}).get('keyframes', [])
            target_frame = None
            
            for kf in keyframes:
                if kf['frame'] == frame_number:
                    target_frame = kf
                    break
            
            if target_frame is None:
                print(f"❌ Frame {frame_number} not found in animation")
                return False
            
            print(f"🔧 Processing frame {frame_number}: {target_frame.get('description', 'unnamed')}")
            
            # Create scene
            scene = Scene(f"{self.name}_frame_{frame_number}")
            
            # Load pieces with frame-specific positions
            pieces_config = self.scene_data.get('pieces', {})
            frame_pieces = target_frame.get('pieces', {})
            
            for piece_name, piece_config in pieces_config.items():
                source_path = piece_config['source']
                
                # Get position for this frame (use frame-specific or default)
                if piece_name in frame_pieces:
                    position = frame_pieces[piece_name].get('position', [0, 0, 0])
                    rotation = frame_pieces[piece_name].get('rotation', [0, 0, 0])
                else:
                    position = [0, 0, 0]
                    rotation = [0, 0, 0]
                
                print(f"   Adding {piece_name}: {source_path} at {position} rotation {rotation}")
                scene.add_piece(piece_name, source_path, position)
                scene.pieces[piece_name]['rotation'] = rotation
            
            # Assemble the scene
            assembled_atoms = scene.assemble()
            
            # Validate geometry
            warnings = validate_scene_geometry(assembled_atoms)
            if warnings:
                print("⚠️  Geometry warnings:")
                for warning in warnings[:5]:  # Show first 5 warnings
                    print(f"   {warning}")
                if len(warnings) > 5:
                    print(f"   ... and {len(warnings) - 5} more warnings")
            
            # Create scene directory structure
            scene_dir = Path(f"data/scenes/{self.name}")
            frame_dir = scene_dir / f"frame_{frame_number}"
            frame_dir.mkdir(parents=True, exist_ok=True)
            
            # Save the frame XYZ
            output_file = frame_dir / f"{self.name}_frame_{frame_number}.xyz"
            scene.save(str(output_file))
            
            # Save frame info
            frame_info = {
                'frame': frame_number,
                'description': target_frame.get('description', 'unnamed'),
                'pieces': frame_pieces,
                'calculation': target_frame.get('calculation', {}),
                'n_atoms': len(assembled_atoms),
                'warnings': warnings
            }
            
            info_file = frame_dir / f"frame_{frame_number}_info.yaml"
            with open(info_file, 'w') as f:
                yaml.dump(frame_info, f, default_flow_style=False)
            
            print(f"✅ Frame saved to {output_file}")
            print(f"📄 Frame info saved to {info_file}")
            print(f"🗂️  Scene dir: {scene_dir}")
            
            return True
            
        except Exception as e:
            print(f"❌ Error generating frame {frame_number}: {e}")
            return False