#!/usr/bin/env python3
"""
Filter NHTBH38 reactions for mechanosynthesis relevance.
Focus on reactions containing our target elements: H, C, O, N, F, Cl (and others later).
"""

import os
from pathlib import Path

# Target elements for mechanosynthesis (starting with basics)
TARGET_ELEMENTS = {'H', 'C', 'O', 'N', 'F', 'Cl'}

# High-relevance reactions for mechanosynthesis
MECHANOSYNTHESIS_REACTIONS = {
    'H_FCH3_HF_CH3': {
        'description': 'H-abstraction from carbon',
        'relevance': 'HIGH',
        'files': ['H.xyz', 'CH3F.xyz', 'HF.xyz', 'CH3.xyz'],
        'reaction': 'H + FCH3 → HF + CH3',
        'barrier_forward': 30.38,
        'barrier_reverse': 57.02
    },
    'H_C2H4_CH3CH2': {
        'description': 'C-H bond formation', 
        'relevance': 'HIGH',
        'files': ['H.xyz', 'C2H4.xyz', 'CH3CH2.xyz'],
        'reaction': 'H + C2H4 → CH3CH2',
        'barrier_forward': 1.72,
        'barrier_reverse': 41.75
    },
    'CH3_C2H4_CH3CH2CH2': {
        'description': 'C-C bond formation',
        'relevance': 'HIGH', 
        'files': ['CH3.xyz', 'C2H4.xyz', 'CH3CH2CH2.xyz'],
        'reaction': 'CH3 + C2H4 → CH3CH2CH2',
        'barrier_forward': 6.85,
        'barrier_reverse': 32.97
    },
    'F_CH3F_FCH3_F': {
        'description': 'SN2 substitution',
        'relevance': 'MEDIUM',
        'files': ['F-.xyz', 'CH3F.xyz', 'CH3F.xyz', 'F-.xyz'], 
        'reaction': 'F- + CH3F → FCH3 + F-',
        'barrier_forward': -0.34,
        'barrier_reverse': -0.34
    },
    'H_CO_HCO': {
        'description': 'Association reaction',
        'relevance': 'MEDIUM',
        'files': ['H.xyz', 'CO.xyz', 'HCO.xyz'],
        'reaction': 'H + CO → HCO', 
        'barrier_forward': 3.17,
        'barrier_reverse': 22.68
    }
}

def check_files_exist():
    """Check which reaction files are available"""
    reactions_dir = Path('reactions')
    available_reactions = {}
    
    for reaction_key, reaction_data in MECHANOSYNTHESIS_REACTIONS.items():
        files_exist = []
        for filename in reaction_data['files']:
            file_path = reactions_dir / filename
            files_exist.append(file_path.exists())
            
        if all(files_exist):
            available_reactions[reaction_key] = reaction_data
            print(f"✓ {reaction_key}: {reaction_data['reaction']}")
        else:
            missing_files = [f for f, exists in zip(reaction_data['files'], files_exist) if not exists]
            print(f"✗ {reaction_key}: Missing {missing_files}")
    
    return available_reactions

def create_subset_metadata(available_reactions):
    """Create metadata for mechanosynthesis subset"""
    subset_info = {
        'total_reactions': len(available_reactions),
        'high_relevance': sum(1 for r in available_reactions.values() if r['relevance'] == 'HIGH'),
        'medium_relevance': sum(1 for r in available_reactions.values() if r['relevance'] == 'MEDIUM'),
        'elements_covered': list(TARGET_ELEMENTS)
    }
    
    print(f"\nMechanosynthesis subset:")
    print(f"Total reactions: {subset_info['total_reactions']}")
    print(f"High relevance: {subset_info['high_relevance']}")
    print(f"Medium relevance: {subset_info['medium_relevance']}")
    print(f"Elements: {', '.join(subset_info['elements_covered'])}")
    
    return subset_info

if __name__ == "__main__":
    print("Filtering NHTBH38 for mechanosynthesis relevance...")
    available = check_files_exist()
    subset_info = create_subset_metadata(available)