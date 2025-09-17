from ase.build import diamond111
from ase.neighborlist import NeighborList
import numpy as np
from scipy.spatial.distance import cdist

def main():
    diamond = diamond111('C', size=(5, 5, 5), vacuum=6.0)


def find_atoms_with_hanging_bonds(atoms, cutoff=1.8, expected_coordination=None):
    """
    Find atoms with hanging bonds (fewer neighbors than expected).
    
    Args:
        atoms: ASE Atoms object
        cutoff: Distance cutoff for neighbor detection
        expected_coordination: Dict mapping atomic numbers to expected coordination numbers.
                             If None, uses default values for common elements.
    """
    if expected_coordination is None:
        # Default coordination numbers for common elements
        expected_coordination = {
            1: 1,   # H
            6: 4,   # C
            14: 4,  # Si
            32: 4,  # Ge
            7: 3,   # N
            8: 2,   # O
            15: 3,  # P
            16: 2,  # S
        }
    
    print("Total atoms: ", len(atoms))
    
    cutoffs = [cutoff/2] * len(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    
    hanging_atoms = []
    for i, atom in enumerate(atoms):
        atomic_number = atom.number
        if atomic_number in expected_coordination:
            indices, offsets = nl.get_neighbors(i)
            num_neighbors = len(indices)
            expected = expected_coordination[atomic_number]
            if num_neighbors < expected:
                hanging_atoms.append(i)
    
    print(f"Найдено атомов с висячими связями: {len(hanging_atoms)}")
    return hanging_atoms


def get_bond_vectors(atoms, atom_index, cutoff=1.8):
    """
    Получить векторы существующих связей для атома.
    
    Args:
        atoms: ASE Atoms object
        atom_index: Индекс атома
        cutoff: Расстояние для определения связей
        
    Returns:
        list: Список векторов связей (направленных от атома к соседям)
    """
    cutoffs = [cutoff/2] * len(atoms)
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)
    
    indices, offsets = nl.get_neighbors(atom_index)
    atom_pos = atoms.positions[atom_index]
    
    bond_vectors = []
    for neighbor_idx, offset in zip(indices, offsets):
        neighbor_pos = atoms.positions[neighbor_idx] + np.dot(offset, atoms.get_cell())
        bond_vector = neighbor_pos - atom_pos
        bond_vectors.append(bond_vector)
    
    return bond_vectors


def get_tetrahedral_passivation_vectors(bond_vectors, bond_length=1.09, surface_normal=None):
    """
    Определить оптимальные векторы пассивации на основе тетраэдрической геометрии.
    
    Args:
        bond_vectors: Список существующих векторов связей
        bond_length: Длина связи для пассивации (например, C-H = 1.09 Å)
        surface_normal: Вектор нормали к поверхности (для предотвращения проникновения в кристалл)
        
    Returns:
        list: Список векторов для размещения атомов пассивации
    """
    # Идеальные тетраэдрические направления (нормализованные)
    tetrahedral_dirs = np.array([
        [1, 1, 1],
        [1, -1, -1],
        [-1, 1, -1],
        [-1, -1, 1]
    ])
    tetrahedral_dirs = tetrahedral_dirs / np.linalg.norm(tetrahedral_dirs, axis=1, keepdims=True)
    
    if len(bond_vectors) == 0:
        # Если нет связей, возвращаем все тетраэдрические направления
        valid_dirs = tetrahedral_dirs
        if surface_normal is not None:
            # Фильтруем направления, указывающие наружу от поверхности
            surface_normal = surface_normal / np.linalg.norm(surface_normal)
            valid_dirs = [d for d in tetrahedral_dirs if np.dot(d, surface_normal) > 0.3]  # увеличен порог
        return [d * bond_length for d in valid_dirs[:4]]
    
    # Нормализуем существующие векторы связей
    bond_vectors = np.array(bond_vectors)
    bond_dirs = bond_vectors / np.linalg.norm(bond_vectors, axis=1, keepdims=True)
    
    # Улучшенная нормаль поверхности: средневзвешенная от существующих связей
    if surface_normal is None and len(bond_vectors) > 0:
        # Вычисляем среднее направление существующих связей и инвертируем
        mean_bond_dir = np.mean(bond_dirs, axis=0)
        surface_normal = -mean_bond_dir / np.linalg.norm(mean_bond_dir)
    
    # Найдем неиспользованные тетраэдрические направления
    passivation_vectors = []
    
    # Определяем количество недостающих связей
    expected_bonds = 4  # Для углерода
    needed_passivations = expected_bonds - len(bond_vectors)
    
    # Если пассивация не нужна или уже достаточно связей
    if needed_passivations <= 0:
        return []
    
    # Сортируем тетраэдрические направления по совместимости с нормалью поверхности
    scored_dirs = []
    for tet_dir in tetrahedral_dirs:
        score = 0
        
        # Проверяем направление относительно нормали поверхности
        if surface_normal is not None:
            surface_normal_norm = surface_normal / np.linalg.norm(surface_normal)
            surface_compatibility = np.dot(tet_dir, surface_normal_norm)
            if surface_compatibility <= 0.3:  # увеличен порог с 0.1 до 0.3
                continue
            score += surface_compatibility
        
        # Проверим, не занято ли это направление существующей связью
        # Используем порог угла 60 градусов (π/3 ≈ 1.05 радиан)
        min_angle = np.min([np.arccos(np.clip(np.dot(tet_dir, bond_dir), -1, 1)) 
                           for bond_dir in bond_dirs])
        
        # Если угол больше 60 градусов, направление свободно
        if min_angle > 1.05:  # увеличен с 0.52 (30°) до 1.05 (60°)
            score += min_angle  # добавляем угол в качестве бонуса
            scored_dirs.append((score, tet_dir))
    
    # Сортируем по убыванию счета и берем лучшие направления
    scored_dirs.sort(key=lambda x: x[0], reverse=True)
    
    for score, tet_dir in scored_dirs[:needed_passivations]:
        passivation_vectors.append(tet_dir * bond_length)
    
    return passivation_vectors


def calculate_passivation_positions(atoms, hanging_atom_indices, passivation_element='H', 
                                  cutoff=1.8, bond_length=1.09):
    """
    Рассчитать позиции атомов пассивации для всех атомов с висячими связями.
    Улучшенная версия с проверкой пересечений и ограничением количества пассиваций.
    
    Args:
        atoms: ASE Atoms object
        hanging_atom_indices: Список индексов атомов с висячими связями
        passivation_element: Элемент для пассивации (по умолчанию водород)
        cutoff: Расстояние для определения связей
        bond_length: Длина связи пассивации
        
    Returns:
        dict: Словарь {atom_index: [список_позиций_пассивации]}
    """
    passivation_positions = {}
    all_passivation_positions = []  # Для проверки пересечений
    
    # Определяем границы кристалла для всех осей
    positions = atoms.positions
    center = np.mean(positions, axis=0)
    
    for atom_idx in hanging_atom_indices:
        atom_pos = atoms.positions[atom_idx]
        
        # Получаем векторы существующих связей
        bond_vectors = get_bond_vectors(atoms, atom_idx, cutoff)
        
        # Вычисляем нормаль поверхности как направление от центра кристалла
        surface_normal = atom_pos - center
        if np.linalg.norm(surface_normal) > 1e-6:
            surface_normal = surface_normal / np.linalg.norm(surface_normal)
        else:
            surface_normal = np.array([0, 0, 1])  # По умолчанию
        
        # Определяем векторы пассивации с более строгими ограничениями
        passivation_vectors = get_tetrahedral_passivation_vectors(bond_vectors, bond_length, surface_normal)
        
        # Фильтруем пассивации, чтобы избежать пересечений
        valid_passivations = []
        for vec in passivation_vectors:
            new_pos = atom_pos + vec
            
            # Проверяем, что новая позиция не слишком близко к другим атомам
            too_close = False
            
            # Проверяем расстояние до существующих атомов (кроме родительского углерода)
            for other_idx, other_atom in enumerate(atoms):
                if other_idx == atom_idx:
                    continue
                dist = np.linalg.norm(new_pos - atoms.positions[other_idx])
                min_dist = 1.0 if other_atom.symbol == 'H' else 1.5  # H-H: 1.0Å, C-H: 1.5Å
                if dist < min_dist:
                    too_close = True
                    break
            
            # Проверяем расстояние до уже размещенных пассиваций
            for prev_pos in all_passivation_positions:
                if np.linalg.norm(new_pos - prev_pos) < 1.0:  # H-H минимум 1.0Å
                    too_close = True
                    break
            
            # Проверяем, что водород направлен наружу от кристалла
            to_center = np.linalg.norm(new_pos - center)
            parent_to_center = np.linalg.norm(atom_pos - center)
            if to_center <= parent_to_center:  # Водород должен быть дальше от центра
                too_close = True
            
            if not too_close:
                valid_passivations.append(vec)
                all_passivation_positions.append(new_pos)
        
        # Рассчитываем абсолютные позиции атомов пассивации
        positions_list = [atom_pos + vec for vec in valid_passivations]
        passivation_positions[atom_idx] = positions_list
        
        print(f"Атом {atom_idx}: {len(bond_vectors)} связей, {len(valid_passivations)}/{len(passivation_vectors)} пассиваций")
    
    return passivation_positions


def create_passivated_structure(atoms, hanging_atom_indices, passivation_element='H', 
                               cutoff=1.8, bond_length=1.09):
    """
    Создать новую структуру с добавленными атомами пассивации.
    
    Args:
        atoms: Исходная ASE Atoms структура
        hanging_atom_indices: Список индексов атомов с висячими связями
        passivation_element: Элемент для пассивации
        cutoff: Расстояние для определения связей
        bond_length: Длина связи пассивации
        
    Returns:
        ASE Atoms object: Новая структура с пассивацией
    """
    from ase import Atoms
    
    # Создаем копию исходной структуры
    passivated_atoms = atoms.copy()
    
    # Рассчитываем позиции пассивации
    passivation_positions = calculate_passivation_positions(
        atoms, hanging_atom_indices, passivation_element, cutoff, bond_length
    )
    
    # Добавляем атомы пассивации
    for atom_idx, positions in passivation_positions.items():
        for pos in positions:
            passivated_atoms.append(passivation_element)
            passivated_atoms.positions[-1] = pos
    
    print(f"Исходная структура: {len(atoms)} атомов")
    print(f"Пассивированная структура: {len(passivated_atoms)} атомов")
    print(f"Добавлено атомов пассивации: {len(passivated_atoms) - len(atoms)}")
    
    return passivated_atoms



def build_diamond_with_highlighting(size, highlight=False, output_file='diamond.xyz', passivate_element=None):
    """Build diamond structure with optional highlighting of unbonded atoms using tags"""
    try:
        size_int = int(size)
        if size_int <= 0:
            raise ValueError("Size must be positive")
    except ValueError:
        print(f"Error: Invalid size '{size}'. Please provide a positive integer.")
        return False
    
    # Create diamond structure
    diamond = diamond111('C', size=(size_int, size_int, size_int), vacuum=6.0)
    print(f"Created diamond structure with {len(diamond)} atoms")
    
    # Apply hydrogen passivation if requested
    if passivate_element:
        hanging_atoms = find_atoms_with_hanging_bonds(diamond)
        if hanging_atoms:
            print(f"Found {len(hanging_atoms)} atoms with hanging bonds - passivating with {passivate_element}")
            diamond = create_passivated_structure(diamond, hanging_atoms, passivate_element)
        else:
            print("No hanging bonds found - structure is fully coordinated")
    
    if highlight:
        # Find atoms with hanging bonds
        hanging_atoms = find_atoms_with_hanging_bonds(diamond)
        print(f"Found {len(hanging_atoms)} atoms with hanging bonds")
        
        # Create a copy for highlighting
        highlighted_diamond = diamond.copy()
        
        # Initialize all atoms with tag 0 (fully coordinated)
        highlighted_diamond.set_tags([0] * len(highlighted_diamond))
        
        if hanging_atoms:
            # Mark hanging bond atoms with different tags based on coordination
            cutoffs = [0.9] * len(highlighted_diamond)
            nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
            nl.update(highlighted_diamond)
            
            tags = highlighted_diamond.get_tags()
            for idx in hanging_atoms:
                indices, offsets = nl.get_neighbors(idx)
                num_neighbors = len(indices)
                
                # Tag based on coordination deficiency
                if num_neighbors == 0:
                    tags[idx] = 5  # No bonds
                elif num_neighbors == 1:
                    tags[idx] = 4  # 1 bond (3 missing)
                elif num_neighbors == 2:
                    tags[idx] = 3  # 2 bonds (2 missing)
                elif num_neighbors == 3:
                    tags[idx] = 2  # 3 bonds (1 missing)
                else:
                    tags[idx] = 1  # Other undercoordinated cases
            
            highlighted_diamond.set_tags(tags)
            print(f"Coordination statistics:")
            for coord in range(6):
                count = sum(1 for t in tags if t == coord)
                if count > 0:
                    if coord == 0:
                        print(f"  Tag {coord}: {count} atoms (fully coordinated)")
                    else:
                        print(f"  Tag {coord}: {count} atoms (undercoordinated)")
        
        # Save the highlighted structure
        from ase.io import write
        write(output_file, highlighted_diamond, format='extxyz')
        print(f"Diamond structure saved to '{output_file}' with tags for coordination highlighting")
        print("Tags: 0=fully coordinated, 1-5=undercoordinated (higher tag = fewer bonds)")
        
        # Also create a separate file with original structure
        original_output = output_file.replace('.xyz', '_original.xyz')
        write(original_output, diamond)
        print(f"Original structure saved to '{original_output}'")
    else:
        # Save without highlighting
        from ase.io import write
        write(output_file, diamond)
        print(f"Diamond structure saved to '{output_file}'")
    
    return True


def main():
    diamond = diamond111('C', size=(5, 5, 5), vacuum=6.0)

    atoms_with_hanging_bonds = find_atoms_with_hanging_bonds(diamond)
    print(len(atoms_with_hanging_bonds))

    # Рассчитываем позиции пассивации
    passivation_positions = calculate_passivation_positions(diamond, atoms_with_hanging_bonds)

    # Выводим информацию о пассивации
    total_passivation_atoms = sum(len(positions) for positions in passivation_positions.values())
    print(f"\nОбщее количество атомов пассивации: {total_passivation_atoms}")

    # Пример для первых нескольких атомов
    print("\nПримеры позиций пассивации:")
    for i, (atom_idx, positions) in enumerate(list(passivation_positions.items())[:3]):
        atom_pos = diamond.positions[atom_idx]
        print(f"\nАтом {atom_idx} в позиции {atom_pos}:")
        for j, pos in enumerate(positions):
            print(f"  Пассивация {j}: {pos}")
            distance = np.linalg.norm(pos - atom_pos)
            print(f"    Расстояние: {distance:.3f} Å")

    # Создаем полную пассивированную структуру
    print("\n" + "="*50)
    print("СОЗДАНИЕ ПАССИВИРОВАННОЙ СТРУКТУРЫ")
    print("="*50)

    passivated_diamond = create_passivated_structure(diamond, atoms_with_hanging_bonds)

    # Сохраняем структуру в файл
    from ase.io import write
    output_file = "diamond_passivated.xyz"
    write(output_file, passivated_diamond)
    print(f"\nПассивированная структура сохранена в файл: {output_file}")

    # Проверяем, что в пассивированной структуре нет висячих связей
    print("\nПроверка пассивированной структуры:")
    hanging_after = find_atoms_with_hanging_bonds(passivated_diamond)
    print(f"Атомов с висячими связями после пассивации: {len(hanging_after)}")

    if len(hanging_after) == 0:
        print("✅ Пассивация успешна! Все висячие связи закрыты.")
    else:
        print("⚠️  Остались атомы с висячими связями:")
        for idx in hanging_after:
            atom = passivated_diamond[idx]
            print(f"  Атом {idx}: {atom.symbol} в позиции {atom.position}")




if __name__ == "__main__":
    main()

