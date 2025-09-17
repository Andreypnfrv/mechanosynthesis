Ок. Пересобрал урок под стек ASE + OVITO (OpenMM — как pip-пакет). LAMMPS/ReaxFF убрал. Реактивность закрываем DFTB+ через ASE.

⸻

Урок 2. Молекулярная механика и подготовка структур (ASE + OVITO)

Цель: настроить быстрые MD через OpenMM (не реактивные) и DFTB+ (полуреактивные), научиться выгружать траектории в OVITO и проверять стабильность алмазоподобных кластеров.

1) Теория — коротко
	•	MM/OpenMM: быстрый расчёт сил по параметризованному FF (без электронов). Хорошо для «физики решётки» и больших систем, без разрывов связей.
	•	DFTB+: полуквантовый метод (электронная структура → можно делать эмпирическую реактивность лучше MM). Медленнее MM, но на кластерах/малых ячейках ок.
	•	OVITO: визуализация/анализ траекторий (.xyz, .traj, lammps-dump и т.д.).

⸻

2) Практика

2.1 OpenMM через ASE: графен ~100 атомов, 300 K (20 ps)

Задача — быстро получить устойчивую 2D-структуру и валидную траекторию для OVITO.

pip install ase openmm

# run_graphene_openmm.py
from ase.build import graphene_nanoribbon
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase.io import write
from ase import units
from ase.calculators.openmm import OpenMM
from ase.md import MDLogger

# ~10x10 «ячейка» даёт ~100–200 атомов — ок для демо
atoms = graphene_nanoribbon(10, 10, type='armchair', vacuum=6.0, saturated=False)
atoms.center(vacuum=6.0)

calc = OpenMM(atoms, forcefield='amber99sb.xml')  # простой FF для демо
atoms.set_calculator(calc)

MaxwellBoltzmannDistribution(atoms, temperature_K=300)
dyn = Langevin(atoms, 1.0 * units.fs, temperature_K=300, friction=0.01)

traj = Trajectory('graphene_md.traj', 'w', atoms)
dyn.attach(traj.write, interval=50)
dyn.attach(MDLogger(atoms, dyn, 'graphene_md.log', stress=False, peratom=False), interval=50)

dyn.run(20000)  # 20 ps

# Экспорт в OVITO-совместимый XYZ
from ase.io.trajectory import Trajectory as Tread
imgs = Tread('graphene_md.traj')
write('graphene_md.xyz', imgs)
print('Wrote graphene_md.xyz')

Открыть в OVITO: File → Load → graphene_md.xyz (или двойной клик по файлу).

⸻

2.2 DFTB+ через ASE: алмазоподобный кластер C≈30, 300 K (демо 10–20 ps)

Реактивность/устойчивость кластера — через полуквантовую MD. Для длинных траекторий используем домашку.

Требования: установлен DFTB+ и sk-набор (например, mio-1-1 или 3ob). Укажи путь в SKDIR.

# run_c30_dftb_md.py
import os, numpy as np
from ase.build import bulk
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.io.trajectory import Trajectory
from ase.io import write
from ase.calculators.dftb import Dftb
from ase.md import MDLogger

# Сферический кластер из алмаза
a = 3.567
dia = bulk('C', 'diamond', a=a) * (4, 4, 4)
pos = dia.positions.copy()
center = dia.get_center_of_mass()
r = 5.6  # подбери радиус ~С30 (± пару атомов — не критично)
mask = np.linalg.norm(pos - center, axis=1) < r
cluster = dia[mask]
cluster.center(vacuum=8.0)

SKDIR = os.environ.get('SKDIR', '/path/to/skf/mio-1-1')  # ВАЖНО: поправь путь

calc = Dftb(
    atoms=cluster,
    label='c30_md',
    Hamiltonian_SCC='Yes',
    Hamiltonian_SCCTolerance=1e-6,
    Hamiltonian_MaxSCCIterations=150,
    Hamiltonian_SlaterKosterFiles_Prefix=SKDIR + '/',
    Hamiltonian_SlaterKosterFiles_Separator='-',
    Hamiltonian_SlaterKosterFiles_Suffix='.skf',
    Hamiltonian_MaxAngularMomentum_='',
    Hamiltonian_MaxAngularMomentum_C='p'
)
cluster.set_calculator(calc)

MaxwellBoltzmannDistribution(cluster, temperature_K=300)
dyn = Langevin(cluster, 0.5 * units.fs, temperature_K=300, friction=0.01)

traj = Trajectory('c30_md.traj', 'w', cluster)
dyn.attach(traj.write, interval=50)
dyn.attach(MDLogger(cluster, dyn, 'c30_md.log', stress=False, peratom=False), interval=50)

dyn.run(20000)  # 10 ps демо (для 200 ps — см. домашку)

from ase.io.trajectory import Trajectory as Tread
imgs = Tread('c30_md.traj')
write('c30_md.xyz', imgs)
print('Wrote c30_md.xyz')

В OVITO: File → Load → c30_md.xyz → добавить модуль Compute bond properties/Cluster analysis при необходимости.

⸻

3) Домашнее задание
	1.	OpenMM/ASE: графен ~100–200 атомов, 300 K, 20–50 ps. Сохранить graphene_md.xyz.
	2.	DFTB+/ASE: алмазоподобный кластер C≈30, 300 K, 200 ps (0.5 fs → 400 000 шагов).
	•	Критерии стабильности: сохранение sp³-координации (соседей ~4), отсутствие длительных sp²-фрагментов; RMSD < ~1.0 Å относительно стартовой (грубая эвристика).
	3.	Выгрузить обе траектории в OVITO. Сделать 3 скрина: старт/середина/финал для каждого кейса.

⸻

4) Полезные заметки
	•	Для DFTB+ не экономь на SCCTolerance и MaxSCCIterations.
	•	Если кластер «раздувает», подними трение friction до 0.05 и/или уменьши шаг до 0.25 fs.
	•	Для графена термостат Ланжевена + 1 fs — норм.
	•	OVITO любит EXTXYZ: write('file.xyz', imgs, format='extxyz') — метаданные полезны.

⸻

Хочешь, сделаю тебе готовый мини-репозиторий: make run-graphene, make run-c30, автоконверсия в .xyz, и ovito запускается на нужный файл одной командой. Это сэкономит время и уберёт ручные шаги.