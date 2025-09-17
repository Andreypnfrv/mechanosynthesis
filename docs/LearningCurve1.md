### Этап 1. Основы квантовой химии на Mac M1

Цель: освоить установку и запуск DFT-(ORCA) и DFTB+-расчётов на простых системах (H₂, Si₂), сравнить скорости и результаты.

⸻

1. Установка окружения
	1.	Miniforge (ARM64)

# в терминале
curl -LO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh  # принять лицензию, установить в ~/miniforge3
source ~/miniforge3/bin/activate


	2.	Создать conda-среду

conda create -n msynth python=3.10 -y
conda activate msynth


	3.	Установить Python-пакеты

conda install -c conda-forge ase orca dftbplus matplotlib -y
pip install py3Dmol


	4.	ORCA (если не из conda)
	•	Скачать Mac ARM-бинарник с orca.eu и распаковать в ~/Programs/orca/
	•	В ~/.zshrc добавить:

export PATH=$HOME/Programs/orca:$PATH
export OSCAR_KEY="ваш_ключ_ORCA"


	•	source ~/.zshrc → orca --version

⸻

2. DFT-расчёт H₂
	1.	Создать файл h2.inp

! PBE0 def2-SVP Opt TightSCF
%pal nprocs 4 end
* xyz 0 1
H    0.000000   0.000000   0.000000
H    0.000000   0.000000   0.740000
*


	2.	Запустить и измерить время

time orca h2.inp > h2.out
grep "FINAL SINGLE POINT ENERGY" h2.out


	3.	Сохранить результат
	•	Энергия: строка FINAL SINGLE POINT ENERGY
	•	Время: из time

⸻

3. DFTB+-расчёт H₂
	1.	Пример входника dftb_h2.in

Geometry = GenFormat {
  <<< "h2.xyz"
}
Hamiltonian = DFTB {
  SCC = Yes
  SlaterKosterFiles = Type2FileNames {
    Prefix = "/opt/miniforge3/share/dftbplus/slakos/"
    Separator = "-"
    Suffix = ".skf"
  }
}
Driver = ConjugateGradient {
  MaxSteps = 1000
}


	2.	Подготовить h2.xyz (в простом XYZ-формате).
	3.	Запустить

time dftb+ > dftb_h2.log
grep "Total SCC" dftb_h2.log


	4.	Записать энергию и время.

⸻

4. Повторить на Si₂
	•	Сменить систему: si2.inp и si2.xyz
	•	Повторить шаги 2 и 3: DFT (ORCA) и DFTB+.

⸻

5. Сравнение и отчёт
	•	Собрать в файл results.md:

| Система | Метод | Энергия (Eh) | Время (s) |
|---------|-------|--------------|-----------|
| H2      | DFT   | …            | …         |
| H2      | DFTB+ | …            | …         |
| Si2     | DFT   | …            | …         |
| Si2     | DFTB+ | …            | …         |


	•	Построить простой график времени vs системы:

import matplotlib.pyplot as plt
# заполнить списки times = [..], labels=[..]
plt.plot(labels, times, marker='o')
plt.ylabel("Time (s)")
plt.show()



⸻

Как использовать меня дальше
Соберите результаты (h2.out, dftb_h2.log, si2.out, dftb_si2.log) в одну папку и пришлите её путь — я автоматически сгенерирую вам Snakemake-пайплайн и Jupyter-отчёт, который за вас выполнит все расчёты и визуализацию за одну команду.