# Cloud DFT Calculation Infrastructure

Автоматизированная инфраструктура для DFT расчетов в Google Cloud с автоматическим скачиванием результатов и остановкой серверов.

## Структура

```
lib/cloud/
├── terraform/           # Terraform конфигурация
│   ├── main.tf         # Основная инфраструктура
│   ├── variables.tf    # Переменные
│   └── outputs.tf      # Выходные значения
├── containers/         # Docker контейнеры
│   ├── dftb/          # DFTB+ контейнер
│   │   ├── Dockerfile
│   │   └── calculate.py
│   └── orca/          # Orca контейнер
│       ├── Dockerfile
│       └── calculate.py
├── scripts/           # Скрипты оркестрации
│   ├── orchestrator.py # Cloud Function
│   └── requirements.txt
├── cloud_runner.py    # Главный скрипт запуска
├── cloudbuild.yaml    # Конфигурация сборки
├── setup.sh          # Скрипт установки
└── README.md         # Эта документация
```

## Быстрый старт

### 1. Установка

```bash
# Запустите скрипт установки
./setup.sh your-project-id us-central1
```

### 2. Подготовка контейнеров

```bash
# Скопируйте ваши бинарники
cp /path/to/dftb+ containers/dftb/
cp -r /path/to/slakos containers/dftb/slakos/
cp /path/to/orca containers/orca/
cp -r /path/to/basis_sets containers/orca/basis_sets/

# Пересоберите контейнеры
gcloud builds submit --config=cloudbuild.yaml .
```

### 3. Запуск расчета

```bash
python3 cloud_runner.py \
  --project your-project-id \
  --bucket your-project-id-dft-calculations \
  --function-url https://your-function-url \
  --input molecule.xyz \
  --calculator dftb \
  --method relax \
  --cpu-count 128 \
  --cleanup
```

## Поддерживаемые методы

### DFTB+
- `relax`: Геометрическая оптимизация
- `hessian`: Расчет матрицы Гессе и частот
- `single`: Одноточечный расчет

### Orca
- `opt-B3LYP-def2-SVP`: Оптимизация с B3LYP/def2-SVP
- `hessian-B3LYP-def2-SVP`: Частоты с B3LYP/def2-SVP  
- `single-M06-2X-def2-TZVP`: Одноточечный с M06-2X/def2-TZVP

## Автоматическое управление ресурсами

Скрипт автоматически:
1. 📤 Загружает входной файл в Cloud Storage
2. 🚀 Запускает Batch Job с preemptible инстансами
3. ⏳ Ждет завершения расчета
4. 📥 Скачивает все результаты
5. 🗑️ Удаляет временные файлы и останавливает сервер (с флагом `--cleanup`)

## Стоимость

Примерные цены для 100 атомов:
- **DFTB+ релаксация**: $1-3 за расчет
- **DFT оптимизация**: $5-15 за расчет  
- **DFT гессиан**: $15-50 за расчет

*С preemptible скидками до 80%*

## Мониторинг

```bash
# Посмотреть активные задачи
gcloud batch jobs list --location=us-central1

# Посмотреть логи задачи
gcloud batch jobs describe JOB_ID --location=us-central1

# Посмотреть файлы в bucket
gsutil ls gs://your-bucket-name/
```

## Очистка ресурсов

```bash
# Удалить всю инфраструктуру
cd terraform
terraform destroy -auto-approve
```

## Устранение неполадок

### Job Failed
```bash
# Посмотрите логи задачи
gcloud logging read "resource.type=batch_task" --limit=50
```

### Container Build Failed  
```bash
# Проверьте права доступа к бинарникам
ls -la containers/*/
```

### Download Failed
```bash
# Проверьте содержимое bucket
gsutil ls -la gs://your-bucket-name/outputs/
```