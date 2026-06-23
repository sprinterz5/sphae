# Установка sphae (GPU) на Windows 11 — пошаговый гайд

> **Для кого:** Windows 11, NVIDIA GPU, хочешь запускать sphae с GPU через WSL2.
> Если GPU нет — всё то же самое, только на шаге 9 добавляешь `--no-use-gpu`.

---

## Шаг 1 — Включить WSL2

В PowerShell (от администратора):

```powershell
wsl --install
```

Перезагрузить компьютер. После перезагрузки откроется терминал Ubuntu — придумать имя пользователя и пароль.

Если Ubuntu уже установлена, просто открыть её из меню Пуск.

---

## Шаг 2 — Настроить DNS (если интернет не работает в WSL)

В терминале Ubuntu:

```bash
sudo bash -c 'echo -e "nameserver 8.8.8.8\nnameserver 1.1.1.1" > /etc/resolv.conf'
sudo bash -c 'echo -e "[network]\ngenerateResolvConf = false" >> /etc/wsl.conf'
```

Проверить: `ping -c 2 google.com` — должен получить ответ.

---

## Шаг 3 — Установить Miniforge (conda)

```bash
curl -LO https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
~/miniforge3/bin/conda init bash
exec bash
```

Проверить: `conda --version` — должна показать версию.

Установить mamba (быстрый solver):

```bash
conda install -n base -c conda-forge mamba -y
```

---

## Шаг 4 — Настроить conda (один раз)

```bash
conda config --set channel_priority flexible
```

> **Почему flexible, не strict?** Strict ломает pharokka — bioconda не может разрешить biopython и phanotate в строгом режиме приоритетов.

---

## Шаг 5 — Зарегистрировать CUDA библиотеки WSL2 (один раз)

NVIDIA-драйвер кладёт `libcuda.so.1` в `/usr/lib/wsl/lib/`, но динамический линкер его
там не ищет по умолчанию. Без этого шага PyTorch не видит GPU:

```bash
echo /usr/lib/wsl/lib | sudo tee /etc/ld.so.conf.d/wsl.conf
sudo ldconfig
```

Проверить:

```bash
ldconfig -p | grep libcuda
# ожидается: libcuda.so.1 => /usr/lib/wsl/lib/libcuda.so.1
```

---

## Шаг 6 — Установить sphae с GPU-форка

```bash
pip install git+https://github.com/sprinterz5/sphae.git@gpu-support
```

> **Важно:** устанавливать именно с этого форка, не с PyPI (`pip install sphae`).
> PyPI-версия использует старые env-файлы без CUDA и с битыми зависимостями.

Проверить: `sphae --version` — должна показать `1.5.5+gpu`.

---

## Шаг 7 — Создать рабочую директорию

```bash
mkdir -p ~/sphae-work && cd ~/sphae-work
```

> **Почему не `/mnt/c/` или Desktop?** Snakemake пишет служебные файлы в текущую папку.
> NTFS-диски (всё что под `/mnt/`) не поддерживают Linux-права доступа → PermissionError.
> Базы данных при этом **можно** хранить на D: (передаются через `--db_dir`).

---

## Шаг 8 — Скачать базы данных

Определить куда сохранять. Если есть второй диск D:

```bash
DB_DIR=/mnt/d/sphae-databases
```

Если только C: (но лучше D: — базы занимают ~15 ГБ):

```bash
DB_DIR=~/sphae-databases
```

Запустить установку:

```bash
cd ~/sphae-work
sphae install --db_dir $DB_DIR --threads 8 --conda-frontend mamba
```

Это займёт **30–90 минут** (скачивает ~15 ГБ). Параллельно скачивает 5–6 баз.

### Если упадёт `download_medaka_models`

Это нормально на WSL2. medaka нужен только для Nanopore-ридов, а его старый PyTorch
несовместим с ядром WSL2. Исправить нельзя, но и не нужно:

```bash
mkdir -p $DB_DIR/medaka_models
touch $DB_DIR/medaka_models/medaka.flag
```

После этого перезапустить `sphae install ...` — он продолжит с того места где остановился
(уже скачанное не трогает), и пропустит medaka.

---

## Шаг 9 — Запустить аннотацию генома фага

Подготовить папку с геномом (важно: `--genome` принимает **папку**, не файл):

```bash
mkdir -p ~/sphae-work/genomes
# скопировать .fasta/.fa/.fna файл(ы) в эту папку
cp /путь/к/фагу.fasta ~/sphae-work/genomes/
```

Запустить:

```bash
cd ~/sphae-work

sphae annotate \
    --genome genomes \
    --output results \
    --db_dir $DB_DIR \
    --threads 8 \
    --use-gpu
```

Без GPU:

```bash
sphae annotate ... --no-use-gpu
```

---

## Что куда скачивается

| База | Размер | Нужна для |
|------|--------|-----------|
| Pfam35.0 | ~280 МБ | аннотация доменов |
| pharokka_db | ~280 МБ | основная аннотация фага |
| checkv-db-v1.5 | ~656 МБ | оценка полноты генома |
| phold (ProstT5 + Foldseek) | ~5.7 ГБ | структурная аннотация белков (GPU) |
| phynteny models | ~100 МБ | предсказание функций генов |
| medaka models | пропускаем | только Nanopore, несовм. с WSL2 |

---

## Частые ошибки

### `PermissionError: Operation not permitted: '.../.snakemake/...'`
Запущено из `/mnt/c/` или Desktop. Перейти в `~/sphae-work`:
```bash
cd ~/sphae-work && sphae install ...
```

### `Genome samples: []` (ничего не аннотируется)
Передан файл вместо папки. `--genome` должен быть папкой с `.fasta` файлами:
```bash
mkdir genomes && mv my_phage.fasta genomes/
sphae annotate --genome genomes ...
```

### `No available GPU was found` / `Using device: cpu`
PyTorch не видит GPU. Проверить:
```bash
ldconfig -p | grep libcuda   # должно найти /usr/lib/wsl/lib/libcuda.so.1
```
Если пусто — выполнить Шаг 5. Если уже выполнен — проверить `nvidia-smi` из WSL2:
```bash
nvidia-smi   # должен показать GPU
```

### `Foldseek not found. Please reinstall phold.`
Старый conda env без foldseek. Удалить env и перезапустить:
```bash
# найти путь к env
ls ~/miniforge3/lib/python*/site-packages/sphae/workflow/conda/
# удалить env с phold
rm -rf ~/miniforge3/lib/python*/site-packages/sphae/workflow/conda/HASH_*
sphae annotate ...   # snakemake пересоберёт env
```

### `ModuleNotFoundError: No module named 'requests'`
Установлена PyPI-версия sphae. Переустановить с форка:
```bash
pip install --force-reinstall git+https://github.com/sprinterz5/sphae.git@gpu-support
sphae install ...
```

### `ValueError: ...torch.load...CVE-2025-32434`
Та же причина — старый phold env из PyPI. То же решение: переустановить с форка и
удалить старый conda env:
```bash
pip install --force-reinstall git+https://github.com/sprinterz5/sphae.git@gpu-support
rm -rf ~/miniforge3/lib/python*/site-packages/sphae/workflow/conda/HASH_*
sphae annotate ...
```

### `LibMambaUnsatisfiableError: ...biopython excluded by strict repo priority`
```bash
conda config --set channel_priority flexible
```

### `ImportError: libtorch_cpu.so: cannot enable executable stack`
Это medaka на WSL2. Создать fake-флаг (см. Шаг 8).

---

## Системные требования

| | Минимум | Рекомендуется |
|--|---------|---------------|
| ОС | Windows 10 21H2 | Windows 11 |
| RAM | 16 ГБ | 32 ГБ |
| Диск | 20 ГБ свободно | 50 ГБ на отдельном диске |
| GPU | любая NVIDIA (driver ≥ 525) | RTX серия, 8+ ГБ VRAM |
| CUDA | 12.x через driver | — |

GPU необязателен — без него phold работает на CPU, просто медленнее (~10–20× дольше).
