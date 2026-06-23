# sphae-gpu — Phold GPU support (fork changes)

This fork enables **GPU execution for `phold` and `phynteny_transformer`** in sphae.
Upstream sphae always calls `phold` with `--cpu` and installs CPU-only PyTorch.
We fix both issues:

1. The `--cpu` flag is no longer hard-coded; it is injected only when `use_gpu: False`.
2. Both `phold.yaml` and `phynteny.yaml` install CUDA-enabled PyTorch via pip from
   the official PyTorch wheel index (`download.pytorch.org/whl/cu124`), bypassing
   conda's solver which consistently picks the CPU-only build from conda-forge.

Tested with NVIDIA RTX 4060 Laptop (8 GB VRAM), driver 610.47, CUDA UMD 13.3,
WSL2 kernel 6.6.x (Microsoft).

---

## Behaviour

- **Default is GPU.** `use_gpu: True` is the default in this fork. Pass `--no-use-gpu`
  (or set `use_gpu: False` in config.yaml) to fall back to CPU.
- CLI override per run:
  - `sphae run   --input <dir> --use-gpu ...`
  - `sphae run   --input <dir> --no-use-gpu ...`
  - `sphae annotate --genome <dir> --use-gpu ...`

---

## What changed

### Conda environments

| File | Change |
|------|--------|
| `sphae/workflow/envs/phold.yaml` | Removed `pytorch`/`nvidia` channels. `torch==2.6.0+cu124` installed via `pip --extra-index-url https://download.pytorch.org/whl/cu124` (guarantees CUDA build). Added `foldseek` as conda dep (pip install of phold doesn't include the binary). `phold>=1.2.5` via pip. |
| `sphae/workflow/envs/phynteny.yaml` | Same pip CUDA torch treatment; `requests` added as conda dep (required by `install_models`); `phynteny_transformer>=0.1.3` via pip. |

### Pipeline rules & CLI

| File | Change |
|------|--------|
| `sphae/config/config.yaml` | Added `args.use_gpu: True` (default GPU on). |
| `sphae/workflow/rules/1.preflight.smk` | Globals `USE_GPU`, `PHOLD_CPU_FLAG` (`""` on GPU, `--cpu` on CPU), `PHOLD_GPU_RESOURCE` (`1`/`0`). |
| `sphae/workflow/rules/1.preflight-annot.smk` | Same globals for the `annotate` workflow. |
| `sphae/workflow/rules/8.phold.smk` | `phold_run_paired`, `phold_run_longreads`: `--cpu` → `{params.cpu}`, `params.cpu`, `resources.gpu`. |
| `sphae/workflow/rules/789.genome.annot.smk` | `phold_run_genome`: same treatment. |
| `sphae/workflow/rules/789.prot.annot.smk` | `phold_run_protein`: same treatment. |
| `sphae/__main__.py` | `--use-gpu/--no-use-gpu` added to `run` and `annotate`; `resolve_db_paths()` maps `--db_dir` to individual db paths; `--conda-prefix` added to `install` command. |
| `pyproject.toml`, `sphae.VERSION` | Version → `1.5.5+gpu`; added `package-data` so test genomes are included. |

---

## How the GPU flag flows

```
sphae annotate --genome genomes/ --use-gpu
  → merge_config['args']['use_gpu'] = True
    → 1.preflight-annot.smk: PHOLD_CPU_FLAG = ""
      → 789.genome.annot.smk: phold predict ... {params.cpu} ...
        → becomes: phold predict ...  ...   (empty string → GPU mode)
```

When `--no-use-gpu`:
```
PHOLD_CPU_FLAG = "--cpu"  →  phold predict ... --cpu ...   (CPU mode)
```

---

## Why conda channels fail to install CUDA torch (and how we fix it)

The `pytorch` + `nvidia` conda channel approach (`pytorch>=2.6.0` + `pytorch-cuda=12.4`)
**does not reliably work** with `channel_priority flexible`:

- conda-forge also has `pytorch` packages (CPU-only builds).
- With flexible priority, conda's solver can choose the conda-forge CPU build even
  when `pytorch` channel is listed first — especially if a newer version exists on
  conda-forge (`2.7.1` appeared while `pytorch-cuda=12.4` only covers ≤ `2.6.x`).
- Pinning to `pytorch=2.6.0` in conda still doesn't guarantee the CUDA build because
  conda-forge also has `pytorch=2.6.0`.

**Fix:** install torch entirely via pip from the official PyTorch CUDA wheel index.
This is also the approach recommended by pytorch.org for non-conda environments:

```yaml
- pip:
    - --extra-index-url https://download.pytorch.org/whl/cu124
    - torch==2.6.0+cu124
```

The `+cu124` local version tag is unambiguous — pip fetches exactly the CUDA 12.4 wheel
and does not fall back to a CPU build.

---

## Why foldseek must be a conda dep

`phold compare` requires the `foldseek` binary on PATH. When phold is installed via pip
(not bioconda), foldseek is **not** installed as a side effect. It must be listed
explicitly as a conda dependency so bioconda provides the binary:

```yaml
dependencies:
    - foldseek    # ← required by phold compare; pip install of phold omits this
    - pip:
        - phold>=1.2.5
```

---

## Sanity checks

```bash
# all phold rules use the placeholder, not hard-coded --cpu
grep -rn "params.cpu\|PHOLD_CPU_FLAG" sphae/workflow/

# CLI exposes the flag
sphae annotate --help | grep -- --use-gpu

# verify GPU is visible inside the phold env (find env path first)
PHOLD_ENV=$(ls /home/$USER/miniforge3/lib/python*/site-packages/sphae/workflow/conda/ \
    | grep -v "^$" | head -1)
conda run -p "$PHOLD_ENV" \
    python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
# expected: 2.6.0 12.4 True
```

---

## CUDA version note

`torch==2.6.0+cu124` bundles its own CUDA 12.4 runtime. It requires:
- NVIDIA driver ≥ 525.60 (CUDA 12.x backward compat)
- `libcuda.so.1` visible to the dynamic linker

On WSL2, `libcuda.so.1` lives in `/usr/lib/wsl/lib/` which is **not** in ldconfig by
default. Register it once:

```bash
echo /usr/lib/wsl/lib | sudo tee /etc/ld.so.conf.d/wsl.conf
sudo ldconfig
```

---

## WSL2 known issues and workarounds

These are one-time setup steps on Windows Subsystem for Linux 2. Once done they stay fixed.

### 1. Always run from `~/sphae-work`, never from `/mnt/c/` or `/mnt/d/`

Snakemake writes `.snakemake/` to the current directory. NTFS mounts (`/mnt/c/`, `/mnt/d/`)
don't support Linux `chmod` — Snakemake crashes with `PermissionError: Operation not permitted`.

```bash
mkdir -p ~/sphae-work && cd ~/sphae-work
sphae install --db_dir /mnt/d/sphae-databases ...
```

Databases can live on `/mnt/d/` (just passed as `--db_dir`). The *working directory* must
be on the ext4 filesystem (`/home/...`).

### 2. Install from this fork, not PyPI

The upstream PyPI package has the old env yamls (CPU-only torch, missing foldseek dep,
missing requests). Always install from the `gpu-support` branch:

```bash
pip install --force-reinstall git+https://github.com/sprinterz5/sphae.git@gpu-support
```

After reinstalling, snakemake detects changed yaml content by hash and rebuilds affected
conda envs automatically on the next run.

### 3. conda channel_priority must be `flexible`

`strict` breaks pharokka installation (biopython excluded, phanotate unsatisfiable).

```bash
conda config --set channel_priority flexible
```

This is persistent — set once, stays across sessions.

### 4. Register CUDA libs for WSL2 (one-time)

`libcuda.so.1` ships with the NVIDIA Windows driver at `/usr/lib/wsl/lib/` but is not
registered in ldconfig by default. Without this, `torch.cuda.is_available()` returns False
even with a CUDA-enabled torch build:

```bash
echo /usr/lib/wsl/lib | sudo tee /etc/ld.so.conf.d/wsl.conf
sudo ldconfig
```

### 5. medaka always fails on WSL2 (Nanopore only — safe to skip)

medaka uses an old PyTorch that requires an executable stack (`PT_GNU_STACK`).
WSL2's kernel (6.6.x Microsoft) disables this security exception → `libtorch_cpu.so: cannot
enable executable stack`. There is no fix short of a custom kernel build.

**medaka is only used for Nanopore long-read polishing.** For Illumina short reads it is
never invoked. Work around it by creating the expected flag file manually:

```bash
mkdir -p /mnt/d/sphae-databases/medaka_models
touch /mnt/d/sphae-databases/medaka_models/medaka.flag
```

Snakemake sees the output file → skips `download_medaka_models` → install completes.

### 6. phold CVE-2025-32434 (torch.load safety check)

`transformers` ≥ 4.51 blocks `torch.load` unless PyTorch ≥ 2.6. Old phold envs (torch 2.1)
crash with `ValueError: Due to a serious vulnerability issue in torch.load...`.

Fixed in this fork: torch is pinned to `2.6.0+cu124` via pip. If you see this error in an
existing env, delete the conda env directory and rerun `sphae annotate` to rebuild it.

### 7. phynteny `ModuleNotFoundError: No module named 'requests'`

The `install_models` script imports `requests` at module level, but upstream `phynteny.yaml`
didn't list it as a dependency. Fixed in this fork: `requests` is now a conda dep.

### 8. `--genome` must be a directory, not a file

`sphae annotate` globs `*.fasta` / `*.fa` / `*.fna` **inside** the path you pass to
`--genome`. Passing a single file path results in `Genome samples: []`.

```bash
# Wrong
sphae annotate --genome my_phage.fasta ...

# Correct
mkdir -p genomes && mv my_phage.fasta genomes/
sphae annotate --genome genomes ...
```

---

## Dependency version rationale

| Package | Why this version |
|---------|-----------------|
| `torch==2.6.0+cu124` | Exact CUDA 12.4 pip wheel; satisfies CVE-2025-32434 (≥2.6); `+cu124` prevents pip falling back to CPU build |
| `foldseek` (conda) | Required binary for `phold compare`; not included when phold installed via pip |
| `phold>=1.2.5` | Older phold had `download_requests()` signature bug (too few positional args) |
| `phynteny_transformer>=0.1.3` | Minimum version supporting PyTorch 2.x |
