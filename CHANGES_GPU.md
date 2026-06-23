# sphae-gpu — Phold GPU support (fork changes)

This fork enables **GPU execution for `phold` and `phynteny_transformer`** in sphae.
Upstream sphae always calls `phold` with `--cpu` and installs CPU-only PyTorch.
We fix both issues:

1. The `--cpu` flag is no longer hard-coded; it is injected only when `use_gpu: False`.
2. Both `phold.yaml` and `phynteny.yaml` now install CUDA-enabled PyTorch via the
   official `pytorch` + `nvidia` conda channels before the bioinformatics tools are
   added, preventing bioconda from overriding the GPU build with a CPU one.

Tested with NVIDIA RTX 4060 Laptop (8 GB VRAM), driver 610.47, CUDA UMD 13.3.

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
| `sphae/workflow/envs/phold.yaml` | Added `pytorch` + `nvidia` channels; `pytorch>=2.6.0` + `pytorch-cuda=12.6` as conda deps (≥2.6 required by CVE-2025-32434); moved `phold>=1.2.5` to `pip:` to prevent bioconda overriding CUDA torch. |
| `sphae/workflow/envs/phynteny.yaml` | Same CUDA torch treatment; added `requests` to conda deps (required by `install_models` script); `phynteny_transformer>=0.1.3` via pip. |

### Pipeline rules & CLI

| File | Change |
|------|--------|
| `sphae/config/config.yaml` | Added `args.use_gpu: True` (default GPU on). |
| `sphae/workflow/rules/1.preflight.smk` | Globals `USE_GPU`, `PHOLD_CPU_FLAG` (`""` on GPU, `--cpu` on CPU), `PHOLD_GPU_RESOURCE` (`1`/`0`). |
| `sphae/workflow/rules/1.preflight-annot.smk` | Same globals for the `annotate` workflow. |
| `sphae/workflow/rules/8.phold.smk` | `phold_run_paired`, `phold_run_longreads`: `--cpu` → `{params.cpu}`, `params.cpu`, `resources.gpu`. |
| `sphae/workflow/rules/789.genome.annot.smk` | `phold_run_genome`: same treatment. |
| `sphae/workflow/rules/789.prot.annot.smk` | `phold_run_protein`: same treatment. |
| `sphae/__main__.py` | `--use-gpu/--no-use-gpu` added to `run` and `annotate`; passed into `config['args']['use_gpu']`. |
| `pyproject.toml`, `sphae.VERSION` | Version → `1.5.5+gpu`. |

---

## How the GPU flag flows

```
sphae run --use-gpu
  → merge_config['args']['use_gpu'] = True
    → 1.preflight.smk: PHOLD_CPU_FLAG = ""
      → 8.phold.smk: phold predict ... {params.cpu} ...
        → becomes: phold predict ...  ...   (empty string → GPU mode)
```

When `--no-use-gpu`:
```
PHOLD_CPU_FLAG = "--cpu"  →  phold predict ... --cpu ...   (CPU mode)
```

---

## Why bioconda would break GPU torch (and how we fix it)

When `phold` is listed as a conda dep, bioconda resolves it and pulls in `pytorch` from
its own channel — which is the **CPU-only** build. Conda then has two `pytorch` specs
that conflict or silently chooses the CPU one.

Fixing it:

1. Declare `pytorch>=2.1.0` + `pytorch-cuda=12.4` from the `pytorch` / `nvidia`
   channels **first** in the deps list.
2. Install `phold` (and `phynteny_transformer`) via `pip:` so they inherit the already-
   resolved CUDA torch rather than triggering bioconda's solver.

---

## Sanity checks

```bash
# all phold rules use the placeholder, not hard-coded --cpu
grep -rn "params.cpu\|PHOLD_CPU_FLAG" sphae/workflow/

# CLI exposes the flag
sphae run --help | grep -- --use-gpu

# verify GPU is visible inside the phold env
conda run -n phold python -c "import torch; print(torch.cuda.is_available())"
# expected: True
```

---

## CUDA version note

`pytorch-cuda=12.6` works with any NVIDIA driver ≥ 525.60 (CUDA 12.x backward compat).
The RTX 4060 laptop (driver 610.47, CUDA UMD 13.3) easily satisfies this.

---

## WSL2 known issues and workarounds

These are one-time setup steps on Windows Subsystem for Linux 2. Once done they stay fixed.

### 1. Always run from `~/sphae-work`, never from `/mnt/c/`

Snakemake writes `.snakemake/` to the current directory. NTFS mounts (`/mnt/c/`, `/mnt/d/`)
don't support Linux `chmod` — Snakemake crashes with `PermissionError: Operation not permitted`.

```bash
mkdir -p ~/sphae-work && cd ~/sphae-work
sphae install --db_dir /mnt/d/sphae-databases ...
```

Databases can live on `/mnt/d/` (just passed as `--db_dir`). The *working directory* must
be on the ext4 filesystem (`/home/...`).

### 2. Install from this fork, not PyPI

The upstream PyPI package has the old `phold.yaml` (torch 2.1, phold 1.1) and `phynteny.yaml`
(no `requests`). Always install from the `gpu-support` branch:

```bash
pip install --force-reinstall git+https://github.com/sprinterz5/sphae.git@gpu-support
```

After reinstalling, delete stale conda envs so Snakemake rebuilds them from the new yamls:

```bash
rm -rf ~/sphae-work/.snakemake/conda/1367cd8cb4d8bb15a870f8137840eb73_*
rm -rf ~/sphae-work/.snakemake/conda/c4c7c12328bdd72e35b132f6b156a281_*
```

### 3. conda channel_priority must be `flexible`

`strict` breaks pharokka installation (biopython excluded, phanotate unsatisfiable).

```bash
conda config --set channel_priority flexible
```

This is persistent — set once, stays across sessions.

### 4. medaka always fails on WSL2 (Nanopore only — safe to skip)

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

### 5. phold CVE-2025-32434 (torch.load safety check)

`transformers` ≥ 4.51 blocks `torch.load` unless PyTorch ≥ 2.6. Old phold envs (torch 2.1)
crash with `ValueError: Due to a serious vulnerability issue in torch.load...`.

Fixed in this fork: `phold.yaml` pins `pytorch>=2.6.0`. If you see this error in an existing
env, delete the env dir (see §2 above) and rerun `sphae install`.

### 6. phynteny `ModuleNotFoundError: No module named 'requests'`

The `install_models` script imports `requests` at module level, but upstream `phynteny.yaml`
didn't list it as a dependency. Fixed in this fork: `requests` is now a conda dep.

---

## Dependency version rationale

| Package | Why this version |
|---------|-----------------|
| `pytorch>=2.6.0` | CVE-2025-32434: `transformers` blocks `torch.load` on older versions |
| `pytorch-cuda=12.6` | Matches driver 610.47 (CUDA UMD 13.3); backward compat ≥ driver 525 |
| `phold>=1.2.5` | Older phold had `download_requests()` signature bug (too few positional args) |
| `phynteny_transformer>=0.1.3` | Minimum version supporting PyTorch 2.x |
