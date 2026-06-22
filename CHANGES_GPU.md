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
| `sphae/workflow/envs/phold.yaml` | Added `pytorch` + `nvidia` channels; added `pytorch>=2.1.0` and `pytorch-cuda=12.4` as conda deps; moved `phold>=1.1.0` to `pip:` section to prevent bioconda overriding CUDA torch. |
| `sphae/workflow/envs/phynteny.yaml` | Same treatment: CUDA torch via conda, `phynteny_transformer>=0.1.3` via pip. |

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

`pytorch-cuda=12.4` works with any NVIDIA driver ≥ 525.60 (CUDA 12.x backward compat).
The RTX 4060 laptop (driver 610.47, CUDA UMD 13.3) easily satisfies this.
If you need a newer CUDA wheel later, replace `12.4` with `12.6` (available from
pytorch nightly) — the rest of the setup stays the same.
