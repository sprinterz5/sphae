# sphae — GPU fork

> **This is a GPU-enabled fork of [sphae](https://github.com/linsalrob/sphae).**
> Upstream sphae always runs phold in CPU mode. This fork adds `--use-gpu` / `--no-use-gpu`
> flags and installs CUDA-enabled PyTorch so phold runs on your NVIDIA GPU.
>
> Tested: RTX 4060 Laptop · driver 610.47 · CUDA 12.x · WSL2 · Windows 11

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/403889262.svg)](https://zenodo.org/doi/10.5281/zenodo.8365088)

---

## What's new in this fork

| Feature | Details |
|---------|---------|
| `--use-gpu` / `--no-use-gpu` | Controls whether phold uses GPU. Default: GPU on. |
| CUDA PyTorch in phold env | `torch==2.6.0+cu124` installed via pip wheel index — conda channels unreliably pick CPU-only builds. |
| CUDA PyTorch in phynteny env | Same treatment for phynteny_transformer. |
| `foldseek` as conda dep | pip install of phold omits the foldseek binary; it's now listed explicitly. |
| `--db_dir` for annotate | `sphae annotate` now accepts `--db_dir` and resolves all database paths automatically. |
| Version | `1.5.5+gpu` |

**Speed comparison (lambda phage, 88 proteins, RTX 4060 Laptop):**

| Step | CPU | GPU |
|------|-----|-----|
| phold predict | ~9 hours | ~3 minutes |
| Full pipeline | ~9+ hours | ~11 minutes |

For full technical details see [CHANGES_GPU.md](CHANGES_GPU.md).
For Windows/WSL2 setup from scratch see [WSL2_SETUP.md](WSL2_SETUP.md).

---

## Install

**Do not** use `pip install sphae` — that installs the upstream PyPI version without GPU support.

```bash
pip install git+https://github.com/sprinterz5/sphae.git@gpu-support
```

Verify:

```bash
sphae --version
# expected: sphae 1.5.5+gpu
```

---

## Install databases

```bash
sphae install --db_dir /path/to/databases --threads 8 --conda-frontend mamba
```

This downloads ~15 GB and takes 30–90 minutes. Databases:

| Database | Size | Used by |
|----------|------|---------|
| pharokka_db | ~280 MB | pharokka annotation |
| Pfam35.0 | ~280 MB | viral_verify contig classification |
| checkv-db-v1.5 | ~656 MB | completeness check |
| phold (ProstT5 + Foldseek) | ~5.7 GB | structural protein annotation (GPU) |
| phynteny models | ~100 MB | gene function prediction |
| medaka models | — | Nanopore only; fails on WSL2, safe to skip |

If `download_medaka_models` fails on WSL2, create a stub and rerun install:

```bash
mkdir -p /path/to/databases/medaka_models
touch /path/to/databases/medaka_models/medaka.flag
sphae install --db_dir /path/to/databases ...
```

---

## Running

### Annotate assembled genomes (most common use case)

`--genome` expects a **directory** containing `.fasta` / `.fa` / `.fna` files, not a single file.

```bash
mkdir -p genomes
cp my_phage.fasta genomes/

# GPU (default)
sphae annotate \
    --genome genomes \
    --output results \
    --db_dir /path/to/databases \
    --threads 8

# CPU fallback
sphae annotate \
    --genome genomes \
    --output results \
    --db_dir /path/to/databases \
    --threads 8 \
    --no-use-gpu
```

Multiple phages: put all `.fasta` files in the same directory — sphae processes them all in one run.

### Full pipeline (QC → assembly → annotation)

```bash
# Illumina paired-end reads
sphae run --input reads_dir --output results --use-gpu

# Nanopore
sphae run --input reads_dir --sequencing longread --output results --use-gpu

# Nanopore without medaka polishing (recommended on WSL2)
sphae run --input reads_dir --sequencing longread --output results --use-gpu --no_medaka
```

### Protein-only annotation

```bash
sphae annotate --protein protein_dir --output results --db_dir /path/to/databases
```

Note: phynteny is skipped in protein-only mode (no GBK file is generated).

---

## Output

### `sphae annotate` → `results/final-annotate/<sample>/`

| File | Description |
|------|-------------|
| `<sample>.gbk` | Full genome annotation (GenBank format, phynteny output) |
| `<sample>_summary.txt` | Summary: length, GC, completeness, taxonomy, gene hits |
| `<sample>_summary.functions` | Per-gene function table (pharokka + phold + phynteny) |
| `<sample>_phold_amr.tsv` | AMR gene hits (CARD database) |
| `<sample>_phold_vfdb.tsv` | Virulence factor hits |
| `<sample>_phold_acr.tsv` | Anti-CRISPR hits |
| `<sample>_phold_defense.tsv` | Defense mechanism hits |
| `plots/` | Circular genome visualization (phold plot) |

### `sphae run` → `results/RESULTS/<sample>/`

Same annotation outputs plus:

| File | Description |
|------|-------------|
| `<sample>.fasta` | Assembled genome (reoriented to terminase if found) |
| `trees/all_terL.nwk` | Terminase large subunit phylogeny |
| `trees/all_portal.nwk` | Portal protein phylogeny |
| `<sample>_phageterm/` | PhageTerm results (paired-end only) |

---

## Windows / WSL2 setup

Full step-by-step guide (DNS, conda, CUDA registration, database paths, common errors):
→ **[WSL2_SETUP.md](WSL2_SETUP.md)**

Quick checklist:
1. Run from `~/sphae-work`, never from `/mnt/c/` (NTFS breaks snakemake)
2. `conda config --set channel_priority flexible` (strict breaks pharokka)
3. Register CUDA libs once: `echo /usr/lib/wsl/lib | sudo tee /etc/ld.so.conf.d/wsl.conf && sudo ldconfig`
4. Install from this fork, not PyPI

---

## FAQ

**"Genome samples: []" — nothing to annotate**
You passed a file path to `--genome`. It must be a directory:
```bash
mkdir genomes && mv my_phage.fasta genomes/
sphae annotate --genome genomes ...
```

**"Foldseek not found. Please reinstall phold."**
Old conda env from a previous install. Delete it and rerun:
```bash
rm -rf ~/miniforge3/lib/python*/site-packages/sphae/workflow/conda/HASH_*
sphae annotate ...
```

**"No available GPU was found" / "Using device: cpu"**
PyTorch can't find the CUDA driver. Check:
```bash
ldconfig -p | grep libcuda   # should return /usr/lib/wsl/lib/libcuda.so.1
nvidia-smi                   # should show your GPU
```
If `ldconfig` returns nothing — run step 3 of the WSL2 checklist above.

**"ValueError: torch.load ... CVE-2025-32434"**
Old phold env with torch < 2.6. Delete old envs and reinstall:
```bash
pip install --force-reinstall git+https://github.com/sprinterz5/sphae.git@gpu-support
rm -rf ~/miniforge3/lib/python*/site-packages/sphae/workflow/conda/HASH_*
sphae annotate ...
```

**"PermissionError: Operation not permitted: '.../.snakemake/...'"**
Running from `/mnt/c/` or Desktop. Move to the Linux filesystem:
```bash
mkdir -p ~/sphae-work && cd ~/sphae-work
```

**"Failed during assembly"**
Assembly couldn't produce contigs — usually low coverage. Check logs at
`results/PROCESSING/assembly/flye/<sample>/assembly_info.txt`.

**"Genome includes multiple contigs, fragmented"**
Assembly produced short fragments instead of a complete genome.
Verify at `results/PROCESSING/assembly/flye/<sample>-assembly-stats_flye.csv`.
Options: resequence for better coverage, or try a different assembler.

**How do I visualize gene annotations?**
Use [Clinker](https://github.com/gamcil/clinker) on the GenBank files from `results/final-annotate/`.
For reoriented genomes, run [dnaapler](https://github.com/gbouras13/dnaapler) first, then `sphae annotate` again on the reoriented fastas.

**How do I change subsampling depth?**
```bash
sphae config   # copies config.yaml to current directory
# edit: bases: 10000000 → your value
sphae run --input ... --config config.yaml
```

---

## Upstream project

This fork is based on [sphae](https://github.com/linsalrob/sphae) by the Edwards Lab.

**Cite sphae:** https://doi.org/10.1093/bioadv/vbaf004

Issues with the GPU fork → open an issue in this repo.
Issues with core sphae functionality → open an issue in the [upstream repo](https://github.com/linsalrob/sphae/issues).
