# OligoMiner2

**A Python toolkit for oligonucleotide probe design.**

OligoMiner2 is a library of composable parts for designing oligo probes: thermodynamic mining,
specificity analysis, duplex-stability models, and probe assembly. Each part is usable on its own,
so a probe design is written as a program against the toolkit rather than configured through one
fixed program.

The [end-to-end pipeline](pipeline.md) is one such program — the conventional FISH probe design
flow, assembled from these parts, and the place to start if it is the flow you want.

```{warning}
This package is under active development. Until the first stable release (`v1.0.0`) the API may
change without notice, so pin to a specific version if you depend on it.
```

## Install

```bash
pip install oligominer
```

That is everything except PyTorch. Mining, specificity, probe assembly and three of the four
duplex-stability models work with it. Only the duplex-BiLSTM needs torch, so it is an extra:

```bash
pip install "oligominer[torch]"
```

### CPU or GPU

There is no separate GPU extra, because there is nothing for it to install. On Linux the default
PyTorch wheel on PyPI already depends on the full CUDA stack, so `oligominer[torch]` is a
GPU-capable install and the BiLSTM uses a GPU when one is present. The GPU is the default, not an
upgrade.

What that costs is size: the CUDA stack runs to several GB. If you have no GPU, or want a small
install, take torch from PyTorch's CPU index first and then the extra is already satisfied:

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
pip install "oligominer[torch]"
```

This cannot be an extra — an extra names packages, and cannot say which index to take one from.

## What is in the toolkit

| Subsystem | What it does |
|---|---|
| `thermodynamics` | Melting temperature and formamide correction, candidate mining from sequence, exact NUPACK pDup |
| `specificity` | Alignment handling, duplex reconstruction from a reference, duplex stability, k-mer frequency |
| `models` | Four registered duplex-stability models selected by name, and retraining for two of them |
| `probe_design` | Probe sets, padlock and split architectures, domain layouts, exclusions, scoring, splitting, I/O |
| `bioinformatics` | Sequence and annotation I/O, transcriptome construction |

## Quick start

```python
from oligominer import mine_fasta

# mine candidate probes with default parameters
df = mine_fasta('genome.fa')
```

Mining is the entry point most designs begin from, but nothing requires the rest of the pipeline
to follow it — the resulting frame is ordinary tabular data.

```{toctree}
:maxdepth: 2
:hidden:

pipeline
api
changelog
```

## Where to go next

- [The pipeline](pipeline.md) — the end-to-end flow: mining, alignment, prediction, and what each
  stage decides.
- [API reference](api.md) — the toolkit, subsystem by subsystem, generated from the source.
- [Changelog](changelog.md) — what changed, newest first.
