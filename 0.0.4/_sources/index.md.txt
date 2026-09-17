# OligoMiner2

**Genome-scale oligonucleotide probe design for DNA and RNA FISH.**

OligoMiner2 takes a target genome or transcriptome and returns candidate probe sequences that are
thermodynamically optimized and filtered for specificity.

```{warning}
This package is under active development. Until the first stable release (`v1.0.0`) the API may
change without notice, so pin to a specific version if you depend on it.
```

## Install

```bash
pip install oligominer
```

The sequence model is optional, because torch is large and mining, screening and scoring with the
tree models do not need it:

```bash
pip install "oligominer[bilstm]"
```

## Quick start

```python
from oligominer import mine_fasta

# mine candidate probes with default parameters
df = mine_fasta('genome.fa')
```

```{toctree}
:maxdepth: 2
:hidden:

pipeline
api
changelog
```

## Where to go next

- [The pipeline](pipeline.md) — mining, alignment and prediction, and what each stage decides.
- [API reference](api.md) — every public module, generated from the source.
- [Changelog](changelog.md) — what changed, newest first.
