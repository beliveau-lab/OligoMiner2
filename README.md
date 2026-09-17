# OligoMiner2

<!-- Badges and links are ABSOLUTE URLs: PyPI renders this file standalone, where a relative
     path resolves against pypi.org and breaks. -->
[![PyPI](https://img.shields.io/pypi/v/oligominer.svg)](https://pypi.org/project/oligominer/)
[![Python versions](https://img.shields.io/pypi/pyversions/oligominer.svg)](https://pypi.org/project/oligominer/)
[![License](https://img.shields.io/pypi/l/oligominer.svg)](https://github.com/beliveau-lab/OligoMiner2/blob/main/LICENSE)
[![CI](https://github.com/beliveau-lab/OligoMiner2/actions/workflows/ci.yml/badge.svg?branch=dev)](https://github.com/beliveau-lab/OligoMiner2/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/badge/docs-oligominer.org-blue.svg)](https://oligominer.org/)

> [!WARNING]
> This package is under active development.
> Until the first stable release (`v1.0.0`), the API may change without notice.
> Expect breaking changes in minor releases, and pin your dependency to a specific version if you use it in production.

**A Python toolkit for oligonucleotide probe design.**

OligoMiner2 is a library of composable parts for designing oligo probes:
thermodynamic mining, specificity analysis, duplex-stability models, and probe
assembly. Each part is usable on its own, so a probe design is written as a
program against the toolkit rather than configured through one fixed program.

The end-to-end pipeline below is one such program: the conventional FISH probe
design flow, assembled from these parts, and the place to start if it is the
flow you want.

## What is in the toolkit

| Subsystem | What it does |
|---|---|
| `thermodynamics` | Melting temperature and formamide correction, candidate mining from sequence, exact NUPACK pDup |
| `specificity` | Alignment handling, duplex reconstruction from a reference, duplex stability, k-mer frequency |
| `models` | Four registered duplex-stability models selected by name, and retraining for two of them |
| `probe_design` | Probe sets, padlock and split architectures, domain layouts, exclusions, scoring, splitting, I/O |
| `bioinformatics` | Sequence and annotation I/O, transcriptome construction |

## The pipeline

OligoMiner2 takes a target genome or transcriptome and returns candidate probe
sequences that are thermodynamically optimized and filtered for specificity.

1. **Mine** candidate probes from FASTA sequences, filtering by melting
   temperature, GC content, length, homopolymer runs, sequence entropy,
   soft-masked repeats, and prohibited subsequences.
2. **Align** candidates against a reference genome with Bowtie2 and
   reconstruct each alignment's duplex from the reference.
3. **Predict** a binding probability for every reconstructed duplex, screening
   with a model and verifying the ones that could change a decision with exact
   NUPACK pDup. K-mer frequency analysis runs as a parallel track.

### Targets

Probes can be mined against a genome, a transcriptome built from a genome and
its annotation, or individual transcripts. Transcript-aware modes cover
exon-exon junctions, regions that discriminate between isoforms, and introns.

### Probe architectures

Beyond conventional single-oligo probes, the package designs padlock probes
with split homology arms, and split architectures where one functional element
spans two oligos binding adjacent sites -- HCR 3.0 split-initiator pairs and
split-FISH bridge pairs. Both members of a pair are tracked as one targeting
unit, so filtering cannot leave a half-probe in an order.

Synthetic domains -- barcodes, primer sites, amplifier initiators -- are
appended through a domain layout, which places a sequence at a named slot
rather than at an end.

### Models

Duplex stability is predicted by one of four registered models, selected by
name: a condition-aware physics-feature gradient-boosted model, a duplex
BiLSTM, the PaintSHOP model, and the OligoMiner v1 classifier. The first two
can be retrained on new data through the package.

### Scale

K-mer counting dispatches between jellyfish and a pure-numpy backend. Alignment
and duplex construction stream rather than loading a genome's alignments into
memory, and every entry point resolves its core count from the scheduler's
grant rather than from the machine size.

## Quick start

```python
from oligominer import mine_fasta

# mine candidate probes with default parameters
df = mine_fasta('genome.fa')
```


## Documentation

Full documentation including API reference and example notebooks is available at the
[OligoMiner2 docs site](https://oligominer.org/).


## License

OligoMiner2 is open-source software. See [LICENSE](LICENSE) for details.
