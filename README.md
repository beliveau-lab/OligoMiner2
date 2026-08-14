# OligoMiner2


> [!WARNING]
> This package is under active development.
> Until the first stable release (`v1.0.0`), the API may change without notice.
> Expect breaking changes in minor releases, and pin your dependency to a specific version if you use it in production.

**Genome-scale oligonucleotide probe design for DNA and RNA FISH.**

OligoMiner2 is a Python package for designing oligonucleotide probes used in
fluorescence *in situ* hybridization (FISH) experiments. It takes a target
genome or transcriptome as input and returns candidate probe sequences that are
thermodynamically optimized and filtered for specificity.

## The pipeline

1. **Mine** candidate probes from FASTA sequences, filtering by melting
   temperature, GC content, length, homopolymer runs, sequence entropy,
   soft-masked repeats, and prohibited subsequences.
2. **Align** candidates against a reference genome with Bowtie2 and
   reconstruct each alignment's duplex from the reference.
3. **Score** specificity with k-mer frequency analysis and duplex stability
   predictions, optionally verifying the alignments that matter with exact
   NUPACK pDup.

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
rather than at an end, and screened for mutual orthogonality.

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
