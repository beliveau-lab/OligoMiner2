# OligoMiner2

**Genome-scale oligonucleotide probe design for DNA and RNA FISH.**

OligoMiner2 is a Python package for designing oligonucleotide probes used in
fluorescence *in situ* hybridization (FISH) experiments. It takes a target
genome or transcriptome as input and returns candidate probe sequences that are
thermodynamically optimized and filtered for specificity.

## The pipeline

1. **Mine** candidate probes from FASTA sequences, filtering by melting
   temperature, GC content, length, homopolymer runs, and prohibited
   subsequences.
2. **Align** candidates against a reference genome with Bowtie2 to identify
   off-target binding sites.
3. **Score** specificity using k-mer frequency analysis (Jellyfish) and
   thermodynamic duplex stability predictions (NUPACK / XGBoost).

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
