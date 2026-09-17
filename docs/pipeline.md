# The pipeline

1. **Mine** candidate probes from FASTA sequences, filtering by melting temperature, GC content,
   length, homopolymer runs, sequence entropy, soft-masked repeats, and prohibited subsequences.
2. **Align** candidates against a reference genome with Bowtie2 and reconstruct each alignment's
   duplex from the reference.
3. **Predict** a binding probability for every reconstructed duplex, screening with a model and
   verifying the ones that could change a decision with exact NUPACK pDup. K-mer frequency
   analysis runs as a parallel track.

## Targets

Probes can be mined against a genome, a transcriptome built from a genome and its annotation, or
individual transcripts. Transcript-aware modes cover exon-exon junctions, regions that discriminate
between isoforms, and introns.

## Probe architectures

Beyond conventional single-oligo probes, the package designs padlock probes with split homology
arms, and split architectures where one functional element spans two oligos binding adjacent
sites — HCR 3.0 split-initiator pairs and split-FISH bridge pairs. Both members of a pair are
tracked as one targeting unit, so filtering cannot leave a half-probe in an order.

Synthetic domains — barcodes, primer sites, amplifier initiators — are appended through a domain
layout, which places a sequence at a named slot rather than at an end.

## Models

Duplex stability is predicted by one of four registered models, selected by name: a
condition-aware physics-feature gradient-boosted model, a duplex BiLSTM, the PaintSHOP model, and
the OligoMiner v1 classifier. The first two can be retrained on new data through the package.

```{note}
The BiLSTM needs torch, which is not installed by default. Install `oligominer[bilstm]` to use it;
every other model works without it.
```

## Scale

K-mer counting dispatches between jellyfish and a pure-numpy backend. Alignment and duplex
construction stream rather than loading a genome's alignments into memory, and every entry point
resolves its core count from the scheduler's grant rather than from the machine size.

## External programs

Some stages shell out to programs that are not Python packages and are installed separately:

| Program | Used by | Notes |
|---|---|---|
| Bowtie2 | alignment | required for the specificity track |
| jellyfish | k-mer counting | optional; a pure-numpy backend is used when it is absent |
| NUPACK | exact pDup verification | licence-gated, so it is installed manually rather than from PyPI |
