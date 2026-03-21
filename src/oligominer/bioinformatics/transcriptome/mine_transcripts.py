"""
RNA transcript probe mining.

Convenience functions for mining probes from transcript features. Each
function extracts the relevant genomic sequences (exons, introns, or
spliced transcript), mines them, and returns a ProbeSet ready for the
standard alignment/scoring pipeline.

Two mining strategies are supported:

  - **Per-interval mining** (mine_exons, mine_introns, mine_flattened_gene):
    each exon or intron is mined independently as a genomic interval.
    Probe seq_id encodes the genomic locus so coordinates can be mapped
    back via local_to_genomic().

  - **Spliced transcript mining** (mine_spliced_transcript): exon
    sequences are concatenated into a virtual mRNA and mined as a single
    sequence. This allows the miner to consider Tm and spacing across
    exon-exon junctions. Probe coordinates are relative to the spliced
    transcript, not the genome.
"""

import os
import tempfile

from oligominer.bioinformatics.file_io import write_fasta
from oligominer.probe_design.probe_set import ProbeSet
from oligominer.thermodynamics.mining import mine_sequence

from .transcript_seq import (
    get_exon_seqs,
    get_intron_seqs,
    get_spliced_seq,
    get_flattened_seqs,
)


# ---------------------------------------------------------------------------
# Per-interval mining
# ---------------------------------------------------------------------------

def mine_exons(gtf_df, fasta_path, transcript_id=None, gene_id=None,
               cores=1, **mining_params):
    """
    Mine probes from exon sequences of a transcript or gene.

    Each exon is mined independently. Probe seq_id values encode the
    genomic locus (e.g. 'chrI:1807-2169(-)') so that coordinates can
    be converted to genomic positions via local_to_genomic().

    Provide either transcript_id for a single isoform or gene_id for
    all exons across all isoforms.

    Args:
        gtf_df (pandas.DataFrame): parsed, exon-filtered GTF with
            parsed attributes (see gtf_io.parse_attributes).
        fasta_path (str or pyfaidx.Fasta): path to the genome FASTA file,
            or a pre-loaded pyfaidx.Fasta object.
        transcript_id (str or None): mine exons for this transcript.
        gene_id (str or None): mine exons for this gene.
        cores (int): number of CPU cores for parallel mining.
        **mining_params: forwarded to mine_sequence() (min_length,
            max_length, min_tm, max_tm, etc.).

    Returns:
        probe_set (ProbeSet): mined probes from all exons.
    """
    seqs = get_exon_seqs(
        gtf_df, fasta_path,
        transcript_id=transcript_id, gene_id=gene_id,
    )
    probe_set = _mine_seq_dict(seqs, cores=cores, **mining_params)

    # success
    return probe_set


def mine_introns(gtf_df, fasta_path, transcript_id, cores=1,
                 **mining_params):
    """
    Mine probes from intron sequences of a transcript.

    Introns are derived from gaps between consecutive exons. Useful for
    designing probes that detect nascent (unspliced) pre-mRNA, a common
    approach for visualizing active transcription sites.

    Args:
        gtf_df (pandas.DataFrame): parsed, exon-filtered GTF with
            parsed attributes.
        fasta_path (str or pyfaidx.Fasta): path to the genome FASTA file,
            or a pre-loaded pyfaidx.Fasta object.
        transcript_id (str): the transcript whose introns to mine.
        cores (int): number of CPU cores for parallel mining.
        **mining_params: forwarded to mine_sequence().

    Returns:
        probe_set (ProbeSet): mined probes from all introns.
    """
    seqs = get_intron_seqs(gtf_df, fasta_path, transcript_id)
    probe_set = _mine_seq_dict(seqs, cores=cores, **mining_params)

    # success
    return probe_set


def mine_flattened_gene(flat_df, fasta_path, gene_id, cores=1,
                        **mining_params):
    """
    Mine probes from the flattened (pan-isoform) exonic segments of a gene.

    Uses the output of flatten_isoforms() to mine probes from regions
    shared across the maximum number of transcript isoforms. This
    maximizes the chance that probes will detect all isoforms of a gene.

    Args:
        flat_df (pandas.DataFrame): output of flatten_isoforms().
        fasta_path (str or pyfaidx.Fasta): path to the genome FASTA file,
            or a pre-loaded pyfaidx.Fasta object.
        gene_id (str): the gene to mine.
        cores (int): number of CPU cores for parallel mining.
        **mining_params: forwarded to mine_sequence().

    Returns:
        probe_set (ProbeSet): mined probes from flattened segments.
    """
    seqs = get_flattened_seqs(flat_df, fasta_path, gene_id)
    probe_set = _mine_seq_dict(seqs, cores=cores, **mining_params)

    # success
    return probe_set


# ---------------------------------------------------------------------------
# Spliced transcript mining
# ---------------------------------------------------------------------------

def mine_spliced_transcript(gtf_df, fasta_path, transcript_id, cores=1,
                            **mining_params):
    """
    Mine probes from a spliced (in silico) transcript sequence.

    Concatenates exon sequences into a virtual mRNA and mines it as a
    single continuous sequence. This allows the miner to consider Tm
    and spacing constraints across exon-exon junctions.

    Probe coordinates in the returned ProbeSet are relative to the
    spliced transcript (not the genome). The seq_id is set to the
    transcript_id.

    Args:
        gtf_df (pandas.DataFrame): parsed, exon-filtered GTF with
            parsed attributes.
        fasta_path (str or pyfaidx.Fasta): path to the genome FASTA file,
            or a pre-loaded pyfaidx.Fasta object.
        transcript_id (str): the transcript to build and mine.
        cores (int): number of CPU cores for parallel mining.
        **mining_params: forwarded to mine_sequence().

    Returns:
        probe_set (ProbeSet): mined probes from the spliced transcript.
    """
    spliced_seq = get_spliced_seq(gtf_df, fasta_path, transcript_id)

    probes = mine_sequence(
        spliced_seq,
        seq_id=transcript_id,
        cores=cores,
        **mining_params,
    )
    probe_set = ProbeSet(probes)

    # success
    return probe_set


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _mine_seq_dict(seqs, cores=1, **mining_params):
    """
    Mine probes from a dict of {label: sequence} and return a ProbeSet.

    Each sequence is mined independently using the label as seq_id.
    Results are aggregated into a single ProbeSet.

    Args:
        seqs (dict): {feature_label: dna_sequence} pairs.
        cores (int): number of CPU cores for parallel mining.
        **mining_params: forwarded to mine_sequence().

    Returns:
        probe_set (ProbeSet): mined probes from all sequences.
    """
    all_probes = []
    for label, seq in seqs.items():
        if not seq:
            continue
        probes = mine_sequence(seq, seq_id=label, cores=cores, **mining_params)
        all_probes.extend(probes)

    probe_set = ProbeSet(all_probes)

    # success
    return probe_set
