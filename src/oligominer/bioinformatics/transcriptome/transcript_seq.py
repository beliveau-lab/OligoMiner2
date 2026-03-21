"""
Transcript sequence extraction primitives.

Given a parsed GTF annotation DataFrame and a genome FASTA, extracts
genomic sequences for transcript features: exons, introns, spliced
transcripts, and flattened gene models. All functions are strand-aware
and return sequences oriented 5' to 3' on the transcript.

These primitives are the building blocks for RNA probe mining. Each
function returns a dict of {feature_id: sequence} suitable for
write_fasta / split_fasta or direct use with mine_sequence.
"""

import numpy as np
from pyfaidx import Fasta

from oligominer.bioinformatics.file_io import load_fasta
from oligominer.utils.seq_utils import rev_comp


# ---------------------------------------------------------------------------
# Exon sequences
# ---------------------------------------------------------------------------

def get_exon_seqs(gtf_df, fasta, transcript_id=None, gene_id=None):
    """
    Extract individual exon sequences from the genome.

    Returns one sequence per exon, keyed by a genomic coordinate string
    that encodes the locus and strand. Sequences are returned on the
    transcript strand (reverse-complemented for minus-strand genes).

    Provide either transcript_id to get exons for a single isoform, or
    gene_id to get all exons across all isoforms of a gene.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF with columns seqid, start,
            end, strand, transcript_id, gene_id. May contain mixed feature
            types (exons are selected automatically).
        fasta (str or pyfaidx.Fasta): path to the genome FASTA file, or
            a pre-loaded pyfaidx.Fasta object.
        transcript_id (str or None): select exons for this transcript.
        gene_id (str or None): select exons for this gene.

    Returns:
        exon_seqs (dict): {exon_label: sequence} where exon_label is
            formatted as 'seqid:start-end(strand)'.
    """
    exons = _select_features(gtf_df, transcript_id=transcript_id, gene_id=gene_id)
    fasta = _resolve_fasta(fasta)

    exon_seqs = {}
    for _, row in exons.iterrows():
        label = _interval_label(row['seqid'], row['start'], row['end'], row['strand'])
        seq = _extract_seq(fasta, row['seqid'], row['start'], row['end'], row['strand'])
        exon_seqs[label] = seq

    # success
    return exon_seqs


# ---------------------------------------------------------------------------
# Intron sequences
# ---------------------------------------------------------------------------

def get_intron_seqs(gtf_df, fasta, transcript_id):
    """
    Extract intron sequences for a specific transcript isoform.

    Introns are derived from the gaps between consecutive exons of the
    transcript, sorted by genomic position. For single-exon transcripts,
    an empty dict is returned. Sequences are returned on the transcript
    strand.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF with columns seqid, start,
            end, strand, transcript_id. May contain mixed feature types
            (exons are selected automatically).
        fasta (str or pyfaidx.Fasta): path to the genome FASTA file, or
            a pre-loaded pyfaidx.Fasta object.
        transcript_id (str): the transcript to extract introns for.

    Returns:
        intron_seqs (dict): {intron_label: sequence} where intron_label
            is formatted as 'seqid:start-end(strand)'.
    """
    exons = _select_features(gtf_df, transcript_id=transcript_id)
    if len(exons) < 2:
        return {}

    fasta = _resolve_fasta(fasta)
    chrom = exons.iloc[0]['seqid']
    strand = exons.iloc[0]['strand']

    # sort exons by genomic position to find gaps
    exons_sorted = exons.sort_values('start')
    ends = exons_sorted['end'].values
    starts = exons_sorted['start'].values

    intron_seqs = {}
    for i in range(len(exons_sorted) - 1):
        intron_start = int(ends[i])
        intron_end = int(starts[i + 1])
        if intron_end <= intron_start:
            continue
        label = _interval_label(chrom, intron_start, intron_end, strand)
        seq = _extract_seq(fasta, chrom, intron_start, intron_end, strand)
        intron_seqs[label] = seq

    # success
    return intron_seqs


# ---------------------------------------------------------------------------
# Spliced transcript sequence
# ---------------------------------------------------------------------------

def get_spliced_seq(gtf_df, fasta, transcript_id):
    """
    Construct the spliced mRNA sequence for a transcript by concatenating
    its exon sequences in transcript order (5' to 3').

    For plus-strand genes, exons are concatenated in ascending genomic
    order. For minus-strand genes, exons are concatenated in descending
    genomic order, and each exon is reverse-complemented before joining.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF with columns seqid, start,
            end, strand, transcript_id. May contain mixed feature types
            (exons are selected automatically).
        fasta (str or pyfaidx.Fasta): path to the genome FASTA file, or
            a pre-loaded pyfaidx.Fasta object.
        transcript_id (str): the transcript to build.

    Returns:
        spliced_seq (str): the full spliced transcript sequence, 5' to 3'.
    """
    exons = _select_features(gtf_df, transcript_id=transcript_id)
    fasta = _resolve_fasta(fasta)
    strand = exons.iloc[0]['strand']

    # sort exons in transcript order
    ascending = (strand == '+')
    exons_sorted = exons.sort_values('start', ascending=ascending)

    # concatenate exon sequences
    parts = []
    for _, row in exons_sorted.iterrows():
        seq = _extract_seq(fasta, row['seqid'], row['start'], row['end'], row['strand'])
        parts.append(seq)

    spliced_seq = ''.join(parts)

    # success
    return spliced_seq


# ---------------------------------------------------------------------------
# Flattened gene-model sequences
# ---------------------------------------------------------------------------

def get_flattened_seqs(flat_df, fasta, gene_id):
    """
    Extract sequences for the flattened (pan-isoform) exonic segments
    of a gene.

    Takes the output of flatten_isoforms() and extracts the corresponding
    genomic sequences. Each segment represents a region shared across the
    maximum number of isoforms.

    Args:
        flat_df (pandas.DataFrame): output of flatten_isoforms() with
            columns seqid, start, end, gene_id, strand.
        fasta (str or pyfaidx.Fasta): path to the genome FASTA file, or
            a pre-loaded pyfaidx.Fasta object.
        gene_id (str): the gene to extract flattened segments for.

    Returns:
        segment_seqs (dict): {segment_label: sequence} where segment_label
            is formatted as 'seqid:start-end(strand)'.
    """
    segments = flat_df[flat_df['gene_id'] == gene_id]
    fasta = _resolve_fasta(fasta)

    segment_seqs = {}
    for _, row in segments.iterrows():
        label = _interval_label(row['seqid'], row['start'], row['end'], row['strand'])
        seq = _extract_seq(fasta, row['seqid'], row['start'], row['end'], row['strand'])
        segment_seqs[label] = seq

    # success
    return segment_seqs


# ---------------------------------------------------------------------------
# Coordinate parsing
# ---------------------------------------------------------------------------

def parse_interval_label(label):
    """
    Parse a genomic interval label back into its components.

    Accepts labels in the format 'seqid:start-end(strand)' as produced
    by get_exon_seqs, get_intron_seqs, and get_flattened_seqs.

    Args:
        label (str): a genomic interval label.

    Returns:
        seqid (str): the chromosome or sequence identifier.
        start (int): the start coordinate (0-based).
        end (int): the end coordinate.
        strand (str): '+' or '-'.
    """
    # split 'chrI:100-200(+)' into components
    coord_part, strand = label.rstrip(')').rsplit('(', 1)
    seqid, span = coord_part.rsplit(':', 1)
    start, end = span.split('-')

    # success
    return seqid, int(start), int(end), strand


def local_to_genomic(seq_id, local_start, local_stop):
    """
    Convert probe coordinates local to a mined interval back to absolute
    genomic coordinates.

    The seq_id is expected to be a genomic interval label as produced by
    get_exon_seqs etc. (e.g. 'chrI:1807-2169(-)').

    For plus-strand intervals, genomic_start = interval_start + local_start.
    For minus-strand intervals, the coordinates are mirrored so that
    genomic coordinates remain in standard ascending order.

    Args:
        seq_id (str): the interval label used as seq_id during mining.
        local_start (int): 0-based start within the mined sequence.
        local_stop (int): 0-based stop within the mined sequence.

    Returns:
        seqid (str): chromosome name.
        genomic_start (int): absolute genomic start coordinate.
        genomic_stop (int): absolute genomic stop coordinate.
        strand (str): '+' or '-'.
    """
    seqid, interval_start, interval_end, strand = parse_interval_label(seq_id)

    if strand == '+':
        genomic_start = interval_start + local_start
        genomic_stop = interval_start + local_stop
    else:
        # mining runs on the rev_comp, so coordinates are mirrored
        interval_len = interval_end - interval_start
        genomic_start = interval_start + (interval_len - local_stop)
        genomic_stop = interval_start + (interval_len - local_start)

    # success
    return seqid, genomic_start, genomic_stop, strand


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _resolve_fasta(fasta):
    """
    Accept a file path or a pre-loaded pyfaidx.Fasta and return a Fasta object.

    This avoids redundant index loads when the caller already has a Fasta
    open (e.g. when calling multiple extraction functions on the same genome).

    Args:
        fasta (str or pyfaidx.Fasta): path to a genome FASTA file, or an
            already-loaded pyfaidx.Fasta object.

    Returns:
        fasta (pyfaidx.Fasta): the loaded FASTA.
    """
    if isinstance(fasta, Fasta):
        return fasta

    # success
    return load_fasta(fasta)


def _select_features(gtf_df, transcript_id=None, gene_id=None):
    """
    Select exon records from a GTF DataFrame by transcript or gene.

    Filters to exon-type records only. If the input GTF contains mixed
    feature types (exon, CDS, transcript, etc.), only exons are returned.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF with transcript_id and
            gene_id columns. May optionally have a 'type' column for
            feature-type filtering.
        transcript_id (str or None): filter to this transcript.
        gene_id (str or None): filter to this gene.

    Returns:
        features (pandas.DataFrame): the selected exon records.

    Raises:
        ValueError: if neither transcript_id nor gene_id is provided,
            or if the selection is empty.
    """
    # filter to exons if the type column is present
    if 'type' in gtf_df.columns:
        gtf_df = gtf_df[gtf_df['type'] == 'exon']

    if transcript_id is not None:
        features = gtf_df[gtf_df['transcript_id'] == transcript_id]
    elif gene_id is not None:
        features = gtf_df[gtf_df['gene_id'] == gene_id]
    else:
        raise ValueError("provide either transcript_id or gene_id")

    if features.empty:
        query = transcript_id or gene_id
        raise ValueError(f"no exon records found for '{query}'")

    # success
    return features


def _extract_seq(fasta, seqid, start, end, strand):
    """
    Extract a genomic sequence, reverse-complementing for minus strand.

    Args:
        fasta (pyfaidx.Fasta): loaded genome FASTA.
        seqid (str): chromosome name.
        start (int): 0-based start coordinate.
        end (int): end coordinate.
        strand (str): '+' or '-'.

    Returns:
        seq (str): the extracted sequence on the transcript strand.
    """
    seq = str(fasta[seqid][int(start):int(end)])
    if strand == '-':
        seq = rev_comp(seq)

    # success
    return seq


def _interval_label(seqid, start, end, strand):
    """
    Format a genomic interval as a string label.

    Args:
        seqid (str): chromosome name.
        start (int): start coordinate.
        end (int): end coordinate.
        strand (str): '+' or '-'.

    Returns:
        label (str): formatted as 'seqid:start-end(strand)'.
    """
    label = f'{seqid}:{int(start)}-{int(end)}({strand})'

    # success
    return label
