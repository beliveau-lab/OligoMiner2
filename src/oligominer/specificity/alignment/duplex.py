"""
# Genomic duplex reconstruction

Recovers the reference sequence under each alignment interval, which is what a
duplex-stability model scores a probe against.

The sequence is looked up in the genome rather than inferred from the CIGAR. A
CIGAR records how many bases matched, mismatched, inserted or deleted, not which
bases the reference holds: an 'M' operation covers both matches and mismatches, so
the reference base under a mismatch is absent from the alignment record.

Each chromosome is read once and held as a string, so every interval on it is a
slice. Intervals are grouped by chromosome, so a genome-scale table reads each
chromosome exactly once.

For a minus-strand alignment the derived sequence is the reverse complement of the
reference span.
"""

import numpy as np
import pandas as pd

from oligominer.utils.seq_utils import rev_comp

BED_COLUMNS = ['align_seqid', 'align_start', 'align_stop', 'seqid',
               'align_score', 'align_strand', 'align_cigar']


def load_reference(fasta_path):
    """
    Open a reference FASTA for repeated interval lookups.

    Args:
        fasta_path (str): path to the reference FASTA. A .fai index is created
            beside it if one is missing.

    Returns:
        fasta (pyfaidx.Fasta): the opened reference, with case preserved.
    """
    from oligominer.bioinformatics.file_io import load_fasta

    # success
    return load_fasta(str(fasta_path))


def chrom_sizes(fasta):
    """
    Return a chromosome name to length mapping.

    Args:
        fasta (pyfaidx.Fasta): an opened reference.

    Returns:
        sizes (dict): chromosome name mapped to its length in bases.
    """
    # success
    return {name: len(fasta[name]) for name in fasta.keys()}


def clamp_intervals(align_df, sizes):
    """
    Clamp alignment intervals to chromosome boundaries.

    An alignment near a chromosome end can produce a span running past the end of
    the sequence. Fetching such a span either errors or silently returns a short
    sequence, so intervals are clamped and the affected rows are flagged.

    Args:
        align_df (pandas.DataFrame): alignments with align_seqid, align_start and
            align_stop columns.
        sizes (dict): chromosome name mapped to length.

    Returns:
        out (pandas.DataFrame): a copy with clamped coordinates and a was_clamped
            boolean column.

    Raises:
        KeyError: if an alignment names a chromosome absent from the reference.
    """
    out = align_df.copy()
    limits = out['align_seqid'].map(sizes)

    if limits.isna().any():
        missing = sorted(out.loc[limits.isna(), 'align_seqid'].unique())[:5]
        raise KeyError(
            f'alignments reference chromosomes absent from the FASTA: {missing}')

    start = out['align_start'].to_numpy(dtype=np.int64)
    stop = out['align_stop'].to_numpy(dtype=np.int64)
    limit = limits.to_numpy(dtype=np.int64)

    new_start = np.clip(start, 0, limit)
    new_stop = np.clip(stop, 0, limit)

    out['was_clamped'] = (new_start != start) | (new_stop != stop)
    out['align_start'] = new_start
    out['align_stop'] = new_stop

    # success
    return out


def fetch_derived_seqs(align_df, fasta, to_upper=True):
    """
    Look up the reference sequence under every alignment interval.

    Args:
        align_df (pandas.DataFrame): clamped intervals with align_seqid,
            align_start, align_stop and align_strand columns.
        fasta (pyfaidx.Fasta): the opened reference.
        to_upper (bool): upper-case the returned sequences.

    Returns:
        seqs (pandas.Series): the derived sequence per row, indexed like align_df.
    """
    seqs = pd.Series('', index=align_df.index, dtype=object)

    starts = align_df['align_start'].to_numpy(dtype=np.int64)
    stops = align_df['align_stop'].to_numpy(dtype=np.int64)
    if 'align_strand' in align_df:
        strands = align_df['align_strand'].to_numpy()
    else:
        strands = np.full(len(align_df), '+')

    for chrom, block in align_df.groupby('align_seqid', sort=False):
        # materialize the chromosome once so each interval is a string slice
        chrom_seq = str(fasta[chrom])
        positions = align_df.index.get_indexer(block.index)
        values = []
        for i in positions:
            piece = chrom_seq[starts[i]:stops[i]]
            if to_upper:
                piece = piece.upper()
            if strands[i] == '-':
                piece = rev_comp(piece.upper())
            values.append(piece)
        seqs.loc[block.index] = values

    # success
    return seqs


def reconstruct(align_df, fasta_path, to_upper=True):
    """
    Add the reference sequence under each alignment to an alignment table.

    Args:
        align_df (pandas.DataFrame): alignment records carrying BED_COLUMNS.
        fasta_path (str): path to the reference FASTA.
        to_upper (bool): upper-case the derived sequences.

    Returns:
        out (pandas.DataFrame): a copy with derived_seq and was_clamped added.
    """
    fasta = load_reference(fasta_path)
    out = clamp_intervals(align_df, chrom_sizes(fasta))
    out['derived_seq'] = fetch_derived_seqs(out, fasta, to_upper=to_upper)

    # success
    return out
