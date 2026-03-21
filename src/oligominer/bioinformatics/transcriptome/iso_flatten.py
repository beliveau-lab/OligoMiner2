"""
Isoform flattening utilities.

Collapses transcript isoforms for each gene into the segments shared
across the maximum number of isoforms. This produces a simplified
annotation where each gene is represented by non-redundant exonic
regions, useful for probe design that targets all isoforms.
"""

import numpy as np
import pandas as pd


def flatten_isoforms(df):
    """
    Collapse transcript isoforms to maximally shared exonic segments.

    Groups exon records by gene_id and, for each gene, identifies the
    genomic intervals that are covered by the largest number of
    isoforms. Single-isoform genes are passed through unchanged.

    Args:
        df (pandas.DataFrame): exon-level annotation records with at
            least seqid, start, end, gene_id, score, strand, and
            transcript_id_full columns.

    Returns:
        flat_df (pandas.DataFrame): flattened annotation with columns
            seqid, start, end, gene_id, score, strand, coverage.
    """
    # group exons by gene
    grouped = df.groupby('gene_id').aggregate({
        'seqid': 'first',
        'start': lambda x: tuple(x),
        'end': lambda x: tuple(x),
        'score': 'first',
        'strand': lambda x: tuple(x),
        'transcript_id_full': lambda x: tuple(x),
    }).reset_index(drop=False)

    # flatten each gene's exons to shared segments
    segment_lists = grouped.apply(_flatten_gene_exons, axis=1)
    merged_segments = [seg for seg_list in segment_lists.values for seg in seg_list]

    flat_df = pd.DataFrame(
        merged_segments,
        columns=['seqid', 'start', 'end', 'gene_id', 'score', 'strand', 'coverage'],
    )

    # success
    return flat_df


def _flatten_gene_exons(row):
    """
    Flatten exon intervals for a single gene to maximally shared segments.

    For genes with overlapping exons across isoforms, finds the segments
    covered by the maximum number of isoforms. For genes with no overlap
    (max coverage == 1), the original exon coordinates are returned.

    Args:
        row (pandas.Series): aggregated gene record with tuple-valued
            start, end, and strand fields.

    Returns:
        segments (list): list of tuples, each containing
            (seqid, start, end, gene_id, score, strand, coverage).
    """
    # convert input coordinates to numpy arrays
    starts = np.array(row.start, dtype=int)
    ends = np.array(row.end, dtype=int)

    num_exons = starts.size
    min_start = np.min(starts)
    max_end = np.max(ends)

    # zero-shift coordinates for coverage calculation
    starts_norm = starts - min_start
    ends_norm = ends - min_start

    # build per-base coverage array
    coverage = np.zeros(max_end - min_start, dtype=int)
    for i in range(num_exons):
        coverage[starts_norm[i]:ends_norm[i]] += 1
    max_coverage = np.max(coverage)

    # determine output intervals
    output_starts = starts
    output_ends = ends
    if max_coverage > 1:

        # find start and end coordinates of maximally shared segments
        mask = (coverage == max_coverage).astype(int)
        shared_starts = np.where(np.append([mask[0]], np.diff(mask)) == 1)[0]
        shared_ends = np.where(np.diff(mask) == -1)[0] + 1
        if shared_ends.size < shared_starts.size:
            shared_ends = np.append(shared_ends, np.max(ends_norm))

        # restore absolute coordinates
        if shared_starts.size > 0:
            output_starts = shared_starts + min_start
            output_ends = shared_ends + min_start

    # format output records
    segments = []
    for i in range(output_starts.size):
        strand = row.strand[i] if i < len(row.strand) else row.strand[0]
        segment = (
            row.seqid,
            int(output_starts[i]),
            int(output_ends[i]),
            row.gene_id,
            row.score,
            strand,
            int(max_coverage),
        )
        segments.append(segment)

    # success
    return segments
