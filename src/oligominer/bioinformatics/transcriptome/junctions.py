"""
# Exon-exon junctions and isoform-discriminating regions

Two capabilities that need the transcript rather than the genome.

A probe spanning an exon-exon junction only binds spliced mRNA, because the two
halves of its target are not adjacent in the genome. `junction_windows` returns
the sequence around each junction of a transcript, together with the offset of
the junction within that window, so a probe can be required to straddle it.

`discriminating_regions` is the counterpart of isoform flattening. Flattening
finds the exonic sequence shared across a gene's isoforms, which is what a probe
set targeting the whole gene wants. Distinguishing one isoform from its siblings
needs the opposite: the transcript intervals that no sibling covers.
"""

import pandas as pd

from .transcript_seq import _resolve_fasta, _select_features, get_spliced_seq

# bases of each flanking exon included in a junction window by default
DEFAULT_FLANK = 40

JUNCTION_COLUMNS = ['transcript_id', 'junction_index', 'junction_offset',
                    'window_start', 'window_stop', 'seq']

DISCRIMINATING_COLUMNS = ['gene_id', 'transcript_id', 'seqid', 'start', 'end',
                          'strand', 'length']


def exon_order(gtf_df, transcript_id):
    """
    Return a transcript's exons in transcript order, 5' to 3'.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF carrying exon records.
        transcript_id (str): the transcript to select.

    Returns:
        exons (pandas.DataFrame): the exons, sorted 5' to 3' along the transcript.
    """
    exons = _select_features(gtf_df, transcript_id=transcript_id)
    ascending = exons.iloc[0]['strand'] == '+'

    # success
    return exons.sort_values('start', ascending=ascending).reset_index(drop=True)


def junction_offsets(gtf_df, transcript_id):
    """
    Return where each exon-exon junction falls in the spliced transcript.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF carrying exon records.
        transcript_id (str): the transcript to select.

    Returns:
        offsets (list): the transcript coordinate of each junction, one per
            adjacent exon pair. Empty for a single-exon transcript.
    """
    exons = exon_order(gtf_df, transcript_id)
    lengths = (exons['end'] - exons['start'] + 1).tolist()

    offsets = []
    running = 0
    for length in lengths[:-1]:
        running += length
        offsets.append(running)

    # success
    return offsets


def junction_windows(gtf_df, fasta, transcript_id, flank=DEFAULT_FLANK):
    """
    Return the spliced sequence around each exon-exon junction.

    A window is the last `flank` bases of one exon followed by the first `flank`
    bases of the next, taken from the spliced transcript so the two halves are
    adjacent. A probe mined from this window is junction-spanning only if it
    covers the junction offset, which `spans_junction` checks.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF carrying exon records.
        fasta (str or pyfaidx.Fasta): the genome FASTA, or a loaded Fasta.
        transcript_id (str): the transcript to select.
        flank (int): bases of each flanking exon to include.

    Returns:
        windows (pandas.DataFrame): one row per junction with JUNCTION_COLUMNS.
            junction_offset is the offset within the window, not the transcript.
    """
    spliced = get_spliced_seq(gtf_df, _resolve_fasta(fasta), transcript_id)
    offsets = junction_offsets(gtf_df, transcript_id)

    rows = []
    for index, offset in enumerate(offsets):
        start = max(0, offset - flank)
        stop = min(len(spliced), offset + flank)
        rows.append({
            'transcript_id': transcript_id,
            'junction_index': index,
            'junction_offset': offset - start,
            'window_start': start,
            'window_stop': stop,
            'seq': spliced[start:stop],
        })

    # success
    return pd.DataFrame(rows, columns=JUNCTION_COLUMNS)


def spans_junction(start, stop, junction_offset, min_overhang=1):
    """
    Report whether a probe covers a junction with enough sequence on both sides.

    Args:
        start (int): probe start within the window.
        stop (int): probe stop within the window, exclusive.
        junction_offset (int): the junction's offset within the window.
        min_overhang (int): bases the probe must place on each side.

    Returns:
        spans (bool): True when the probe straddles the junction.
    """
    # success
    return (start <= junction_offset - min_overhang
            and stop >= junction_offset + min_overhang)


def junction_probes(probes, junction_offset, min_overhang=1):
    """
    Keep only the probes that straddle a junction.

    Args:
        probes (list): (seq_id, start, stop, probe_seq, tm) tuples mined from a
            junction window.
        junction_offset (int): the junction's offset within that window.
        min_overhang (int): bases the probe must place on each side.

    Returns:
        kept (list): the probes that span the junction.
    """
    # success
    return [p for p in probes
            if spans_junction(p[1], p[2], junction_offset, min_overhang)]


def discriminating_regions(gtf_df, gene_id=None, min_length=1):
    """
    Return the genomic intervals unique to each transcript of a gene.

    An interval is discriminating when it is exonic in one transcript and in no
    other transcript of the same gene. A transcript whose exons are all shared
    has no discriminating region and is reported with none, which is the
    honest answer to an isoform-resolved request for such a transcript.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF carrying exon records with
            gene_id and transcript_id.
        gene_id (str, optional): restrict to one gene. All genes when None.
        min_length (int): discard discriminating intervals shorter than this.

    Returns:
        regions (pandas.DataFrame): one row per discriminating interval, with
            DISCRIMINATING_COLUMNS.
    """
    exons = gtf_df
    if 'feature' in exons.columns:
        exons = exons[exons['feature'] == 'exon']
    if gene_id is not None:
        exons = exons[exons['gene_id'] == gene_id]

    rows = []
    for gene, gene_exons in exons.groupby('gene_id', sort=False):
        by_transcript = {
            transcript: _interval_set(block)
            for transcript, block in gene_exons.groupby('transcript_id', sort=False)
        }
        if len(by_transcript) < 2:
            continue

        for transcript, intervals in by_transcript.items():
            others = set()
            for sibling, sibling_intervals in by_transcript.items():
                if sibling != transcript:
                    others |= sibling_intervals

            unique = intervals - others
            if not unique:
                continue

            seqid = gene_exons['seqid'].iloc[0]
            strand = gene_exons['strand'].iloc[0]
            for start, end in _merge_positions(unique):
                length = end - start + 1
                if length >= min_length:
                    rows.append({
                        'gene_id': gene, 'transcript_id': transcript,
                        'seqid': seqid, 'start': start, 'end': end,
                        'strand': strand, 'length': length,
                    })

    # success
    return pd.DataFrame(rows, columns=DISCRIMINATING_COLUMNS)


def _interval_set(block):
    """
    Return every genomic position covered by a set of exons.

    Args:
        block (pandas.DataFrame): exon records with start and end.

    Returns:
        positions (set): the covered genomic positions.
    """
    positions = set()
    for start, end in zip(block['start'], block['end']):
        positions.update(range(int(start), int(end) + 1))

    # success
    return positions


def _merge_positions(positions):
    """
    Collapse a set of positions into contiguous intervals.

    Args:
        positions (set): genomic positions.

    Returns:
        intervals (list): (start, end) tuples, inclusive, sorted ascending.
    """
    intervals = []
    ordered = sorted(positions)
    if not ordered:
        return intervals

    start = previous = ordered[0]
    for position in ordered[1:]:
        if position == previous + 1:
            previous = position
        else:
            intervals.append((start, previous))
            start = previous = position
    intervals.append((start, previous))

    # success
    return intervals
