"""
# Interval exclusion

Removes probes overlapping regions a user wants avoided, named by a BED file.

One stage serves every annotation of that shape: segmental duplications, common
SNPs, blacklists, centromeres, or any interval set a user supplies. Rather than a
separate feature per annotation type, the caller names the file.

Overlap is computed by sorting the exclusion intervals per chromosome and binary
searching, so a genome-scale probe table costs one sort per chromosome rather
than a scan per probe.
"""

import numpy as np
import pandas as pd

BED_COLUMNS = ['chrom', 'start', 'stop']


def read_bed(path):
    """
    Read the first three columns of a BED file.

    Lines beginning with 'track', 'browser' or '#' are skipped.

    Args:
        path (str): path to the BED file.

    Returns:
        intervals (pandas.DataFrame): chrom, start and stop columns.
    """
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith(('#', 'track', 'browser')):
                continue
            fields = line.split('\t')
            if len(fields) < 3:
                fields = line.split()
            if len(fields) < 3:
                continue
            rows.append((fields[0], int(fields[1]), int(fields[2])))

    # success
    return pd.DataFrame(rows, columns=BED_COLUMNS)


def _merge_intervals(starts, stops):
    """
    Merge overlapping intervals on one chromosome.

    Args:
        starts (numpy.ndarray): interval starts.
        stops (numpy.ndarray): interval stops.

    Returns:
        merged_starts (numpy.ndarray): the merged starts, sorted ascending.
        merged_stops (numpy.ndarray): the merged stops.
    """
    order = np.argsort(starts, kind='stable')
    starts, stops = starts[order], stops[order]

    merged_starts, merged_stops = [], []
    for start, stop in zip(starts, stops):
        if merged_stops and start <= merged_stops[-1]:
            merged_stops[-1] = max(merged_stops[-1], stop)
        else:
            merged_starts.append(start)
            merged_stops.append(stop)

    # success
    return np.array(merged_starts, dtype=np.int64), np.array(merged_stops, dtype=np.int64)


def overlaps_intervals(probe_df, intervals, chrom_col='seq_id',
                       start_col='start', stop_col='stop'):
    """
    Return which probes overlap any of the given intervals.

    Args:
        probe_df (pandas.DataFrame): probes with chromosome, start and stop.
        intervals (pandas.DataFrame): exclusion intervals with BED_COLUMNS.
        chrom_col (str): probe column naming the chromosome.
        start_col (str): probe column holding the start coordinate.
        stop_col (str): probe column holding the stop coordinate.

    Returns:
        overlapping (numpy.ndarray): boolean, one entry per probe.
    """
    overlapping = np.zeros(len(probe_df), dtype=bool)
    if len(intervals) == 0 or len(probe_df) == 0:
        return overlapping

    probe_starts = probe_df[start_col].to_numpy(dtype=np.int64)
    probe_stops = probe_df[stop_col].to_numpy(dtype=np.int64)

    for chrom, block in intervals.groupby('chrom', sort=False):
        on_chrom = (probe_df[chrom_col].to_numpy() == chrom)
        if not on_chrom.any():
            continue

        starts, stops = _merge_intervals(
            block['start'].to_numpy(dtype=np.int64),
            block['stop'].to_numpy(dtype=np.int64))

        positions = np.flatnonzero(on_chrom)
        # the interval that could contain a probe's start is the last one
        # beginning at or before it, so one binary search decides each probe
        candidate = np.searchsorted(starts, probe_starts[positions], side='right') - 1

        hit = np.zeros(len(positions), dtype=bool)
        valid = candidate >= 0
        hit[valid] = probe_starts[positions][valid] < stops[candidate[valid]]

        # a probe can also start before an interval and run into it
        following = candidate + 1
        in_range = following < len(starts)
        hit[in_range] |= probe_stops[positions][in_range] > starts[following[in_range]]

        overlapping[positions] = hit

    # success
    return overlapping


def exclude_intervals(probe_df, bed_path, chrom_col='seq_id',
                      start_col='start', stop_col='stop'):
    """
    Drop probes overlapping any interval in a BED file.

    Args:
        probe_df (pandas.DataFrame): probes with chromosome, start and stop.
        bed_path (str): path to the BED file naming regions to avoid.
        chrom_col (str): probe column naming the chromosome.
        start_col (str): probe column holding the start coordinate.
        stop_col (str): probe column holding the stop coordinate.

    Returns:
        kept (pandas.DataFrame): the probes that overlap nothing.
        n_dropped (int): how many probes were removed.
    """
    intervals = read_bed(bed_path)
    overlapping = overlaps_intervals(
        probe_df, intervals,
        chrom_col=chrom_col, start_col=start_col, stop_col=stop_col)

    kept = probe_df[~overlapping].copy()

    # success
    return kept, int(overlapping.sum())
