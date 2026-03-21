"""
# Probe CSV I/O

Read and write probe and alignment DataFrames as CSV files.
"""

import pandas as pd

from oligominer.utils import check_input_exists, check_dir_exists
from oligominer.thermodynamics.mining import PROBE_COLUMNS

# expected column order for alignment CSV files
ALIGN_COLUMNS = [
    'align_seqid', 'align_start', 'align_stop', 'seqid',
    'align_score', 'align_strand', 'align_cigar', 'derived_seq'
]


def read_probe_csv(path):
    """
    Read a probe CSV file into a DataFrame.

    Args:
        path (str): path to the probe CSV file.

    Returns:
        probe_df (pandas.DataFrame): the probe data.
    """
    check_input_exists(path)
    probe_df = pd.read_csv(path)

    # success
    return probe_df


def write_probe_csv(probe_df, path):
    """
    Write a probe DataFrame to a CSV file.

    Args:
        probe_df (pandas.DataFrame): the probe data.
        path (str): output file path.
    """
    check_dir_exists(path, parent_dir=True, create=True)
    probe_df.to_csv(path, index=False, float_format="%.2f")


def read_align_csv(path):
    """
    Read an alignment CSV file into a DataFrame.

    Args:
        path (str): path to the alignment CSV file.

    Returns:
        align_df (pandas.DataFrame): the alignment data.
    """
    check_input_exists(path)
    align_df = pd.read_csv(path)

    # success
    return align_df


def write_align_csv(align_df, path):
    """
    Write an alignment DataFrame to a CSV file.

    Args:
        align_df (pandas.DataFrame): the alignment data.
        path (str): output file path.
    """
    check_dir_exists(path, parent_dir=True, create=True)
    align_df.to_csv(path, index=False)
