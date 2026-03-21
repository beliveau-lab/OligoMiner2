"""
# BED I/O

Utilities for reading and converting BED format data into pandas DataFrames
for downstream analysis.
"""

from io import StringIO

import pandas as pd

def bed_to_df(bed_data):
    """
    Convert BED data to a pandas DataFrame.

    Args:
        bed_data (str): BED data as a string.

    Returns:
        df (pandas.DataFrame): the converted DataFrame.
    """
    df = pd.read_csv(StringIO(bed_data), sep='\t', header=None)

    # success
    return df
