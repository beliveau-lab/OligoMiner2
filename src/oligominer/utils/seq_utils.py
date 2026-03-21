"""
Sequence utility functions.

General-purpose helpers for working with DNA sequences, including reverse
complement, GC content calculation, and numeric clamping.
"""

import numpy as np

COMPLEMENT = str.maketrans('ACGTacgt', 'TGCAtgca')


def rev_comp(seq):
    """
    Return the reverse complement of a DNA sequence.

    Args:
        seq (str): the input DNA sequence.

    Returns:
        rc (str): the reverse complement.
    """
    rc = seq.translate(COMPLEMENT)[::-1]

    # success
    return rc



def calc_gc(seq, as_percent=False):
    """
    Calculate GC content as a decimal (0.0 - 1.0) or percentage (0-100).
    
    Args:
        seq (str): the input DNA sequence.
        as_percent (bool): if True, return GC content as a percentage (0-100).

    Returns:
        gc_content (float): GC content as a decimal (0.0 - 1.0) or percentage (0-100).
    """
    # ensure uppercase for counting
    seq = seq.upper()
    
    # count G and C, then compute decimal GC content
    gc_count = seq.count('G') + seq.count('C')
    gc_content = (gc_count / len(seq)) if len(seq) > 0 else 0.0
    if as_percent:
        gc_content *= 100.0

    # success
    return gc_content


def clamp(values, lo, hi):
    """
    Clamp array or scalar values to [lo, hi] range.

    Args:
        values (numpy.ndarray or float): the values to clamp.
        lo (float): lower bound.
        hi (float): upper bound.

    Returns:
        clamped (numpy.ndarray or float): values clamped to [lo, hi].
    """
    clamped = np.clip(values, lo, hi)

    # success
    return clamped
