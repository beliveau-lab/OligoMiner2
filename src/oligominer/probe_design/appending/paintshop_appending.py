"""
PaintSHOP-specific probe appending functions.

Ports the PaintSHOP Shiny app appending interface (R) to Python.
The probe anatomy after full PaintSHOP appending is:

    5'--[Outer]--[Inner]--[Homology]--[Inner]--[Outer]--3'

with optional bridge sequences and SABER concatemers.

Includes SABER concatemer appending and MERFISH barcode-encoded bridge
appending. These functions build on the generic appending primitives in
``appending.py``.
"""

import numpy as np
import pandas as pd

from oligominer.utils.exceptions import InvalidInputError
from .appending import append_sequences


# ------------------------------------------------------------------
# SABER handler
# ------------------------------------------------------------------

def append_saber(probes, sequences, scheme, target_column=None,
                 n_per_target=None, ranges=None):
    """
    Append SABER concatemer sequences to probes.

    SABER sequences are always appended to the 3' end in forward
    orientation.

    Args:
        probes (pandas.DataFrame): probe data with a ``sequence``
            column.
        sequences (pandas.DataFrame): SABER sequences with ``id``
            and ``seq`` columns.
        scheme (str): appending scheme (``"same"``, ``"unique"``,
            ``"multiple"``, or ``"custom"``).
        target_column (str, optional): column for grouping by target.
        n_per_target (int, optional): sequences per target for
            ``"multiple"`` scheme.
        ranges (list of str, optional): ranges for ``"custom"`` scheme.

    Returns:
        result (pandas.DataFrame): probes with SABER sequences
            appended to the 3' end.
        entries (pandas.Series): tracking strings.
    """
    result, entries = append_sequences(
        probes, sequences, scheme,
        target_column=target_column,
        n_per_target=n_per_target,
        ranges=ranges,
        left=False, rc=False,
    )

    # success
    return result, entries


# ------------------------------------------------------------------
# MERFISH barcode appending
# ------------------------------------------------------------------

def collect_indices(barcode):
    """
    Convert a 16-bit MHD4 barcode string to bridge indices.

    Args:
        barcode (str): a 16-character binary string, e.g.
            ``"1000000010001001"``.

    Returns:
        indices (list of int): 1-based positions where the barcode
            character is ``"1"``. For example, ``[1, 9, 13, 16]``.
    """
    indices = [i + 1 for i, c in enumerate(barcode) if c == "1"]

    # success
    return indices


def add_bridges(probes, bridges, indices, seed=None):
    """
    Append MERFISH bridges to probes for a single target.

    Each probe gets 3 of the 4 bridges indicated by *indices*. One
    bridge is dropped randomly per probe. If the dropped bridge is in
    the first half (index 1 or 2), the 5' side gets 2 bridges and the
    3' side gets 1. Otherwise the 3' side gets 2 and the 5' side
    gets 1. Bridge sequences are directly concatenated (no TTT linker).

    Args:
        probes (pandas.DataFrame): probes for a single target, with a
            ``sequence`` column.
        bridges (pandas.DataFrame): the 16 MERFISH bridges with ``id``
            and ``seq`` columns.
        indices (list of int): four 1-based bridge indices from the
            barcode.
        seed (int, optional): random seed for reproducibility.

    Returns:
        result (pandas.DataFrame): probes with bridges appended.
    """
    result = probes.copy()
    n_probes = len(result)

    rng = np.random.default_rng(seed)

    # choose which of the 4 bridges to drop for each probe (0-3)
    dropped = rng.integers(0, 4, size=n_probes)

    for i, idx in enumerate(result.index):
        # remove the dropped bridge
        curr_indices = [indices[j] for j in range(4) if j != dropped[i]]

        # shuffle the remaining 3 bridges
        probe_rng = np.random.default_rng(seed=i if seed is not None else None)
        curr_indices = list(probe_rng.permutation(curr_indices))

        # determine heavy side based on which bridge was dropped
        # if dropped bridge was index 0 or 1 (first half), 3' is heavy
        # if dropped bridge was index 2 or 3 (second half), 5' is heavy
        heavy_five = dropped[i] >= 2

        seq = result.at[idx, 'sequence']

        if heavy_five:
            # 5' gets 2 bridges, 3' gets 1
            b1 = bridges['seq'].iloc[curr_indices[0] - 1]
            b2 = bridges['seq'].iloc[curr_indices[1] - 1]
            b3 = bridges['seq'].iloc[curr_indices[2] - 1]
            seq = b2 + b1 + seq + b3
        else:
            # 5' gets 1 bridge, 3' gets 2
            b1 = bridges['seq'].iloc[curr_indices[0] - 1]
            b2 = bridges['seq'].iloc[curr_indices[1] - 1]
            b3 = bridges['seq'].iloc[curr_indices[2] - 1]
            seq = b1 + seq + b2 + b3

        result.at[idx, 'sequence'] = seq

    # success
    return result


def append_barcodes(probes, bridges, barcodes, target_column="refseq"):
    """
    Append MERFISH barcode-encoded bridges to a probe set.

    Each unique target is assigned a barcode from *barcodes*. The
    barcode determines which 4 of 16 bridges are used for that target.
    Each probe then gets 3 of those 4 bridges appended to its 5' and
    3' ends.

    Args:
        probes (pandas.DataFrame): probe data with ``sequence`` and
            *target_column* columns.
        bridges (pandas.DataFrame): the 16 MERFISH bridges with ``id``
            and ``seq`` columns.
        barcodes (pandas.DataFrame): barcodes with a ``barcode`` column
            containing 16-character binary strings. Must have at least
            as many rows as unique targets.
        target_column (str): column name used to group probes by target.

    Returns:
        result (pandas.DataFrame): probes with bridges appended.

    Raises:
        ValueError: if there are fewer barcodes than unique targets.
    """
    unique_targets = probes[target_column].unique()

    if len(barcodes) < len(unique_targets):
        raise InvalidInputError(
            f"Not enough barcodes ({len(barcodes)}) for "
            f"{len(unique_targets)} unique targets."
        )

    result_parts = []
    for i, target in enumerate(unique_targets):
        target_probes = probes[probes[target_column] == target].copy()
        barcode_str = barcodes['barcode'].iloc[i]
        indices = collect_indices(barcode_str)
        target_result = add_bridges(target_probes, bridges, indices, seed=i)
        result_parts.append(target_result)

    result = pd.concat(result_parts)

    # success
    return result
