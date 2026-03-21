"""
Generic sequence appending primitives.

Functions for appending synthetic sequences to probe homology regions
using various assignment schemes: same, unique, multiple, or custom.

Each core function returns a tuple of ``(result_df, entries)`` where
*entries* is a pandas Series mapping each probe row index to a
``"{id}_{seq}"`` tracking string.

Appending schemes
-----------------
- **same**: append the same sequence to every probe
- **unique**: one unique sequence per target (gene / genomic region)
- **multiple**: N sequences per target in round-robin across probes
- **custom**: one sequence per user-defined probe index range
"""

import pandas as pd

from oligominer.utils.seq_utils import rev_comp
from oligominer.utils.exceptions import InvalidInputError
from .config import LINKER


# ------------------------------------------------------------------
# internal helpers
# ------------------------------------------------------------------

def _maybe_rc(seq, rc):
    """Optionally reverse-complement a sequence."""
    if rc:
        return rev_comp(seq)
    return seq


def _parse_ranges(range_strings):
    """
    Convert a list of ``"start-stop"`` strings to (start, stop) tuples.

    Indices are 1-based and inclusive.

    Args:
        range_strings (list of str): e.g. ``["1-5", "6-10"]``.

    Returns:
        ranges (list of tuple): e.g. ``[(1, 5), (6, 10)]``.
    """
    ranges = []
    for s in range_strings:
        parts = s.strip().split("-")
        ranges.append((int(parts[0]), int(parts[1])))

    # success
    return ranges


def _join(left_seq, right_seq, linker):
    """Concatenate two sequences with a linker."""
    return f"{left_seq}{linker}{right_seq}"


# ------------------------------------------------------------------
# core appending functions
# ------------------------------------------------------------------

def append_same(probes, sequences, left=True, rc=False, linker=LINKER):
    """
    Append the same sequence to every probe.

    Args:
        probes (pandas.DataFrame): probe data with a ``sequence`` column.
        sequences (pandas.DataFrame): appending sequences with ``id``
            and ``seq`` columns. Only the first row is used.
        left (bool): if True, prepend to the 5' end. If False, append
            to the 3' end.
        rc (bool): if True, reverse-complement the sequence before
            appending.
        linker (str): separator inserted between the appended sequence
            and the probe.

    Returns:
        result (pandas.DataFrame): copy of *probes* with modified
            ``sequence`` column.
        entries (pandas.Series): tracking strings indexed like *result*.
    """
    result = probes.copy()

    seq = _maybe_rc(sequences['seq'].iloc[0], rc)
    seq_id = sequences['id'].iloc[0]
    entry = f"{seq_id}_{sequences['seq'].iloc[0]}"

    if left:
        result['sequence'] = result['sequence'].apply(
            lambda s: _join(seq, s, linker)
        )
    else:
        result['sequence'] = result['sequence'].apply(
            lambda s: _join(s, seq, linker)
        )

    entries = pd.Series(entry, index=result.index)

    # success
    return result, entries


def append_unique(probes, sequences, target_column, left=True, rc=False,
                  linker=LINKER):
    """
    Append one unique sequence per target.

    Each unique value in *target_column* is assigned a different
    sequence from *sequences*. All probes sharing the same target get
    the same appended sequence.

    Args:
        probes (pandas.DataFrame): probe data with ``sequence`` and
            *target_column* columns.
        sequences (pandas.DataFrame): appending sequences with ``id``
            and ``seq`` columns. Must have at least as many rows as
            unique targets.
        target_column (str): column name used to group probes by target.
        left (bool): if True, prepend to the 5' end.
        rc (bool): if True, reverse-complement sequences before
            appending.
        linker (str): separator inserted between the appended sequence
            and the probe.

    Returns:
        result (pandas.DataFrame): copy of *probes* with modified
            ``sequence`` column.
        entries (pandas.Series): tracking strings indexed like *result*.

    Raises:
        InvalidInputError: if there are fewer sequences than unique targets.
    """
    unique_targets = probes[target_column].unique()

    if len(sequences) < len(unique_targets):
        raise InvalidInputError(
            f"Not enough sequences ({len(sequences)}) for "
            f"{len(unique_targets)} unique targets."
        )

    result = probes.copy()
    entries = pd.Series("", index=result.index)

    for i, target in enumerate(unique_targets):
        mask = result[target_column] == target
        seq = _maybe_rc(sequences['seq'].iloc[i], rc)
        seq_id = sequences['id'].iloc[i]

        if left:
            result.loc[mask, 'sequence'] = result.loc[mask, 'sequence'].apply(
                lambda s: _join(seq, s, linker)
            )
        else:
            result.loc[mask, 'sequence'] = result.loc[mask, 'sequence'].apply(
                lambda s: _join(s, seq, linker)
            )

        entries.loc[mask] = f"{seq_id}_{sequences['seq'].iloc[i]}"

    # success
    return result, entries


def append_multiple(probes, sequences, n_per_target, target_column,
                    left=True, rc=False, linker=LINKER):
    """
    Append N sequences per target in round-robin fashion.

    Each target is assigned *n_per_target* consecutive sequences from
    *sequences*. Within a target's probes, sequences are cycled so
    that every probe gets one of the N assigned sequences.

    Args:
        probes (pandas.DataFrame): probe data with ``sequence`` and
            *target_column* columns.
        sequences (pandas.DataFrame): appending sequences with ``id``
            and ``seq`` columns.
        n_per_target (int): number of distinct sequences per target.
        target_column (str): column name used to group probes by target.
        left (bool): if True, prepend to the 5' end.
        rc (bool): if True, reverse-complement sequences before
            appending.
        linker (str): separator inserted between the appended sequence
            and the probe.

    Returns:
        result (pandas.DataFrame): copy of *probes* with modified
            ``sequence`` column.
        entries (pandas.Series): tracking strings indexed like *result*.

    Raises:
        InvalidInputError: if there are not enough sequences for all targets.
    """
    unique_targets = probes[target_column].unique()
    needed = len(unique_targets) * n_per_target

    if len(sequences) < needed:
        raise InvalidInputError(
            f"Not enough sequences ({len(sequences)}) for "
            f"{len(unique_targets)} targets x {n_per_target} per target "
            f"({needed} needed)."
        )

    result = probes.copy()
    entries = pd.Series("", index=result.index)

    seq_offset = 0

    for target in unique_targets:
        mask = result[target_column] == target
        target_indices = result.index[mask]

        for j, idx in enumerate(target_indices):
            # cycle through this target's N sequences
            seq_pos = seq_offset + (j % n_per_target)
            seq = _maybe_rc(sequences['seq'].iloc[seq_pos], rc)
            seq_id = sequences['id'].iloc[seq_pos]
            probe_seq = result.at[idx, 'sequence']

            if left:
                result.at[idx, 'sequence'] = _join(seq, probe_seq, linker)
            else:
                result.at[idx, 'sequence'] = _join(probe_seq, seq, linker)

            entries.at[idx] = (
                f"{seq_id}_{sequences['seq'].iloc[seq_pos]}"
            )

        seq_offset += n_per_target

    # success
    return result, entries


def append_custom(probes, sequences, ranges, left=True, rc=False,
                  linker=LINKER):
    """
    Append one unique sequence per custom probe index range.

    Each range in *ranges* receives a different sequence from
    *sequences*. Ranges are 1-based and inclusive.

    Args:
        probes (pandas.DataFrame): probe data with a ``sequence``
            column.
        sequences (pandas.DataFrame): appending sequences with ``id``
            and ``seq`` columns. Must have at least as many rows as
            ranges.
        ranges (list of str): range strings like ``["1-5", "6-10"]``.
        left (bool): if True, prepend to the 5' end.
        rc (bool): if True, reverse-complement sequences before
            appending.
        linker (str): separator inserted between the appended sequence
            and the probe.

    Returns:
        result (pandas.DataFrame): copy of *probes* with modified
            ``sequence`` column.
        entries (pandas.Series): tracking strings indexed like *result*.

    Raises:
        InvalidInputError: if there are fewer sequences than ranges.
    """
    parsed = _parse_ranges(ranges)

    if len(sequences) < len(parsed):
        raise InvalidInputError(
            f"Not enough sequences ({len(sequences)}) for "
            f"{len(parsed)} custom ranges."
        )

    result = probes.copy()
    entries = pd.Series("", index=result.index)

    for i, (start, stop) in enumerate(parsed):
        # convert 1-based inclusive to 0-based iloc slice
        idx_slice = slice(start - 1, stop)
        seq = _maybe_rc(sequences['seq'].iloc[i], rc)
        seq_id = sequences['id'].iloc[i]
        iloc_indices = result.index[idx_slice]

        if left:
            result.loc[iloc_indices, 'sequence'] = (
                result.loc[iloc_indices, 'sequence'].apply(
                    lambda s: _join(seq, s, linker)
                )
            )
        else:
            result.loc[iloc_indices, 'sequence'] = (
                result.loc[iloc_indices, 'sequence'].apply(
                    lambda s: _join(s, seq, linker)
                )
            )

        entries.loc[iloc_indices] = (
            f"{seq_id}_{sequences['seq'].iloc[i]}"
        )

    # success
    return result, entries


# ------------------------------------------------------------------
# dispatcher
# ------------------------------------------------------------------

def append_sequences(probes, sequences, scheme, target_column=None,
                     n_per_target=None, ranges=None, left=True, rc=False,
                     linker=LINKER):
    """
    Append sequences to probes using the specified scheme.

    This is a convenience dispatcher that routes to the appropriate
    appending function based on *scheme*.

    Args:
        probes (pandas.DataFrame): probe data with a ``sequence``
            column.
        sequences (pandas.DataFrame): appending sequences with ``id``
            and ``seq`` columns.
        scheme (str): one of ``"same"``, ``"unique"``, ``"multiple"``,
            or ``"custom"``.
        target_column (str, optional): column name for grouping probes
            by target. Required for ``"unique"`` and ``"multiple"``.
        n_per_target (int, optional): number of distinct sequences per
            target. Required for ``"multiple"``.
        ranges (list of str, optional): range strings for ``"custom"``.
        left (bool): if True, prepend to the 5' end.
        rc (bool): if True, reverse-complement sequences.
        linker (str): separator inserted between the appended sequence
            and the probe.

    Returns:
        result (pandas.DataFrame): probes with modified ``sequence``
            column.
        entries (pandas.Series): tracking strings.
    """
    if scheme == "same":
        result, entries = append_same(
            probes, sequences, left=left, rc=rc, linker=linker
        )
    elif scheme == "unique":
        result, entries = append_unique(
            probes, sequences, target_column, left=left, rc=rc, linker=linker
        )
    elif scheme == "multiple":
        result, entries = append_multiple(
            probes, sequences, n_per_target, target_column,
            left=left, rc=rc, linker=linker
        )
    elif scheme == "custom":
        result, entries = append_custom(
            probes, sequences, ranges, left=left, rc=rc, linker=linker
        )
    else:
        raise InvalidInputError(
            f"Unknown appending scheme: {scheme!r}. "
            f"Must be one of 'same', 'unique', 'multiple', 'custom'."
        )

    # success
    return result, entries


# ------------------------------------------------------------------
# appending table assembly
# ------------------------------------------------------------------

def build_appending_table(probes, entries_dict):
    """
    Assemble a table recording which sequences were appended to each probe.

    The table has one row per probe and one column per appending step.

    Args:
        probes (pandas.DataFrame): the original probe DataFrame (used
            for the index).
        entries_dict (dict): mapping of label strings to
            ``pandas.Series`` of tracking entries, as returned by the
            core appending functions.

    Returns:
        table (pandas.DataFrame): appending table with columns for each
            appending step.
    """
    table = pd.DataFrame(index=probes.index)
    for label, entries in entries_dict.items():
        table[label] = entries

    # success
    return table
