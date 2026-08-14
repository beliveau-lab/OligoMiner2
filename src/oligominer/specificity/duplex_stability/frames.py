"""
# Duplex frames

Shapes a merged probe and alignment table into the frame the duplex-stability
encoders read.

The encoders consume aligned columns rather than raw sequences. Each row needs:

- ``target_seq``: the reference sequence at the alignment site
- ``ops``: the CIGAR expanded to one character per alignment column
- ``probe_aln`` and ``target_aln``: the gapped alignment strings
- ``label_celsius`` and ``label_sodium``: the condition the physics features are
  computed at

The physics model recomputes its features at each row's own temperature and salt,
so those two must be per-row columns. Passing ``celsius=None`` keeps whatever the
frame already carries, which is what a multi-temperature corpus requires.
"""

import numpy as np
import pandas as pd

# alignment operations that consume a base from both the probe and the target.
# S is included because the extracted target window spans the full probe
# footprint, so a soft-clipped flank lines up positionally rather than gapping
_CONSUMES_BOTH = '=XMS'

# consumes the probe only, leaving a gap in the target
_CONSUMES_PROBE = 'I'

# consumes the target only, leaving a gap in the probe
_CONSUMES_TARGET = 'DN'

DEFAULT_CELSIUS = 47.0
DEFAULT_SODIUM = 0.39


def expand_cigar(cigar):
    """
    Expand a CIGAR string into one operation character per alignment column.

    '36M' becomes 36 'M' characters.

    Args:
        cigar (str): a SAM CIGAR string.

    Returns:
        ops (str): one operation character per alignment column.
    """
    ops = []
    count = ''
    for char in str(cigar):
        if char.isdigit():
            count += char
        else:
            ops.append(char * int(count or 1))
            count = ''

    # success
    return ''.join(ops)


def build_aln(probe, target, ops):
    """
    Build the gapped alignment strings for one duplex.

    Args:
        probe (str): the probe sequence.
        target (str): the reference sequence spanning the probe's footprint.
        ops (str): the expanded CIGAR, one character per column.

    Returns:
        probe_aln (str or None): the gapped probe, or None if either sequence ran
            short against the CIGAR.
        target_aln (str or None): the gapped target, or None in the same case.
    """
    probe_out, target_out = [], []
    probe_i = target_i = 0

    for op in ops:
        if op in _CONSUMES_BOTH:
            if probe_i >= len(probe) or target_i >= len(target):
                return None, None
            probe_out.append(probe[probe_i])
            target_out.append(target[target_i])
            probe_i += 1
            target_i += 1
        elif op == _CONSUMES_PROBE:
            if probe_i >= len(probe):
                return None, None
            probe_out.append(probe[probe_i])
            target_out.append('-')
            probe_i += 1
        elif op in _CONSUMES_TARGET:
            if target_i >= len(target):
                return None, None
            probe_out.append('-')
            target_out.append(target[target_i])
            target_i += 1

    # success
    return ''.join(probe_out), ''.join(target_out)


def build_duplex_frame(merged_df, celsius=DEFAULT_CELSIUS, sodium=DEFAULT_SODIUM,
                       drop_malformed=True):
    """
    Shape a merged probe and alignment table into an encoder-ready frame.

    Rows whose sequences run short against their CIGAR cannot be aligned and are
    dropped by default, because a truncated alignment scores without error. The
    number dropped is recorded in the returned frame's attrs under
    'n_dropped_malformed'.

    Args:
        merged_df (pandas.DataFrame): merged probe and alignment table carrying
            probe_seq, derived_seq and align_cigar.
        celsius (float, optional): temperature the physics features are computed
            at. None keeps each row's existing label_celsius, which a
            multi-temperature corpus requires.
        sodium (float, optional): sodium molarity. None keeps each row's existing
            label_sodium.
        drop_malformed (bool): drop rows whose alignment could not be built.

    Returns:
        df (pandas.DataFrame): the encoder-ready frame.

    Raises:
        KeyError: if celsius is None and the frame carries no label_celsius.
    """
    df = merged_df.copy()

    if df.empty:
        for column in ('target_seq', 'ops', 'probe_aln', 'target_aln'):
            df[column] = pd.Series(dtype=object)
        df['align_score'] = pd.Series(dtype=float)
        df['length'] = pd.Series(dtype=int)
        df.attrs['n_dropped_malformed'] = 0
        return df

    if 'target_seq' not in df:
        df['target_seq'] = df['derived_seq']
    if 'ops' not in df:
        df['ops'] = df['align_cigar'].map(expand_cigar)

    probe_aln, target_aln = [], []
    for probe, target, ops in zip(df['probe_seq'], df['target_seq'], df['ops']):
        aligned_probe, aligned_target = build_aln(str(probe), str(target), str(ops))
        probe_aln.append(aligned_probe)
        target_aln.append(aligned_target)

    df['probe_aln'] = probe_aln
    df['target_aln'] = target_aln

    n_before = len(df)
    malformed = df['probe_aln'].isna() | df['target_aln'].isna()
    if drop_malformed and malformed.any():
        df = df[~malformed].copy()
    df.attrs['n_dropped_malformed'] = int(n_before - len(df))

    # alignment scores arrive from BED as strings and can be the literal "NA"
    if 'align_score' in df:
        df['align_score'] = pd.to_numeric(df['align_score'], errors='coerce').fillna(0.0)
    else:
        df['align_score'] = 0.0

    df = _apply_condition(df, celsius, sodium)

    if 'length' not in df:
        df['length'] = df['probe_seq'].str.len()

    _check_ops_distinguish_matches(df)

    # success
    return df


def _check_ops_distinguish_matches(df):
    """
    Verify the CIGARs distinguish matches from mismatches.

    The thermodynamic features are built from runs of '=' operations. A CIGAR
    written with 'M', which bowtie2 emits without --xeq, marks matches and
    mismatches alike, so no helix is found and every physics feature computes as
    zero without raising.

    Args:
        df (pandas.DataFrame): the frame being built, carrying an ops column.

    Returns:
        ok (bool): True when at least one row distinguishes matches.

    Raises:
        ValueError: if no row's operations contain '=' or 'X'.
    """
    if len(df) == 0:
        return True

    ops = df['ops'].astype(str)
    if not ops.str.contains('[=X]', regex=True).any():
        raise ValueError(
            "no alignment in this frame uses '=' or 'X' operations, so matches "
            "and mismatches cannot be told apart and every thermodynamic feature "
            "would compute as zero. Align with xeq=True, which is the default for "
            "oligominer's bowtie2 wrapper."
        )

    # success
    return True


def _apply_condition(df, celsius, sodium):
    """
    Set the per-row condition columns the physics encoder reads.

    Args:
        df (pandas.DataFrame): the frame being built.
        celsius (float or None): temperature, or None to keep the frame's own.
        sodium (float or None): sodium molarity, or None to keep the frame's own.

    Returns:
        df (pandas.DataFrame): the frame with label_celsius and label_sodium set.

    Raises:
        KeyError: if celsius is None and no label_celsius column exists.
    """
    if celsius is not None:
        if 'label_celsius' in df.columns:
            df['label_celsius'] = df['label_celsius'].fillna(float(celsius))
        else:
            df['label_celsius'] = float(celsius)
    elif 'label_celsius' not in df.columns:
        raise KeyError(
            'celsius=None requires the frame to carry a label_celsius column '
            'per row, and this one does not')

    if sodium is not None:
        if 'label_sodium' in df.columns:
            df['label_sodium'] = df['label_sodium'].fillna(float(sodium))
        else:
            df['label_sodium'] = float(sodium)

    # success
    return df
