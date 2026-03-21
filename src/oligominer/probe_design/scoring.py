"""
# Probe Scoring

Functions for computing on-target and off-target scores from duplex
stability predictions.

The PaintSHOP pipeline scores probes by:

- **on_target_score**: the duplex prediction at the probe's intended
  genomic binding site.
- **off_target_score**: the sum of duplex predictions at all other
  alignment locations.

These scores help rank probes by specificity: ideal probes have high
on-target scores and low off-target scores.
"""

import pandas as pd


def label_on_target(merged_df):
    """
    Add an ``on_target`` boolean column to a merged duplex DataFrame.

    A row is on-target when the alignment location matches the probe's
    origin coordinates (same sequence ID and start position).

    Args:
        merged_df (pandas.DataFrame): merged probe-alignment table with
            columns ``seq_id``, ``start``, ``align_seqid``, and
            ``align_start``.

    Returns:
        merged_df (pandas.DataFrame): the input DataFrame with an added
            ``on_target`` column.
    """
    merged_df = merged_df.copy()
    merged_df['on_target'] = (
        (merged_df['seq_id'] == merged_df['align_seqid']) &
        (merged_df['start'] == merged_df['align_start'])
    )

    # success
    return merged_df


def score_probes(merged_df, pred_column="duplex_pred"):
    """
    Compute on-target and off-target scores for each probe.

    Groups the merged duplex table by probe ``seqid`` and aggregates
    the duplex predictions into per-probe scores.

    Args:
        merged_df (pandas.DataFrame): merged probe-alignment table with
            columns ``seqid``, ``on_target``, and *pred_column*. The
            ``on_target`` column should be boolean (see
            :func:`label_on_target`).
        pred_column (str): name of the column containing duplex
            predictions.

    Returns:
        scores_df (pandas.DataFrame): one row per probe with columns
            ``seqid``, ``on_target_score``, and ``off_target_score``.
    """
    # on-target: sum of predictions where on_target is True
    on_target = (
        merged_df[merged_df['on_target']]
        .groupby('seqid')[pred_column]
        .sum()
        .reset_index(name='on_target_score')
    )

    # off-target: sum of predictions where on_target is False
    off_target = (
        merged_df[~merged_df['on_target']]
        .groupby('seqid')[pred_column]
        .sum()
        .reset_index(name='off_target_score')
    )

    # combine, filling missing values with 0
    scores_df = on_target.set_index('seqid').join(
        off_target.set_index('seqid'), how='outer'
    ).fillna(0).reset_index()

    # success
    return scores_df
