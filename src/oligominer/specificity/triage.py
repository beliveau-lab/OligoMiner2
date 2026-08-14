"""
# Two-stage duplex triage

Screens every alignment with a fast model, then computes exact NUPACK pDup for
only the ones that could matter.

Exact pDup on every alignment of a genome-wide run is the accurate answer and an
expensive one. Most alignments out of a permissive aligner are not credible
duplexes at the hybridization condition, and a model rejects those far more
cheaply than the physics does. Running the physics only on the survivors keeps
the exact answer where it changes a decision.

## Two strand conventions meet here

A duplex frame stores ``derived_seq`` on the same strand as the probe, so an
on-target row has ``derived_seq == probe_seq``. That is the frame the models were
fitted in: a match is a column where the two agree.

NUPACK is given the molecules that actually hybridize, so it needs the reverse
complement of the target. Passing the same-strand sequence instead asks what a
probe does against a copy of itself: measured on a real probe at its own locus,
that returns pDup 0.000276 where the duplex is 0.999948.

## The two columns are never merged

A triaged table carries the model score and the exact value in separate columns,
plus a column saying which rows were verified. Writing an estimate and a
measurement into one column makes them indistinguishable downstream, and the rows
that were never verified are exactly the ones a reader would most want to know
about.
"""

import numpy as np
import pandas as pd

from oligominer.models import load
from oligominer.models.registry import DEFAULT_MODEL

# rows scoring at or above this are sent to the physics
DEFAULT_THRESHOLD = 0.05

# column names the triage writes
MODEL_COLUMN = 'pdup_model'
EXACT_COLUMN = 'pdup_exact'
SOURCE_COLUMN = 'pdup_source'
FINAL_COLUMN = 'pdup'


def triage(frame, model_name=DEFAULT_MODEL, threshold=DEFAULT_THRESHOLD,
           verify=True, nupack_model=None, max_verify=None):
    """
    Score duplexes with a model and verify the credible ones with NUPACK.

    Args:
        frame (pandas.DataFrame): a duplex frame carrying the aligned columns,
            as built by duplex_stability.frames.build_duplex_frame.
        model_name (str): the screening model, from oligominer.models.
        threshold (float): rows scoring at or above this are verified.
        verify (bool): run the physics on the survivors. When False the model
            score is reported alone, which is what a run without NUPACK does.
        nupack_model (nupack.Model, optional): the thermodynamic model. Defaults
            to the standard FISH condition.
        max_verify (int, optional): cap on rows sent to the physics. The
            highest-scoring rows are verified when the cap binds.

    Returns:
        out (pandas.DataFrame): a copy carrying pdup_model, pdup_exact,
            pdup_source and pdup. pdup_exact is NaN where nothing was verified,
            and pdup_source names the origin of each pdup value.

    Raises:
        ValueError: if the screening model does not emit pDup, since a decision
            score has no threshold in pDup units.
    """
    model = load(model_name)
    if not model.outputs_pdup:
        raise ValueError(
            f'{model_name} emits a decision score rather than pDup, so it cannot '
            f'be thresholded at {threshold} in pDup units. Use a model whose '
            f'outputs_pdup is True.')

    out = frame.copy()
    out[MODEL_COLUMN] = model.predict(frame)
    out[EXACT_COLUMN] = np.nan

    selected = _select_for_verification(out[MODEL_COLUMN], threshold, max_verify)

    if verify and selected.any():
        out.loc[selected, EXACT_COLUMN] = _exact_pdup(
            out.loc[selected], nupack_model=nupack_model)

    verified = out[EXACT_COLUMN].notna()
    out[SOURCE_COLUMN] = np.where(verified, 'nupack', model_name)
    out[FINAL_COLUMN] = np.where(verified, out[EXACT_COLUMN], out[MODEL_COLUMN])

    out.attrs['triage'] = {
        'model': model_name,
        'threshold': threshold,
        'n_rows': int(len(out)),
        'n_selected': int(selected.sum()),
        'n_verified': int(verified.sum()),
        'fraction_verified': float(verified.mean()) if len(out) else 0.0,
    }

    # success
    return out


def _select_for_verification(scores, threshold, max_verify):
    """
    Choose which rows go to the physics.

    Args:
        scores (pandas.Series): the model's pDup estimates.
        threshold (float): minimum score to be verified.
        max_verify (int or None): cap on the number verified.

    Returns:
        selected (pandas.Series): boolean, True for rows to verify.
    """
    selected = scores >= threshold

    if max_verify is not None and selected.sum() > max_verify:
        # keep the highest-scoring rows, which are the ones whose exact value is
        # most likely to change a decision
        cutoff = scores[selected].nlargest(max_verify).min()
        selected = selected & (scores >= cutoff)

    # success
    return selected


def _exact_pdup(frame, nupack_model=None):
    """
    Compute exact pDup for a frame, grouped by probe.

    derived_seq is stored on the same strand as the probe, so it is reverse
    complemented here to give NUPACK the strand that actually hybridizes.

    Args:
        frame (pandas.DataFrame): rows carrying probe_seq and derived_seq.
        nupack_model (nupack.Model, optional): the thermodynamic model.

    Returns:
        values (numpy.ndarray): exact pDup per row.
    """
    from oligominer.thermodynamics.nupack import calc_pdup_one_to_many
    from oligominer.utils.seq_utils import rev_comp

    values = np.zeros(len(frame), dtype=np.float64)
    positions = pd.Series(np.arange(len(frame)), index=frame.index)

    for probe, block in frame.groupby('probe_seq', sort=False):
        rows = positions.loc[block.index].to_numpy()
        targets = [rev_comp(str(t).upper()) for t in block['derived_seq']]
        values[rows] = calc_pdup_one_to_many(str(probe), targets,
                                             model=nupack_model)

    # success
    return values


def triage_summary(out):
    """
    Return how much work the triage avoided.

    Args:
        out (pandas.DataFrame): a frame returned by triage().

    Returns:
        summary (dict): the triage record, plus the speedup implied by the
            fraction of rows that reached the physics.
    """
    summary = dict(out.attrs.get('triage', {}))
    fraction = summary.get('fraction_verified')
    if fraction:
        summary['physics_calls_avoided'] = summary['n_rows'] - summary['n_verified']
        summary['implied_speedup'] = 1.0 / fraction

    # success
    return summary
