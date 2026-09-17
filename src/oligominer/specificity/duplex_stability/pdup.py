"""
# Two-stage pDup prediction

Screens every alignment with a fast model, then computes exact NUPACK pDup for
only the ones that could matter.

This is the duplex-stability stage of specificity analysis: alignment proposes
candidate duplexes, duplex reconstruction builds the frame, and this module puts
a binding probability on every row of it.

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

The returned table carries the model score and the exact value in separate columns,
plus a column saying which rows were verified. Writing an estimate and a
measurement into one column makes them indistinguishable downstream, and the rows
that were never verified are exactly the ones a reader would most want to know
about.
"""

import numpy as np
import pandas as pd

from oligominer.models.registry import DEFAULT_MODEL

# rows scoring at or above this are sent to the physics
DEFAULT_VERIFY_ABOVE = 0.05

# column names this stage writes
MODEL_COLUMN = 'pdup_model'
EXACT_COLUMN = 'pdup_exact'
SOURCE_COLUMN = 'pdup_source'
FINAL_COLUMN = 'pdup'

# where the per-run record is stashed on the returned frame
ATTRS_KEY = 'pdup_prediction'


def predict_pdup(frame, model=DEFAULT_MODEL, verify_above=DEFAULT_VERIFY_ABOVE,
                 verify=True, nupack_model=None, max_verify=None):
    """
    Score duplexes with a model and verify the credible ones with NUPACK.

    Args:
        frame (pandas.DataFrame): a duplex frame carrying the aligned columns,
            as built by duplex_stability.frames.build_duplex_frame.
        model (str): the screening model, from oligominer.models.
        verify_above (float): rows scoring at or above this are sent to the
            physics. This is a verification cutoff, not a binder threshold.
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
            score cannot be compared against a cutoff in pDup units.
    """
    from oligominer.models import load

    loaded = load(model)
    if not loaded.outputs_pdup:
        raise ValueError(
            f'{model} emits a decision score rather than pDup, so it cannot be '
            f'compared against verify_above={verify_above} in pDup units. Use a '
            f'model whose outputs_pdup is True.')

    out = frame.copy()
    out[MODEL_COLUMN] = loaded.predict(frame)
    out[EXACT_COLUMN] = np.nan

    verifiable = _verifiable(out)
    selected = _select_for_verification(out[MODEL_COLUMN], verify_above, max_verify)
    selected = selected & verifiable

    if verify and selected.any():
        out.loc[selected, EXACT_COLUMN] = _exact_pdup(
            out.loc[selected], nupack_model=nupack_model)

    verified = out[EXACT_COLUMN].notna()
    out[SOURCE_COLUMN] = np.where(verified, 'nupack', model)
    out[FINAL_COLUMN] = np.where(verified, out[EXACT_COLUMN], out[MODEL_COLUMN])

    out.attrs[ATTRS_KEY] = {
        'model': model,
        'verify_above': verify_above,
        'n_rows': int(len(out)),
        'n_selected': int(selected.sum()),
        'n_verified': int(verified.sum()),
        'n_unverifiable': int((~verifiable).sum()),
        'fraction_verified': float(verified.mean()) if len(out) else 0.0,
    }

    # success
    return out


def _verifiable(frame):
    """
    Report which rows the physics can accept.

    A target read from a genome can contain N where the assembly has a gap, and
    NUPACK rejects any character outside its alphabet. Such a row keeps its model
    score and is reported as unverified, rather than failing the run.

    Args:
        frame (pandas.DataFrame): rows carrying probe_seq and derived_seq.

    Returns:
        verifiable (pandas.Series): boolean, True where both sequences are ACGT.
    """
    pure = r'^[ACGTacgt]+$'
    verifiable = (frame['probe_seq'].astype(str).str.match(pure)
                  & frame['derived_seq'].astype(str).str.match(pure))

    # success
    return verifiable.fillna(False)


def _select_for_verification(scores, verify_above, max_verify):
    """
    Choose which rows go to the physics.

    Args:
        scores (pandas.Series): the model's pDup estimates.
        verify_above (float): minimum score to be verified.
        max_verify (int or None): cap on the number verified.

    Returns:
        selected (pandas.Series): boolean, True for rows to verify.
    """
    selected = scores >= verify_above

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


def pdup_summary(out):
    """
    Return how much physics the two-stage prediction avoided.

    Args:
        out (pandas.DataFrame): a frame returned by predict_pdup().

    Returns:
        summary (dict): the prediction record, plus the speedup implied by the
            fraction of rows that reached the physics.
    """
    summary = dict(out.attrs.get(ATTRS_KEY, {}))
    fraction = summary.get('fraction_verified')
    if fraction:
        summary['physics_calls_avoided'] = summary['n_rows'] - summary['n_verified']
        summary['implied_speedup'] = 1.0 / fraction

    # success
    return summary
