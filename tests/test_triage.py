"""Tests for two-stage duplex triage.

The screen decides which alignments reach the physics. The property that matters
is not the speedup but that an estimate and a measurement never end up
indistinguishable in the same column.
"""

import random

import numpy as np
import pandas as pd
import pytest

from oligominer.specificity.duplex_stability.frames import build_duplex_frame
from oligominer.specificity.triage import (
    EXACT_COLUMN,
    FINAL_COLUMN,
    MODEL_COLUMN,
    SOURCE_COLUMN,
    triage,
    triage_summary,
)
from oligominer.utils.seq_utils import rev_comp

nupack = pytest.importorskip('nupack')


def cigar_from(probe, target):
    """Build an xeq CIGAR from two equal-length sequences."""
    ops = ['=' if a == b else 'X' for a, b in zip(probe, target)]
    out, run, current = [], 0, ops[0]
    for op in ops:
        if op == current:
            run += 1
        else:
            out.append(f'{run}{current}')
            current, run = op, 1
    out.append(f'{run}{current}')
    return ''.join(out)


@pytest.fixture(scope='module')
def frame():
    """Duplexes spanning near-perfect to heavily mismatched."""
    random.seed(6)
    rows = []
    for n_mismatch in (0, 1, 2, 4, 8, 12):
        for _ in range(4):
            probe = ''.join(random.choice('ACGT') for _ in range(32))
            target = list(probe)
            for j in random.sample(range(32), n_mismatch):
                target[j] = random.choice([b for b in 'ACGT' if b != probe[j]])
            target = ''.join(target)
            rows.append({'probe_seq': probe, 'derived_seq': target,
                         'align_cigar': cigar_from(probe, target),
                         'align_score': -3 * n_mismatch})

    # success
    return build_duplex_frame(pd.DataFrame(rows), celsius=69.5, sodium=0.39)


class TestColumns:

    def test_the_model_score_is_kept_separately(self, frame):
        out = triage(frame, threshold=0.5)
        assert out[MODEL_COLUMN].notna().all()

    def test_unverified_rows_have_no_exact_value(self, frame):
        out = triage(frame, threshold=0.5)
        unverified = out[SOURCE_COLUMN] != 'nupack'
        assert out.loc[unverified, EXACT_COLUMN].isna().all()

    def test_the_source_column_names_the_origin_of_every_value(self, frame):
        out = triage(frame, threshold=0.5)
        assert set(out[SOURCE_COLUMN]) <= {'nupack', 'physics-xgb'}

    def test_the_final_column_takes_the_exact_value_where_it_exists(self, frame):
        out = triage(frame, threshold=0.5)
        verified = out[EXACT_COLUMN].notna()
        assert np.allclose(out.loc[verified, FINAL_COLUMN],
                           out.loc[verified, EXACT_COLUMN])

    def test_the_final_column_falls_back_to_the_model(self, frame):
        out = triage(frame, threshold=0.5)
        unverified = out[EXACT_COLUMN].isna()
        assert np.allclose(out.loc[unverified, FINAL_COLUMN],
                           out.loc[unverified, MODEL_COLUMN])


class TestSelection:

    def test_a_low_threshold_verifies_more_than_a_high_one(self, frame):
        low = triage(frame, threshold=0.001).attrs['triage']['n_verified']
        high = triage(frame, threshold=0.9).attrs['triage']['n_verified']
        assert low >= high

    def test_a_threshold_above_every_score_verifies_nothing(self, frame):
        out = triage(frame, threshold=1.1)
        assert out[EXACT_COLUMN].isna().all()
        assert (out[SOURCE_COLUMN] != 'nupack').all()

    def test_verification_can_be_capped(self, frame):
        out = triage(frame, threshold=0.0, max_verify=3)
        assert out.attrs['triage']['n_verified'] <= 3

    def test_the_cap_keeps_the_highest_scoring_rows(self, frame):
        out = triage(frame, threshold=0.0, max_verify=3)
        verified = out[out[EXACT_COLUMN].notna()]
        unverified = out[out[EXACT_COLUMN].isna()]
        if len(verified) and len(unverified):
            assert verified[MODEL_COLUMN].min() >= unverified[MODEL_COLUMN].max()

    def test_verification_can_be_skipped_entirely(self, frame):
        """A run without NUPACK reports the model score alone."""
        out = triage(frame, threshold=0.0, verify=False)
        assert out[EXACT_COLUMN].isna().all()
        assert np.allclose(out[FINAL_COLUMN], out[MODEL_COLUMN])


class TestStrandConvention:
    """derived_seq is same-strand; NUPACK needs the strand that hybridizes."""

    def test_a_perfect_on_target_duplex_scores_near_one(self):
        probe = 'CGACAATGCACGACAGAGGAAGCAGAACAGA'
        on_target = build_duplex_frame(
            pd.DataFrame([{'probe_seq': probe, 'derived_seq': probe,
                           'align_cigar': f'{len(probe)}=', 'align_score': 0}]),
            celsius=69.5, sodium=0.39)

        out = triage(on_target, threshold=0.0)
        assert out[EXACT_COLUMN].iloc[0] > 0.9


class TestExactness:

    def test_verified_values_match_calc_pdup(self, frame):
        from oligominer.thermodynamics.nupack import calc_pdup

        out = triage(frame, threshold=0.0)
        verified = out[out[EXACT_COLUMN].notna()].head(6)
        for _, row in verified.iterrows():
            # derived_seq is same-strand, so the duplex is against its complement
            reference = calc_pdup(str(row['probe_seq']),
                                  rev_comp(str(row['derived_seq'])),
                                  conc_a=1e-6, conc_b=1e-12)
            assert abs(row[EXACT_COLUMN] - reference) < 1e-9

    def test_rows_stay_with_their_own_duplex_when_grouped_by_probe(self, frame):
        """Grouping by probe must not permute exact values onto other rows."""
        from oligominer.thermodynamics.nupack import calc_pdup

        duplicated = pd.concat([frame, frame.head(4)], ignore_index=True)
        duplicated = build_duplex_frame(
            duplicated[['probe_seq', 'derived_seq', 'align_cigar', 'align_score']],
            celsius=69.5, sodium=0.39)

        out = triage(duplicated, threshold=0.0)
        for _, row in out[out[EXACT_COLUMN].notna()].head(8).iterrows():
            reference = calc_pdup(str(row['probe_seq']),
                                  rev_comp(str(row['derived_seq'])),
                                  conc_a=1e-6, conc_b=1e-12)
            assert abs(row[EXACT_COLUMN] - reference) < 1e-9


class TestGuards:

    def test_a_decision_score_model_is_refused(self, frame):
        with pytest.raises(ValueError, match='decision score'):
            triage(frame, model_name='om1-lda', threshold=0.5)

    def test_the_summary_reports_what_was_avoided(self, frame):
        out = triage(frame, threshold=0.5)
        summary = triage_summary(out)
        assert summary['n_rows'] == len(frame)
        assert summary['n_verified'] <= summary['n_rows']
        if summary.get('fraction_verified'):
            assert summary['physics_calls_avoided'] == (
                summary['n_rows'] - summary['n_verified'])
