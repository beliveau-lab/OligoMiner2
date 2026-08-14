"""Tests for the model registry, loading and scoring.

Covers the registry contract, the loading traps that would otherwise return
plausible wrong numbers, that the vectorized encoder is bit-identical to the
reference one, and that each model responds to condition the way its registry
entry declares.
"""

import random

import numpy as np
import pandas as pd
import pytest

from oligominer import models
from oligominer.specificity.duplex_stability import features_fast
from oligominer.specificity.duplex_stability.frames import (
    build_aln,
    build_duplex_frame,
    expand_cigar,
)
from oligominer.specificity.duplex_stability.l4t import features as l4t_features


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

    # success
    return ''.join(out)


@pytest.fixture
def duplex_df():
    """Thirty 36-mer duplexes carrying zero to three mismatches each."""
    random.seed(1)
    rows = []
    for i in range(30):
        probe = ''.join(random.choice('ACGT') for _ in range(36))
        target = list(probe)
        for _ in range(random.randint(0, 3)):
            j = random.randrange(36)
            target[j] = random.choice([b for b in 'ACGT' if b != probe[j]])
        target = ''.join(target)
        rows.append({'seqid': f'p{i}', 'probe_seq': probe, 'derived_seq': target,
                     'align_cigar': cigar_from(probe, target), 'align_score': -5})

    # success
    return pd.DataFrame(rows)


class TestRegistry:

    def test_four_models_are_registered(self):
        assert models.available() == ['duplex-BiLSTM', 'om1-lda', 'physics-xgb', 'ps-xgb']

    def test_every_artifact_is_present(self):
        assert models.missing_artifacts() == []

    def test_om2_and_baseline_families_are_disjoint_and_complete(self):
        from oligominer.models.registry import BASELINE_MODELS, OM2_MODELS
        assert set(OM2_MODELS) == {'physics-xgb', 'duplex-BiLSTM'}
        assert set(BASELINE_MODELS) == {'ps-xgb', 'om1-lda'}
        assert set(OM2_MODELS) & set(BASELINE_MODELS) == set()

    def test_unknown_model_raises_and_names_the_valid_keys(self):
        with pytest.raises(KeyError, match='registered models'):
            models.spec('physics-xgboost')

    def test_default_model_is_registered(self):
        from oligominer.models.registry import DEFAULT_MODEL
        assert DEFAULT_MODEL in models.REGISTRY


class TestCigarAndAlignment:

    def test_expand_cigar_repeats_each_operation(self):
        assert expand_cigar('3=1X2=') == '===X=='

    def test_expand_cigar_handles_a_single_block(self):
        assert expand_cigar('36M') == 'M' * 36

    def test_build_aln_on_a_clean_match(self):
        probe_aln, target_aln = build_aln('ACGT', 'ACGT', '====')
        assert probe_aln == 'ACGT'
        assert target_aln == 'ACGT'

    def test_insertion_gaps_the_target(self):
        probe_aln, target_aln = build_aln('ACGT', 'ACT', '==I=')
        assert probe_aln == 'ACGT'
        assert target_aln == 'AC-T'

    def test_deletion_gaps_the_probe(self):
        probe_aln, target_aln = build_aln('ACT', 'ACGT', '==D=')
        assert probe_aln == 'AC-T'
        assert target_aln == 'ACGT'

    def test_a_short_sequence_returns_none_rather_than_truncating(self):
        assert build_aln('AC', 'ACGT', '====') == (None, None)


class TestFrameBuilding:

    def test_frame_carries_the_columns_the_encoders_read(self, duplex_df):
        frame = build_duplex_frame(duplex_df)
        for column in ('probe_aln', 'target_aln', 'ops', 'label_celsius', 'label_sodium'):
            assert column in frame.columns

    def test_condition_is_recorded_per_row(self, duplex_df):
        frame = build_duplex_frame(duplex_df, celsius=69.5, sodium=0.39)
        assert (frame['label_celsius'] == 69.5).all()
        assert (frame['label_sodium'] == 0.39).all()

    def test_celsius_none_requires_the_frame_to_carry_it(self, duplex_df):
        with pytest.raises(KeyError, match='label_celsius'):
            build_duplex_frame(duplex_df, celsius=None)

    def test_celsius_none_preserves_a_per_row_condition(self, duplex_df):
        duplex_df = duplex_df.copy()
        duplex_df['label_celsius'] = np.linspace(20, 80, len(duplex_df))
        frame = build_duplex_frame(duplex_df, celsius=None)
        assert frame['label_celsius'].nunique() == len(duplex_df)

    def test_malformed_rows_are_dropped_and_counted(self, duplex_df):
        broken = duplex_df.copy()
        broken.loc[0, 'align_cigar'] = '80='
        frame = build_duplex_frame(broken)
        assert frame.attrs['n_dropped_malformed'] == 1
        assert len(frame) == len(duplex_df) - 1

    def test_na_alignment_scores_become_zero(self, duplex_df):
        duplex_df = duplex_df.copy()
        duplex_df['align_score'] = 'NA'
        frame = build_duplex_frame(duplex_df)
        assert (frame['align_score'] == 0.0).all()

    def test_m_only_cigars_raise_rather_than_scoring_as_zero(self, duplex_df):
        """Without =/X the physics finds no helix and every feature computes as zero."""
        duplex_df = duplex_df.copy()
        duplex_df['align_cigar'] = '36M'
        with pytest.raises(ValueError, match="'=' or 'X'"):
            build_duplex_frame(duplex_df)


class TestEncoderParity:

    def test_fast_encoder_is_bit_identical_to_the_reference(self, duplex_df):
        frame = build_duplex_frame(duplex_df, celsius=47.0, sodium=0.39)
        reference = l4t_features.enc_om2(frame)
        fast = features_fast.enc_om2_fast(frame)
        assert list(reference.columns) == list(fast.columns)
        assert np.array_equal(reference.values, fast.values)

    def test_the_encoding_is_103_features(self, duplex_df):
        frame = build_duplex_frame(duplex_df)
        assert features_fast.enc_om2_fast(frame).shape[1] == 103

    def test_paintshop_features_are_a_strict_subset(self, duplex_df):
        """physics-xgb contains the PaintSHOP features verbatim, ps_ prefixed."""
        frame = build_duplex_frame(duplex_df)
        ps = l4t_features.enc_paintshop37(frame)
        om2 = features_fast.enc_om2_fast(frame)
        assert {f'ps_{c}' for c in ps.columns} <= set(om2.columns)
        assert om2.shape[1] == ps.shape[1] + 66


class TestLoadingAndScoring:

    @pytest.mark.parametrize('name', ['physics-xgb', 'ps-xgb', 'om1-lda', 'duplex-BiLSTM'])
    def test_every_model_loads_and_scores(self, name, duplex_df):
        frame = build_duplex_frame(duplex_df)
        values = models.load(name).predict(frame)
        assert len(values) == len(frame)
        assert np.isfinite(values).all()

    @pytest.mark.parametrize('name', ['physics-xgb', 'ps-xgb', 'duplex-BiLSTM'])
    def test_pdup_models_return_probabilities(self, name, duplex_df):
        frame = build_duplex_frame(duplex_df)
        values = models.load(name).predict(frame)
        assert values.min() >= 0.0
        assert values.max() <= 1.0

    def test_lda_returns_a_ranking_not_hard_class_labels(self, duplex_df):
        """decision_function is the ordering; predict would collapse to 0/1."""
        frame = build_duplex_frame(duplex_df)
        values = models.load('om1-lda').predict(frame)
        assert len(np.unique(values)) > 2
        assert not models.load('om1-lda').outputs_pdup

    def test_xgb_output_is_a_probability_not_a_logit(self, duplex_df):
        """A model fitted with a logit link needs the inverse sigmoid applied."""
        frame = build_duplex_frame(duplex_df, celsius=47.0)
        model = models.load('physics-xgb')
        raw = model.model.predict(
            features_fast.enc_om2_fast(frame).values.astype(np.float32)
        )
        assert raw.min() < 0.0 or raw.max() > 1.0
        assert model.predict(frame).max() <= 1.0

    def test_feature_names_come_from_the_encoder(self, duplex_df):
        frame = build_duplex_frame(duplex_df)
        names = models.load('physics-xgb').feature_names(frame)
        assert len(names) == 103

    def test_bilstm_reports_no_feature_names(self, duplex_df):
        frame = build_duplex_frame(duplex_df)
        assert models.load('duplex-BiLSTM').feature_names(frame) == []

    def test_load_all_returns_every_model(self, duplex_df):
        loaded = models.load_all()
        assert set(loaded) == set(models.available())

    def test_load_all_can_exclude_baselines(self):
        loaded = models.load_all(include_baselines=False)
        assert set(loaded) == {'physics-xgb', 'duplex-BiLSTM'}


class TestEmptyFrame:
    """Every probe on a chromosome can be dropped upstream, leaving no rows."""

    @pytest.fixture
    def empty_frame(self):
        empty = pd.DataFrame(
            columns=['probe_seq', 'derived_seq', 'align_cigar', 'align_score'])
        return build_duplex_frame(empty, celsius=69.5)

    @pytest.mark.parametrize('name', ['physics-xgb', 'ps-xgb', 'om1-lda'])
    def test_predicting_an_empty_frame_returns_nothing(self, name, empty_frame):
        values = models.load(name).predict(empty_frame)
        assert len(values) == 0

    def test_the_encoder_returns_the_full_width_on_an_empty_frame(self, empty_frame):
        assert features_fast.enc_om2_fast(empty_frame).shape == (0, 103)

    def test_an_empty_frame_does_not_trip_the_cigar_guard(self, empty_frame):
        assert len(empty_frame) == 0


class TestConditionResponse:
    """Each model must respond to temperature as its registry entry declares."""

    def _mean_by_temperature(self, name, duplex_df, temperatures):
        model = models.load(name)
        return [model.predict(build_duplex_frame(duplex_df, celsius=t, sodium=0.39)).mean()
                for t in temperatures]

    def test_physics_xgb_melts_as_temperature_rises(self, duplex_df):
        values = self._mean_by_temperature('physics-xgb', duplex_df, [17, 37, 57, 87])
        assert values == sorted(values, reverse=True)
        assert values[0] - values[-1] > 0.5

    def test_ps_xgb_is_condition_blind(self, duplex_df):
        values = self._mean_by_temperature('ps-xgb', duplex_df, [17, 47, 87])
        assert len(set(np.round(values, 10))) == 1

    def test_bilstm_is_condition_blind(self, duplex_df):
        values = self._mean_by_temperature('duplex-BiLSTM', duplex_df, [17, 47, 87])
        assert len(set(np.round(values, 10))) == 1

    def test_the_declared_condition_awareness_matches_behaviour(self, duplex_df):
        for name in models.available():
            declared = models.spec(name)['condition_aware']
            values = self._mean_by_temperature(name, duplex_df, [17, 87])
            moved = abs(values[0] - values[1]) > 1e-9
            assert moved == declared, f'{name} declares {declared} but moved={moved}'


class TestEmptyDuplexFrame:
    """A chromosome can lose every alignment upstream."""

    def test_an_empty_table_builds_an_empty_frame(self):
        empty = pd.DataFrame(columns=['probe_seq', 'derived_seq',
                                      'align_cigar'])
        out = build_duplex_frame(empty)

        assert out.empty
        assert 'probe_aln' in out.columns
        assert out.attrs['n_dropped_malformed'] == 0

    def test_an_empty_frame_does_not_trip_the_cigar_guard(self):
        # the guard raises when no row distinguishes matches, which no row can
        # do when there are no rows
        empty = pd.DataFrame(columns=['probe_seq', 'derived_seq',
                                      'align_cigar'])
        assert build_duplex_frame(empty).empty

    def test_an_empty_frame_predicts_nothing(self):
        empty = pd.DataFrame(columns=['probe_seq', 'derived_seq',
                                      'align_cigar'])
        from oligominer.models import load

        assert len(load('physics-xgb').predict(build_duplex_frame(empty))) == 0
