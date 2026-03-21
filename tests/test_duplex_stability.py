"""Tests for duplex stability prediction methods.

Covers three approaches to estimating duplex formation probability:
  1. NUPACK pDup (direct thermodynamic simulation)
  2. PaintSHOP XGBoost model (fast approximation)
  3. Legacy OligoMiner v1 LDA model (simple classifier)
"""

import numpy as np
import pandas as pd
import pytest

from oligominer.utils.exceptions import ConfigurationError


# ---------------------------------------------------------------------------
# test sequences and fixtures
# ---------------------------------------------------------------------------

PROBE_SEQ = "ATCGATCGATCGATCGATCGATCG"
DERIVED_SEQ = "CGATCGATCGATCGATCGATCGAT"
ALIGN_SCORE = -10

@pytest.fixture
def single_pair_df():
    """DataFrame with a single probe-target pair for batch prediction."""
    return pd.DataFrame({
        'probe_seq': [PROBE_SEQ],
        'derived_seq': [DERIVED_SEQ],
        'align_score': [ALIGN_SCORE],
    })


@pytest.fixture
def multi_pair_df():
    """DataFrame with several probe-target pairs for batch prediction."""
    return pd.DataFrame({
        'probe_seq': [
            "ATCGATCGATCGATCGATCGATCG",
            "GGCCGGCCGGCCGGCCGGCCGGCC",
            "ATATATATATATATATATATATATAT",
            "GCGCGCGCGCGCGCGCGCGCGCGC",
            "AACCTTGGAACCTTGGAACCTTGG",
        ],
        'derived_seq': [
            "CGATCGATCGATCGATCGATCGAT",
            "GGCCGGCCGGCCGGCCGGCCGGCC",
            "ATATATATATATATATATATATATAT",
            "GCGCGCGCGCGCGCGCGCGCGCGC",
            "CCAAGGTTCCAAGGTTCCAAGGTT",
        ],
        'align_score': [-10, -5, -20, -3, -15],
    })


# ---------------------------------------------------------------------------
# NUPACK pDup
# ---------------------------------------------------------------------------

class TestNupackPdup:

    @pytest.fixture(autouse=True)
    def _skip_if_no_nupack(self):
        pytest.importorskip("nupack")

    def test_perfect_complement_high_pdup(self):
        """A probe paired with its perfect complement should have high pDup."""
        from oligominer.thermodynamics.nupack import calc_pdup

        pdup = calc_pdup(PROBE_SEQ)
        assert 0.0 <= pdup <= 1.0
        assert pdup > 0.5

    def test_mismatched_pair_lower_pdup(self):
        """A mismatched pair should have lower pDup than a perfect match."""
        from oligominer.thermodynamics.nupack import calc_pdup

        pdup_perfect = calc_pdup(PROBE_SEQ)
        pdup_mismatch = calc_pdup(PROBE_SEQ, seq_b="AAAAAAAAAAAAAAAAAAAAAAAAA")
        assert pdup_mismatch < pdup_perfect

    def test_pdup_in_unit_range(self):
        """pDup should always be between 0 and 1."""
        from oligominer.thermodynamics.nupack import calc_pdup

        pdup = calc_pdup(PROBE_SEQ, seq_b=DERIVED_SEQ)
        assert 0.0 <= pdup <= 1.0

    def test_custom_concentrations(self):
        """Changing strand concentrations should affect pDup."""
        from oligominer.thermodynamics.nupack import calc_pdup

        pdup_default = calc_pdup(PROBE_SEQ)
        pdup_high_b = calc_pdup(PROBE_SEQ, conc_b=1e-6)
        # different concentrations should give a different result
        assert pdup_default != pytest.approx(pdup_high_b, abs=1e-6)


# ---------------------------------------------------------------------------
# PaintSHOP XGBoost model
# ---------------------------------------------------------------------------

class TestPaintshopXgboost:

    @pytest.fixture(autouse=True)
    def _skip_if_no_xgboost(self):
        pytest.importorskip("xgboost")

    def test_load_model_valid_temperatures(self):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            load_model, AVAILABLE_TEMPERATURES,
        )

        for temp in AVAILABLE_TEMPERATURES:
            model = load_model(temp)
            assert model is not None

    def test_load_model_invalid_temperature(self):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            load_model,
        )

        with pytest.raises(ConfigurationError, match="No model available"):
            load_model(99)

    def test_load_model_caching(self):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            load_model,
        )

        model_a = load_model(37)
        model_b = load_model(37)
        assert model_a is model_b

    def test_compute_features_shape(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            compute_features, FEATURE_COLUMNS,
        )

        features = compute_features(multi_pair_df)
        assert features.shape == (len(multi_pair_df), len(FEATURE_COLUMNS))

    def test_compute_features_columns(self, single_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            compute_features, FEATURE_COLUMNS,
        )

        features = compute_features(single_pair_df)
        assert list(features.columns) == list(FEATURE_COLUMNS)

    def test_predict_batch_shape(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            predict_duplex_batch,
        )

        preds = predict_duplex_batch(multi_pair_df, temperature=37)
        assert preds.shape == (len(multi_pair_df),)

    def test_predict_batch_normalized_range(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            predict_duplex_batch,
        )

        preds = predict_duplex_batch(multi_pair_df, temperature=37, normalize=True)
        assert np.all(preds >= 0.0)
        assert np.all(preds <= 1.0)

    def test_predict_batch_unnormalized_range(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            predict_duplex_batch,
        )

        preds = predict_duplex_batch(multi_pair_df, temperature=37, normalize=False)
        assert np.all(preds >= 0.0)
        assert np.all(preds <= 100.0)

    def test_predict_single(self):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            predict_duplex,
        )

        pred = predict_duplex(PROBE_SEQ, DERIVED_SEQ, ALIGN_SCORE, temperature=37)
        assert isinstance(pred, float)
        assert 0.0 <= pred <= 1.0

    def test_predict_single_matches_batch(self, single_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            predict_duplex, predict_duplex_batch,
        )

        single = predict_duplex(PROBE_SEQ, DERIVED_SEQ, ALIGN_SCORE, temperature=42)
        batch = predict_duplex_batch(single_pair_df, temperature=42)
        assert single == pytest.approx(float(batch[0]), abs=1e-6)

    def test_different_temperatures_give_different_results(self, single_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            predict_duplex_batch,
        )

        pred_37 = predict_duplex_batch(single_pair_df, temperature=37)[0]
        pred_60 = predict_duplex_batch(single_pair_df, temperature=60)[0]
        assert pred_37 != pytest.approx(pred_60, abs=1e-6)


# ---------------------------------------------------------------------------
# Legacy OligoMiner v1 LDA model
# ---------------------------------------------------------------------------

class TestLegacyLda:

    @pytest.fixture(autouse=True)
    def _skip_if_no_sklearn(self):
        pytest.importorskip("sklearn")

    def test_load_model_valid_temperatures(self):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            load_model, AVAILABLE_TEMPERATURES,
        )

        for temp in AVAILABLE_TEMPERATURES:
            model = load_model(temp)
            assert model is not None

    def test_load_model_invalid_temperature(self):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            load_model,
        )

        with pytest.raises(ConfigurationError, match="No model available"):
            load_model(99)

    def test_load_model_caching(self):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            load_model,
        )

        model_a = load_model(42)
        model_b = load_model(42)
        assert model_a is model_b

    def test_compute_features_shape(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            compute_features, FEATURE_COLUMNS,
        )

        features = compute_features(multi_pair_df)
        assert features.shape == (len(multi_pair_df), len(FEATURE_COLUMNS))

    def test_compute_features_columns(self, single_pair_df):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            compute_features, FEATURE_COLUMNS,
        )

        features = compute_features(single_pair_df)
        assert list(features.columns) == list(FEATURE_COLUMNS)

    def test_compute_features_gc_is_percentage(self, single_pair_df):
        """GC content should be 0-100 (percentage) to match BioPython GC()."""
        from oligominer.specificity.duplex_stability.legacy_lda import (
            compute_features,
        )

        features = compute_features(single_pair_df)
        gc = features['probe_gc'].iloc[0]
        assert 0.0 <= gc <= 100.0
        # ATCGATCG... is 50% GC
        assert gc == pytest.approx(50.0, abs=1.0)

    def test_predict_batch_shape(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex_batch,
        )

        preds = predict_duplex_batch(multi_pair_df, temperature=42)
        assert preds.shape == (len(multi_pair_df),)

    def test_predict_batch_normalized_range(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex_batch,
        )

        preds = predict_duplex_batch(multi_pair_df, temperature=42, normalize=True)
        assert np.all(preds >= 0.0)
        assert np.all(preds <= 1.0)

    def test_predict_batch_unnormalized_range(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex_batch,
        )

        preds = predict_duplex_batch(multi_pair_df, temperature=42, normalize=False)
        assert np.all(preds >= 0.0)
        assert np.all(preds <= 100.0)

    def test_predict_single(self):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex,
        )

        pred = predict_duplex(PROBE_SEQ, ALIGN_SCORE, temperature=42)
        assert isinstance(pred, float)
        assert 0.0 <= pred <= 1.0

    def test_predict_single_matches_batch(self, single_pair_df):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex, predict_duplex_batch,
        )

        single = predict_duplex(PROBE_SEQ, ALIGN_SCORE, temperature=42)
        batch = predict_duplex_batch(single_pair_df, temperature=42)
        assert single == pytest.approx(float(batch[0]), abs=1e-6)

    def test_different_temperatures_give_different_results(self, single_pair_df):
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex_batch,
        )

        pred_32 = predict_duplex_batch(single_pair_df, temperature=32)[0]
        pred_57 = predict_duplex_batch(single_pair_df, temperature=57)[0]
        assert pred_32 != pytest.approx(pred_57, abs=1e-6)

    def test_high_gc_higher_offtarget_at_high_temp(self):
        """At high temperatures, high-GC probes should have higher off-target
        probability than low-GC probes (GC pairs are more stable)."""
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex_batch,
        )

        df = pd.DataFrame({
            'probe_seq': [
                "ATATATATATATATATATATATATAT",  # low GC
                "GCGCGCGCGCGCGCGCGCGCGCGCG",  # high GC
            ],
            'align_score': [-10, -10],
        })

        preds = predict_duplex_batch(df, temperature=57)
        assert preds[1] > preds[0]


# ---------------------------------------------------------------------------
# cross-model consistency
# ---------------------------------------------------------------------------

class TestCrossModelConsistency:
    """Verify that XGBoost and LDA models have consistent interfaces and
    produce outputs on the same scale, even if the values differ."""

    @pytest.fixture(autouse=True)
    def _skip_if_missing_deps(self):
        pytest.importorskip("xgboost")
        pytest.importorskip("sklearn")

    def test_both_models_return_same_shape(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            predict_duplex_batch as xgb_predict,
        )
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex_batch as lda_predict,
        )

        xgb_preds = xgb_predict(multi_pair_df, temperature=37)
        lda_preds = lda_predict(multi_pair_df, temperature=37)
        assert xgb_preds.shape == lda_preds.shape

    def test_both_models_in_unit_range(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            predict_duplex_batch as xgb_predict,
        )
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex_batch as lda_predict,
        )

        xgb_preds = xgb_predict(multi_pair_df, temperature=37, normalize=True)
        lda_preds = lda_predict(multi_pair_df, temperature=37, normalize=True)

        for preds in (xgb_preds, lda_preds):
            assert np.all(preds >= 0.0)
            assert np.all(preds <= 1.0)

    def test_both_models_in_percent_range(self, multi_pair_df):
        from oligominer.specificity.duplex_stability.paintshop_xgboost import (
            predict_duplex_batch as xgb_predict,
        )
        from oligominer.specificity.duplex_stability.legacy_lda import (
            predict_duplex_batch as lda_predict,
        )

        xgb_preds = xgb_predict(multi_pair_df, temperature=37, normalize=False)
        lda_preds = lda_predict(multi_pair_df, temperature=37, normalize=False)

        for preds in (xgb_preds, lda_preds):
            assert np.all(preds >= 0.0)
            assert np.all(preds <= 100.0)
