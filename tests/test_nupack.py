"""Tests for the NUPACK thermodynamics subpackage.

Covers: pDup calculation, competitive pDup, and input validation.

All tests are skipped when NUPACK is not installed.
"""

import pytest

nupack = pytest.importorskip("nupack", reason="NUPACK is not installed")

from oligominer.thermodynamics.nupack.pdup import calc_pdup, calc_competitive_pdup


# a short probe / target pair for fast tests
PROBE_A = "ATCGATCGATCGATCGATCG"
PROBE_B = "GCTAGCTAGCTAGCTAGCTA"
TARGET = "CGATCGATCGATCGATCGAT"


# ---------------------------------------------------------------------------
# calc_pdup
# ---------------------------------------------------------------------------

class TestCalcPdup:

    def test_returns_float(self):
        result = calc_pdup(PROBE_A)
        assert isinstance(result, float)

    def test_pdup_between_zero_and_one(self):
        result = calc_pdup(PROBE_A)
        assert 0.0 <= result <= 1.0

    def test_perfect_complement_high_pdup(self):
        """A probe with its exact reverse complement should have high pDup."""
        result = calc_pdup(PROBE_A, conc_a=1e-6, conc_b=1e-6)
        assert result > 0.1

    def test_explicit_seq_b(self):
        result = calc_pdup(PROBE_A, seq_b=TARGET)
        assert 0.0 <= result <= 1.0

    def test_empty_seq_a_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            calc_pdup("")

    def test_empty_seq_b_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            calc_pdup(PROBE_A, seq_b="")


# ---------------------------------------------------------------------------
# calc_competitive_pdup
# ---------------------------------------------------------------------------

class TestCalcCompetitivePdup:

    def test_returns_two_floats(self):
        at_pdup, bt_pdup = calc_competitive_pdup(PROBE_A, PROBE_B, TARGET)
        assert isinstance(at_pdup, float)
        assert isinstance(bt_pdup, float)

    def test_values_between_zero_and_one(self):
        at_pdup, bt_pdup = calc_competitive_pdup(PROBE_A, PROBE_B, TARGET)
        assert 0.0 <= at_pdup <= 1.0
        assert 0.0 <= bt_pdup <= 1.0

    def test_empty_seq_a_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            calc_competitive_pdup("", PROBE_B, TARGET)

    def test_empty_seq_b_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            calc_competitive_pdup(PROBE_A, "", TARGET)

    def test_empty_target_raises(self):
        with pytest.raises(ValueError, match="non-empty"):
            calc_competitive_pdup(PROBE_A, PROBE_B, "")
