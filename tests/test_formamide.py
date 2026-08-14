"""Tests for the two formamide conversions and the default NUPACK model.

Formamide enters a calculation in two opposite directions: it lowers a melting
temperature, and it raises the formamide-free temperature that represents a
hybridization condition. Using one where the other belongs simulates a different
experiment without raising.
"""

import pytest

from oligominer.thermodynamics import effective_hyb_temperature, formamide_correction


class TestFormamideCorrection:
    """Lowering a Tm to account for formamide in the buffer."""

    def test_tm_is_depressed(self):
        assert formamide_correction(80.0, 50, 0.65) == pytest.approx(80.0 - 32.5)

    def test_zero_formamide_is_a_no_op(self):
        assert formamide_correction(37.0, 0) == 37.0

    def test_the_factor_scales_the_depression(self):
        assert formamide_correction(37.0, 50, 0.5) == pytest.approx(37.0 - 25.0)

    def test_more_formamide_lowers_tm_further(self):
        assert formamide_correction(80.0, 60) < formamide_correction(80.0, 30)


class TestEffectiveHybTemperature:
    """Raising a hybridization temperature to its formamide-free equivalent."""

    def test_the_standard_fish_condition_is_69_5(self):
        assert effective_hyb_temperature(37.0, 50, 0.65) == pytest.approx(69.5)

    def test_zero_formamide_is_a_no_op(self):
        assert effective_hyb_temperature(37.0, 0) == 37.0

    def test_more_formamide_raises_the_equivalent_temperature(self):
        assert effective_hyb_temperature(37.0, 60) > effective_hyb_temperature(37.0, 30)

    def test_the_two_conversions_are_inverses(self):
        assert formamide_correction(
            effective_hyb_temperature(37.0, 50, 0.65), 50, 0.65
        ) == pytest.approx(37.0)

    def test_the_two_conversions_move_in_opposite_directions(self):
        assert effective_hyb_temperature(37.0, 50) > 37.0
        assert formamide_correction(37.0, 50) < 37.0


class TestDefaultNupackModel:

    def test_model_runs_at_the_effective_temperature(self):
        pytest.importorskip('nupack')
        from oligominer.thermodynamics.nupack.config import DEFAULT_NUPACK_MODEL

        celsius = DEFAULT_NUPACK_MODEL.temperature - 273.15
        assert celsius == pytest.approx(69.5, abs=0.01)

    def test_pdup_discriminates_at_the_default_condition(self):
        """At a depressed rather than effective temperature every duplex saturates."""
        pytest.importorskip('nupack')
        from oligominer.thermodynamics.nupack import calc_pdup
        from oligominer.utils.seq_utils import rev_comp

        probe = 'GGATCACAGTCTACACTGCTCACTCCAACC'
        perfect = rev_comp(probe)
        mismatched = list(perfect)
        for j in range(0, len(mismatched), 5):
            mismatched[j] = 'A' if mismatched[j] != 'A' else 'C'

        strong = calc_pdup(probe, perfect)
        weak = calc_pdup(probe, ''.join(mismatched))

        assert strong > 0.9
        assert weak < 0.01
        assert strong / weak > 100
