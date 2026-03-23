"""Tests for the thermodynamics subpackage.

Covers: integer encoding, Tm grid calculation, formamide correction,
and nearest-neighbor parameter tables.
"""

import numpy as np
import pytest

from oligominer.thermodynamics.mining.int_encoding import seq_to_8bit, DNA_ASCII_LUT
from oligominer.thermodynamics.mining.calc_tm_2d import (
    get_dinuc_grid, get_dS_grid, get_dH_grid, get_tm_grid,
)
from oligominer.thermodynamics.mining.config import GET_DEFAULT_MINING_CONFIG
from oligominer.thermodynamics import formamide_correction


# ---------------------------------------------------------------------------
# integer encoding
# ---------------------------------------------------------------------------

class TestIntEncoding:

    def test_acgt_encoding(self):
        arr = seq_to_8bit("ACGT")
        assert list(arr) == [0, 1, 2, 3]

    def test_lowercase_encoding(self):
        arr = seq_to_8bit("acgt")
        assert list(arr) == [0, 1, 2, 3]

    def test_n_encoding(self):
        arr = seq_to_8bit("N")
        assert arr[0] == 4

    def test_mixed_sequence(self):
        arr = seq_to_8bit("ACNGT")
        assert list(arr) == [0, 1, 4, 2, 3]

    def test_output_dtype(self):
        arr = seq_to_8bit("ATCG")
        assert arr.dtype == np.uint8

    def test_empty_sequence(self):
        arr = seq_to_8bit("")
        assert len(arr) == 0


# ---------------------------------------------------------------------------
# Tm grid calculation
# ---------------------------------------------------------------------------

class TestTmGrid:

    @pytest.fixture
    def config(self):
        config = GET_DEFAULT_MINING_CONFIG()
        config['min_length'] = 18
        config['max_length'] = 22
        return config

    @pytest.fixture
    def nuc_array(self):
        # a 100-base sequence with no Ns
        seq = "ATCGATCGATCGATCGATCG" * 5
        return seq_to_8bit(seq)

    def test_dinuc_grid_shape(self, nuc_array):
        grid = get_dinuc_grid(nuc_array, max_length=22)
        n_positions = len(nuc_array) - 22 + 1
        assert grid.shape[0] == n_positions
        assert grid.shape[1] == 2  # 5' base, 3' base
        assert grid.shape[2] == 21  # max_length - 1

    def test_tm_grid_shape(self, nuc_array, config):
        tm = get_tm_grid(nuc_array, config)
        n_positions = len(nuc_array) - config['max_length'] + 1
        n_lengths = config['max_length'] - config['min_length'] + 1
        assert tm.shape == (n_positions, n_lengths)

    def test_tm_values_are_reasonable(self, nuc_array, config):
        config['pct_formamide'] = 0
        tm = get_tm_grid(nuc_array, config)
        # Tm should be in a reasonable range for short DNA oligos
        assert tm.min() > 0
        assert tm.max() < 120

    def test_longer_probes_have_higher_tm(self, nuc_array, config):
        """For most sequences, longer probes should have higher Tm."""
        config['pct_formamide'] = 0
        tm = get_tm_grid(nuc_array, config)
        # check that column means increase with length (roughly)
        col_means = tm.mean(axis=0)
        # not strictly monotonic for all sequences, but generally true
        assert col_means[-1] > col_means[0]

    def test_salt_correction_affects_ds(self):
        seq = "ATCGATCGATCGATCGATCG" * 5
        nuc_array = seq_to_8bit(seq)
        grid = get_dinuc_grid(nuc_array, max_length=22)

        ds_low = get_dS_grid(grid, 18, 22, Na=50)
        ds_high = get_dS_grid(grid, 18, 22, Na=500)
        assert not np.allclose(ds_low, ds_high)


# ---------------------------------------------------------------------------
# formamide correction
# ---------------------------------------------------------------------------

class TestFormamideCorrection:

    def test_basic_correction(self):
        result = formamide_correction(37.0, 50, 0.65)
        expected = 37.0 - (50 * 0.65)
        assert result == pytest.approx(expected)

    def test_zero_formamide(self):
        result = formamide_correction(37.0, 0)
        assert result == 37.0

    def test_custom_factor(self):
        result = formamide_correction(37.0, 50, 0.5)
        expected = 37.0 - (50 * 0.5)
        assert result == pytest.approx(expected)
