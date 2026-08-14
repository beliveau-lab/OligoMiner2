"""Tests for soft-mask aware mining.

Soft-masked (repeat and low-complexity) regions are marked by lowercase in genome
FASTA files. Excluding them requires both that the encoder distinguishes case and
that the miner does not upper-case the sequence before encoding it, so these
tests cover the encoder, the miner and the casing of the probes they produce.
"""

import warnings

import numpy as np
import pytest

from oligominer.thermodynamics.mining import mine_sequence
from oligominer.thermodynamics.mining.int_encoding import (
    DNA_ASCII_LUT,
    SOFTMASK_ASCII_LUT,
    seq_to_8bit,
)

# 60 nt with the middle 20 soft-masked
UPPER = 'ACGTACGTACGTACGTACGT'
LOWER = 'acgtacgtacgtacgtacgt'
MIXED = UPPER + LOWER + UPPER

MINING = dict(seq_id='t', min_length=20, max_length=20, min_tm=0, max_tm=200,
              min_gc=0, max_gc=100, max_homopolymer=None)


class TestEncoding:

    def test_default_lut_folds_case(self):
        assert np.array_equal(seq_to_8bit('acgt'), seq_to_8bit('ACGT'))

    def test_softmask_lut_sends_lowercase_to_the_ambiguous_code(self):
        encoded = seq_to_8bit('acgt', mask_soft=True)
        assert np.array_equal(encoded, np.full(4, 4, dtype=np.uint8))

    def test_softmask_lut_leaves_uppercase_alone(self):
        assert np.array_equal(
            seq_to_8bit('ACGT', mask_soft=True), seq_to_8bit('ACGT')
        )

    def test_n_still_encodes_as_ambiguous_under_both_tables(self):
        assert seq_to_8bit('N')[0] == 4
        assert seq_to_8bit('N', mask_soft=True)[0] == 4

    def test_the_two_tables_differ_only_in_lowercase_bases(self):
        differing = np.nonzero(DNA_ASCII_LUT != SOFTMASK_ASCII_LUT)[0]
        assert sorted(chr(i) for i in differing) == ['a', 'c', 'g', 't']


class TestMiningHonorsTheMask:

    def test_masked_region_is_mined_when_the_flag_is_off(self):
        probes = mine_sequence(MIXED, mask_soft=False, **MINING)
        starts = {start for _, start, _, _, _ in probes}
        # a probe starting inside the masked middle third
        assert any(20 <= s < 40 for s in starts)

    def test_masked_region_is_excluded_when_the_flag_is_on(self):
        probes = mine_sequence(MIXED, mask_soft=True, **MINING)
        for _, start, stop, _, _ in probes:
            assert stop <= 20 or start >= 40, (
                f"probe {start}-{stop} overlaps the soft-masked region 20-40"
            )

    def test_masking_reduces_the_probe_count(self):
        unmasked = mine_sequence(MIXED, mask_soft=False, **MINING)
        masked = mine_sequence(MIXED, mask_soft=True, **MINING)
        assert len(masked) < len(unmasked)

    def test_a_fully_masked_sequence_yields_nothing(self):
        probes = mine_sequence(LOWER * 3, mask_soft=True, **MINING)
        assert probes == []

    def test_an_unmasked_sequence_is_unaffected_by_the_flag(self):
        """The flag must be a no-op on input that carries no soft-masking."""
        seq = UPPER * 3
        with pytest.warns(RuntimeWarning):
            masked = mine_sequence(seq, mask_soft=True, **MINING)
        assert masked == mine_sequence(seq, mask_soft=False, **MINING)


class TestUnmaskedReferenceWarning:
    """Asking to mask a sequence that carries no lowercase warns.

    Not every assembly ships a soft-masked FASTA, and on one that does not, the
    filter returns every probe. The warning distinguishes that from a filter
    that ran and found nothing to remove.
    """

    def test_warns_when_asked_to_mask_an_unmasked_sequence(self):
        with pytest.warns(RuntimeWarning, match="no lowercase"):
            mine_sequence(UPPER * 3, mask_soft=True, **MINING)

    def test_the_warning_names_the_sequence(self):
        params = dict(MINING)
        params['seq_id'] = 'chr19'
        with pytest.warns(RuntimeWarning, match="chr19"):
            mine_sequence(UPPER * 3, mask_soft=True, **params)

    def test_no_warning_when_the_sequence_is_masked(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            mine_sequence(MIXED, mask_soft=True, **MINING)

    def test_no_warning_when_masking_is_off(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            mine_sequence(UPPER * 3, mask_soft=False, **MINING)


class TestOutputCasing:
    """Encoding from the raw sequence must not leak lowercase into the output."""

    def test_probe_sequences_are_upper_case_even_from_mixed_case_input(self):
        probes = mine_sequence(MIXED, mask_soft=False, **MINING)
        assert probes, "expected probes from the mixed-case sequence"
        for _, _, _, probe_seq, _ in probes:
            assert probe_seq == probe_seq.upper()

    def test_lower_case_input_gives_the_same_probes_as_upper_case(self):
        """With masking off, case must not change the result at all."""
        lower = mine_sequence(LOWER * 3, mask_soft=False, **MINING)
        upper = mine_sequence(UPPER * 3, mask_soft=False, **MINING)
        assert lower == upper
