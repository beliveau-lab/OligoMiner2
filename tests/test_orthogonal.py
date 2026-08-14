"""Tests for orthogonal sequence design.

The load-bearing claim is that nominating candidate cross-hybridizing pairs with
a k-mer join is exact rather than approximate: a shared complementary run of
length k or more is exactly a shared k-mer between one sequence and the reverse
complement of the other. These tests assert that against a brute-force screen.
"""

import random

import pandas as pd
import pytest

from oligominer.probe_design import (
    evict,
    nominate,
    screen_bruteforce,
    verify_against_screen,
    worst_pairs,
)
from oligominer.utils.seq_utils import rev_comp


@pytest.fixture
def barcodes():
    """Sixty random 25-mers."""
    random.seed(4)

    # success
    return [''.join(random.choice('ACGT') for _ in range(25)) for _ in range(60)]


class TestExactness:
    """The k-mer join must return the brute-force pair set, not approximate it."""

    def test_the_two_screens_agree(self, barcodes):
        report = verify_against_screen(barcodes, min_run=7)
        assert report['equal'], (
            f"missed {report['missed_by_nominator']}, "
            f"extra {report['extra_in_nominator']}")

    def test_nothing_is_missed_at_several_run_lengths(self, barcodes):
        for min_run in (6, 7, 8, 10):
            report = verify_against_screen(barcodes, min_run=min_run)
            assert report['missed_by_nominator'] == []

    def test_nothing_spurious_is_nominated(self, barcodes):
        for min_run in (6, 7, 8, 10):
            report = verify_against_screen(barcodes, min_run=min_run)
            assert report['extra_in_nominator'] == []

    def test_a_planted_complementary_run_is_found(self):
        # flanks must not be complementary to each other, or they extend the run
        run = 'ACGTTGCAG'
        a = 'CCCC' + run + 'CCCC'
        b = 'CCCC' + rev_comp(run) + 'CCCC'
        pairs = nominate([a, b], min_run=len(run))
        assert len(pairs) == 1
        assert set(pairs.iloc[0][['i', 'j']]) == {0, 1}

    def test_a_run_shorter_than_the_threshold_is_not_nominated(self):
        run = 'ACGTTG'
        a = 'CCCC' + run + 'CCCC'
        b = 'CCCC' + rev_comp(run) + 'CCCC'
        assert len(nominate([a, b], min_run=len(run) + 1)) == 0

    def test_a_longer_run_is_still_found_at_a_shorter_threshold(self):
        run = 'ACGTTGCAGGAT'
        a = 'CCCC' + run + 'CCCC'
        b = 'CCCC' + rev_comp(run) + 'CCCC'
        assert len(nominate([a, b], min_run=7)) == 1


class TestPairTable:

    def test_indices_are_ordered_within_a_pair(self, barcodes):
        pairs = nominate(barcodes, min_run=7)
        assert (pairs['i'] < pairs['j']).all()

    def test_each_pair_appears_once(self, barcodes):
        pairs = nominate(barcodes, min_run=7)
        assert not pairs.duplicated(subset=['i', 'j']).any()

    def test_a_sequence_is_never_paired_with_itself(self, barcodes):
        pairs = nominate(barcodes, min_run=7)
        assert (pairs['i'] != pairs['j']).all()

    def test_the_shared_kmer_is_recorded(self, barcodes):
        pairs = nominate(barcodes, min_run=7)
        assert (pairs['shared_kmer'].str.len() == 7).all()

    def test_an_empty_set_nominates_nothing(self):
        assert len(nominate([], min_run=7)) == 0

    def test_a_single_sequence_nominates_nothing(self):
        assert len(nominate(['ACGTACGTACGTACGTACGTACGTA'], min_run=7)) == 0


class TestOrdering:

    def test_pairs_are_ordered_worst_first(self, barcodes):
        pairs = nominate(barcodes, min_run=7)
        scores = list(range(len(pairs)))
        ordered = worst_pairs(pairs, scores)
        assert list(ordered['score']) == sorted(scores, reverse=True)

    def test_the_worst_n_can_be_taken(self, barcodes):
        pairs = nominate(barcodes, min_run=7)
        ordered = worst_pairs(pairs, range(len(pairs)), n=5)
        assert len(ordered) == 5


class TestEviction:

    def test_eviction_leaves_no_nominated_pair(self, barcodes):
        pairs = nominate(barcodes, min_run=7)
        kept, _ = evict(barcodes, pairs)
        survivors = [barcodes[i] for i in kept]
        assert len(nominate(survivors, min_run=7)) == 0

    def test_every_sequence_is_kept_or_dropped(self, barcodes):
        pairs = nominate(barcodes, min_run=7)
        kept, dropped = evict(barcodes, pairs)
        assert sorted(kept + dropped) == list(range(len(barcodes)))

    def test_nothing_is_dropped_when_nothing_is_nominated(self, barcodes):
        empty = pd.DataFrame(columns=['i', 'j', 'shared_kmer', 'n_shared'])
        kept, dropped = evict(barcodes, empty)
        assert dropped == []
        assert len(kept) == len(barcodes)

    def test_the_kept_side_can_be_chosen(self):
        run = 'ACGTTGCAG'
        seqs = ['TTTT' + run + 'TTTT', 'AAAA' + rev_comp(run) + 'AAAA']
        pairs = nominate(seqs, min_run=len(run))
        assert evict(seqs, pairs, keep_first=True)[1] == [1]
        assert evict(seqs, pairs, keep_first=False)[1] == [0]


class TestBruteForceReference:

    def test_the_reference_finds_a_planted_run(self):
        run = 'ACGTTGCAGG'
        seqs = ['CCCC' + run + 'CCCC', 'CCCC' + rev_comp(run) + 'CCCC']
        assert screen_bruteforce(seqs, min_run=len(run)) == {(0, 1)}

    def test_the_reference_rejects_a_short_run(self):
        run = 'ACGTT'
        seqs = ['CCCC' + run + 'CCCC', 'CCCC' + rev_comp(run) + 'CCCC']
        assert screen_bruteforce(seqs, min_run=len(run) + 2) == set()

    def test_self_complementary_flanks_extend_a_run(self):
        """TTTT and AAAA pair with each other, which is why flanks are chosen
        not to be complementary in the tests above."""
        run = 'ACGTT'
        seqs = ['TTTT' + run + 'TTTT', 'AAAA' + rev_comp(run) + 'AAAA']
        assert screen_bruteforce(seqs, min_run=len(run) + 2) == {(0, 1)}
