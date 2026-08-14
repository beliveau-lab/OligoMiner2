"""Tests for querying a Jellyfish index directly."""

import shutil
import subprocess

import pytest

from oligominer.specificity.kmers.jellyfish_query import (
    _get_kmers, calc_max_kmer, calc_max_kmer_multi, jellyfish_query,
)

K = 12

needs_jellyfish = pytest.mark.skipif(shutil.which('jellyfish') is None,
                                     reason='jellyfish is not installed')

# a 12-mer placed three times, and one placed once
REPEATED = 'ACGTACGTACGT'
UNIQUE = 'TTGGCCAATTGA'


@pytest.fixture
def index(tmp_path):
    """A Jellyfish index over a sequence with a known repeat structure."""
    spacer = 'GGGGGGGGGGGGGGGG'
    genome = tmp_path / 'genome.fa'
    genome.write_text('>chr1\n'
                      + spacer.join([REPEATED, REPEATED, REPEATED, UNIQUE])
                      + '\n')

    path = tmp_path / 'index.jf'
    subprocess.run(['jellyfish', 'count', '-m', str(K), '-s', '10000',
                    '-o', str(path), str(genome)], check=True)
    return str(path)


class TestGetKmers:

    def test_a_sequence_yields_one_kmer_per_offset(self):
        assert len(_get_kmers('A' * 20, K)) == 20 - K + 1

    def test_the_kmers_are_the_sliding_windows(self):
        assert _get_kmers('ACGTA', 3) == ['ACG', 'CGT', 'GTA']

    def test_a_sequence_of_exactly_k_yields_one(self):
        assert _get_kmers('A' * K, K) == ['A' * K]

    def test_a_sequence_shorter_than_k_yields_none(self):
        assert _get_kmers('ACGT', K) == []


@needs_jellyfish
class TestCalcMaxKmer:

    def test_a_repeated_sequence_scores_its_repeat_count(self, index):
        assert calc_max_kmer(index, REPEATED, K) == 3

    def test_a_unique_sequence_scores_one(self, index):
        assert calc_max_kmer(index, UNIQUE, K) == 1

    def test_the_maximum_over_a_mixed_sequence_wins(self, index):
        # a sequence containing both reports the repeated one
        assert calc_max_kmer(index, REPEATED + UNIQUE, K) == 3

    def test_a_sequence_shorter_than_k_scores_zero(self, index):
        # it holds no k-mer at all; the backend dispatcher reports 0 for this
        # case, and this helper must not disagree with it
        assert calc_max_kmer(index, 'ACGT', K) == 0

    def test_it_agrees_with_the_backend_dispatcher(self, index):
        from oligominer.specificity.kmers import max_kmer

        for seq in (REPEATED, UNIQUE, REPEATED + UNIQUE, 'ACGT'):
            assert calc_max_kmer(index, seq, K) == int(
                max_kmer(index, [seq], k=K, backend='jellyfish')[0])


@needs_jellyfish
class TestCalcMaxKmerMulti:

    def test_one_count_per_sequence_in_order(self, index):
        assert calc_max_kmer_multi(index, [REPEATED, UNIQUE], K) == [3, 1]

    def test_no_sequences_yields_no_counts(self, index):
        assert calc_max_kmer_multi(index, [], K) == []


@needs_jellyfish
class TestJellyfishQuery:

    def test_each_queried_mer_comes_back(self, index):
        out = jellyfish_query(index, mers=[REPEATED, UNIQUE])
        lines = [line for line in out.strip().splitlines() if line]

        assert len(lines) == 2

    def test_the_count_is_the_second_field(self, index):
        out = jellyfish_query(index, mers=[REPEATED])
        assert int(out.split()[1]) == 3
