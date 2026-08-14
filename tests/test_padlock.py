"""Tests for padlock (split-homology) probe mining.

A padlock is one oligo whose two ends carry the target homology and whose middle
is a non-hybridizing backbone. Probe and target are antiparallel, so the probe's
5' arm pairs with the downstream target segment and its 3' arm with the upstream
one. The identity that makes this work, and that these tests assert, is that the
two arms concatenated in probe order are the reverse complement of one contiguous
target window.
"""

import random

import pytest

from oligominer.probe_design import (
    PADLOCK_COLUMNS,
    arm_params,
    check_identity,
    mine_padlock_sequence,
    padlocks_to_df,
)
from oligominer.utils.seq_utils import rev_comp


@pytest.fixture(scope='module')
def target():
    """A random 20 kb target sequence."""
    random.seed(2)

    # success
    return ''.join(random.choice('ACGT') for _ in range(20000))


class TestArmParams:

    def test_defaults_are_complete(self):
        params = arm_params()
        for key in ('min_length', 'max_length', 'min_tm', 'max_tm', 'min_gc', 'max_gc'):
            assert key in params

    def test_overrides_are_applied(self):
        assert arm_params(min_length=18)['min_length'] == 18

    def test_an_unknown_parameter_raises(self):
        with pytest.raises(ValueError, match='unknown arm parameter'):
            arm_params(min_lenght=18)


class TestMining:

    def test_padlocks_are_found(self, target):
        assert len(mine_padlock_sequence(target, seq_id='chr1')) > 0

    def test_rows_match_the_declared_columns(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1'))
        assert list(df.columns) == PADLOCK_COLUMNS

    def test_a_negative_gap_raises(self, target):
        with pytest.raises(ValueError, match='gap must be'):
            mine_padlock_sequence(target, seq_id='chr1', gap=-1)

    def test_arm_length_bounds_are_honored(self, target):
        df = padlocks_to_df(mine_padlock_sequence(
            target, seq_id='chr1',
            arm5=arm_params(min_length=18, max_length=20),
            arm3=arm_params(min_length=18, max_length=20)))
        assert df['arm5_len'].between(18, 20).all()
        assert df['arm3_len'].between(18, 20).all()

    def test_footprints_do_not_overlap_by_default(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1'))
        spans = sorted(zip(df['start'], df['stop']))
        assert all(a[1] <= b[0] for a, b in zip(spans, spans[1:]))


class TestIdentity:
    """The two arms in probe order are the reverse complement of the target window."""

    def test_gapless_padlocks_satisfy_the_identity(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1', gap=0))
        assert len(df) > 0
        assert all(check_identity(row, target) for _, row in df.iterrows())

    def test_gapped_padlocks_satisfy_the_identity(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1', gap=5))
        assert len(df) > 0
        assert all(check_identity(row, target) for _, row in df.iterrows())

    def test_probe_seq_is_the_revcomp_of_the_window_when_gapless(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1', gap=0))
        row = df.iloc[0]
        assert row['probe_seq'] == rev_comp(target[row['start']:row['stop']])

    def test_probe_seq_is_the_two_arms_concatenated(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1', gap=0))
        row = df.iloc[0]
        assert row['probe_seq'] == row['arm_5p'] + row['arm_3p']

    def test_the_junction_is_where_the_arms_meet(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1', gap=0))
        row = df.iloc[0]
        assert row['junction'] == len(row['arm_5p'])

    def test_the_gap_sequence_is_recorded(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1', gap=5))
        row = df.iloc[0]
        assert len(row['gap_seq']) == 5
        assert row['gap_len'] == 5

    def test_a_gapless_padlock_records_no_gap(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1', gap=0))
        assert (df['gap_len'] == 0).all()


class TestDownstreamCompatibility:
    """probe_seq must be usable by stages that know nothing about padlocks."""

    def test_probe_seq_is_plain_acgt(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1'))
        assert df['probe_seq'].str.fullmatch('[ACGT]+').all()

    def test_probe_length_equals_the_two_arms(self, target):
        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1'))
        assert (df['probe_seq'].str.len() == df['arm5_len'] + df['arm3_len']).all()

    def test_padlocks_can_be_screened_for_kmers(self, target, tmp_path):
        """A padlock's homology goes through k-mer screening as one sequence."""
        from oligominer.specificity.kmers import build_index, max_kmer

        fasta = tmp_path / 'ref.fa'
        fasta.write_text('>chr1\n' + '\n'.join(
            target[i:i + 60] for i in range(0, len(target), 60)) + '\n')
        index = tmp_path / 'ref.npz'
        build_index(str(fasta), str(index), k=12, backend='numpy', min_count=2)

        df = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1'))
        counts = max_kmer(index, df['probe_seq'].tolist()[:20], k=12, backend='numpy')
        assert len(counts) == 20
