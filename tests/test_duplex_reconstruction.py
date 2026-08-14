"""Tests for genomic duplex reconstruction.

Covers interval clamping, strand-aware sequence lookup, and agreement with the
bedtools path the lookup replaces. The bedtools comparison skips when the binary
is not installed.
"""

import random
import shutil

import pandas as pd
import pytest

from oligominer.specificity.alignment.duplex import (
    BED_COLUMNS,
    chrom_sizes,
    clamp_intervals,
    fetch_derived_seqs,
    load_reference,
    reconstruct,
)
from oligominer.utils.seq_utils import rev_comp

needs_bedtools = pytest.mark.skipif(
    shutil.which('bedtools') is None, reason='bedtools not installed'
)


@pytest.fixture
def genome(tmp_path):
    """A two-chromosome reference with known sequence."""
    random.seed(7)
    chr1 = ''.join(random.choice('ACGT') for _ in range(1000))
    chr2 = ''.join(random.choice('ACGT') for _ in range(500))

    path = tmp_path / 'genome.fa'
    with open(path, 'w') as handle:
        for name, seq in (('chr1', chr1), ('chr2', chr2)):
            handle.write(f'>{name}\n')
            for i in range(0, len(seq), 60):
                handle.write(seq[i:i + 60] + '\n')

    # success
    return path, {'chr1': chr1, 'chr2': chr2}


def _align_row(seqid, start, stop, strand, name='p1'):
    """Build one alignment record."""
    return {
        'align_seqid': seqid, 'align_start': start, 'align_stop': stop,
        'seqid': name, 'align_score': 0, 'align_strand': strand,
        'align_cigar': f'{stop - start}M',
    }


class TestChromSizes:

    def test_sizes_match_the_reference(self, genome):
        path, seqs = genome
        sizes = chrom_sizes(load_reference(path))
        assert sizes == {'chr1': 1000, 'chr2': 500}


class TestClamping:

    def test_interval_past_the_end_is_clamped(self, genome):
        path, _ = genome
        df = pd.DataFrame([_align_row('chr2', 480, 600, '+')])
        out = clamp_intervals(df, chrom_sizes(load_reference(path)))
        assert out['align_stop'].iloc[0] == 500
        assert bool(out['was_clamped'].iloc[0])

    def test_negative_start_is_clamped(self, genome):
        path, _ = genome
        df = pd.DataFrame([_align_row('chr1', -10, 50, '+')])
        out = clamp_intervals(df, chrom_sizes(load_reference(path)))
        assert out['align_start'].iloc[0] == 0
        assert bool(out['was_clamped'].iloc[0])

    def test_interval_inside_the_chromosome_is_untouched(self, genome):
        path, _ = genome
        df = pd.DataFrame([_align_row('chr1', 100, 140, '+')])
        out = clamp_intervals(df, chrom_sizes(load_reference(path)))
        assert out['align_start'].iloc[0] == 100
        assert out['align_stop'].iloc[0] == 140
        assert not bool(out['was_clamped'].iloc[0])

    def test_unknown_chromosome_raises(self, genome):
        path, _ = genome
        df = pd.DataFrame([_align_row('chrZ', 0, 40, '+')])
        with pytest.raises(KeyError, match='absent from the FASTA'):
            clamp_intervals(df, chrom_sizes(load_reference(path)))


class TestDerivedSequence:

    def test_plus_strand_returns_the_reference_span(self, genome):
        path, seqs = genome
        df = pd.DataFrame([_align_row('chr1', 100, 140, '+')])
        out = reconstruct(df, path)
        assert out['derived_seq'].iloc[0] == seqs['chr1'][100:140]

    def test_minus_strand_returns_the_reverse_complement(self, genome):
        path, seqs = genome
        df = pd.DataFrame([_align_row('chr1', 100, 140, '-')])
        out = reconstruct(df, path)
        assert out['derived_seq'].iloc[0] == rev_comp(seqs['chr1'][100:140])

    def test_the_two_strands_differ(self, genome):
        """Strand handling is only meaningful if it changes the answer."""
        path, _ = genome
        df = pd.DataFrame([_align_row('chr1', 100, 140, '+'),
                           _align_row('chr1', 100, 140, '-')])
        out = reconstruct(df, path)
        assert out['derived_seq'].iloc[0] != out['derived_seq'].iloc[1]

    def test_each_row_gets_its_own_sequence(self, genome):
        path, seqs = genome
        rows = [_align_row('chr1', 0, 40, '+'),
                _align_row('chr2', 10, 50, '+'),
                _align_row('chr1', 900, 940, '-')]
        out = reconstruct(pd.DataFrame(rows), path)
        assert out['derived_seq'].iloc[0] == seqs['chr1'][0:40]
        assert out['derived_seq'].iloc[1] == seqs['chr2'][10:50]
        assert out['derived_seq'].iloc[2] == rev_comp(seqs['chr1'][900:940])

    def test_interleaved_chromosomes_stay_aligned_to_their_rows(self, genome):
        """Grouping by chromosome must not permute the results."""
        path, seqs = genome
        rows = []
        expected = []
        for i in range(20):
            chrom = 'chr1' if i % 2 == 0 else 'chr2'
            start = i * 10
            rows.append(_align_row(chrom, start, start + 30, '+'))
            expected.append(seqs[chrom][start:start + 30])
        out = reconstruct(pd.DataFrame(rows), path)
        assert list(out['derived_seq']) == expected

    def test_non_default_index_is_preserved(self, genome):
        path, seqs = genome
        df = pd.DataFrame([_align_row('chr1', 0, 40, '+'),
                           _align_row('chr2', 0, 40, '+')],
                          index=['x', 'y'])
        out = reconstruct(df, path)
        assert out.loc['x', 'derived_seq'] == seqs['chr1'][0:40]
        assert out.loc['y', 'derived_seq'] == seqs['chr2'][0:40]

    def test_derived_sequence_length_matches_the_interval(self, genome):
        path, _ = genome
        df = pd.DataFrame([_align_row('chr1', 100, 137, '+')])
        out = reconstruct(df, path)
        assert len(out['derived_seq'].iloc[0]) == 37

    def test_clamped_interval_returns_the_truncated_span(self, genome):
        path, seqs = genome
        df = pd.DataFrame([_align_row('chr2', 480, 600, '+')])
        out = reconstruct(df, path)
        assert out['derived_seq'].iloc[0] == seqs['chr2'][480:500]

    def test_soft_masked_reference_is_upper_cased_by_default(self, tmp_path):
        path = tmp_path / 'masked.fa'
        path.write_text('>chr1\n' + 'acgt' * 25 + '\n')
        df = pd.DataFrame([_align_row('chr1', 0, 20, '+')])
        out = reconstruct(df, path)
        assert out['derived_seq'].iloc[0] == 'ACGT' * 5

    def test_case_can_be_preserved(self, tmp_path):
        path = tmp_path / 'masked.fa'
        path.write_text('>chr1\n' + 'acgt' * 25 + '\n')
        df = pd.DataFrame([_align_row('chr1', 0, 20, '+')])
        out = reconstruct(df, path, to_upper=False)
        assert out['derived_seq'].iloc[0] == 'acgt' * 5


@needs_bedtools
class TestBedtoolsParity:
    """The lookup must reproduce what bedtools slop + getfasta -s returned."""

    def test_derived_sequences_match_bedtools(self, genome):
        from oligominer.specificity.alignment import get_fasta, trim_bed_coords

        path, _ = genome
        random.seed(11)
        rows = []
        for _ in range(200):
            chrom = random.choice(['chr1', 'chr2'])
            limit = 1000 if chrom == 'chr1' else 500
            start = random.randint(0, limit - 50)
            rows.append(_align_row(chrom, start, start + random.randint(20, 45),
                                   random.choice('+-')))
        align_df = pd.DataFrame(rows)

        fast = reconstruct(align_df, path)['derived_seq']

        bed_data = '\n'.join(
            '\t'.join(str(align_df.iloc[i][c]) for c in BED_COLUMNS)
            for i in range(len(align_df))
        ) + '\n'
        trimmed = trim_bed_coords(bed_data=bed_data, fasta_path=str(path))
        slow = [s.upper() for s in
                get_fasta(bed_data=trimmed, fasta_path=str(path)).strip().split('\n')]

        mismatches = [i for i, (a, b) in enumerate(zip(fast, slow)) if a != b]
        assert mismatches == [], (
            f'{len(mismatches)}/{len(align_df)} rows differ, first at {mismatches[:3]}'
        )
