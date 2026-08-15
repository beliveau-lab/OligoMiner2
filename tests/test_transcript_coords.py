"""
Tests for transcript sequence extraction and coordinate mapping.

GTF coordinates are 1-based and inclusive of both ends; probe coordinates are
0-based and half-open. Getting that conversion wrong shifts every probe mined
from a transcript by one base without any error, so the round trip is asserted
directly against the reference sequence.
"""

import pandas as pd
import pytest

from oligominer.bioinformatics.transcriptome.transcript_seq import (
    get_exon_seqs, get_intron_seqs, get_spliced_seq, local_to_genomic,
    parse_interval_label,
)
from oligominer.thermodynamics.mining import mine_sequence
from oligominer.utils.seq_utils import rev_comp


@pytest.fixture
def genome(tmp_path):
    """A single-chromosome FASTA with a known sequence."""
    import random
    random.seed(7)
    sequence = ''.join(random.choice('ACGT') for _ in range(3000))

    path = tmp_path / 'genome.fa'
    with open(path, 'w') as handle:
        handle.write('>chr1\n')
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i:i + 60] + '\n')

    return str(path), sequence


def gtf(rows):
    """Build a GTF frame from (start, end, strand, transcript) tuples."""
    return pd.DataFrame([
        {'seqid': 'chr1', 'type': 'exon', 'start': start, 'end': end,
         'strand': strand, 'transcript_id': transcript, 'gene_id': 'G1',
         'score': '.'}
        for start, end, strand, transcript in rows
    ])


class TestParseIntervalLabel:

    def test_round_trips_a_label(self):
        assert parse_interval_label('chrI:1807-2169(-)') == \
            ('chrI', 1807, 2169, '-')

    def test_a_contig_name_containing_a_colon_survives(self):
        seqid, start, end, strand = parse_interval_label('HLA:A*01:100-200(+)')
        assert seqid == 'HLA:A*01'
        assert (start, end, strand) == (100, 200, '+')


class TestExtractionIsOneBasedInclusive:

    def test_an_exon_spans_end_minus_start_plus_one_bases(self, genome):
        path, sequence = genome
        # GTF 101-200 is 100 bases, not 99
        seqs = get_exon_seqs(gtf([(101, 200, '+', 'T1')]), path,
                             transcript_id='T1')

        assert len(next(iter(seqs.values()))) == 100

    def test_a_plus_strand_exon_matches_the_reference(self, genome):
        path, sequence = genome
        seqs = get_exon_seqs(gtf([(101, 200, '+', 'T1')]), path,
                             transcript_id='T1')

        # 1-based inclusive 101..200 is python slice [100:200]
        assert next(iter(seqs.values())) == sequence[100:200]

    def test_a_minus_strand_exon_is_reverse_complemented(self, genome):
        path, sequence = genome
        seqs = get_exon_seqs(gtf([(101, 200, '-', 'T1')]), path,
                             transcript_id='T1')

        assert next(iter(seqs.values())) == rev_comp(sequence[100:200])

    def test_an_intron_lies_between_the_exons(self, genome):
        path, sequence = genome
        seqs = get_intron_seqs(gtf([(101, 200, '+', 'T1'),
                                    (301, 400, '+', 'T1')]), path, 'T1')

        assert next(iter(seqs.values())) == sequence[200:300]

    def test_an_intron_excludes_the_flanking_exon_bases(self, genome):
        # exons 101-200 and 301-400 leave an intron of exactly 201-300; taking
        # the exon boundaries themselves includes one exonic base at each end
        path, sequence = genome
        seqs = get_intron_seqs(gtf([(101, 200, '+', 'T1'),
                                    (301, 400, '+', 'T1')]), path, 'T1')

        label, seq = next(iter(seqs.items()))
        assert label == 'chr1:201-300(+)'
        assert len(seq) == 100

    def test_adjacent_exons_yield_no_intron(self, genome):
        path, _ = genome
        seqs = get_intron_seqs(gtf([(101, 200, '+', 'T1'),
                                    (201, 300, '+', 'T1')]), path, 'T1')

        assert seqs == {}

    def test_a_single_base_intron_is_one_base(self, genome):
        path, sequence = genome
        seqs = get_intron_seqs(gtf([(101, 200, '+', 'T1'),
                                    (202, 300, '+', 'T1')]), path, 'T1')

        assert len(next(iter(seqs.values()))) == 1

    def test_a_minus_strand_intron_is_reverse_complemented(self, genome):
        path, sequence = genome
        seqs = get_intron_seqs(gtf([(101, 200, '-', 'T1'),
                                    (301, 400, '-', 'T1')]), path, 'T1')

        assert next(iter(seqs.values())) == rev_comp(sequence[200:300])

    def test_a_spliced_transcript_concatenates_its_exons(self, genome):
        path, sequence = genome
        spliced = get_spliced_seq(gtf([(101, 200, '+', 'T1'),
                                       (301, 400, '+', 'T1')]), path, 'T1')

        assert spliced == sequence[100:200] + sequence[300:400]


class TestLocalToGenomic:

    def test_a_plus_strand_offset_is_added_to_the_interval_start(self):
        # 1-based 101 is 0-based 100, so local 0 is genomic 100
        assert local_to_genomic('chr1:101-200(+)', 0, 30) == \
            ('chr1', 100, 130, '+')

    def test_a_plus_strand_probe_at_the_interval_end(self):
        assert local_to_genomic('chr1:101-200(+)', 70, 100) == \
            ('chr1', 170, 200, '+')

    def test_a_minus_strand_offset_is_mirrored(self):
        # mining runs on the reverse complement, so local 0 is the interval's
        # genomic END, and coordinates must still come back ascending
        assert local_to_genomic('chr1:101-200(-)', 0, 30) == \
            ('chr1', 170, 200, '-')

    def test_a_minus_strand_probe_at_the_far_end(self):
        assert local_to_genomic('chr1:101-200(-)', 70, 100) == \
            ('chr1', 100, 130, '-')

    def test_the_interval_length_is_inclusive(self):
        # a 1-base interval must map local [0,1) onto exactly one base
        assert local_to_genomic('chr1:500-500(+)', 0, 1) == \
            ('chr1', 499, 500, '+')
        assert local_to_genomic('chr1:500-500(-)', 0, 1) == \
            ('chr1', 499, 500, '-')

    def test_coordinates_stay_ascending_on_both_strands(self):
        for strand in ('+', '-'):
            _, start, stop, _ = local_to_genomic(f'chr1:101-200({strand})',
                                                 10, 40)
            assert start < stop


class TestMinedProbesMapBackToTheReference:
    """The property that a coordinate error would break silently."""

    @pytest.mark.parametrize('strand', ['+', '-'])
    def test_a_probe_mined_from_an_interval_is_found_at_its_coordinates(
            self, genome, strand):
        path, sequence = genome
        label = f'chr1:101-900({strand})'

        seqs = get_exon_seqs(gtf([(101, 900, strand, 'T1')]), path,
                             transcript_id='T1')
        interval_seq = next(iter(seqs.values()))

        probes = mine_sequence(interval_seq, seq_id=label, min_length=30,
                               max_length=37, min_tm=0, max_tm=100,
                               min_gc=0, max_gc=100)
        assert probes, 'no probes mined from the interval'

        for _, local_start, local_stop, probe_seq, _ in probes[:20]:
            _, start, stop, _ = local_to_genomic(label, local_start,
                                                 local_stop)
            reference = sequence[start:stop]
            expected = reference if strand == '+' else rev_comp(reference)
            assert probe_seq == expected


class TestFeatureSelectorContract:
    """
    transcript_id and gene_id select different record sets. Honouring one when
    both are given returns features the caller did not ask for.
    """

    def test_giving_both_selectors_is_rejected(self):
        from oligominer.bioinformatics.transcriptome.transcript_seq import (
            _select_features,
        )

        with pytest.raises(ValueError, match='both'):
            _select_features(gtf([(101, 200, '+', 'T1')]),
                             transcript_id='T1', gene_id='G1')

    def test_giving_neither_is_rejected(self):
        from oligominer.bioinformatics.transcriptome.transcript_seq import (
            _select_features,
        )

        with pytest.raises(ValueError):
            _select_features(gtf([(101, 200, '+', 'T1')]))

    def test_one_selector_is_accepted(self):
        from oligominer.bioinformatics.transcriptome.transcript_seq import (
            _select_features,
        )

        assert len(_select_features(gtf([(101, 200, '+', 'T1')]),
                                    transcript_id='T1')) == 1
