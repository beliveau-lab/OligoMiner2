"""Tests for turning alignment records into a duplex-ready frame."""

import shutil
import subprocess

import pytest

from oligominer.specificity.alignment.process_alignments import (
    process_alignments,
)
from oligominer.utils.exceptions import InvalidInputError

SAM_HEADER = '@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:60\n'

# one forward and one reverse alignment of two probes
SAM_RECORDS = (
    'probe1\t0\tchr1\t1\t42\t8=\t*\t0\t0\tACGTACGT\t~~~~~~~~\tAS:i:16\n'
    'probe2\t16\tchr1\t21\t42\t8=\t*\t0\t0\tTTTTGGGG\t~~~~~~~~\tAS:i:14\n'
)

SAM_DATA = SAM_HEADER + SAM_RECORDS

EXPECTED_COLUMNS = ['align_seqid', 'align_start', 'align_stop', 'seqid',
                    'align_score', 'align_strand', 'align_cigar']


@pytest.fixture
def reference(tmp_path):
    """A short reference the alignments sit inside."""
    path = tmp_path / 'ref.fa'
    path.write_text('>chr1\n' + 'ACGTACGT' + 'A' * 12 + 'ccccgggg'
                    + 'T' * 32 + '\n')
    return str(path)


class TestProcessAlignments:

    def test_requires_exactly_one_input(self):
        with pytest.raises(InvalidInputError):
            process_alignments()

        with pytest.raises(InvalidInputError):
            process_alignments(sam_data=SAM_DATA, bam_path='x.bam')

    def test_header_lines_are_not_alignments(self):
        # a SAM read from disk carries its header; the spec forbids a read name
        # starting with '@', so header lines are safe to skip and must be
        with_header = process_alignments(sam_data=SAM_DATA)
        without_header = process_alignments(sam_data=SAM_RECORDS)

        assert len(with_header) == len(without_header) == 2

    def test_one_row_per_alignment(self):
        df = process_alignments(sam_data=SAM_DATA)
        assert len(df) == 2

    def test_columns_are_named(self):
        df = process_alignments(sam_data=SAM_DATA)
        assert list(df.columns)[:len(EXPECTED_COLUMNS)] == EXPECTED_COLUMNS

    def test_the_probe_name_is_carried_through(self):
        df = process_alignments(sam_data=SAM_DATA)
        assert set(df['seqid']) == {'probe1', 'probe2'}

    def test_strand_comes_from_the_sam_flag(self):
        df = process_alignments(sam_data=SAM_DATA).set_index('seqid')
        assert df.loc['probe1', 'align_strand'] == '+'
        assert df.loc['probe2', 'align_strand'] == '-'

    def test_no_derived_sequence_without_a_reference(self):
        assert 'derived_seq' not in process_alignments(sam_data=SAM_DATA)

    def test_a_reference_adds_the_derived_sequence(self, reference):
        df = process_alignments(sam_data=SAM_DATA, ref_fasta=reference)

        assert 'derived_seq' in df
        assert df['derived_seq'].str.len().tolist() == [8, 8]

    def test_derived_sequences_are_uppercased_by_default(self, reference):
        # the reference is soft-masked where probe2 aligns
        df = process_alignments(sam_data=SAM_DATA,
                                ref_fasta=reference).set_index('seqid')
        assert df.loc['probe2', 'derived_seq'].isupper()

    def test_case_is_preserved_when_asked(self, reference):
        df = process_alignments(sam_data=SAM_DATA, ref_fasta=reference,
                                to_upper=False).set_index('seqid')
        assert not df.loc['probe2', 'derived_seq'].isupper()


@pytest.mark.skipif(shutil.which('samtools') is None,
                    reason='samtools is not installed')
class TestFromBam:

    def test_a_bam_gives_the_same_frame_as_its_sam(self, tmp_path):
        sam = tmp_path / 'aln.sam'
        sam.write_text(SAM_DATA)

        bam = tmp_path / 'aln.bam'
        with open(bam, 'wb') as handle:
            subprocess.run(['samtools', 'view', '-bS', str(sam)],
                           stdout=handle, check=True)

        from_sam = process_alignments(sam_data=SAM_DATA)
        from_bam = process_alignments(bam_path=str(bam))

        assert from_bam.equals(from_sam)
