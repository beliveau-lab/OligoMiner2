"""Tests for FASTQ formatting and SAM/BAM loading."""

import shutil
import subprocess

import pytest

from oligominer.bioinformatics.file_io.fastq_io import seqs_to_fastq
from oligominer.bioinformatics.file_io.sam_bam_io import (
    load_bam_file, load_sam_file,
)

SAM_TEXT = (
    '@HD\tVN:1.0\tSO:unsorted\n'
    '@SQ\tSN:chr1\tLN:100\n'
    'probe1\t0\tchr1\t1\t42\t4=\t*\t0\t0\tACGT\t~~~~\n'
    'probe2\t16\tchr1\t10\t42\t4=\t*\t0\t0\tTGCA\t~~~~\n'
)


class TestSeqsToFastq:

    def test_each_sequence_becomes_a_four_line_record(self):
        text = seqs_to_fastq(['ACGT', 'TTTT'])
        assert len(text.splitlines()) == 8

    def test_the_quality_line_matches_the_sequence_length(self):
        lines = seqs_to_fastq(['ACGTACGT']).splitlines()
        assert lines[2] == '+'
        assert len(lines[3]) == len(lines[1]) == 8

    def test_identifiers_are_used_when_given(self):
        text = seqs_to_fastq(['ACGT', 'TTTT'], ['chr1:1-4', 'chr1:10-14'])
        assert text.startswith('@chr1:1-4\n')
        assert '@chr1:10-14\n' in text

    def test_identifiers_default_to_a_positional_name(self):
        text = seqs_to_fastq(['ACGT', 'TTTT'])
        assert '@seq_0\n' in text
        assert '@seq_1\n' in text

    def test_fewer_identifiers_than_sequences_truncates(self):
        # zip stops at the shorter list, so a short id list silently drops
        # sequences rather than mislabelling them
        text = seqs_to_fastq(['ACGT', 'TTTT'], ['only-one'])
        assert len(text.splitlines()) == 4

    def test_no_sequences_yields_an_empty_string(self):
        assert seqs_to_fastq([]) == ''


class TestLoadSam:

    def test_reads_the_file_verbatim(self, tmp_path):
        path = tmp_path / 'aln.sam'
        path.write_text(SAM_TEXT)

        assert load_sam_file(str(path)) == SAM_TEXT

    def test_a_missing_file_raises(self, tmp_path):
        with pytest.raises(Exception):
            load_sam_file(str(tmp_path / 'absent.sam'))


@pytest.mark.skipif(shutil.which('samtools') is None,
                    reason='samtools is not installed')
class TestLoadBam:

    @pytest.fixture
    def bam(self, tmp_path):
        sam = tmp_path / 'aln.sam'
        sam.write_text(SAM_TEXT)

        path = tmp_path / 'aln.bam'
        with open(path, 'wb') as handle:
            subprocess.run(['samtools', 'view', '-bS', str(sam)],
                           stdout=handle, check=True)
        return str(path)

    def test_returns_the_alignment_records(self, bam):
        text = load_bam_file(bam)

        assert 'probe1' in text
        assert 'probe2' in text

    def test_the_header_is_not_included(self, bam):
        # samtools view without -h emits records only
        assert '@SQ' not in load_bam_file(bam)

    def test_one_line_per_record(self, bam):
        assert len([line for line in load_bam_file(bam).splitlines()
                    if line.strip()]) == 2
