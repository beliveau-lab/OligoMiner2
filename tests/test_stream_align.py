"""Tests for streaming alignment.

Covers FASTQ writing, the shared bowtie2 command builder, and that the streaming
path produces byte-identical BED to the buffered path it replaces. Tests needing
an aligner skip when bowtie2 is not installed.
"""

import random
import shutil

import pandas as pd
import pytest

from oligominer.specificity.alignment import (
    align_to_bed,
    bowtie_build,
    build_bowtie2_cmd,
    write_fastq,
)

needs_bowtie2 = pytest.mark.skipif(
    shutil.which('bowtie2') is None, reason='bowtie2 not installed'
)


@pytest.fixture
def genome(tmp_path):
    """A single-chromosome reference."""
    random.seed(3)
    seq = ''.join(random.choice('ACGT') for _ in range(20000))
    path = tmp_path / 'ref.fa'
    with open(path, 'w') as handle:
        handle.write('>chr1\n')
        for i in range(0, len(seq), 60):
            handle.write(seq[i:i + 60] + '\n')

    # success
    return path, seq


@pytest.fixture
def probes(genome):
    """Probes taken from the reference, so they align."""
    _, seq = genome
    rows = [{'seqid': f'p{i}', 'probe_seq': seq[i * 400:i * 400 + 36]}
            for i in range(20)]

    # success
    return pd.DataFrame(rows)


@pytest.fixture
def bt2_index(tmp_path, genome):
    """A bowtie2 index over the reference."""
    path, _ = genome
    prefix = str(tmp_path / 'idx')
    bowtie_build(str(path), prefix)

    # success
    return prefix


class TestWriteFastq:

    def test_record_count(self, tmp_path):
        out = tmp_path / 'p.fastq'
        n = write_fastq(['ACGT', 'TTTT'], ['a', 'b'], out)
        assert n == 2

    def test_record_format(self, tmp_path):
        out = tmp_path / 'p.fastq'
        write_fastq(['ACGTA'], ['probe1'], out)
        assert out.read_text() == '@probe1\nACGTA\n+\nIIIII\n'

    def test_quality_string_matches_sequence_length(self, tmp_path):
        out = tmp_path / 'p.fastq'
        write_fastq(['ACGTACGTAC'], ['x'], out)
        lines = out.read_text().splitlines()
        assert len(lines[3]) == len(lines[1])


class TestCommandBuilder:

    def test_index_and_input_are_present(self):
        cmd = build_bowtie2_cmd('/tmp/idx', input_file='/tmp/p.fastq')
        assert cmd[0] == 'bowtie2'
        assert '-x' in cmd
        assert '-U' in cmd

    def test_reporting_depth_is_passed(self):
        cmd = build_bowtie2_cmd('/tmp/idx', input_file='/tmp/p.fastq', k=100)
        assert cmd[cmd.index('-k') + 1] == '100'

    def test_preset_overrides_seed_parameters(self):
        cmd = build_bowtie2_cmd('/tmp/idx', input_file='/tmp/p.fastq',
                                preset={'D': 20, 'R': 3, 'N': 1, 'L': 20})
        assert cmd[cmd.index('-N') + 1] == '1'
        assert cmd[cmd.index('-L') + 1] == '20'

    def test_no_unal_is_emitted(self):
        cmd = build_bowtie2_cmd('/tmp/idx', input_file='/tmp/p.fastq', no_unal=True)
        assert '--no-unal' in cmd

    def test_threads_flag_present_when_above_one(self):
        cmd = build_bowtie2_cmd('/tmp/idx', input_file='/tmp/p.fastq', threads=4)
        assert cmd[cmd.index('-p') + 1] == '4'


@needs_bowtie2
class TestStreamingAlignment:

    def test_bed_is_written(self, probes, bt2_index, tmp_path):
        out = tmp_path / 'out.bed'
        info = align_to_bed(probes, bt2_index, out, threads=1, k=10)
        assert out.exists()
        assert info['n_reads'] == len(probes)
        assert info['n_rows'] > 0

    def test_every_probe_aligns_somewhere(self, probes, bt2_index, tmp_path):
        out = tmp_path / 'out.bed'
        align_to_bed(probes, bt2_index, out, threads=1, k=10)
        aligned = {line.split('\t')[3] for line in out.read_text().splitlines()}
        assert aligned == set(probes['seqid'])

    def test_bed_rows_have_the_expected_field_count(self, probes, bt2_index, tmp_path):
        out = tmp_path / 'out.bed'
        align_to_bed(probes, bt2_index, out, threads=1, k=10)
        for line in out.read_text().splitlines():
            assert len(line.split('\t')) == 7

    def test_coordinates_are_non_negative(self, probes, bt2_index, tmp_path):
        out = tmp_path / 'out.bed'
        align_to_bed(probes, bt2_index, out, threads=1, k=10)
        for line in out.read_text().splitlines():
            fields = line.split('\t')
            assert int(fields[1]) >= 0
            assert int(fields[2]) > int(fields[1])

    def test_keep_fastq_retains_the_intermediate(self, probes, bt2_index, tmp_path):
        out = tmp_path / 'out.bed'
        info = align_to_bed(probes, bt2_index, out, threads=1, k=10, keep_fastq=True)
        assert info['fastq'] is not None
        from pathlib import Path
        assert Path(info['fastq']).exists()

    def test_a_bad_index_raises_rather_than_writing_an_empty_bed(self, probes, tmp_path):
        from oligominer.utils.exceptions import ExternalCommandFailed
        with pytest.raises(ExternalCommandFailed):
            align_to_bed(probes, str(tmp_path / 'absent_index'),
                         tmp_path / 'out.bed', threads=1)

    def test_streaming_matches_the_buffered_path(self, probes, bt2_index, tmp_path):
        """The streamed BED must be identical to what the buffered path produces."""
        from oligominer.specificity.alignment import bam_to_bed, bowtie_align

        out = tmp_path / 'streamed.bed'
        align_to_bed(probes, bt2_index, out, threads=1, k=10, reorder=True)
        streamed = sorted(line for line in out.read_text().splitlines() if line)

        fastq = tmp_path / 'buffered.fastq'
        write_fastq(probes['probe_seq'], probes['seqid'], fastq)
        sam = bowtie_align(bt2_index, input_file=str(fastq), threads=1, k=10,
                           no_unal=True, reorder=True)
        buffered = sorted(
            line for line in bam_to_bed(bam_data=sam).strip().split('\n') if line
        )

        assert streamed == buffered
