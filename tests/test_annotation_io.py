"""Tests for writing, splitting and merging annotation files."""

import os

import pandas as pd
import pytest

from oligominer.bioinformatics.file_io.bed_io import bed_to_df
from oligominer.bioinformatics.file_io.gtf_io import (
    merge_annotation_beds, split_gtf, write_bed, write_gtf,
)
from oligominer.bioinformatics.file_io.exceptions import EmptyExportError


@pytest.fixture
def annotation():
    """A parsed annotation spanning two chromosomes."""
    return pd.DataFrame([
        {'seqid': 'chr1', 'source': 'test', 'type': 'exon', 'start': 100,
         'end': 200, 'score': '.', 'strand': '+', 'phase': '.',
         'transcript_id': 'T1', 'transcript_id_full': 'T1.2',
         'gene_id': 'G1'},
        {'seqid': 'chr1', 'source': 'test', 'type': 'exon', 'start': 300,
         'end': 400, 'score': '.', 'strand': '+', 'phase': '.',
         'transcript_id': 'T1', 'transcript_id_full': 'T1.2',
         'gene_id': 'G1'},
        {'seqid': 'chr2', 'source': 'test', 'type': 'exon', 'start': 500,
         'end': 600, 'score': '.', 'strand': '-', 'phase': '.',
         'transcript_id': 'T2', 'transcript_id_full': 'T2.1',
         'gene_id': 'G2'},
    ])


class TestWriteGtf:

    def test_writes_every_record_without_a_header(self, annotation, tmp_path):
        path = tmp_path / 'out.gtf'
        write_gtf(annotation, str(path))

        written = pd.read_csv(path, sep='\t', header=None)
        assert len(written) == 3

    def test_creates_the_parent_directory(self, annotation, tmp_path):
        path = tmp_path / 'nested' / 'out.gtf'
        write_gtf(annotation, str(path))

        assert path.exists()

    def test_an_empty_annotation_raises_rather_than_writing_nothing(
            self, tmp_path):
        with pytest.raises(EmptyExportError):
            write_gtf(pd.DataFrame(), str(tmp_path / 'out.gtf'))


class TestWriteBed:

    def test_writes_the_bed_columns_in_order(self, annotation, tmp_path):
        path = tmp_path / 'out.bed'
        write_bed(annotation, str(path))

        written = pd.read_csv(path, sep='\t', header=None)
        assert written[0].tolist() == ['chr1', 'chr1', 'chr2']
        assert written[1].tolist() == [100, 300, 500]
        assert written[3].tolist() == ['T1', 'T1', 'T2']

    def test_the_source_and_phase_columns_are_not_written(self, annotation,
                                                          tmp_path):
        path = tmp_path / 'out.bed'
        write_bed(annotation, str(path))

        written = pd.read_csv(path, sep='\t', header=None)
        assert len(written.columns) == 8

    def test_missing_optional_columns_are_skipped(self, annotation, tmp_path):
        path = tmp_path / 'out.bed'
        write_bed(annotation.drop(columns=['transcript_id_full']), str(path))

        written = pd.read_csv(path, sep='\t', header=None)
        assert len(written.columns) == 7

    def test_an_empty_annotation_raises(self, tmp_path):
        with pytest.raises(EmptyExportError):
            write_bed(pd.DataFrame(), str(tmp_path / 'out.bed'))


class TestSplitGtf:

    def test_one_file_per_chromosome(self, annotation, tmp_path):
        paths = split_gtf(annotation, str(tmp_path / 'split'))

        assert len(paths) == 2
        assert {os.path.basename(p) for p in paths} == {
            'chr1_filtered_gtf.tsv', 'chr2_filtered_gtf.tsv'}

    def test_a_file_holds_only_its_own_chromosome(self, annotation, tmp_path):
        split_gtf(annotation, str(tmp_path / 'split'))

        written = pd.read_csv(tmp_path / 'split' / 'chr1_filtered_gtf.tsv',
                              sep='\t')
        assert set(written['seqid']) == {'chr1'}
        assert len(written) == 2

    def test_the_suffix_is_configurable(self, annotation, tmp_path):
        paths = split_gtf(annotation, str(tmp_path / 'split'), suffix='.tsv')

        assert all(p.endswith('.tsv') for p in paths)
        assert os.path.basename(paths[0]) == 'chr1.tsv'

    def test_an_empty_annotation_raises(self, tmp_path):
        with pytest.raises(EmptyExportError):
            split_gtf(pd.DataFrame(), str(tmp_path / 'split'))


class TestMergeAnnotationBeds:

    def test_concatenates_every_bed_in_the_directory(self, tmp_path):
        source = tmp_path / 'beds'
        source.mkdir()
        (source / 'chr1.bed').write_text('chr1\t100\t200\n')
        (source / 'chr2.bed').write_text('chr2\t500\t600\n')

        out = tmp_path / 'merged.bed'
        merge_annotation_beds(str(source), str(out))

        assert len(out.read_text().strip().splitlines()) == 2

    def test_files_of_other_extensions_are_ignored(self, tmp_path):
        source = tmp_path / 'beds'
        source.mkdir()
        (source / 'chr1.bed').write_text('chr1\t100\t200\n')
        (source / 'notes.txt').write_text('ignore me\n')

        out = tmp_path / 'merged.bed'
        merge_annotation_beds(str(source), str(out))

        assert 'ignore me' not in out.read_text()

    def test_an_empty_directory_raises(self, tmp_path):
        source = tmp_path / 'beds'
        source.mkdir()

        with pytest.raises(EmptyExportError):
            merge_annotation_beds(str(source), str(tmp_path / 'merged.bed'))


class TestBedToDf:

    def test_parses_tab_separated_bed_text(self):
        df = bed_to_df('chr1\t100\t200\tprobe1\nchr2\t300\t400\tprobe2\n')

        assert len(df) == 2
        assert df[0].tolist() == ['chr1', 'chr2']
        assert df[1].tolist() == [100, 300]

    def test_columns_are_positional_not_named(self):
        # BED has no header, so callers index by position
        df = bed_to_df('chr1\t100\t200\n')
        assert list(df.columns) == [0, 1, 2]
