"""Tests for the RNA capabilities: junctions, discriminating regions, transcriptome.

A probe spanning an exon-exon junction has no contiguous genomic target, so it
can only be designed against the spliced transcript and only screened against a
transcriptome. These tests use a synthetic two-isoform gene where the exon
structure, and therefore every expected coordinate, is known exactly.
"""

import pandas as pd
import pytest

from oligominer.bioinformatics.transcriptome import (
    build_transcriptome,
    discriminating_regions,
    exon_order,
    junction_offsets,
    junction_probes,
    junction_windows,
    spans_junction,
    transcript_ids,
    transcript_lengths,
)


@pytest.fixture
def genome(tmp_path):
    """A single random chromosome.

    Random rather than periodic: a junction window joins two sequences that are
    not adjacent in the genome, and in a periodic genome that joined sequence
    would occur by chance, so the test that it does not could not fail.
    """
    import random
    random.seed(11)
    seq = ''.join(random.choice('ACGT') for _ in range(2000))
    path = tmp_path / 'genome.fa'
    with open(path, 'w') as handle:
        handle.write('>chr1\n')
        for i in range(0, len(seq), 60):
            handle.write(seq[i:i + 60] + '\n')

    # success
    return path, seq


@pytest.fixture
def annotation():
    """A gene with two isoforms sharing exon 1 and differing after it.

    tx1: exons 101-200, 301-400, 501-600
    tx2: exons 101-200,            501-600, 701-800
    So 301-400 is unique to tx1 and 701-800 is unique to tx2.
    """
    rows = []
    for transcript, spans in (('tx1', [(101, 200), (301, 400), (501, 600)]),
                              ('tx2', [(101, 200), (501, 600), (701, 800)])):
        for start, end in spans:
            rows.append({'seqid': 'chr1', 'feature': 'exon', 'start': start,
                         'end': end, 'strand': '+', 'gene_id': 'geneA',
                         'transcript_id': transcript})

    # success
    return pd.DataFrame(rows)


class TestExonOrder:

    def test_plus_strand_exons_are_ascending(self, annotation):
        exons = exon_order(annotation, 'tx1')
        assert list(exons['start']) == [101, 301, 501]

    def test_minus_strand_exons_are_descending(self, annotation):
        minus = annotation.copy()
        minus['strand'] = '-'
        exons = exon_order(minus, 'tx1')
        assert list(exons['start']) == [501, 301, 101]


class TestJunctionOffsets:

    def test_offsets_fall_at_cumulative_exon_lengths(self, annotation):
        # exons are 100 bases each, so junctions sit at 100 and 200
        assert junction_offsets(annotation, 'tx1') == [100, 200]

    def test_a_transcript_has_one_fewer_junction_than_exons(self, annotation):
        exons = exon_order(annotation, 'tx1')
        assert len(junction_offsets(annotation, 'tx1')) == len(exons) - 1

    def test_a_single_exon_transcript_has_no_junctions(self):
        single = pd.DataFrame([{'seqid': 'chr1', 'feature': 'exon', 'start': 1,
                                'end': 100, 'strand': '+', 'gene_id': 'g',
                                'transcript_id': 'only'}])
        assert junction_offsets(single, 'only') == []


class TestJunctionWindows:

    def test_one_window_per_junction(self, genome, annotation):
        path, _ = genome
        windows = junction_windows(annotation, str(path), 'tx1', flank=20)
        assert len(windows) == 2

    def test_a_window_is_two_flanks_wide(self, genome, annotation):
        path, _ = genome
        windows = junction_windows(annotation, str(path), 'tx1', flank=20)
        assert (windows['seq'].str.len() == 40).all()

    def test_the_junction_sits_in_the_middle_of_the_window(self, genome, annotation):
        path, _ = genome
        windows = junction_windows(annotation, str(path), 'tx1', flank=20)
        assert (windows['junction_offset'] == 20).all()

    def test_the_window_joins_sequence_that_is_not_genomically_adjacent(
            self, genome, annotation):
        """This is what makes a junction probe transcript-specific."""
        path, raw = genome
        windows = junction_windows(annotation, str(path), 'tx1', flank=20)
        window = windows.iloc[0]['seq']

        # the window's two halves come from exon 1 and exon 2
        left, right = window[:20], window[20:]
        assert left == raw[180:200]     # last 20 bases of exon 1 (101-200, 1-based)
        assert right == raw[300:320]    # first 20 bases of exon 2 (301-400)
        # and the joined sequence does not occur in the genome
        assert window not in raw

    def test_a_flank_larger_than_the_exon_is_clipped(self, genome, annotation):
        path, _ = genome
        windows = junction_windows(annotation, str(path), 'tx1', flank=500)
        assert (windows['seq'].str.len() > 0).all()


class TestSpansJunction:

    def test_a_probe_covering_the_junction_spans_it(self):
        assert spans_junction(10, 30, junction_offset=20)

    def test_a_probe_entirely_left_of_the_junction_does_not(self):
        assert not spans_junction(0, 19, junction_offset=20)

    def test_a_probe_entirely_right_of_the_junction_does_not(self):
        assert not spans_junction(21, 40, junction_offset=20)

    def test_a_minimum_overhang_is_enforced(self):
        assert spans_junction(15, 25, junction_offset=20, min_overhang=5)
        assert not spans_junction(18, 25, junction_offset=20, min_overhang=5)

    def test_probes_are_filtered_to_those_that_span(self):
        probes = [('w', 0, 15, 'A' * 15, 42.0),      # left only
                  ('w', 10, 30, 'A' * 20, 42.0),     # spans
                  ('w', 25, 40, 'A' * 15, 42.0)]     # right only
        kept = junction_probes(probes, junction_offset=20)
        assert [p[1] for p in kept] == [10]


class TestDiscriminatingRegions:

    def test_the_unique_exon_of_each_isoform_is_found(self, annotation):
        regions = discriminating_regions(annotation, gene_id='geneA')
        by_transcript = {row['transcript_id']: (row['start'], row['end'])
                         for _, row in regions.iterrows()}
        assert by_transcript['tx1'] == (301, 400)
        assert by_transcript['tx2'] == (701, 800)

    def test_shared_exons_are_not_reported(self, annotation):
        regions = discriminating_regions(annotation, gene_id='geneA')
        starts = set(regions['start'])
        assert 101 not in starts     # shared exon 1
        assert 501 not in starts     # shared exon 3

    def test_a_single_isoform_gene_has_no_discriminating_regions(self):
        single = pd.DataFrame([{'seqid': 'chr1', 'feature': 'exon', 'start': 1,
                                'end': 100, 'strand': '+', 'gene_id': 'g',
                                'transcript_id': 'only'}])
        assert len(discriminating_regions(single)) == 0

    def test_short_regions_can_be_filtered_out(self, annotation):
        assert len(discriminating_regions(annotation, min_length=500)) == 0

    def test_the_reported_length_matches_the_interval(self, annotation):
        regions = discriminating_regions(annotation, gene_id='geneA')
        assert (regions['length'] == regions['end'] - regions['start'] + 1).all()


class TestBuildTranscriptome:

    def test_every_transcript_is_written(self, genome, annotation, tmp_path):
        path, _ = genome
        out = tmp_path / 'tx.fa'
        info = build_transcriptome(annotation, str(path), str(out))
        assert info['n_written'] == 2

    def test_records_are_named_by_transcript_id(self, genome, annotation, tmp_path):
        path, _ = genome
        out = tmp_path / 'tx.fa'
        build_transcriptome(annotation, str(path), str(out))
        names = [line[1:].strip() for line in out.read_text().splitlines()
                 if line.startswith('>')]
        assert names == ['tx1', 'tx2']

    def test_sequence_length_matches_the_spliced_length(self, genome, annotation,
                                                        tmp_path):
        path, _ = genome
        out = tmp_path / 'tx.fa'
        build_transcriptome(annotation, str(path), str(out))

        lengths, current = {}, None
        for line in out.read_text().splitlines():
            if line.startswith('>'):
                current = line[1:].strip()
                lengths[current] = 0
            else:
                lengths[current] += len(line.strip())

        assert lengths == transcript_lengths(annotation)

    def test_the_transcriptome_can_be_indexed_and_searched(self, genome,
                                                           annotation, tmp_path):
        """A junction probe must align to the transcriptome and not the genome."""
        import shutil
        import subprocess

        if shutil.which('bowtie2-build') is None:
            pytest.skip('bowtie2 not installed')

        path, _ = genome
        out = tmp_path / 'tx.fa'
        build_transcriptome(annotation, str(path), str(out))

        index = str(tmp_path / 'txidx')
        subprocess.run(['bowtie2-build', '-q', str(out), index], check=True)

        window = junction_windows(annotation, str(path), 'tx1', flank=20).iloc[0]
        probe = window['seq'][5:35]

        aligned = subprocess.run(
            ['bowtie2', '-x', index, '-c', '-U', probe, '--no-unal', '--xeq'],
            capture_output=True, text=True)
        hits = [line for line in aligned.stdout.splitlines()
                if line and not line.startswith('@')]
        assert hits, 'junction probe did not align to the transcriptome'

    def test_transcript_ids_are_listed(self, annotation):
        assert transcript_ids(annotation) == ['tx1', 'tx2']

    def test_transcript_ids_can_be_restricted_to_a_gene(self, annotation):
        assert transcript_ids(annotation, gene_id='geneA') == ['tx1', 'tx2']
