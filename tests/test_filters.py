"""Tests for the low-complexity filter and interval exclusion.

Covers the dinucleotide entropy filter across sequence classes, and BED-based
exclusion including the overlap cases that a naive containment check misses.
"""

import numpy as np
import pandas as pd
import pytest

from oligominer.probe_design import exclude_intervals, overlaps_intervals, read_bed
from oligominer.thermodynamics.mining import mine_sequence

MINING = dict(seq_id='t', min_length=30, max_length=30, min_tm=0, max_tm=200,
              min_gc=0, max_gc=100, max_homopolymer=None)

RANDOM_SEQ = 'ACGTTGCAAGCTTAGCATCGGATCAGTCAGGCATTAGCCAGTAC' * 2
DINUC_REPEAT = 'AT' * 25
HOMOPOLYMER = 'A' * 45


class TestEntropyFilter:

    def test_random_sequence_survives(self):
        assert len(mine_sequence(RANDOM_SEQ, min_entropy=0.5, **MINING)) > 0

    def test_a_dinucleotide_repeat_is_rejected(self):
        assert mine_sequence(DINUC_REPEAT, min_entropy=0.5, **MINING) == []

    def test_a_homopolymer_is_rejected(self):
        assert mine_sequence(HOMOPOLYMER, min_entropy=0.5, **MINING) == []

    def test_the_filter_is_off_by_default(self):
        assert len(mine_sequence(DINUC_REPEAT, **MINING)) > 0

    def test_a_threshold_of_zero_rejects_nothing(self):
        assert (len(mine_sequence(DINUC_REPEAT, min_entropy=0.0, **MINING))
                == len(mine_sequence(DINUC_REPEAT, **MINING)))

    def test_raising_the_threshold_is_monotonic(self):
        counts = [len(mine_sequence(RANDOM_SEQ, min_entropy=e, **MINING))
                  for e in (0.0, 0.3, 0.6, 0.9, 1.0)]
        assert counts == sorted(counts, reverse=True)

    def test_a_threshold_above_one_rejects_everything(self):
        assert mine_sequence(RANDOM_SEQ, min_entropy=1.01, **MINING) == []

    def test_random_scores_higher_than_a_repeat(self):
        """The filter must separate the two classes, not merely reject both."""
        for threshold in (0.4, 0.5, 0.6):
            assert len(mine_sequence(RANDOM_SEQ, min_entropy=threshold, **MINING)) > 0
            assert mine_sequence(DINUC_REPEAT, min_entropy=threshold, **MINING) == []


class TestReadBed:

    def test_three_columns_are_read(self, tmp_path):
        path = tmp_path / 'x.bed'
        path.write_text('chr1\t100\t200\nchr2\t50\t75\n')
        intervals = read_bed(path)
        assert len(intervals) == 2
        assert list(intervals.columns) == ['chrom', 'start', 'stop']

    def test_extra_columns_are_ignored(self, tmp_path):
        path = tmp_path / 'x.bed'
        path.write_text('chr1\t100\t200\tname\t0\t+\n')
        assert read_bed(path)['stop'].iloc[0] == 200

    def test_headers_and_comments_are_skipped(self, tmp_path):
        path = tmp_path / 'x.bed'
        path.write_text('track name=x\n# a comment\nchr1\t100\t200\n')
        assert len(read_bed(path)) == 1

    def test_an_empty_file_gives_an_empty_frame(self, tmp_path):
        path = tmp_path / 'x.bed'
        path.write_text('')
        assert len(read_bed(path)) == 0


class TestOverlap:

    def _probes(self, spans, chrom='chr1'):
        return pd.DataFrame([{'seq_id': chrom, 'start': s, 'stop': e}
                             for s, e in spans])

    def _intervals(self, spans, chrom='chr1'):
        return pd.DataFrame([{'chrom': chrom, 'start': s, 'stop': e}
                             for s, e in spans])

    def test_a_probe_inside_an_interval_overlaps(self):
        hit = overlaps_intervals(self._probes([(120, 150)]), self._intervals([(100, 200)]))
        assert hit.tolist() == [True]

    def test_a_probe_before_an_interval_does_not(self):
        hit = overlaps_intervals(self._probes([(10, 50)]), self._intervals([(100, 200)]))
        assert hit.tolist() == [False]

    def test_a_probe_after_an_interval_does_not(self):
        hit = overlaps_intervals(self._probes([(300, 350)]), self._intervals([(100, 200)]))
        assert hit.tolist() == [False]

    def test_a_probe_straddling_the_start_overlaps(self):
        """The probe begins outside the interval and runs into it."""
        hit = overlaps_intervals(self._probes([(80, 120)]), self._intervals([(100, 200)]))
        assert hit.tolist() == [True]

    def test_a_probe_straddling_the_end_overlaps(self):
        hit = overlaps_intervals(self._probes([(180, 220)]), self._intervals([(100, 200)]))
        assert hit.tolist() == [True]

    def test_a_probe_containing_an_interval_overlaps(self):
        hit = overlaps_intervals(self._probes([(50, 300)]), self._intervals([(100, 200)]))
        assert hit.tolist() == [True]

    def test_touching_at_the_boundary_does_not_overlap(self):
        """BED is half-open, so a probe ending where an interval starts is clear."""
        hit = overlaps_intervals(self._probes([(50, 100)]), self._intervals([(100, 200)]))
        assert hit.tolist() == [False]

    def test_a_different_chromosome_does_not_overlap(self):
        probes = self._probes([(120, 150)], chrom='chr2')
        hit = overlaps_intervals(probes, self._intervals([(100, 200)], chrom='chr1'))
        assert hit.tolist() == [False]

    def test_overlapping_exclusion_intervals_are_merged(self):
        intervals = self._intervals([(100, 200), (150, 300)])
        hit = overlaps_intervals(self._probes([(250, 260)]), intervals)
        assert hit.tolist() == [True]

    def test_unsorted_intervals_are_handled(self):
        intervals = self._intervals([(500, 600), (100, 200), (300, 400)])
        probes = self._probes([(150, 160), (250, 260), (350, 360), (550, 560)])
        assert overlaps_intervals(probes, intervals).tolist() == [True, False, True, True]

    def test_no_intervals_means_no_overlap(self):
        hit = overlaps_intervals(self._probes([(1, 2)]),
                                 pd.DataFrame(columns=['chrom', 'start', 'stop']))
        assert hit.tolist() == [False]


class TestExcludeIntervals:

    def test_overlapping_probes_are_dropped_and_counted(self, tmp_path):
        path = tmp_path / 'x.bed'
        path.write_text('chr1\t100\t200\n')
        probes = pd.DataFrame([
            {'seq_id': 'chr1', 'start': 10, 'stop': 40},
            {'seq_id': 'chr1', 'start': 120, 'stop': 150},
            {'seq_id': 'chr1', 'start': 500, 'stop': 530},
        ])
        kept, n_dropped = exclude_intervals(probes, path)
        assert n_dropped == 1
        assert list(kept['start']) == [10, 500]

    def test_nothing_is_dropped_when_nothing_overlaps(self, tmp_path):
        path = tmp_path / 'x.bed'
        path.write_text('chr9\t100\t200\n')
        probes = pd.DataFrame([{'seq_id': 'chr1', 'start': 120, 'stop': 150}])
        kept, n_dropped = exclude_intervals(probes, path)
        assert n_dropped == 0
        assert len(kept) == 1

    def test_other_columns_are_preserved(self, tmp_path):
        path = tmp_path / 'x.bed'
        path.write_text('chr1\t100\t200\n')
        probes = pd.DataFrame([{'seq_id': 'chr1', 'start': 10, 'stop': 40,
                                'probe_seq': 'ACGT', 'tm': 42.0}])
        kept, _ = exclude_intervals(probes, path)
        assert kept['probe_seq'].iloc[0] == 'ACGT'
        assert kept['tm'].iloc[0] == 42.0
