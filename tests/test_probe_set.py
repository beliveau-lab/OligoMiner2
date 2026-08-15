"""Tests for the ProbeSet class.

Covers: construction from tuples and DataFrames, from_fasta, from_csv,
export methods (to_csv, to_bed, to_fastq), CSV round-trip, and the
seqid column auto-generation.
"""

import os

import shutil

import pandas as pd
import pytest

from oligominer import ProbeSet, mine_fasta
from oligominer.thermodynamics.mining import probes_to_df
from oligominer.utils.exceptions import PipelineStateError


# ---------------------------------------------------------------------------
# construction
# ---------------------------------------------------------------------------

class TestProbeSetConstruction:

    def test_from_tuples(self, example_probes):
        ps = ProbeSet(example_probes)
        assert isinstance(ps.df, pd.DataFrame)
        assert len(ps) == len(example_probes)
        assert 'seqid' in ps.df.columns

    def test_from_dataframe(self, example_probes):
        df = probes_to_df(example_probes)
        ps = ProbeSet(df)
        assert len(ps) == len(df)
        assert 'seqid' in ps.df.columns

    def test_from_dataframe_with_seqid(self, example_probes):
        """If seqid column already exists, it should not be overwritten."""
        df = probes_to_df(example_probes)
        df['seqid'] = 'custom_id'
        ps = ProbeSet(df)
        assert ps.df['seqid'].iloc[0] == 'custom_id'

    def test_df_is_a_copy(self, example_probes):
        """Modifying the original data should not affect the ProbeSet."""
        df = probes_to_df(example_probes)
        ps = ProbeSet(df)
        df['probe_seq'] = 'AAAA'
        assert ps.df['probe_seq'].iloc[0] != 'AAAA'

    def test_from_fasta(self, example_fasta_path):
        ps = ProbeSet.from_fasta(example_fasta_path)
        assert len(ps) > 0
        assert 'seqid' in ps.df.columns

    def test_empty_tuples(self):
        ps = ProbeSet([])
        assert len(ps) == 0
        assert isinstance(ps.df, pd.DataFrame)

    def test_repr(self, example_probe_set):
        r = repr(example_probe_set)
        assert 'ProbeSet' in r
        assert 'probes=' in r


# ---------------------------------------------------------------------------
# export
# ---------------------------------------------------------------------------

class TestProbeSetExport:

    def test_to_csv(self, example_probe_set, tmp_path):
        path = str(tmp_path / "probes.csv")
        example_probe_set.to_csv(path)
        assert os.path.exists(path)
        df = pd.read_csv(path)
        assert len(df) == len(example_probe_set)

    def test_to_bed(self, example_probe_set, tmp_path):
        path = str(tmp_path / "probes.bed")
        example_probe_set.to_bed(path)
        assert os.path.exists(path)
        with open(path) as f:
            lines = f.readlines()
        assert len(lines) == len(example_probe_set)

    def test_to_fastq(self, example_probe_set, tmp_path):
        path = str(tmp_path / "probes.fastq")
        example_probe_set.to_fastq(path)
        assert os.path.exists(path)
        with open(path) as f:
            content = f.read()
        assert content.startswith('@')


# ---------------------------------------------------------------------------
# CSV round-trip
# ---------------------------------------------------------------------------

class TestProbeSetRoundTrip:

    def test_csv_round_trip(self, example_probe_set, tmp_path):
        path = str(tmp_path / "rt.csv")
        example_probe_set.to_csv(path)
        ps2 = ProbeSet.from_csv(path)
        assert len(ps2) == len(example_probe_set)
        assert list(ps2.df.columns) == list(example_probe_set.df.columns)

    def test_round_trip_preserves_data(self, example_probe_set, tmp_path):
        path = str(tmp_path / "rt.csv")
        example_probe_set.to_csv(path)
        ps2 = ProbeSet.from_csv(path)
        # probe sequences should be identical
        assert list(ps2.df['probe_seq']) == list(example_probe_set.df['probe_seq'])


# ---------------------------------------------------------------------------
# initial state
# ---------------------------------------------------------------------------

class TestProbeSetInitialState:

    def test_align_df_is_none(self, example_probe_set):
        assert example_probe_set.align_df is None

    def test_merged_df_is_none(self, example_probe_set):
        assert example_probe_set.merged_df is None

    def test_score_df_is_none(self, example_probe_set):
        assert example_probe_set.score_df is None

    def test_merge_without_align_raises(self, example_probe_set):
        with pytest.raises(PipelineStateError, match="align"):
            example_probe_set.merge()


class TestKmerProvenance:
    """
    A forward-only index counts a probe's k-mers on one strand; a canonical one
    counts both. The same probe gets a different answer, so two probe sets
    built against indexes of different canonicality are not comparable and the
    manifest has to say which was used.
    """

    @pytest.fixture
    def reference(self, tmp_path):
        path = tmp_path / 'ref.fa'
        path.write_text('>chr1\n' + 'ACGTACGTACGTACGTACGT' * 6 + '\n')
        return str(path)

    def _probe_set(self):
        return ProbeSet(pd.DataFrame({
            'seq_id': ['chr1', 'chr1'],
            'start': [0, 20],
            'stop': [20, 40],
            'probe_seq': ['ACGTACGTACGTACGTACGT'] * 2,
            'tm': [45.0, 45.0],
        }))

    def test_a_forward_only_index_is_recorded_as_such(self, reference,
                                                       tmp_path):
        from oligominer.specificity.kmers import build_index

        index = tmp_path / 'forward.npz'
        build_index(reference, str(index), k=18, backend='numpy',
                    min_count=1, canonical=False)

        probes = self._probe_set()
        probes.compute_max_kmer(str(index), k=18)

        stage = [s for s in probes.manifest['stages']
                 if s['stage'] == 'max_kmer'][-1]
        assert stage['params']['is_canonical'] is False

    @pytest.mark.skipif(shutil.which('jellyfish') is None,
                        reason='jellyfish is required for a canonical index')
    def test_a_canonical_index_is_recorded_as_such(self, reference, tmp_path):
        # the numpy backend counts one strand and refuses canonical outright,
        # so a canonical index has to come from jellyfish
        from oligominer.specificity.kmers import build_index

        index = tmp_path / 'canonical.jf'
        build_index(reference, str(index), k=18, backend='jellyfish',
                    min_count=1, canonical=True)

        probes = self._probe_set()
        probes.compute_max_kmer(str(index), k=18)

        stage = [s for s in probes.manifest['stages']
                 if s['stage'] == 'max_kmer'][-1]
        assert stage['params']['is_canonical'] is True

    def test_the_numpy_backend_refuses_to_pretend_it_is_canonical(self,
                                                                  reference,
                                                                  tmp_path):
        # forward-only counts labelled canonical would be a wrong number with a
        # correct-looking label
        from oligominer.specificity.kmers import build_index
        from oligominer.specificity.kmers.exceptions import KmerIndexError

        with pytest.raises(KmerIndexError, match='forward strand'):
            build_index(reference, str(tmp_path / 'x.npz'), k=18,
                        backend='numpy', min_count=1, canonical=True)

    def test_an_index_without_metadata_records_it_as_unknown(self, reference,
                                                             tmp_path):
        # a jellyfish index built outside the package has no sidecar, and
        # guessing its canonicality would be worse than recording that we do
        # not know
        from oligominer.specificity.kmers import build_index, sidecar_path

        index = tmp_path / 'bare.npz'
        build_index(reference, str(index), k=18, backend='numpy', min_count=1)
        sidecar_path(str(index)).unlink()

        probes = self._probe_set()
        probes.compute_max_kmer(str(index), k=18)

        stage = [s for s in probes.manifest['stages']
                 if s['stage'] == 'max_kmer'][-1]
        assert stage['params']['is_canonical'] is None
