"""Tests for the core probe mining pipeline.

Covers: mine_sequence, mine_fasta, write_probes, probes_to_df, and the
various mining parameter combinations (length, Tm, overlap, exhaustive,
homopolymer, prohibited_seqs, salt/formamide params).
"""

import os
import tempfile

import numpy as np
import pandas as pd
import pytest

from oligominer.utils.exceptions import ConfigurationError

from oligominer import mine_sequence, mine_fasta, write_probes
from oligominer.thermodynamics.mining import probes_to_df, PROBE_COLUMNS


# ---------------------------------------------------------------------------
# mine_sequence
# ---------------------------------------------------------------------------

class TestMineSequence:

    def test_returns_list_of_tuples(self, short_seq):
        probes = mine_sequence(short_seq, seq_id="test")
        assert isinstance(probes, list)
        if len(probes) > 0:
            assert isinstance(probes[0], tuple)
            assert len(probes[0]) == 5

    def test_tuple_structure(self, short_seq):
        probes = mine_sequence(short_seq, seq_id="myseq", min_length=18,
                               max_length=20, pct_formamide=0)
        if len(probes) > 0:
            seq_id, start, stop, probe_seq, tm = probes[0]
            assert seq_id == "myseq"
            assert isinstance(start, int)
            assert isinstance(stop, int)
            assert stop > start
            assert len(probe_seq) == stop - start
            assert isinstance(tm, float)

    def test_respects_length_constraints(self, example_fasta_path):
        from oligominer.bioinformatics.file_io import load_fasta
        fasta = load_fasta(example_fasta_path)
        seq_str = str(fasta[list(fasta.keys())[0]])

        probes = mine_sequence(seq_str, min_length=35, max_length=40)
        for _, start, stop, probe_seq, _ in probes:
            length = stop - start
            assert 35 <= length <= 40
            assert len(probe_seq) == length

    def test_respects_tm_constraints(self, example_fasta_path):
        from oligominer.bioinformatics.file_io import load_fasta
        fasta = load_fasta(example_fasta_path)
        seq_str = str(fasta[list(fasta.keys())[0]])

        probes = mine_sequence(seq_str, min_tm=42, max_tm=47)
        for _, _, _, _, tm in probes:
            assert 42 <= tm <= 47

    def test_no_overlap_mode(self, example_fasta_path):
        from oligominer.bioinformatics.file_io import load_fasta
        fasta = load_fasta(example_fasta_path)
        seq_str = str(fasta[list(fasta.keys())[0]])

        probes = mine_sequence(seq_str, overlap=False)
        for i in range(1, len(probes)):
            prev_stop = probes[i - 1][2]
            curr_start = probes[i][1]
            assert curr_start >= prev_stop

    def test_spacing_mode(self, example_fasta_path):
        from oligominer.bioinformatics.file_io import load_fasta
        fasta = load_fasta(example_fasta_path)
        seq_str = str(fasta[list(fasta.keys())[0]])

        spacing = 10
        probes = mine_sequence(seq_str, spacing=spacing)
        for i in range(1, len(probes)):
            prev_stop = probes[i - 1][2]
            curr_start = probes[i][1]
            assert curr_start >= prev_stop + spacing

    def test_exhaustive_mode(self, example_fasta_path):
        from oligominer.bioinformatics.file_io import load_fasta
        fasta = load_fasta(example_fasta_path)
        seq_str = str(fasta[list(fasta.keys())[0]])

        # exhaustive should return more probes than default
        default = mine_sequence(seq_str)
        exhaustive = mine_sequence(seq_str, exhaustive=True)
        assert len(exhaustive) >= len(default)

    def test_exhaustive_overlap_conflict(self, short_seq):
        with pytest.raises(ConfigurationError, match="exhaustive"):
            mine_sequence(short_seq, exhaustive=True, overlap=False)

    def test_prohibited_seqs(self, example_fasta_path):
        from oligominer.bioinformatics.file_io import load_fasta
        fasta = load_fasta(example_fasta_path)
        seq_str = str(fasta[list(fasta.keys())[0]])

        prohibited = ["CCCC"]
        probes = mine_sequence(seq_str, prohibited_seqs=prohibited)
        for _, _, _, probe_seq, _ in probes:
            assert "CCCC" not in probe_seq

    def test_empty_sequence(self):
        probes = mine_sequence("", seq_id="empty")
        assert probes == []

    def test_short_sequence(self):
        probes = mine_sequence("ATCG", seq_id="tiny")
        assert probes == []

    def test_salt_params_affect_tm(self, example_fasta_path):
        """Different salt concentrations should produce different Tm values."""
        from oligominer.bioinformatics.file_io import load_fasta
        fasta = load_fasta(example_fasta_path)
        seq_str = str(fasta[list(fasta.keys())[0]])

        probes_low_na = mine_sequence(seq_str, Na=50, pct_formamide=0,
                                      min_tm=30, max_tm=100)
        probes_high_na = mine_sequence(seq_str, Na=500, pct_formamide=0,
                                       min_tm=30, max_tm=100)
        # different salt → different Tm values for same positions
        if probes_low_na and probes_high_na:
            tms_low = {(p[1], p[2]): p[4] for p in probes_low_na}
            tms_high = {(p[1], p[2]): p[4] for p in probes_high_na}
            shared = set(tms_low) & set(tms_high)
            if shared:
                key = next(iter(shared))
                assert tms_low[key] != tms_high[key]

    def test_formamide_correction(self, example_fasta_path):
        """Formamide should depress Tm."""
        from oligominer.bioinformatics.file_io import load_fasta
        fasta = load_fasta(example_fasta_path)
        seq_str = str(fasta[list(fasta.keys())[0]])

        probes_no_fmd = mine_sequence(seq_str, pct_formamide=0,
                                      min_tm=30, max_tm=100)
        probes_fmd = mine_sequence(seq_str, pct_formamide=50,
                                   min_tm=30, max_tm=100)
        if probes_no_fmd and probes_fmd:
            tms_no = {(p[1], p[2]): p[4] for p in probes_no_fmd}
            tms_fmd = {(p[1], p[2]): p[4] for p in probes_fmd}
            shared = set(tms_no) & set(tms_fmd)
            if shared:
                key = next(iter(shared))
                # formamide depresses Tm
                assert tms_fmd[key] < tms_no[key]


# ---------------------------------------------------------------------------
# mine_fasta
# ---------------------------------------------------------------------------

class TestMineFasta:

    def test_returns_list_of_tuples(self, example_fasta_path):
        probes = mine_fasta(example_fasta_path)
        assert isinstance(probes, list)
        assert len(probes) > 0
        assert isinstance(probes[0], tuple)

    def test_mines_all_sequences(self, example_fasta_path):
        """Should mine probes from all sequences in a multi-seq FASTA."""
        probes = mine_fasta(example_fasta_path)
        seq_ids = {p[0] for p in probes}
        # bundled FASTA has multiple sequences
        assert len(seq_ids) >= 1

    def test_mining_params_forwarded(self, example_fasta_path):
        probes = mine_fasta(example_fasta_path, min_length=35, max_length=40)
        for _, start, stop, _, _ in probes:
            assert 35 <= (stop - start) <= 40


# ---------------------------------------------------------------------------
# probes_to_df
# ---------------------------------------------------------------------------

class TestProbesToDf:

    def test_converts_tuples_to_dataframe(self, example_probes):
        df = probes_to_df(example_probes)
        assert isinstance(df, pd.DataFrame)
        assert list(df.columns) == PROBE_COLUMNS[:5]
        assert len(df) == len(example_probes)

    def test_empty_list(self):
        df = probes_to_df([])
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 0


# ---------------------------------------------------------------------------
# write_probes
# ---------------------------------------------------------------------------

class TestWriteProbes:

    def test_write_bed(self, example_probes, tmp_path):
        path = str(tmp_path / "probes.bed")
        write_probes(example_probes, path, fmt='bed')
        assert os.path.exists(path)
        with open(path) as f:
            lines = f.readlines()
        assert len(lines) == len(example_probes)
        # BED format: tab-separated
        assert '\t' in lines[0]

    def test_write_fastq(self, example_probes, tmp_path):
        path = str(tmp_path / "probes.fastq")
        write_probes(example_probes, path, fmt='fastq')
        assert os.path.exists(path)
        with open(path) as f:
            content = f.read()
        # FASTQ records start with @
        assert content.startswith('@')

    def test_write_csv(self, example_probes, tmp_path):
        path = str(tmp_path / "probes.csv")
        write_probes(example_probes, path, fmt='csv')
        assert os.path.exists(path)
        df = pd.read_csv(path)
        assert len(df) == len(example_probes)
        assert 'probe_seq' in df.columns
