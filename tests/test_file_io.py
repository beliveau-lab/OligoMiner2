"""Tests for the bioinformatics file I/O module.

Covers: FASTA loading/writing/splitting/filtering, FASTQ generation,
sequence classification, and chrom sizes.
"""

import os
import tempfile

import pytest

from oligominer.bioinformatics.file_io import (
    load_fasta,
    write_fasta,
    split_fasta,
    seqs_to_fasta,
    seqs_to_fastq,
    filter_seq_ids,
    filter_seqs,
    merge_fastas,
    classify_seq_ids,
)
from oligominer.bioinformatics.file_io.chrom_sizes import (
    get_chrom_sizes,
    get_or_create_fai,
)


# ---------------------------------------------------------------------------
# FASTA loading
# ---------------------------------------------------------------------------

class TestLoadFasta:

    def test_load(self, example_fasta_path):
        fasta = load_fasta(example_fasta_path)
        assert len(fasta.keys()) > 0

    def test_sequence_access(self, example_fasta_path):
        fasta = load_fasta(example_fasta_path)
        seq_id = list(fasta.keys())[0]
        seq_str = str(fasta[seq_id])
        assert len(seq_str) > 0
        assert set(seq_str).issubset(set("ACGTN"))


# ---------------------------------------------------------------------------
# FASTA writing
# ---------------------------------------------------------------------------

class TestWriteFasta:

    def test_write_and_reload(self, example_fasta_path, tmp_path):
        fasta = load_fasta(example_fasta_path)
        out_path = str(tmp_path / "out.fa")
        write_fasta(fasta, out_path)
        assert os.path.exists(out_path)

        reloaded = load_fasta(out_path)
        assert set(reloaded.keys()) == set(fasta.keys())

    def test_split_fasta(self, example_fasta_path, tmp_path):
        fasta = load_fasta(example_fasta_path)
        # convert to {id: str} dict as expected by split_fasta
        seqs = {k: str(fasta[k]) for k in fasta.keys()}
        out_dir = str(tmp_path / "split")
        os.makedirs(out_dir)
        paths = split_fasta(seqs, out_dir)
        assert len(paths) == len(seqs)
        for p in paths:
            assert os.path.exists(p)

    def test_merge_fastas(self, example_fasta_path, tmp_path):
        fasta = load_fasta(example_fasta_path)
        seqs = {k: str(fasta[k]) for k in fasta.keys()}

        # split then merge
        split_dir = str(tmp_path / "split")
        os.makedirs(split_dir)
        split_fasta(seqs, split_dir)

        merged_path = str(tmp_path / "merged.fa")
        merge_fastas(split_dir, merged_path)
        assert os.path.exists(merged_path)

        merged = load_fasta(merged_path)
        assert set(merged.keys()) == set(fasta.keys())


# ---------------------------------------------------------------------------
# FASTA string generation
# ---------------------------------------------------------------------------

class TestSeqsToFasta:

    def test_from_list(self):
        seqs = ["ATCG", "GCTA"]
        result = seqs_to_fasta(seqs)
        assert ">seq_0" in result
        assert "ATCG" in result
        assert ">seq_1" in result

    def test_with_custom_ids(self):
        seqs = ["ATCG"]
        ids = ["myseq"]
        result = seqs_to_fasta(seqs, seq_id_list=ids)
        assert ">myseq" in result


# ---------------------------------------------------------------------------
# FASTQ string generation
# ---------------------------------------------------------------------------

class TestSeqsToFastq:

    def test_basic(self):
        seqs = ["ATCGATCG"]
        ids = ["read1"]
        result = seqs_to_fastq(seqs, seq_id_list=ids)
        assert "@read1" in result
        assert "ATCGATCG" in result
        # quality line should be present
        assert "+" in result


# ---------------------------------------------------------------------------
# filtering
# ---------------------------------------------------------------------------

class TestFiltering:

    def test_filter_seq_ids_include(self, example_fasta_path):
        fasta = load_fasta(example_fasta_path)
        first_id = list(fasta.keys())[0]
        filtered = filter_seq_ids(fasta, incl_str=first_id)
        assert first_id in filtered

    def test_filter_seq_ids_exclude(self, example_fasta_path):
        fasta = load_fasta(example_fasta_path)
        first_id = list(fasta.keys())[0]
        filtered = filter_seq_ids(fasta, excl_str=first_id)
        assert first_id not in filtered

    def test_filter_seqs(self, example_fasta_path):
        fasta = load_fasta(example_fasta_path)
        first_id = list(fasta.keys())[0]
        filtered = filter_seqs(fasta, incl_str=first_id)
        assert len(list(filtered.keys())) >= 1


# ---------------------------------------------------------------------------
# classification
# ---------------------------------------------------------------------------

class TestClassify:

    def test_classify_seq_ids(self, example_fasta_path):
        fasta = load_fasta(example_fasta_path)
        df = classify_seq_ids(fasta)
        assert 'category' in df.columns
        assert len(df) == len(fasta.keys())


# ---------------------------------------------------------------------------
# chrom sizes
# ---------------------------------------------------------------------------

class TestChromSizes:

    def test_get_chrom_sizes(self, example_fasta_path):
        sizes = get_chrom_sizes(example_fasta_path)
        assert isinstance(sizes, dict)
        assert len(sizes) > 0
        for name, size in sizes.items():
            assert isinstance(size, int)
            assert size > 0

    def test_get_or_create_fai(self, example_fasta_path):
        fai_path = get_or_create_fai(example_fasta_path)
        assert os.path.exists(fai_path)
