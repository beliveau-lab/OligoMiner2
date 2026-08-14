"""Tests for k-mer index construction, query and backend dispatch.

Covers the numpy index end to end, the dispatch rules, the metadata guards that
reject a mismatched query, and agreement between the two backends on the same
reference. The Jellyfish tests skip when the binary is not installed.
"""

import json
import random

import numpy as np
import pytest

from oligominer.specificity.kmers import (
    KmerIndex,
    build_index,
    have_jellyfish,
    max_kmer,
    read_metadata,
    resolve_backend,
)
from oligominer.specificity.kmers.exceptions import KmerIndexError
from oligominer.specificity.kmers.numpy_index import (
    SENTINEL,
    encode_sequence,
    kmers_to_uint64,
)

K = 8

needs_jellyfish = pytest.mark.skipif(
    have_jellyfish() is None, reason='jellyfish binary not installed'
)


@pytest.fixture
def reference(tmp_path):
    """A small reference FASTA with a sequence repeated a known number of times."""
    random.seed(0)
    filler = ''.join(random.choice('ACGT') for _ in range(2000))
    repeat = 'ACGTACGTACGTACGTACGT'

    # the repeat unit appears 5 times, separated by unique filler
    parts = []
    for i in range(5):
        parts.append(filler[i * 300:(i + 1) * 300])
        parts.append(repeat)
    seq = ''.join(parts)

    path = tmp_path / 'ref.fa'
    path.write_text('>chr1\n' + '\n'.join(seq[i:i + 60] for i in range(0, len(seq), 60)) + '\n')

    # success
    return path, seq, repeat


@pytest.fixture
def numpy_index(tmp_path, reference):
    """A built numpy index over the small reference."""
    fasta, _, _ = reference
    out = tmp_path / 'ref.npz'
    build_index(str(fasta), str(out), k=K, backend='numpy', min_count=2)

    # success
    return out


class TestEncoding:

    def test_bases_encode_to_their_codes(self):
        assert list(encode_sequence('ACGT')) == [0, 1, 2, 3]

    def test_encoding_folds_case(self):
        assert list(encode_sequence('acgt')) == list(encode_sequence('ACGT'))

    def test_non_acgt_encodes_to_255(self):
        assert encode_sequence('N')[0] == 255

    def test_kmer_hashes_are_distinct_for_distinct_kmers(self):
        hashes = kmers_to_uint64(encode_sequence('ACGTACGA'), 4)
        assert len(set(hashes.tolist())) == len(hashes)

    def test_identical_kmers_hash_identically(self):
        a = kmers_to_uint64(encode_sequence('ACGTAC'), 4)
        b = kmers_to_uint64(encode_sequence('TTACGTACTT'), 4)
        assert a[0] in b.tolist()

    def test_window_containing_n_is_sentinelled(self):
        hashes = kmers_to_uint64(encode_sequence('ACNGTACG'), 4)
        assert hashes[0] == SENTINEL
        assert hashes[-1] != SENTINEL

    def test_sequence_shorter_than_k_yields_no_kmers(self):
        assert len(kmers_to_uint64(encode_sequence('ACG'), 8)) == 0


class TestNumpyIndex:

    def test_index_builds_and_reloads_identically(self, numpy_index):
        index = KmerIndex.load(str(numpy_index))
        assert index.k == K
        assert len(index) > 0
        assert np.array_equal(np.sort(index.kmers), index.kmers)

    def test_repeated_kmer_reports_its_true_count(self, numpy_index, reference):
        _, _, repeat = reference
        counts = max_kmer(numpy_index, [repeat], k=K, backend='numpy')
        assert counts[0] >= 5

    def test_unique_sequence_reports_a_low_count(self, numpy_index, reference):
        _, seq, _ = reference
        unique_window = seq[100:100 + 40]
        counts = max_kmer(numpy_index, [unique_window], k=K, backend='numpy')
        assert counts[0] <= 2

    def test_probe_shorter_than_k_reports_zero(self, numpy_index):
        counts = max_kmer(numpy_index, ['ACG'], k=K, backend='numpy')
        assert counts[0] == 0

    def test_probe_of_all_n_reports_zero(self, numpy_index):
        counts = max_kmer(numpy_index, ['N' * 40], k=K, backend='numpy')
        assert counts[0] == 0

    def test_one_count_per_input_sequence(self, numpy_index, reference):
        _, seq, repeat = reference
        seqs = [repeat, seq[100:140], 'ACG', 'N' * 30, seq[500:540]]
        counts = max_kmer(numpy_index, seqs, k=K, backend='numpy')
        assert len(counts) == len(seqs)

    def test_batching_does_not_change_the_answer(self, numpy_index, reference):
        _, seq, repeat = reference
        seqs = [seq[i:i + 40] for i in range(0, 800, 20)] + [repeat]
        index = KmerIndex.load(str(numpy_index))
        assert np.array_equal(
            index.query_probes(seqs, batch=3), index.query_probes(seqs, batch=10_000)
        )

    def test_counts_do_not_saturate_below_the_dtype_maximum(self, tmp_path):
        """A k-mer repeated more than 255 times reports its true count."""
        unit = 'ACGTACGTAC'
        seq = ''.join(unit for _ in range(400))
        fasta = tmp_path / 'rep.fa'
        fasta.write_text('>c\n' + seq + '\n')
        out = tmp_path / 'rep.npz'
        build_index(str(fasta), str(out), k=10, backend='numpy', min_count=2)

        counts = max_kmer(out, [unit], k=10, backend='numpy')
        assert counts[0] > 255


class TestDispatch:

    def test_npz_suffix_selects_numpy(self):
        assert resolve_backend('index.npz') == 'numpy'

    def test_explicit_numpy_is_honored_for_any_suffix(self):
        assert resolve_backend('index.jf', backend='numpy') == 'numpy'

    def test_unknown_backend_raises(self):
        with pytest.raises(KmerIndexError, match='unknown backend'):
            resolve_backend('index.npz', backend='bowtie')

    @pytest.mark.skipif(have_jellyfish() is not None,
                        reason='jellyfish is installed')
    def test_requesting_jellyfish_without_the_binary_raises(self):
        with pytest.raises(KmerIndexError, match='not.*on PATH'):
            resolve_backend('index.jf', backend='jellyfish')

    @needs_jellyfish
    def test_jf_suffix_selects_jellyfish(self):
        assert resolve_backend('index.jf') == 'jellyfish'

    def test_missing_index_raises(self, tmp_path):
        with pytest.raises(KmerIndexError, match='no such k-mer index'):
            max_kmer(tmp_path / 'absent.npz', ['ACGTACGTACGT'], k=K)


class TestMetadataGuards:

    def test_metadata_is_written_beside_the_index(self, numpy_index):
        meta = read_metadata(numpy_index)
        assert meta['k'] == K
        assert meta['is_canonical'] is False
        assert meta['backend'] == 'numpy'

    def test_querying_at_the_wrong_k_raises(self, numpy_index):
        with pytest.raises(KmerIndexError, match='built at k='):
            max_kmer(numpy_index, ['ACGTACGTACGTACGT'], k=K + 1, backend='numpy')

    def test_canonicality_mismatch_raises(self, numpy_index):
        with pytest.raises(KmerIndexError, match='is_canonical'):
            max_kmer(numpy_index, ['ACGTACGTACGTACGT'], k=K,
                     backend='numpy', expect_canonical=True)

    def test_matching_canonicality_is_accepted(self, numpy_index):
        counts = max_kmer(numpy_index, ['ACGTACGTACGTACGT'], k=K,
                          backend='numpy', expect_canonical=False)
        assert len(counts) == 1

    def test_query_works_without_a_sidecar(self, numpy_index):
        """An index with no sidecar is queryable, just unchecked."""
        sidecar = numpy_index.with_name(numpy_index.name + '.om2meta.json')
        sidecar.unlink()
        counts = max_kmer(numpy_index, ['ACGTACGTACGTACGT'], k=K, backend='numpy')
        assert len(counts) == 1

    def test_numpy_backend_refuses_canonical_builds(self, tmp_path, reference):
        fasta, _, _ = reference
        with pytest.raises(KmerIndexError, match='forward strand'):
            build_index(str(fasta), str(tmp_path / 'c.npz'), k=K,
                        backend='numpy', canonical=True)


@needs_jellyfish
class TestBackendAgreement:
    """Both backends must return the same counts for the same reference and probes."""

    def test_backends_agree_on_every_probe(self, tmp_path, reference):
        fasta, seq, repeat = reference

        npz = tmp_path / 'ref.npz'
        jf = tmp_path / 'ref.jf'
        build_index(str(fasta), str(npz), k=K, backend='numpy', min_count=1)
        build_index(str(fasta), str(jf), k=K, backend='jellyfish')

        probes = [seq[i:i + 40] for i in range(0, len(seq) - 40, 37)] + [repeat]

        from_numpy = max_kmer(npz, probes, k=K, backend='numpy')
        from_jellyfish = max_kmer(jf, probes, k=K, backend='jellyfish')

        mismatches = np.flatnonzero(from_numpy != from_jellyfish)
        assert mismatches.size == 0, (
            f'{mismatches.size}/{len(probes)} probes disagree, '
            f'first at {mismatches[:5].tolist()}: '
            f'numpy={from_numpy[mismatches[:5]].tolist()} '
            f'jellyfish={from_jellyfish[mismatches[:5]].tolist()}'
        )


class TestConcurrentBuildScratch:
    """Builds run concurrently and must not share their bin files."""

    def test_concurrent_builds_sharing_a_tmp_dir_agree_with_a_lone_build(
            self, tmp_path, reference):
        import concurrent.futures

        fasta, _, _ = reference

        alone = tmp_path / 'alone.npz'
        build_index(str(fasta), str(alone), k=K, backend='numpy', min_count=2)
        expected = KmerIndex.load(str(alone))

        shared = tmp_path / 'scratch'
        shared.mkdir()

        def build(name):
            out = tmp_path / f'{name}.npz'
            KmerIndex.build(str(fasta), k=K, min_count=2,
                            tmp_dir=str(shared)).save(str(out))
            return KmerIndex.load(str(out))

        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
            built = list(pool.map(build, ['a', 'b', 'c']))

        # the bin files are opened for append, so a shared directory would have
        # each build's k-mers accumulate into the others
        for index in built:
            assert len(index.kmers) == len(expected.kmers)
            assert index.counts.max() == expected.counts.max()

    def test_the_scratch_directory_is_removed(self, tmp_path, reference):
        fasta, _, _ = reference
        shared = tmp_path / 'scratch'
        shared.mkdir()

        KmerIndex.build(str(fasta), k=K, min_count=2, tmp_dir=str(shared))

        # the shared parent survives; only this build's own directory goes
        assert shared.exists()
        assert list(shared.iterdir()) == []
