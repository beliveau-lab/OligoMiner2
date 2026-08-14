"""
# Pure-numpy k-mer counting

A k-mer count index built and queried with numpy alone, so k-mer screening works
without Jellyfish installed.

Building encodes every k-mer of a reference as a uint64, distributes them to bin
files on disk, then sorts and counts each bin independently. Peak memory is set by
the bin size rather than the genome size. Querying encodes the k-mers of each probe
and binary-searches the sorted index.

Counts are stored as uint32, so an index reports true counts for every k-mer of a
mammalian genome rather than saturating.

K-mers are counted on the forward strand only. A k-mer containing any non-ACGT base
is assigned a sentinel hash and excluded.
"""

import os
import shutil
import tempfile

import numpy as np

from oligominer.utils.cores import resolve_cores
from .exceptions import KmerIndexError

# A=0, C=1, G=2, T=3 in either case; everything else 255
BASE_ENCODE = np.full(256, 255, dtype=np.uint8)
BASE_ENCODE[[ord(b) for b in 'ACGTacgt']] = [0, 1, 2, 3, 0, 1, 2, 3]

# hash assigned to any k-mer window containing a non-ACGT base
SENTINEL = np.uint64(2 ** 63)

# counts are stored in this dtype and saturate at its maximum
COUNT_DTYPE = np.uint32

# bases read per pass when binning a chromosome
CHUNK_BP = 10_000_000

# probes hashed per batch when querying
QUERY_BATCH = 400_000


def encode_sequence(seq):
    """
    Encode a DNA sequence as a uint8 array.

    Args:
        seq (str): the DNA sequence.

    Returns:
        encoded (numpy.ndarray): uint8, A=0, C=1, G=2, T=3, any other base 255.
    """
    encoded = BASE_ENCODE[np.frombuffer(seq.encode('ascii'), dtype=np.uint8)]

    # success
    return encoded


def kmers_to_uint64(encoded, k):
    """
    Convert an encoded sequence into one uint64 hash per k-mer window.

    Accumulates one column at a time so memory stays at two arrays of n_kmers
    rather than an n-by-k matrix.

    Args:
        encoded (numpy.ndarray): uint8 encoded bases.
        k (int): k-mer length.

    Returns:
        hashes (numpy.ndarray): uint64, SENTINEL where the window contained a
            non-ACGT base.
    """
    n = len(encoded)
    if n < k:
        return np.array([], dtype=np.uint64)

    n_kmers = n - k + 1
    hashes = np.zeros(n_kmers, dtype=np.uint64)
    has_invalid = np.zeros(n_kmers, dtype=bool)

    for j in range(k):
        power = np.uint64(4 ** (k - 1 - j))
        col = encoded[j:j + n_kmers]
        has_invalid |= (col == 255)
        hashes += col.astype(np.uint64) * power

    hashes[has_invalid] = SENTINEL

    # success
    return hashes


def _encode_concat(seqs, k):
    """
    Concatenate sequences separated by k-1 N's and encode them in one call.

    The separators guarantee that any k-mer window spanning two sequences contains
    an N and is therefore sentinelled, so no boundary bookkeeping is needed.

    Args:
        seqs (list): the sequences.
        k (int): k-mer length.

    Returns:
        encoded (numpy.ndarray): uint8 encoding of the concatenated buffer.
        starts (numpy.ndarray): int64, each sequence's offset into encoded.
        lengths (numpy.ndarray): int64, each sequence's length.
    """
    sep = 'N' * (k - 1)
    lengths = np.fromiter((len(s) for s in seqs), dtype=np.int64, count=len(seqs))
    blob = sep.join(s.upper() for s in seqs)
    encoded = BASE_ENCODE[np.frombuffer(blob.encode('ascii'), dtype=np.uint8)]

    starts = np.zeros(len(lengths), dtype=np.int64)
    if len(lengths) > 1:
        starts[1:] = (np.cumsum(lengths[:-1])
                      + np.arange(1, len(lengths), dtype=np.int64) * (k - 1))

    # success
    return encoded, starts, lengths


class KmerIndex:
    """
    A sorted k-mer count index supporting maximum-count queries.

    Attributes:
        kmers (numpy.ndarray): uint64 k-mer hashes, sorted ascending.
        counts (numpy.ndarray): count per entry of kmers, parallel to it.
        k (int): the k-mer length the index was built at.
        min_count (int): k-mers rarer than this were not stored.
    """

    def __init__(self, kmers, counts, k, min_count=2):
        """
        Args:
            kmers (numpy.ndarray): sorted uint64 k-mer hashes, no sentinels.
            counts (numpy.ndarray): counts parallel to kmers.
            k (int): k-mer length.
            min_count (int): the minimum count stored when building.
        """
        self.kmers = kmers
        self.counts = counts
        self.k = int(k)
        self.min_count = int(min_count)

    def __len__(self):
        return len(self.kmers)

    def __repr__(self):
        return (f"KmerIndex(k={self.k}, entries={len(self.kmers):,}, "
                f"min_count={self.min_count})")

    @staticmethod
    def build(fasta_path, k=18, min_count=2, n_bins=256, tmp_dir=None):
        """
        Build a k-mer count index from a FASTA file.

        Args:
            fasta_path (str): path to the reference FASTA.
            k (int): k-mer length.
            min_count (int): minimum count to store. K-mers rarer than this are
                omitted, and queries report them as min_count - 1.
            n_bins (int): number of hash bins. More bins means less memory per bin.
            tmp_dir (str, optional): parent directory to create scratch under.
                A system temporary directory is used and removed when None.

        Returns:
            index (KmerIndex): the built index.
        """
        from oligominer.bioinformatics.file_io import load_fasta

        fasta = load_fasta(fasta_path)

        # a private directory per build: the bin files are opened for append,
        # so two builds sharing a directory would append into each other's bins
        if tmp_dir is not None:
            os.makedirs(tmp_dir, exist_ok=True)
        tmp_dir = tempfile.mkdtemp(prefix='om2_kmer_build_', dir=tmp_dir)

        # k-mers are binned on the top bits of the 2k-bit hash space, so the bins
        # partition that space in ascending order and concatenate already sorted
        bin_shift = 2 * k - int(np.log2(n_bins))

        try:
            bin_files = [open(os.path.join(tmp_dir, f'bin_{i:04d}.u64'), 'ab')
                         for i in range(n_bins)]
            try:
                for seq_id in fasta.keys():
                    seq_len = len(fasta[seq_id])
                    if seq_len < k:
                        continue
                    for chunk_start in range(0, seq_len, CHUNK_BP):
                        # overlap by k-1 so windows spanning a chunk edge are counted
                        end = min(chunk_start + CHUNK_BP + k - 1, seq_len)
                        chunk = str(fasta[seq_id][chunk_start:end]).upper()
                        hashes = kmers_to_uint64(encode_sequence(chunk), k)
                        hashes = hashes[hashes != SENTINEL]
                        if len(hashes) == 0:
                            continue
                        _write_bins(hashes, bin_shift, bin_files)
            finally:
                for bin_file in bin_files:
                    bin_file.close()

            result_kmers = []
            result_counts = []
            for b in range(n_bins):
                bin_path = os.path.join(tmp_dir, f'bin_{b:04d}.u64')
                if not os.path.exists(bin_path) or os.path.getsize(bin_path) == 0:
                    continue
                bin_hashes = np.fromfile(bin_path, dtype=np.uint64)
                bin_hashes.sort()
                unique, raw_counts = np.unique(bin_hashes, return_counts=True)
                keep = raw_counts >= min_count
                if keep.any():
                    result_kmers.append(unique[keep])
                    result_counts.append(
                        np.minimum(raw_counts[keep], np.iinfo(COUNT_DTYPE).max
                                   ).astype(COUNT_DTYPE))
                del bin_hashes, unique, raw_counts
        finally:
            shutil.rmtree(tmp_dir, ignore_errors=True)

        if result_kmers:
            all_kmers = np.concatenate(result_kmers)
            all_counts = np.concatenate(result_counts)
        else:
            all_kmers = np.array([], dtype=np.uint64)
            all_counts = np.array([], dtype=COUNT_DTYPE)

        # success
        return KmerIndex(all_kmers, all_counts, k, min_count=min_count)

    def save(self, path):
        """
        Save the index to a compressed .npz file.

        Args:
            path (str): destination path.

        Returns:
            path (str): the path written.
        """
        np.savez_compressed(path, kmers=self.kmers, counts=self.counts,
                            k=np.array([self.k]),
                            min_count=np.array([self.min_count]))

        # success
        return path

    @staticmethod
    def load(path):
        """
        Load an index from a .npz file.

        Args:
            path (str): the index path.

        Returns:
            index (KmerIndex): the loaded index.
        """
        data = np.load(path)
        min_count = int(data['min_count'][0]) if 'min_count' in data else 2

        # success
        return KmerIndex(data['kmers'], data['counts'], int(data['k'][0]),
                         min_count=min_count)

    def query_probes(self, probe_seqs, batch=QUERY_BATCH):
        """
        Return the maximum k-mer count within each probe sequence.

        Probes are hashed in batches, concatenated and encoded in one call per
        batch rather than one call per probe.

        Args:
            probe_seqs (list): the probe sequences.
            batch (int): probes per batch, which bounds peak memory.

        Returns:
            counts (numpy.ndarray): int64, one maximum count per probe. Probes
                shorter than k, or containing no valid k-mer, report 0.
        """
        k = self.k
        probe_seqs = list(probe_seqs)
        out = np.zeros(len(probe_seqs), dtype=np.int64)

        for start in range(0, len(probe_seqs), batch):
            chunk = probe_seqs[start:start + batch]
            encoded, offsets, _ = _encode_concat(chunk, k)
            hashes = kmers_to_uint64(encoded, k)
            if hashes.size == 0:
                continue

            positions = np.flatnonzero(hashes != SENTINEL)
            if positions.size == 0:
                continue
            valid_hashes = hashes[positions]

            # a k-mer window starting inside a probe cannot reach past it without
            # crossing a separator and being sentinelled, so its start position
            # alone identifies the probe it belongs to
            probe_of = np.searchsorted(offsets, positions, side='right') - 1

            counts = self.lookup(valid_hashes)
            np.maximum.at(out, start + probe_of, counts)

        # success
        return out

    def lookup(self, hashes):
        """
        Return the stored count for each k-mer hash.

        A hash absent from the index occurred fewer than min_count times when the
        index was built, so it is reported as min_count - 1.

        Args:
            hashes (numpy.ndarray): uint64 k-mer hashes.

        Returns:
            counts (numpy.ndarray): int64 counts, parallel to hashes.
        """
        counts = np.full(len(hashes), self.min_count - 1, dtype=np.int64)
        if len(self.kmers) == 0:
            return counts

        insert_at = np.searchsorted(self.kmers, hashes)
        in_bounds = insert_at < len(self.kmers)
        hit = np.zeros(len(hashes), dtype=bool)
        hit[in_bounds] = self.kmers[insert_at[in_bounds]] == hashes[in_bounds]
        counts[hit] = self.counts[insert_at[hit]]

        # success
        return counts


def _write_bins(hashes, bin_shift, bin_files):
    """
    Append each hash to the bin file selected by its top bits.

    Args:
        hashes (numpy.ndarray): uint64 hashes with sentinels already removed.
        bin_shift (int): right shift that maps a hash to its bin id.
        bin_files (list): open binary file handles, indexed by bin id.

    Returns:
        n (int): the number of hashes written.
    """
    bin_ids = (hashes >> np.uint64(bin_shift)).astype(np.int32)
    order = np.argsort(bin_ids, kind='stable')
    hashes_sorted = hashes[order]
    bins_sorted = bin_ids[order]

    changes = np.where(np.diff(bins_sorted))[0] + 1
    starts = np.concatenate([[0], changes])
    ends = np.concatenate([changes, [len(bins_sorted)]])

    for i in range(len(starts)):
        hashes_sorted[starts[i]:ends[i]].tofile(bin_files[bins_sorted[starts[i]]])

    # success
    return len(hashes)


def build_numpy_index(fasta_path, output_file, k=18, min_count=2, n_bins=256,
                      tmp_dir=None, cores=None):
    """
    Build a numpy k-mer index from a FASTA file and save it.

    Args:
        fasta_path (str): path to the reference FASTA.
        output_file (str): destination .npz path.
        k (int): k-mer length.
        min_count (int): minimum count to store.
        n_bins (int): number of hash bins.
        tmp_dir (str, optional): directory for temporary bin files.
        cores (int, optional): accepted for interface parity with the Jellyfish
            backend. The build is single-threaded and bounded by disk.

    Returns:
        index_path (str): the path written.
    """
    resolve_cores(cores)
    if k < 1 or 4 ** k >= 2 ** 63:
        raise KmerIndexError(
            f"k={k} is outside the range this index can hash; "
            f"k must satisfy 4**k < 2**63")

    index = KmerIndex.build(fasta_path, k=k, min_count=min_count,
                            n_bins=n_bins, tmp_dir=tmp_dir)
    index.save(output_file)

    # success
    return output_file
