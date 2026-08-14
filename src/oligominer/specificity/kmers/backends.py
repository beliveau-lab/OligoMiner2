"""
# K-mer backend dispatch

One k-mer screening interface served by two interchangeable backends: Jellyfish
when it is installed, and the pure-numpy index otherwise.

    counts = max_kmer(index_path, probe_seqs, k=18)
    counts = max_kmer(index_path, probe_seqs, k=18, backend='jellyfish')

Backends are selected in this order:

1. an explicit ``backend`` argument, which raises if that backend is unusable
2. the index file's suffix, ``.jf`` for Jellyfish and ``.npz`` for numpy
3. Jellyfish if its binary is on PATH, otherwise numpy

Two properties of an index cannot be inferred from the index file itself and are
recorded in a sidecar written beside it: the k it was built at, and whether it
counts canonical k-mers. Querying against a mismatch in either returns plausible
numbers rather than an error, so both are checked and raise.
"""

import json
import shutil
import subprocess
import threading
from pathlib import Path

import numpy as np

from oligominer.utils.cores import resolve_cores
from .exceptions import KmerIndexError
from .numpy_index import KmerIndex, build_numpy_index

JELLYFISH_SUFFIXES = {'.jf'}
NUMPY_SUFFIXES = {'.npz'}

SIDECAR_SUFFIX = '.om2meta.json'


def have_jellyfish():
    """
    Return the path to the jellyfish binary, or None if it is not installed.

    Returns:
        path (str or None): the executable path, or None.
    """
    # success
    return shutil.which('jellyfish')


def resolve_backend(index_path, backend='auto'):
    """
    Decide which backend answers for an index.

    Args:
        index_path (str or pathlib.Path): the k-mer index path.
        backend (str): 'auto', 'jellyfish' or 'numpy'.

    Returns:
        name (str): 'jellyfish' or 'numpy'.

    Raises:
        KmerIndexError: if an explicitly requested backend cannot be used, or the
            backend name is not recognized.
    """
    suffix = Path(index_path).suffix.lower()

    if backend == 'jellyfish':
        if have_jellyfish() is None:
            raise KmerIndexError(
                "backend='jellyfish' was requested but the jellyfish binary is not "
                "on PATH; install it or use backend='numpy'")
        name = 'jellyfish'
    elif backend == 'numpy':
        name = 'numpy'
    elif backend == 'auto':
        if suffix in JELLYFISH_SUFFIXES:
            name = 'jellyfish' if have_jellyfish() else 'numpy'
        elif suffix in NUMPY_SUFFIXES:
            name = 'numpy'
        else:
            name = 'jellyfish' if have_jellyfish() else 'numpy'
    else:
        raise KmerIndexError(
            f"unknown backend {backend!r}; expected 'auto', 'jellyfish' or 'numpy'")

    # success
    return name


def sidecar_path(index_path):
    """
    Return the metadata sidecar path for an index.

    Args:
        index_path (str or pathlib.Path): the index path.

    Returns:
        path (pathlib.Path): the sidecar path.
    """
    # success
    return Path(str(index_path) + SIDECAR_SUFFIX)


def write_metadata(index_path, info):
    """
    Write an index's metadata sidecar.

    Args:
        index_path (str or pathlib.Path): the index path.
        info (dict): the metadata to record.

    Returns:
        path (pathlib.Path): the sidecar written.
    """
    path = sidecar_path(index_path)
    path.write_text(json.dumps(info, indent=1))

    # success
    return path


def read_metadata(index_path):
    """
    Read an index's metadata sidecar.

    Args:
        index_path (str or pathlib.Path): the index path.

    Returns:
        info (dict or None): the sidecar contents, or None when no sidecar exists.
    """
    path = sidecar_path(index_path)
    info = json.loads(path.read_text()) if path.exists() else None

    # success
    return info


def build_index(fasta_path, output_file, k=18, backend='auto', cores=None,
                min_count=2, size=None, canonical=False, verbose=False):
    """
    Build a k-mer index with whichever backend is selected.

    Writes a metadata sidecar recording k and canonicality beside the index.

    Args:
        fasta_path (str): path to the reference FASTA.
        output_file (str): destination index path. Under backend='auto' the
            suffix selects the backend: '.jf' for Jellyfish, '.npz' for numpy.
        k (int): k-mer length.
        backend (str): 'auto', 'jellyfish' or 'numpy'.
        cores (int, optional): worker threads. None resolves the scheduler grant.
        min_count (int): numpy backend only. K-mers rarer than this are not stored.
        size (str, optional): Jellyfish hash size, e.g. '3300M'. Sized from the
            input FASTA when None.
        canonical (bool): count a k-mer and its reverse complement together.
        verbose (bool): print the backend command.

    Returns:
        info (dict): backend, path, k, is_canonical and min_count.

    Raises:
        KmerIndexError: if the backend cannot build the requested index.
    """
    fasta_path, output_file = Path(fasta_path), Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    name = resolve_backend(output_file, backend)
    cores = resolve_cores(cores)

    if name == 'jellyfish':
        from .jellyfish_build import jellyfish_build

        if size is None:
            size = _size_hash_from_fasta(fasta_path)
        jellyfish_build(str(fasta_path), str(output_file), k=k, size=size,
                        cores=cores, canonical=canonical, verbose=verbose)
    else:
        if canonical:
            raise KmerIndexError(
                'the numpy backend counts the forward strand only; '
                'canonical=True requires backend="jellyfish"')
        build_numpy_index(str(fasta_path), str(output_file), k=k,
                          min_count=min_count, cores=cores)

    info = {
        'backend': name,
        'path': str(output_file),
        'k': int(k),
        'is_canonical': bool(canonical),
        'min_count': int(min_count) if name == 'numpy' else None,
    }
    write_metadata(output_file, info)

    # success
    return info


def _size_hash_from_fasta(fasta_path):
    """
    Choose a Jellyfish hash size from the size of the input FASTA.

    Args:
        fasta_path (pathlib.Path): the reference FASTA.

    Returns:
        size (str): a Jellyfish size string such as '512M'.
    """
    n_bases = fasta_path.stat().st_size

    # success
    return f"{max(16, int(n_bases / 1e6) * 2)}M"


def max_kmer(index_path, seqs, k=18, backend='auto', expect_canonical=None,
             verbose=False):
    """
    Return the maximum k-mer count within each sequence.

    A probe whose most abundant k-mer occurs once in the reference has no
    repetitive content at that k, so lower counts indicate higher specificity.

    Args:
        index_path (str or pathlib.Path): a Jellyfish '.jf' or numpy '.npz' index.
        seqs (list): the sequences to query.
        k (int): k-mer length. Must match the k the index was built at.
        backend (str): 'auto', 'jellyfish' or 'numpy'.
        expect_canonical (bool, optional): when set, raise unless the index was
            built with that canonicality.
        verbose (bool): print the backend and index being used.

    Returns:
        counts (numpy.ndarray): int64, one maximum count per input sequence.

    Raises:
        KmerIndexError: if the index is missing, or its recorded k or canonicality
            disagrees with the query.
    """
    index_path = Path(index_path)
    if not index_path.exists():
        raise KmerIndexError(f'no such k-mer index: {index_path}')

    _check_metadata(index_path, k, expect_canonical)

    name = resolve_backend(index_path, backend)
    if verbose:
        print(f'  max_kmer via {name} on {index_path.name} ({len(seqs)} seqs)')

    if name == 'jellyfish':
        counts = _max_kmer_jellyfish(index_path, seqs, k)
    else:
        counts = _max_kmer_numpy(index_path, seqs, k)

    # success
    return counts


def _check_metadata(index_path, k, expect_canonical):
    """
    Verify a query's k and canonicality against the index's recorded metadata.

    Args:
        index_path (pathlib.Path): the index path.
        k (int): the k the caller intends to query at.
        expect_canonical (bool or None): the canonicality the caller expects, or
            None to accept either.

    Returns:
        checked (bool): True when a sidecar was present and agreed, False when
            there was no sidecar to check against.

    Raises:
        KmerIndexError: on a k or canonicality mismatch.
    """
    meta = read_metadata(index_path)
    if meta is None:
        return False

    if int(meta.get('k', k)) != int(k):
        raise KmerIndexError(
            f"index {index_path.name} was built at k={meta['k']} but queried at "
            f"k={k}; the counts would not correspond to the requested k-mers")

    if (expect_canonical is not None
            and bool(meta.get('is_canonical')) != bool(expect_canonical)):
        raise KmerIndexError(
            f"index {index_path.name} has is_canonical={meta.get('is_canonical')} "
            f"but the caller expects {expect_canonical}; the two conventions "
            f"disagree on which probes pass a given cutoff")

    # success
    return True


def _max_kmer_numpy(index_path, seqs, k):
    """
    Query a numpy index for the maximum k-mer count of each sequence.

    Args:
        index_path (pathlib.Path): the '.npz' index.
        seqs (list): the sequences to query.
        k (int): k-mer length.

    Returns:
        counts (numpy.ndarray): int64 maximum counts.

    Raises:
        KmerIndexError: if the index was built at a different k.
    """
    index = KmerIndex.load(str(index_path))
    if index.k != int(k):
        raise KmerIndexError(
            f'numpy index was built at k={index.k} but queried at k={k}')

    # success
    return index.query_probes(list(seqs))


def _max_kmer_jellyfish(index_path, seqs, k):
    """
    Query a Jellyfish index for the maximum k-mer count of each sequence.

    Runs one subprocess for the whole batch and consumes its output line by line,
    so neither the k-mer stream nor the result stream is held in memory.

    `jellyfish query` prints '<kmer> <count>' when k-mers are passed as arguments
    and '<count>' alone when they arrive on stdin under -i, so the count is read
    as the last field of each line.

    Args:
        index_path (pathlib.Path): the '.jf' index.
        seqs (list): the sequences to query.
        k (int): k-mer length.

    Returns:
        counts (numpy.ndarray): int64 maximum counts, 0 for sequences shorter than k.
    """
    seqs = list(seqs)
    n_kmers = np.array([max(0, len(s) - k + 1) for s in seqs], dtype=np.int64)
    counts = np.zeros(len(seqs), dtype=np.int64)
    if n_kmers.sum() == 0:
        return counts

    # -i is required: without it jellyfish query ignores stdin, expects the k-mers
    # as trailing arguments, and exits having printed nothing
    cmd = [have_jellyfish(), 'query', '-i', str(index_path)]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, text=True, bufsize=1 << 20)

    def _feed():
        try:
            for seq in seqs:
                upper = seq.upper()
                for i in range(len(upper) - k + 1):
                    proc.stdin.write(upper[i:i + k] + '\n')
            proc.stdin.close()
        except (BrokenPipeError, ValueError):
            pass

    feeder = threading.Thread(target=_feed, daemon=True)
    feeder.start()

    seq_idx = 0
    remaining = n_kmers[0]
    running_max = 0

    # sequences shorter than k produce no output lines, so skip past them
    while seq_idx < len(seqs) and remaining == 0:
        seq_idx += 1
        remaining = n_kmers[seq_idx] if seq_idx < len(seqs) else 0

    for line in proc.stdout:
        parts = line.split()
        if not parts:
            continue
        try:
            value = int(parts[-1])
        except ValueError:
            continue
        if value > running_max:
            running_max = value
        remaining -= 1
        if remaining == 0:
            counts[seq_idx] = running_max
            running_max = 0
            seq_idx += 1
            while seq_idx < len(seqs) and n_kmers[seq_idx] == 0:
                seq_idx += 1
            if seq_idx >= len(seqs):
                break
            remaining = n_kmers[seq_idx]

    proc.stdout.close()
    proc.wait()
    feeder.join(timeout=5)

    # success
    return counts
