"""
# Jellyfish Query

Functions for querying a Jellyfish kmer count index, including validation,
raw queries, and higher-level kmer count operations on sequences. For building
new indexes, see jellyfish_build.py.
"""

import re
import tempfile

from oligominer.utils.shell_pipeline import run_cmd
from oligominer.utils import get_abs_path, check_input_exists, ensure_executable
from oligominer.bioinformatics.file_io.fasta_io import seqs_to_fasta

from .exceptions import JellyfishIndexError, MissingJellyfishIndexError


def validate_index(index_path, k=None, verbose=False):
    """
    Validate a Jellyfish index file and extract its metadata.

    Checks that the file exists, runs jellyfish info to extract the k value
    and canonical flag, and optionally verifies that k matches an expected
    value.

    Args:
        index_path (str): path to the Jellyfish index file.
        k (int, optional): expected k-mer length. If provided and does not
            match the index, raises JellyfishIndexError.
        verbose (bool): if True, print stdout and stderr to the terminal.

    Returns:
        info (dict): index metadata with keys 'index_path' (str), 'k' (int),
            and 'is_canonical' (bool).

    Raises:
        MissingJellyfishIndexError: if the index file does not exist.
        JellyfishIndexError: if k is provided and does not match the index.
    """
    ensure_executable('jellyfish')
    index_path = get_abs_path(index_path)
    check_input_exists(index_path)

    # get metadata from jellyfish index file
    jellyfish_info = run_cmd(['jellyfish', 'info', index_path], verbose=verbose)

    # extract k value from the jellyfish count command string
    command_str = re.findall('command: (.+)', jellyfish_info).pop()
    idx_k = int(re.findall(r'\s+-m\s+(\d+)', command_str).pop())

    # if a k value was provided, ensure it matches the index
    if (k is not None) and (k != idx_k):
        raise JellyfishIndexError(index_path, idx_k, k)

    is_canonical = 'yes' in re.findall('canonical: (.+)', jellyfish_info)

    info = {
        'index_path': index_path,
        'k': idx_k,
        'is_canonical': is_canonical,
    }

    # success
    return info


def jellyfish_query(index_path, mers=None, fasta_path=None, output=None,
                    load=False, no_load=False, verbose=False):
    """
    Query a Jellyfish database (wrapper for 'jellyfish query').

    Args:
        index_path (str): path to the Jellyfish index file.
        mers (list or None): list of kmer strings to query.
        fasta_path (str or None): path to a FASTA file; queries all kmers
            from all sequences in the file.
        output (str or None): path to output file. If None, returns stdout.
        load (bool): force pre-loading of database file into memory.
        no_load (bool): disable pre-loading of database file into memory.
        verbose (bool): if True, print stdout and stderr to the terminal.

    Returns:
        result (str): query results as a tab-delimited string of kmer/count
            pairs, one per line.
    """
    ensure_executable('jellyfish')
    check_input_exists(index_path)

    cmd = ['jellyfish', 'query']
    if fasta_path:
        cmd.extend(['-s', fasta_path])
    if output:
        cmd.extend(['-o', output])
    if load:
        cmd.append('-l')
    if no_load:
        cmd.append('-L')
    cmd.append(index_path)
    if mers:
        cmd.extend(mers)

    result = run_cmd(cmd, verbose=verbose)

    # success
    return result


def _get_kmers(seq, k):
    """
    Decompose a sequence into overlapping k-mers.

    Args:
        seq (str): the input DNA sequence.
        k (int): k-mer length.

    Returns:
        kmers (list): list of k-mer strings.
    """
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]

    # success
    return kmers


def calc_max_kmer(index_path, seq, k, load=False, no_load=False, verbose=False):
    """
    Find the maximum k-mer count in a given sequence using a Jellyfish database.

    Decomposes the sequence into overlapping k-mers, queries the index for
    each, and returns the highest count.

    Args:
        index_path (str): path to the Jellyfish index file.
        seq (str): the sequence to decompose into k-mers.
        k (int): k-mer length.
        load (bool): force pre-loading of database file into memory.
        no_load (bool): disable pre-loading of database file into memory.
        verbose (bool): if True, print stdout and stderr to the terminal.

    Returns:
        max_count (int): the maximum k-mer count.
    """
    kmers = _get_kmers(seq, k)
    counts = jellyfish_query(index_path, mers=kmers, load=load, no_load=no_load, verbose=verbose)

    # extract counts from the query result
    count_lines = counts.strip().split('\n')
    count_values = [int(line.split()[1]) for line in count_lines]

    max_count = max(count_values)

    # success
    return max_count


def calc_max_kmer_multi(index_path, seqs, k, verbose=False):
    """
    Calculate the maximum k-mer count for each sequence in a list.

    Writes all sequences to a temporary FASTA file, queries the Jellyfish
    index for all k-mers in one batch, then partitions the results by
    sequence to find the per-sequence maximum.

    Args:
        index_path (str): path to the Jellyfish index file.
        seqs (list): a list of nucleotide sequence strings.
        k (int): k-mer length.
        verbose (bool): if True, print stdout and stderr to the terminal.

    Returns:
        max_kmer_values (list): a list containing the maximum k-mer count
            for each input sequence.
    """
    # write all sequences to a temp fasta file for batch query
    fasta_str = seqs_to_fasta(seqs)
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fa') as f:
        f.write(fasta_str)
        f.flush()

        # query the jellyfish index
        result = jellyfish_query(index_path, fasta_path=f.name, verbose=verbose)

    # extract counts from the query result
    count_values = [int(m.group(1)) for m in re.finditer(r'\s(\d+)$', result, re.MULTILINE)]

    # calculate the number of kmers per sequence
    kmer_counts_per_seq = [len(seq) - k + 1 for seq in seqs]

    # partition counts by sequence and find the max for each
    max_kmer_values = []
    count_index = 0
    for count in kmer_counts_per_seq:
        current_counts = count_values[count_index:count_index + count]
        count_index += count
        max_kmer_values.append(max(current_counts))

    # success
    return max_kmer_values
