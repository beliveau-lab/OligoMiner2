"""
# Jellyfish Index Building

Builds a Jellyfish kmer count index from a FASTA/FASTQ file using
jellyfish count. After building, validates the index. See jellyfish_query.py
for querying a built index.
"""

import os

from oligominer.utils.cores import resolve_cores
from oligominer.utils.shell_pipeline import run_cmd
from oligominer.utils import (
    get_abs_path, check_dir_exists, check_input_exists, check_output_exists,
    ensure_executable
)


# smallest hash to ask for, so a tiny input does not get a degenerate table
MIN_HASH_SIZE = 1_000_000


def estimate_hash_size(input_file, floor=MIN_HASH_SIZE):
    """
    Estimate Jellyfish's initial hash size from the input.

    Jellyfish sizes its hash by the number of distinct k-mers it expects, which
    for a genome is close to its length. Asking for a fixed large size makes a
    small genome allocate for a large one; asking for too little makes
    Jellyfish write intermediate files and merge them.

    The estimate is the file's size in bytes, which overshoots the sequence
    length by the FASTA headers and newlines -- an overshoot being the cheaper
    error of the two.

    Args:
        input_file (str): path to the input FASTA/FASTQ.
        floor (int): smallest size to return.

    Returns:
        size (int): the hash size to pass to '-s'.
    """
    # success
    return max(int(floor), os.path.getsize(input_file))


def jellyfish_build(input_file, output_file, k=18, size=None, cores=None,
                    canonical=False, out_counter_len=4, text=False, disk=False,
                    sam=None, min_qual_char=None, min_quality=None,
                    reprobes=None, shell=None, L=None, U=None, verbose=False):
    """
    Build a Jellyfish kmer count index (wrapper for 'jellyfish count').

    Args:
        input_file (str): path to the input FASTA/FASTQ file.
        output_file (str): output path for the kmer count file (.jf).
        k (int): length of k-mers.
        size (str or int, optional): initial hash size. None estimates it from
            the input's size, so a small genome does not allocate for a large
            one.
        cores (int, optional): number of threads. None resolves the batch
            scheduler's granted core count (see oligominer.utils.cores).
        canonical (bool): count canonical representation of k-mers.
        out_counter_len (int): counter field length in bytes.
        text (bool): output in text format.
        disk (bool): disk operation, avoid in-memory doubling.
        sam (str, optional): path to SAM/BAM/CRAM file.
        min_qual_char (str, optional): any base below this quality is changed
            to N.
        min_quality (int, optional): minimum quality; below this becomes N.
        reprobes (int, optional): maximum number of reprobes.
        shell (str, optional): shell to run generator commands.
        L (int, optional): don't output k-mers with count < L.
        U (int, optional): don't output k-mers with count > U.
        verbose (bool): if True, print stdout and stderr to the terminal.

    Returns:
        index_path (str): absolute path to the created index file.
    """
    ensure_executable('jellyfish')
    check_input_exists(input_file)
    index_path = get_abs_path(output_file)
    cores = resolve_cores(cores)
    if size is None:
        size = estimate_hash_size(input_file)

    # create output directory as needed
    check_dir_exists(output_file, parent_dir=True, create=True)

    # build the jellyfish count command
    cmd = [
        'jellyfish', 'count',
        '-m', str(k),
        '-s', str(size),
        '-t', str(cores),
        '-o', str(output_file),
        '--out-counter-len', str(out_counter_len),
        input_file
    ]
    if canonical:
        cmd.append('--canonical')
    if text:
        cmd.append('--text')
    if disk:
        cmd.append('--disk')
    if sam:
        cmd.append(f'--sam={sam}')
    if min_qual_char:
        cmd.extend(['-Q', str(min_qual_char)])
    if min_quality is not None:
        cmd.append(f'--min-quality={min_quality}')
    if reprobes is not None:
        cmd.extend(['-p', str(reprobes)])
    if shell:
        cmd.extend(['-S', shell])
    if L is not None:
        cmd.extend(['-L', str(L)])
    if U is not None:
        cmd.extend(['-U', str(U)])

    run_cmd(cmd, verbose=verbose)
    check_output_exists(index_path)

    # success
    return index_path
