"""
FASTA file I/O utilities.

Loading uses pyfaidx for indexed random access (ideal for large genomes).
Writing, splitting, filtering, and merging use plain Python I/O.
"""

import os
import re

from pyfaidx import Fasta

from oligominer.utils import get_abs_path, get_dir_name, check_dir_exists, check_input_exists, check_output_exists
from .chrom_sizes import get_or_create_fai
from .exceptions import EmptyExportError
from .config import FASTA_EXTENSIONS, merge_files_by_extension

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_fasta(fasta_path):
    """
    Load a FASTA file using pyfaidx for indexed random access.

    Args:
        fasta_path (str): path to the input FASTA file.

    Returns:
        fasta (pyfaidx.Fasta): dict-like object (.keys(), [seq_id] -> sequence).
    """
    get_or_create_fai(fasta_path)
    return Fasta(fasta_path, sequence_always_upper=True)


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def write_fasta(seqs, filepath, line_width=60):
    """
    Write sequences to a multi-FASTA file.

    Accepts a plain dict or a pyfaidx.Fasta object.

    Args:
        seqs (dict or pyfaidx.Fasta): {seqid: seq_str} mapping sequence IDs
            to sequence strings, or a pyfaidx.Fasta object.
        filepath (str): destination file path.
        line_width (int): characters per line (None/0 for no wrapping).
    """
    if not seqs:
        raise EmptyExportError(filepath, 'FASTA')

    # normalize pyfaidx.Fasta to plain dict
    if isinstance(seqs, Fasta):
        seqs = {k: str(seqs[k][:]) for k in seqs.keys()}

    # ensure output directory exists
    check_dir_exists(filepath, parent_dir=True, create=True)

    with open(filepath, 'w') as fh:
        for seqid, seq in seqs.items():
            fh.write(f">{seqid}\n")
            if line_width:
                for i in range(0, len(seq), line_width):
                    fh.write(seq[i:i + line_width] + '\n')
            else:
                fh.write(seq + '\n')


def split_fasta(seqs, target_dir, line_width=60):
    """
    Write each sequence to its own FASTA file (<seqid>.fa).

    Args:
        seqs (dict): {seqid: seq_str} mapping sequence IDs to sequence strings.
        target_dir (str): output directory (created if needed).
        line_width (int): characters per line.

    Returns:
        written_paths (list): file paths written.
    """
    if not seqs:
        raise EmptyExportError(target_dir, 'FASTA')

    check_dir_exists(target_dir, create=True)

    written_paths = []
    for seqid, seq in seqs.items():
        out_path = os.path.join(target_dir, f"{seqid}.fa")
        write_fasta({seqid: seq}, out_path, line_width=line_width)
        written_paths.append(out_path)

    return written_paths


def seqs_to_fasta(seq_list, seq_id_list=None):
    """
    Convert a list of sequences into a FASTA formatted string.

    Args:
        seq_list (list): a list of sequences.
        seq_id_list (list or None): a list of sequence IDs. If None, default
            IDs will be used.

    Returns:
        fasta_str (str): a FASTA formatted string.
    """
    if seq_id_list is None:
        seq_id_list = [f"seq_{i}" for i in range(len(seq_list))]

    lines = [f">{seq_id}\n{seq}\n" for seq_id, seq in zip(seq_id_list, seq_list)]
    fasta_str = "".join(lines)

    # success
    return fasta_str


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def filter_seq_ids(seq_source, incl_str=None, excl_str=None):
    """
    Filter sequence IDs from any dict-like source (pyfaidx.Fasta or dict).

    When both incl_str and excl_str are provided, inclusion is applied
    first, then exclusion.

    Args:
        seq_source (dict-like): any object with .keys() returning sequence IDs.
        incl_str (str): regex — keep IDs that match.
        excl_str (str): regex — drop IDs that match.

    Returns:
        seq_ids (list): filtered sequence IDs.
    """
    seq_ids = list(seq_source.keys())

    if incl_str is not None:
        pattern = re.compile(incl_str)
        seq_ids = [sid for sid in seq_ids if pattern.search(sid)]

    if excl_str is not None:
        pattern = re.compile(excl_str)
        seq_ids = [sid for sid in seq_ids if not pattern.search(sid)]

    return seq_ids


def filter_seqs(seq_source, incl_str=None, excl_str=None):
    """
    Filter sequences from any dict-like source (pyfaidx.Fasta or dict).

    Convenience wrapper around filter_seq_ids that returns the filtered
    sequences as a plain dict, ready for write_fasta or split_fasta.

    Args:
        seq_source (dict-like): any object with .keys() and [] access.
        incl_str (str): regex — keep IDs that match.
        excl_str (str): regex — drop IDs that match.

    Returns:
        seqs (dict): {seqid: seq_str} for the filtered sequence IDs.
    """
    filtered_ids = filter_seq_ids(seq_source, incl_str=incl_str, excl_str=excl_str)

    # extract sequences as plain strings
    seqs = {sid: str(seq_source[sid][:]) for sid in filtered_ids}

    # success
    return seqs


# ---------------------------------------------------------------------------
# Merging
# ---------------------------------------------------------------------------


def merge_fastas(input_dir, output_path):
    """
    Merge all FASTA files in a directory into a single file.

    Only files with recognized FASTA extensions (.fa, .fasta, .fna, .fas)
    are included, sorted alphabetically by filename.

    Args:
        input_dir (str): directory containing input FASTA files.
        output_path (str): path for the merged output file.

    Returns:
        output_path (str): path to the created file.
    """
    input_dir = get_abs_path(input_dir)
    check_dir_exists(input_dir)
    check_dir_exists(get_dir_name(output_path), create=True)

    merged = merge_files_by_extension(input_dir, output_path, FASTA_EXTENSIONS)
    if not merged:
        raise EmptyExportError(output_path, 'FASTA')

    check_output_exists(output_path)

    # success
    return output_path
