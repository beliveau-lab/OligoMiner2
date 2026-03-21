"""
# FASTQ I/O

Utilities for converting sequences to FASTQ format. Generates synthetic
quality scores suitable for downstream alignment tools that require FASTQ input.
"""


def seqs_to_fastq(seq_list, seq_id_list=None):
    """
    Convert a list of sequences into a FASTQ formatted string.

    Args:
        seq_list (list): a list of sequences.
        seq_id_list (list or None): a list of sequence IDs. If None, default
            IDs will be used.

    Returns:
        fastq_str (str): a FASTQ formatted string.
    """
    if seq_id_list is None:
        seq_id_list = [f"seq_{i}" for i in range(len(seq_list))]

    lines = [f"@{seq_id}\n{seq}\n+\n{'~' * len(seq)}\n" for seq_id, seq in zip(seq_id_list, seq_list)]
    fastq_str = "".join(lines)

    # success
    return fastq_str
