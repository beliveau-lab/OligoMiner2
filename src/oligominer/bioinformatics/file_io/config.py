"""
Configuration constants and shared helpers for the file_io sub-package.
"""

import os

# FASTA file extensions to check for when loading FASTA files
FASTA_EXTENSIONS = {'.fa', '.fasta', '.fna', '.fas'}

# GTF/GFF file extensions to check for when loading annotation files
GTF_EXTENSIONS = {'.gtf', '.gff', '.gff3'}

# standard GTF column names (9-column format)
GTF_COLUMNS = [
    'seqid', 'source', 'type', 'start', 'end',
    'score', 'strand', 'phase', 'attributes',
]

# default sequence record classification rules
# maps category name -> regex pattern; evaluated in order, first match wins
# unmatched IDs receive the default_category label (default: 'canonical')
DEFAULT_CLASSIFICATION_RULES = {
    'alt': r'_alt',
    'hap': r'_hap',
    'fix': r'_fix',
    'unlocalized': r'_random',
    'unplaced': r'Un_',
}


# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

def merge_files_by_extension(input_dir, output_path, extensions):
    """
    Concatenate all files matching the given extensions into a single file.

    Files are read in sorted filename order. Each file is written verbatim,
    with a newline appended if the file does not end with one.

    Args:
        input_dir (str): directory containing files to merge.
        output_path (str): path for the merged output file.
        extensions (set): file extensions to include (e.g. {'.fa', '.fasta'}).

    Returns:
        merged_paths (list): paths of the files that were merged.
    """
    merged_paths = []
    for fname in sorted(os.listdir(input_dir)):
        fpath = os.path.join(input_dir, fname)
        if not os.path.isfile(fpath):
            continue
        if os.path.splitext(fname)[1].lower() not in extensions:
            continue
        merged_paths.append(fpath)

    with open(output_path, 'w') as outfile:
        for fpath in merged_paths:
            with open(fpath, 'r') as infile:
                for line in infile:
                    outfile.write(line)
                if not line.endswith('\n'):
                    outfile.write('\n')

    # success
    return merged_paths

