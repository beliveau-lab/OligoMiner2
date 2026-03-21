"""
# Bowtie2 Index Building

Builds a Bowtie2 index from a FASTA file using bowtie2-build. After building,
validates that all expected index files were created. See bowtie_align.py for
aligning reads against a built index.
"""

from oligominer.utils.shell_pipeline import run_cmd
from oligominer.utils import get_abs_path, check_dir_exists, check_output_exists, ensure_executable
from .bowtie_presets import BT2_INDEX_EXTENSIONS


def bowtie_build(input_file, output_file, verbose=False, cores=None):
    """
    Build a Bowtie2 index from a FASTA file (wrapper for 'bowtie2-build').

    Args:
        input_file (str): path to the input FASTA file.
        output_file (str): base path for the output index files.
        verbose (bool): if True, print stdout and stderr of each command
            to the terminal.
        cores (int, optional): number of CPU cores to use for building
            the index.

    Returns:
        index_path (str): absolute path to the created index base.
    """
    ensure_executable('bowtie2-build')
    index_path = get_abs_path(output_file)

    # create the output directory as needed
    check_dir_exists(output_file, parent_dir=True, create=True)

    # run bowtie2-build
    cmd = ['bowtie2-build', input_file, output_file]
    if cores:
        cmd.extend(['--threads', str(cores)])
    run_cmd(cmd, verbose=verbose)

    # check that all index files were created
    for ext in BT2_INDEX_EXTENSIONS:
        check_output_exists(output_file + ext)

    # success
    return index_path
