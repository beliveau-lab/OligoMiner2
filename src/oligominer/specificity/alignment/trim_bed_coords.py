"""
# Trim BED Coordinates

Clamps BED coordinates to chromosome boundaries using bedtools slop, ensuring
that downstream tools like bedtools getfasta do not attempt to access positions
outside the reference genome.
"""

from oligominer.utils import check_input_exists, check_output_exists, run_cmd, get_abs_path, require_one_of, ensure_executable
from oligominer.bioinformatics.file_io import get_or_create_fai

def trim_bed_coords(bed_path=None, bed_data=None, fasta_path=None, fai_path=None, output_file=None, verbose=False):
    """
    Trim BED coordinates to fit within chromosome boundaries.

    Uses bedtools slop with -b 0 to clamp coordinates that extend beyond
    the ends of chromosomes. Requires either a FASTA path (from which a
    .fai index is derived) or a direct .fai path.

    Args:
        bed_path (str, optional): path to a BED file.
        bed_data (str, optional): BED data as a string.
        fasta_path (str, optional): path to a reference FASTA file.
        fai_path (str, optional): path to a .fai genome index file.
        output_file (str, optional): path to write output. If None, result
            is returned as a string.
        verbose (bool): if True, print stdout and stderr of each command
            to the terminal.

    Returns:
        result (str or None): trimmed BED data, or None if output_file
            is provided.

    Raises:
        ValueError: if neither or both of fasta_path and fai_path are provided,
            or if neither or both of bed_path and bed_data are provided.
    """
    ensure_executable('bedtools')

    # get path to fai file
    require_one_of(fasta_path, fai_path, 'fasta_path', 'fai_path')
    if fasta_path is not None:
        # get or create a .fai file path from the provided fasta path
        fai_path = get_or_create_fai(fasta_path)
    else:
        # ensure the target .fai file exists
        check_input_exists(fai_path)
        fai_path = get_abs_path(fai_path)

    # run bedtools slop to trim BED file coords to chromosome boundaries
    require_one_of(bed_path, bed_data, 'bed_path', 'bed_data')
    if bed_path is not None:
        # pass a BED file path into bedtools slop
        result = run_cmd(
            ['bedtools', 'slop', '-i', bed_path, '-g', fai_path, '-b', '0'],
            output_file=output_file, verbose=verbose
        )
    else:
        # pipe BED data into bedtools slop
        result = run_cmd(
            ['bedtools', 'slop', '-g', fai_path, '-b', '0'],
            input_data=bed_data, output_file=output_file, verbose=verbose
        )

    # check if the output file was created
    if output_file is not None:
        check_output_exists(output_file)

    # success
    return result
