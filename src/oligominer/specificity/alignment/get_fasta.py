"""
# Get FASTA

Extracts sequences from a reference genome at positions specified by BED
coordinates, using bedtools getfasta under the hood.
"""

from oligominer.utils import check_input_exists, check_output_exists, require_one_of, ensure_executable
from oligominer.utils.shell_pipeline import run_cmd, ShellPipeline

def get_fasta(bed_path=None, bed_data=None, fasta_path=None, output_file=None, verbose=False, strip_col=True):
    """
    Extract sequences from a reference FASTA at BED-specified coordinates.

    Uses bedtools getfasta to look up sequences. Accepts BED input either as
    a file path or as piped string data. When strip_col is True, pipes through
    awk to return bare sequences without the coordinate column.

    Args:
        bed_path (str, optional): path to a BED file.
        bed_data (str, optional): BED data as a string.
        fasta_path (str): path to the reference FASTA file.
        output_file (str, optional): path to write output. If None, result
            is returned as a string.
        verbose (bool): if True, print stdout and stderr of each command
            to the terminal.
        strip_col (bool): if True, strip the coordinate column from the
            bedtools output, returning only sequences.

    Returns:
        result (str or None): extracted sequences, or None if output_file
            is provided.

    Raises:
        ValueError: if neither or both of bed_path and bed_data are provided.
    """
    ensure_executable('bedtools')

    # ensure the target fasta file exists
    check_input_exists(fasta_path)

    require_one_of(bed_path, bed_data, 'bed_path', 'bed_data')

    # build the bedtools getfasta command
    bedtools_cmd = ['bedtools', 'getfasta', '-tab', '-s', '-fi', fasta_path, '-bed']
    if bed_path is not None:
        bedtools_cmd.append(bed_path)
        input_data = None
    else:
        bedtools_cmd.append('-')
        input_data = bed_data

    # single command when no stripping needed, pipeline when piping through awk
    if strip_col:
        pipeline = ShellPipeline()
        pipeline.add(bedtools_cmd)
        pipeline.add(['awk', '{ print $2}'])
        result = pipeline.run(input_data=input_data, output_file=output_file, verbose=verbose)
    else:
        result = run_cmd(bedtools_cmd, input_data=input_data, output_file=output_file, verbose=verbose)

    # check if the output file was created
    if output_file is not None:
        check_output_exists(output_file)

    # success
    return result
