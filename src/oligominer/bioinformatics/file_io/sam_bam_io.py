"""
# SAM/BAM I/O

Utilities for loading alignment data from SAM and BAM files. BAM loading
uses samtools via run_cmd to convert binary data to text format.
"""

from oligominer.utils.shell_pipeline import run_cmd
from oligominer.utils import check_input_exists, ensure_executable

def load_sam_file(input_file):
    """
    Load SAM data from a file.

    Args:
        input_file (str): path to the input SAM file.

    Returns:
        sam_data (str): the SAM data as a string.
    """
    check_input_exists(input_file)

    with open(input_file, "r") as f:
        sam_data = f.read()

    # success
    return sam_data

def load_bam_file(input_file):
    """
    Load SAM data from a BAM file.

    Args:
        input_file (str): path to the input BAM file.

    Returns:
        sam_data (str): the SAM data as a string.
    """
    ensure_executable('samtools')
    check_input_exists(input_file)

    sam_data = run_cmd(['samtools', 'view', input_file])

    # success
    return sam_data
