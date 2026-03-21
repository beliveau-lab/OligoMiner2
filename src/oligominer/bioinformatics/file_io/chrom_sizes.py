"""
Chromosome sizes utilities.

Provides functions for looking up and printing chromosome sizes from
FASTA index (.fai) files, creating them as needed via pyfaidx.
"""

import os

from pyfaidx import Fasta

from oligominer.utils import get_abs_path, get_dir_name, check_input_exists, check_output_exists
from .exceptions import FastaPermissionError


def print_chrom_sizes(fasta_path):
    """
    Print a genome/chrom.sizes file to stdout for the input fasta file.

    Args:
        fasta_path (str): path to the input fasta file.
    """
    chrom_sizes = get_chrom_sizes(fasta_path)
    for chrom in chrom_sizes:
        print(f'{chrom}\t{chrom_sizes[chrom]}')


def get_chrom_sizes(fasta_path):
    """
    Lookup chromosome sizes for the input fasta file.

    The .fai genome file corresponding to the input fasta file is created
    as needed, and then its contents are loaded and parsed, storing the
    size of each seq record in the fasta file in a python dict.

    Args:
        fasta_path (str): path to the input fasta file.

    Returns:
        chrom_sizes (dict): the size of each sequence record in the input fasta.
    """
    fai_path = get_or_create_fai(fasta_path)

    # parse .fai (tab-separated: name, length, offset, linebases, linewidth)
    chrom_sizes = {}
    with open(fai_path, 'r') as infile:
        for line in infile:
            fields = line.strip().split('\t')
            chrom_sizes[fields[0]] = int(fields[1])

    # sort by size descending
    chrom_sizes = dict(sorted(chrom_sizes.items(), key=lambda x: -x[1]))

    # success
    return chrom_sizes

def get_or_create_fai(fasta_path):
    """
    Returns the path of the .fai genome file corresponding to the input
    fasta, creating it as needed.

    Args:
        fasta_path (str): path to the input fasta file.

    Returns:
        fai_path (str): path to the .fai genome file.
    """
    # ensure input file exists
    fasta_path = get_abs_path(fasta_path)
    check_input_exists(fasta_path)

    # get or create .fai file
    fai_path = fasta_path + '.fai'
    if not os.path.exists(fai_path):
        fai_path = create_fai(fasta_path)
        
    # success    
    return fai_path


def create_fai(fasta_path):
    """
    Creates a .fai genome file using pyfaidx.

    See: https://pypi.org/project/pyfaidx/

    Args:
        fasta_path (str): path to the input fasta file.

    Returns:
        fai_path (str): path to the created .fai genome file.
    """
    # ensure input file exists
    fasta_path = get_abs_path(fasta_path)
    check_input_exists(fasta_path)

    # ensure destination dir for fasta index is writeable
    fasta_dir = get_dir_name(fasta_path)
    if not os.access(fasta_dir, os.W_OK):
        raise FastaPermissionError(fasta_path)

    # create .fai index (pyfaidx writes it as a side effect)
    Fasta(fasta_path)

    # ensure target output file is present
    fai_path = fasta_path + '.fai'
    check_output_exists(fai_path)

    # success
    return fai_path
