"""
# Process Alignments

End-to-end pipeline for processing alignment results. Takes SAM data or a BAM
file, converts to BED, builds a DataFrame of alignment metadata, and optionally
looks up derived sequences from a reference genome.
"""

import pandas as pd

from oligominer.specificity.alignment import bam_to_bed, trim_bed_coords, get_fasta
from oligominer.bioinformatics.file_io.bed_io import bed_to_df
from oligominer.bioinformatics.file_io.sam_bam_io import load_bam_file
from oligominer.utils import require_one_of

def process_alignments(sam_data=None, bam_path=None, ref_fasta=None, to_upper=True):
    """
    Process alignment data into a structured DataFrame.

    Converts SAM/BAM alignment data to BED format, parses it into a DataFrame
    with alignment metadata, and optionally looks up derived sequences from
    a reference genome at the aligned coordinates.

    Args:
        sam_data (str, optional): SAM data as a string.
        bam_path (str, optional): path to an input BAM file.
        ref_fasta (str, optional): path to a reference FASTA file. If provided,
            derived sequences are looked up and added to the output.
        to_upper (bool): if True, convert derived sequences to uppercase.

    Returns:
        align_df (pandas.DataFrame): alignment results with columns
            align_seqid, align_start, align_stop, seqid, align_score,
            align_strand, align_cigar, and optionally derived_seq.

    Raises:
        ValueError: if neither or both of sam_data and bam_path are provided.
    """
    # load bam file as needed
    require_one_of(sam_data, bam_path, 'sam_data', 'bam_path')
    if bam_path is not None:
        sam_data = load_bam_file(bam_path)

    # run om2 custom bamtobed
    bed_data = bam_to_bed(bam_data=sam_data)

    # convert BED to a pandas dataframe
    align_df = bed_to_df(bed_data=bed_data)
    align_df.columns = [
        'align_seqid',
        'align_start',
        'align_stop',
        'seqid',
        'align_score',
        'align_strand',
        'align_cigar',
    ]

    # optionally lookup derived sequences from reference genome
    if ref_fasta is not None:

        # trim bed coordinates to fit reference genome features
        trimmed_bed_data = trim_bed_coords(bed_data=bed_data, fasta_path=ref_fasta)

        # lookup sequences at probe alignment sites in reference genome
        fasta_result = get_fasta(bed_data=trimmed_bed_data, fasta_path=ref_fasta)
        derived_seqs = fasta_result.strip().split('\n')
        align_df['derived_seq'] = derived_seqs

        if to_upper:
            align_df['derived_seq'] = align_df['derived_seq'].str.upper()

    # success
    return align_df
