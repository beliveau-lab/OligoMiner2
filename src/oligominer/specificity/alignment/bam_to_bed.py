"""
# BAM to BED Conversion

Converts SAM/BAM alignment data to BED format using an AWK script that
reconstructs full probe-length genomic coordinates from alignment positions
and extended CIGAR strings.

Standard aligners like Bowtie2 report the aligned portion of a read, which
may be shorter than the original probe when soft-clipping occurs (e.g. the
probe partially overhangs a chromosome boundary or contains a mismatched
flank). The reported POS in SAM is where the aligned portion starts, not
where the full probe starts.

The AWK script below recovers the full probe footprint on the genome:

  1. Start with the SAM POS field (1-based leftmost position of the
     aligned portion).

  2. If the CIGAR begins with a soft-clip (e.g. "5S..."), the first N
     bases of the probe were not aligned. Subtract that count from the
     start position to recover where the full probe would begin.

  3. Sum all operation lengths in the CIGAR to get the total probe length,
     then add that to the (adjusted) start to get the end coordinate.

  4. Convert to 0-based half-open BED coordinates (start-1, end-1).

  5. Extract strand from the SAM FLAG field (bit 0x10 = reverse complement)
     and the alignment score from the AS:i: optional tag.

Output columns: chrom, start, end, read_name, align_score, strand, cigar.
"""

from oligominer.utils.shell_pipeline import run_cmd
from oligominer.utils import require_one_of, check_output_exists, ensure_executable

# AWK script that parses SAM records into BED with full probe coordinates
AWK_SCRIPT = """
{
    OFS = "\\t";
    start = $4;
    end = start;
    cigar = $6;

    if (match(cigar, /^[0-9]+S/)) {
        split(cigar, parts, /[A-Z]/);
        start_adjustment = parts[1];
        start -= start_adjustment;
    }

    n = split(cigar, cig_ops, /[0-9]+/);
    m = split(cigar, cig_lengths, /[MIDNSHPX=]/);

    total_length = 0;
    for (i = 1; i <= n; i++) {
        total_length += cig_lengths[i];
    }
    end = start + total_length;

    strand = (int($2 / 16) % 2) ? "-" : "+";

    align_score = "NA";
    for (i = 12; i <= NF; i++) {
        if (substr($i, 1, 2) == "AS") {
            align_score = substr($i, 6);
            break;
        }
    }

    bed_start = start - 1;
    if (bed_start < 0) bed_start = 0;
    bed_end = end - 1;
    if (bed_end < 0) bed_end = 0;
    print $3, bed_start, bed_end, $1, align_score, strand, cigar;
}
"""

def bam_to_bed(input_file=None, bam_data=None, output_file=None, verbose=False):
    """
    Convert BAM alignments to BED format.

    Args:
        input_file (str, optional): path to the input BAM file.
        bam_data (str, optional): BAM data as a string.
        output_file (str, optional): path to the output BED file.
        verbose (bool): if True, print stdout and stderr of each command
            to the terminal.

    Returns:
        result (str): the BED data as a string.

    Raises:
        ValueError: if neither input_file nor bam_data is provided.
    """
    require_one_of(input_file, bam_data, 'input_file', 'bam_data')

    if input_file is not None:
        ensure_executable('samtools')
        bam_data = run_cmd(['samtools', 'view', input_file])

    result = run_cmd(['awk', AWK_SCRIPT], input_data=bam_data, verbose=verbose)
    if output_file is not None:
        with open(output_file, "w") as f:
            f.write(result)
        check_output_exists(output_file)

    # success
    return result


def sam_to_bed(sam_data):
    """
    Convert SAM alignments to BED format.

    Args:
        sam_data (str): SAM data as a string.

    Returns:
        result (str): the converted BED data.
    """
    result = run_cmd(['awk', AWK_SCRIPT], input_data=sam_data)

    # success
    return result
