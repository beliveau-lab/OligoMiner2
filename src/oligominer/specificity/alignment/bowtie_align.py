"""
# Bowtie2 Alignment

Aligns reads to a reference genome using Bowtie2. Uses run_cmd for single
commands (index validation) and ShellPipeline for multi-step pipes (alignment
with BAM conversion). For building new indexes, see bowtie_build.py. For
preset parameter configurations, see bowtie_presets.py.
"""

from oligominer.utils.shell_pipeline import run_cmd, ShellPipeline
from oligominer.utils import (
    get_abs_path, check_dir_exists, check_output_exists, require_one_of,
    ensure_executable
)
from .exceptions import MissingBowtieIndexError
from .bowtie_presets import BT2_INDEX_EXTENSIONS


def check_index_exists(index_path):
    """
    Check that all expected Bowtie2 index files exist on disk.

    Args:
        index_path (str): base path to the Bowtie2 index.

    Returns:
        result (bool): True if all index files exist.

    Raises:
        MissingOutputFile: if any index file is missing.
    """
    for ext in BT2_INDEX_EXTENSIONS:
        check_output_exists(index_path + ext)

    # success
    return True


def validate_index(index_path, verbose=False):
    """
    Validate a Bowtie2 index using bowtie2-inspect.

    Args:
        index_path (str): base path to the Bowtie2 index.
        verbose (bool): if True, print stdout and stderr to the terminal.

    Returns:
        result (bool): True if the index is valid.

    Raises:
        MissingBowtieIndexError: if the index files do not exist.
    """
    ensure_executable('bowtie2-inspect')
    if not check_index_exists(index_path):
        raise MissingBowtieIndexError(index_path)
    run_cmd(['bowtie2-inspect', '-s', index_path], verbose=verbose)

    # success
    return True


def bowtie_align(index_path, input_file=None, input_data=None,
                 sam_output_file=None, bam_output_file=None,
                 preset=None, D=None, R=None, N=None, L=None, i=None,
                 local=False, no_1mm_upfront=False, nofw=False, norc=False,
                 dpad=None, gbar=None, ignore_quals=False,
                 n_ceil=None, ma=None, mp=None, np=None, rdg=None, rfg=None,
                 score_min=None,
                 k=None, a=False, threads=1, reorder=False, mm=False,
                 fasta_input=False, no_unal=False, no_hd=True, xeq=True,
                 no_sq=False, time=False, verbose=False, bt2_verbose=False):
    """
    Align reads to a reference genome using Bowtie2.

    Args:
        index_path (str): base path to the Bowtie2 index files.
        input_file (str, optional): path to the input FASTA/FASTQ file.
        input_data (str, optional): FASTQ (or FASTA with fasta_input=True) data
            as a string to be piped into Bowtie2.
        sam_output_file (str, optional): path to the output SAM file.
        bam_output_file (str, optional): path to the output BAM file.
        preset (dict, optional): use preset parameters from bowtie_presets
            (e.g. bowtie_presets.VERY_FAST).
        D (int, optional): max consecutive seed extension attempts that can fail
            before Bowtie2 moves on.
        R (int, optional): for reads with repetitive seeds, try this many sets
            of seeds.
        N (int, optional): max mismatches in seed alignment (0 or 1). Default: 0.
        L (int, optional): length of seed substrings to align.
        i (str, optional): interval between seed substrings.
            Format 'S,1,0.50' means f(x)=1+0.5*sqrt(x).
        local (bool): use local alignment mode instead of end-to-end.
        no_1mm_upfront (bool): skip the one-mismatch search before the
            multiseed heuristic.
        nofw (bool): do not align forward (Watson) strand.
        norc (bool): do not align reverse-complement (Crick) strand.
        dpad (int, optional): pad dynamic programming problems by this many
            columns on each side.
        gbar (int, optional): disallow gaps within this many positions of
            read extremes.
        ignore_quals (bool): treat all quality values as high.
        n_ceil (str, optional): function for max number of Ns allowed.
            Format 'L,0,0.15' means f(x)=0+0.15*x.
        ma (int, optional): match bonus in local mode (>0).
        mp (str, optional): max and min mismatch penalties (MX,MN).
        np (int, optional): penalty for ambiguous chars (Ns) in read or ref.
        rdg (str, optional): read gap open,extend penalties (INT1,INT2).
        rfg (str, optional): reference gap open,extend penalties (INT1,INT2).
        score_min (str, optional): min acceptable alignment score w.r.t read
            length (L,0,-0.6).
        k (int, optional): report up to k distinct alignments per read.
        a (bool): report all alignments per read (very slow).
        threads (int): number of parallel search threads. Default: 1.
        reorder (bool): keep SAM output in order of input reads.
        mm (bool): use memory-mapped I/O for index.
        fasta_input (bool): input files are in FASTA format (-f).
        no_unal (bool): suppress SAM records for reads that failed to align.
        no_hd (bool): suppress SAM header lines (starting with @).
        xeq (bool): use '='/'X' instead of 'M' to specify matches/mismatches
            in SAM record.
        no_sq (bool): suppress @SQ SAM header lines.
        time (bool): print wall-clock time for loading index and aligning.
        verbose (bool): if True, print stdout and stderr of each command
            to the terminal.
        bt2_verbose (bool): if True, run bowtie2 in --verbose mode.

    Returns:
        result (str or None): SAM file content if no output file is
            specified, otherwise None.

    Raises:
        ValueError: if neither input_file nor input_data is provided.
        RuntimeError: if the alignment fails.
    """
    require_one_of(input_file, input_data, 'input_file', 'input_data')
    ensure_executable('bowtie2')
    if bam_output_file:
        ensure_executable('samtools')

    index_path = get_abs_path(index_path)
    cmd = ['bowtie2', '-x', index_path]

    if input_file:
        cmd.extend(['-U', input_file])
    elif input_data:
        cmd.extend(['-U', '-'])

    if sam_output_file:
        cmd.extend(['-S', sam_output_file])

    if fasta_input:
        cmd.append('-f')

    # apply preset parameters if specified
    if preset is not None:
        D = preset.get('D', D)
        R = preset.get('R', R)
        N = preset.get('N', N)
        L = preset.get('L', L)
        i = preset.get('i', i)
        local = preset.get('local', local)

    # add seed parameters if specified
    if D is not None:
        cmd.extend(['-D', str(D)])
    if R is not None:
        cmd.extend(['-R', str(R)])
    if N is not None:
        cmd.extend(['-N', str(N)])
    if L is not None:
        cmd.extend(['-L', str(L)])
    if i is not None:
        cmd.extend(['-i', str(i)])

    # add alignment mode flags
    if local:
        cmd.append('--local')

    # add other options
    if no_1mm_upfront:
        cmd.append('--no-1mm-upfront')
    if nofw:
        cmd.append('--nofw')
    if norc:
        cmd.append('--norc')
    if dpad is not None:
        cmd.extend(['--dpad', str(dpad)])
    if gbar is not None:
        cmd.extend(['--gbar', str(gbar)])
    if ignore_quals:
        cmd.append('--ignore-quals')
    if n_ceil is not None:
        cmd.extend(['--n-ceil', str(n_ceil)])

    # add scoring parameters
    if ma is not None:
        cmd.extend(['--ma', str(ma)])
    if mp is not None:
        cmd.extend(['--mp', str(mp)])
    if np is not None:
        cmd.extend(['--np', str(np)])
    if rdg is not None:
        cmd.extend(['--rdg', str(rdg)])
    if rfg is not None:
        cmd.extend(['--rfg', str(rfg)])
    if score_min is not None:
        cmd.extend(['--score-min', str(score_min)])

    # add reporting options
    if k is not None:
        cmd.extend(['-k', str(k)])
    if a:
        cmd.append('-a')

    # add threading options
    if threads > 1:
        cmd.extend(['-p', str(threads)])
    if reorder:
        cmd.append('--reorder')
    if mm:
        cmd.append('--mm')

    # add output format options
    if no_unal:
        cmd.append('--no-unal')
    if no_hd:
        cmd.append('--no-hd')
    if xeq:
        cmd.append('--xeq')
    if no_sq:
        cmd.append('--no-sq')
    if time:
        cmd.append('-t')

    if bt2_verbose:
        cmd.append('--verbose')

    # create the pipeline
    pipeline = ShellPipeline(binary=True)
    pipeline.add(cmd)

    # if BAM output is specified, add samtools command to the pipeline
    if bam_output_file:
        pipeline.add(['samtools', 'view', '-bS'])

    # ensure output directory exists
    if bam_output_file is not None:
        check_dir_exists(bam_output_file, parent_dir=True, create=True)
    elif sam_output_file is not None:
        check_dir_exists(sam_output_file, parent_dir=True, create=True)

    if input_data:
        result = pipeline.run(input_data=input_data.encode('utf-8'), output_file=bam_output_file, verbose=verbose)
    else:
        result = pipeline.run(output_file=bam_output_file, verbose=verbose)

    # TODO make sure that each case makes sense:
    #   * pipe vs file input to alignment (in fastq or fasta format)
    #   * output to file vs return the result (in sam or bam format)

    # if no output file is specified, return the result
    if not sam_output_file and not bam_output_file:
        return result.decode()

    return None
