"""
# Streaming alignment

Aligns probes and writes the resulting BED straight to disk, connecting bowtie2 to
the BED conversion with an OS pipe so neither the SAM nor the BED is ever a Python
object.

    bowtie2 -x index -U probes.fastq | awk <AWK_SCRIPT> > out.bed

The buffered path in bowtie_align() decodes the whole SAM into a string while the
bytes are still referenced, then re-encodes it to feed the converter, so several
copies of the largest object in the pipeline are live at once. Streaming leaves
peak memory set by the bowtie2 index.

Both paths build their bowtie2 argv with build_bowtie2_cmd() and share the AWK
program in bam_to_bed, so the flags and the coordinate arithmetic cannot diverge.
"""

import shutil
import subprocess
import tempfile
from pathlib import Path

from oligominer.utils.cores import resolve_cores
from oligominer.utils import ensure_executable
from oligominer.utils.exceptions import ExternalCommandFailed
from .bam_to_bed import AWK_SCRIPT
from .bowtie_align import build_bowtie2_cmd


def write_fastq(seqs, seq_ids, out_path):
    """
    Write a FASTQ file one record at a time.

    Args:
        seqs (iterable): the sequences.
        seq_ids (iterable): the identifiers, parallel to seqs.
        out_path (str): destination path.

    Returns:
        n (int): the number of records written.
    """
    n = 0
    with open(out_path, 'w') as handle:
        for seq_id, seq in zip(seq_ids, seqs):
            seq = str(seq)
            handle.write(f'@{seq_id}\n{seq}\n+\n{"I" * len(seq)}\n')
            n += 1

    # success
    return n


def align_to_bed(probe_df, bt2_index, out_bed, seq_col='probe_seq',
                 id_col='seqid', threads=None, no_unal=True, keep_fastq=False,
                 **bt2_kwargs):
    """
    Align probes and write BED to disk without holding the SAM in memory.

    Args:
        probe_df (pandas.DataFrame): probes carrying seq_col and id_col.
        bt2_index (str): bowtie2 index prefix.
        out_bed (str): destination BED path.
        seq_col (str): column holding the probe sequence.
        id_col (str): column holding the probe identifier.
        threads (int, optional): bowtie2 threads. None resolves the scheduler grant.
        no_unal (bool): suppress SAM records for reads that did not align. The BED
            conversion discards them regardless, so emitting them is wasted work.
        keep_fastq (bool): keep the intermediate FASTQ instead of deleting it.
        **bt2_kwargs: forwarded to build_bowtie2_cmd (preset, k, and the rest).

    Returns:
        info (dict): n_reads, n_rows, out_bed, fastq and the argv that was run.

    Raises:
        ExternalCommandFailed: if bowtie2 or the BED conversion exits non-zero.
    """
    ensure_executable('bowtie2')
    ensure_executable('awk')

    threads = resolve_cores(threads)
    out_bed = Path(out_bed)
    out_bed.parent.mkdir(parents=True, exist_ok=True)

    tmp_dir = Path(tempfile.mkdtemp(prefix='om2_align_'))
    fastq = tmp_dir / 'probes.fastq'

    try:
        n_reads = write_fastq(probe_df[seq_col], probe_df[id_col], fastq)
        cmd = build_bowtie2_cmd(bt2_index, input_file=str(fastq), threads=threads,
                                no_unal=no_unal, **bt2_kwargs)
        n_rows = _run_stream(cmd, out_bed)
    finally:
        if keep_fastq:
            kept = out_bed.with_suffix('.fastq')
            shutil.move(str(fastq), kept)
            fastq = kept
        shutil.rmtree(tmp_dir, ignore_errors=True)

    info = {
        'n_reads': n_reads,
        'n_rows': n_rows,
        'out_bed': str(out_bed),
        'fastq': str(fastq) if keep_fastq else None,
        'cmd': cmd,
    }

    # success
    return info


def _run_stream(cmd, out_bed):
    """
    Run bowtie2 piped into the BED conversion, writing to a file.

    Args:
        cmd (list): the bowtie2 argv.
        out_bed (pathlib.Path): destination BED path.

    Returns:
        n_rows (int): the number of BED rows written.

    Raises:
        ExternalCommandFailed: if either process exits non-zero.
    """
    with open(out_bed, 'wb') as fout:
        aligner = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
        converter = subprocess.Popen(['awk', AWK_SCRIPT], stdin=aligner.stdout,
                                     stdout=fout, stderr=subprocess.PIPE)

        # close this end so bowtie2 receives SIGPIPE if the converter dies, rather
        # than blocking forever on a full pipe
        aligner.stdout.close()

        converter_err = converter.communicate()[1]
        aligner_err = aligner.stderr.read()
        aligner.stderr.close()
        aligner.wait()

    if aligner.returncode != 0:
        raise ExternalCommandFailed(
            ' '.join(cmd), aligner.returncode,
            stderr=aligner_err.decode(errors='replace'))
    if converter.returncode != 0:
        raise ExternalCommandFailed(
            'awk', converter.returncode,
            stderr=converter_err.decode(errors='replace'))

    with open(out_bed, 'rb') as handle:
        n_rows = sum(1 for _ in handle)

    # success
    return n_rows
