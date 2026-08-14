"""
# Transcriptome construction

Builds a transcriptome FASTA from a reference genome and an annotation, so
probes can be screened against transcripts rather than against genomic sequence.

Screening against the genome answers "where else in the genome could this
bind"; screening against a transcriptome answers "which other transcripts could
this bind", which is the question an RNA FISH experiment asks. A probe spanning
an exon-exon junction has no genomic target at all, so it can only be screened
this way.

Each record is one spliced transcript, named by its transcript id, so an aligner
index built from the output reports transcript identifiers directly.
"""

from .transcript_seq import _resolve_fasta, _select_features, get_spliced_seq

# sequence characters per line in the written FASTA
LINE_WIDTH = 60


def transcript_ids(gtf_df, gene_id=None):
    """
    Return the transcript ids present in an annotation.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF carrying transcript_id.
        gene_id (str, optional): restrict to one gene. All genes when None.

    Returns:
        ids (list): the transcript ids, in first-seen order.
    """
    records = gtf_df
    if 'feature' in records.columns:
        records = records[records['feature'] == 'exon']
    if gene_id is not None:
        records = records[records['gene_id'] == gene_id]

    # success
    return list(dict.fromkeys(records['transcript_id']))


def build_transcriptome(gtf_df, fasta, out_fasta, gene_id=None,
                        min_length=1, skip_errors=True):
    """
    Write a FASTA of spliced transcript sequences.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF carrying exon records with
            transcript_id.
        fasta (str or pyfaidx.Fasta): the genome FASTA, or a loaded Fasta.
        out_fasta (str): destination FASTA path.
        gene_id (str, optional): restrict to one gene. All genes when None.
        min_length (int): skip transcripts shorter than this.
        skip_errors (bool): skip a transcript whose sequence cannot be built,
            rather than failing the whole transcriptome. Skipped ids are
            returned so the omission is visible.

    Returns:
        info (dict): out_fasta, n_written, total_bases and skipped ids.
    """
    genome = _resolve_fasta(fasta)
    ids = transcript_ids(gtf_df, gene_id=gene_id)

    n_written = 0
    total_bases = 0
    skipped = []

    with open(out_fasta, 'w') as handle:
        for transcript in ids:
            try:
                sequence = get_spliced_seq(gtf_df, genome, transcript)
            except Exception:
                if not skip_errors:
                    raise
                skipped.append(transcript)
                continue

            if len(sequence) < min_length:
                skipped.append(transcript)
                continue

            handle.write(f'>{transcript}\n')
            for i in range(0, len(sequence), LINE_WIDTH):
                handle.write(sequence[i:i + LINE_WIDTH] + '\n')

            n_written += 1
            total_bases += len(sequence)

    info = {
        'out_fasta': out_fasta,
        'n_written': n_written,
        'total_bases': total_bases,
        'n_skipped': len(skipped),
        'skipped': skipped,
    }

    # success
    return info


def transcript_lengths(gtf_df, gene_id=None):
    """
    Return each transcript's spliced length without building its sequence.

    Args:
        gtf_df (pandas.DataFrame): parsed GTF carrying exon records.
        gene_id (str, optional): restrict to one gene. All genes when None.

    Returns:
        lengths (dict): transcript id mapped to its spliced length in bases.
    """
    records = gtf_df
    if 'feature' in records.columns:
        records = records[records['feature'] == 'exon']
    if gene_id is not None:
        records = records[records['gene_id'] == gene_id]

    lengths = {}
    for transcript, block in records.groupby('transcript_id', sort=False):
        lengths[transcript] = int((block['end'] - block['start'] + 1).sum())

    # success
    return lengths
