"""
End-to-end test of designing probes against a transcriptome.

The package can mine a genome or a transcriptome built from a genome plus its
annotation. Aligning against the transcriptome rather than the genome is what
makes a probe's off-target set the other *transcripts* it would bind, which is
the question an RNA FISH experiment asks. This exercises that whole path so the
claim rests on a run rather than on the pieces existing.
"""

import random
import shutil
import subprocess

import pandas as pd
import pytest

from oligominer.bioinformatics.transcriptome import build_transcriptome
from oligominer.specificity.alignment.process_alignments import (
    process_alignments,
)
from oligominer.specificity.duplex_stability.frames import build_duplex_frame
from oligominer.thermodynamics.mining import mine_fasta

needs_aligner = pytest.mark.skipif(
    shutil.which('bowtie2') is None or shutil.which('samtools') is None,
    reason='bowtie2 and samtools are required')


@pytest.fixture
def annotated_genome(tmp_path):
    """A genome with two transcripts that share one exon."""
    random.seed(23)
    seq = ''.join(random.choice('ACGT') for _ in range(6000))

    genome = tmp_path / 'genome.fa'
    with open(genome, 'w') as handle:
        handle.write('>chr1\n')
        for i in range(0, len(seq), 60):
            handle.write(seq[i:i + 60] + '\n')

    # T1 and T2 share exon 1001-1600; each has a private second exon
    gtf = pd.DataFrame([
        {'seqid': 'chr1', 'type': 'exon', 'start': 1001, 'end': 1600,
         'strand': '+', 'transcript_id': 'T1', 'gene_id': 'G1', 'score': '.'},
        {'seqid': 'chr1', 'type': 'exon', 'start': 2001, 'end': 2600,
         'strand': '+', 'transcript_id': 'T1', 'gene_id': 'G1', 'score': '.'},
        {'seqid': 'chr1', 'type': 'exon', 'start': 1001, 'end': 1600,
         'strand': '+', 'transcript_id': 'T2', 'gene_id': 'G1', 'score': '.'},
        {'seqid': 'chr1', 'type': 'exon', 'start': 3001, 'end': 3600,
         'strand': '+', 'transcript_id': 'T2', 'gene_id': 'G1', 'score': '.'},
    ])

    return str(genome), gtf, seq


class TestBuildTranscriptome:

    def test_one_record_per_transcript(self, annotated_genome, tmp_path):
        genome, gtf, _ = annotated_genome
        out = tmp_path / 'transcriptome.fa'

        info = build_transcriptome(gtf, genome, str(out))
        assert info['n_written'] == 2

    def test_a_transcript_is_its_exons_spliced(self, annotated_genome,
                                               tmp_path):
        genome, gtf, seq = annotated_genome
        out = tmp_path / 'transcriptome.fa'
        build_transcriptome(gtf, genome, str(out))

        from oligominer.bioinformatics.file_io import load_fasta
        built = load_fasta(str(out))

        # 1-based inclusive 1001-1600 and 2001-2600
        assert str(built['T1']) == seq[1000:1600] + seq[2000:2600]

    def test_the_shared_exon_appears_in_both(self, annotated_genome, tmp_path):
        genome, gtf, seq = annotated_genome
        out = tmp_path / 'transcriptome.fa'
        build_transcriptome(gtf, genome, str(out))

        from oligominer.bioinformatics.file_io import load_fasta
        built = load_fasta(str(out))
        shared = seq[1000:1600]

        assert shared in str(built['T1'])
        assert shared in str(built['T2'])


@needs_aligner
class TestDesignAgainstTheTranscriptome:

    def test_a_probe_in_a_shared_exon_hits_both_transcripts(
            self, annotated_genome, tmp_path):
        genome, gtf, _ = annotated_genome
        transcriptome = tmp_path / 'transcriptome.fa'
        build_transcriptome(gtf, genome, str(transcriptome))

        probes = mine_fasta(str(transcriptome), min_length=30, max_length=37,
                            min_tm=0, max_tm=100, min_gc=0, max_gc=100)
        assert probes, 'no probes mined from the transcriptome'

        index = str(tmp_path / 'tx_idx')
        subprocess.run(['bowtie2-build', '-q', str(transcriptome), index],
                       check=True)

        fastq = tmp_path / 'probes.fastq'
        with open(fastq, 'w') as handle:
            for seq_id, start, stop, probe, _ in probes:
                handle.write(f'@{seq_id}:{start}-{stop}\n'
                             f'{probe}\n+\n{"~" * len(probe)}\n')

        aligned = subprocess.run(
            ['bowtie2', '-x', index, '-U', str(fastq), '--very-sensitive-local',
             '--xeq', '-k', '10', '--no-unal'],
            capture_output=True, text=True, check=True)

        df = process_alignments(sam_data=aligned.stdout,
                                ref_fasta=str(transcriptome))

        # the off-target set is now other transcripts, which is the question an
        # RNA experiment asks
        assert set(df['align_seqid']) <= {'T1', 'T2'}

        hits = df.groupby('seqid')['align_seqid'].nunique()
        assert (hits == 2).any(), 'no probe bound both transcripts'

    def test_the_duplex_frame_builds_from_transcriptome_alignments(
            self, annotated_genome, tmp_path):
        genome, gtf, _ = annotated_genome
        transcriptome = tmp_path / 'transcriptome.fa'
        build_transcriptome(gtf, genome, str(transcriptome))

        probes = mine_fasta(str(transcriptome), min_length=30, max_length=37,
                            min_tm=0, max_tm=100, min_gc=0, max_gc=100)

        index = str(tmp_path / 'tx_idx')
        subprocess.run(['bowtie2-build', '-q', str(transcriptome), index],
                       check=True)

        fastq = tmp_path / 'probes.fastq'
        with open(fastq, 'w') as handle:
            for seq_id, start, stop, probe, _ in probes[:40]:
                handle.write(f'@{seq_id}:{start}-{stop}\n'
                             f'{probe}\n+\n{"~" * len(probe)}\n')

        aligned = subprocess.run(
            ['bowtie2', '-x', index, '-U', str(fastq), '--very-sensitive-local',
             '--xeq', '-k', '10', '--no-unal'],
            capture_output=True, text=True, check=True)

        df = process_alignments(sam_data=aligned.stdout,
                                ref_fasta=str(transcriptome))

        lookup = {f'{seq_id}:{start}-{stop}': probe
                  for seq_id, start, stop, probe, _ in probes}
        df['probe_seq'] = df['seqid'].map(lookup)
        df = df.dropna(subset=['probe_seq'])

        frame = build_duplex_frame(df, celsius=69.5, sodium=0.39)
        assert len(frame) > 0
        assert frame['probe_aln'].notna().all()
