from ..._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    'classify_and_write':     '.classify',
    'classify_seq_ids':       '.classify',
    'filter_gtf':             '.gtf_io',
    'filter_seq_ids':         '.fasta_io',
    'filter_seqs':            '.fasta_io',
    'get_or_create_fai':      '.chrom_sizes',
    'load_bam_file':          '.sam_bam_io',
    'load_fasta':             '.fasta_io',
    'load_gtf':               '.gtf_io',
    'load_sam_file':          '.sam_bam_io',
    'merge_annotation_beds':  '.gtf_io',
    'merge_fastas':           '.fasta_io',
    'parse_attributes':       '.gtf_io',
    'seqs_to_fasta':          '.fasta_io',
    'seqs_to_fastq':          '.fastq_io',
    'split_fasta':            '.fasta_io',
    'split_gtf':              '.gtf_io',
    'write_bed':              '.gtf_io',
    'write_fasta':            '.fasta_io',
    'write_gtf':              '.gtf_io',
}

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS)
