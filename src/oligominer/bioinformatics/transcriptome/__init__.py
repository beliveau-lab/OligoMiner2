from ..._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    'build_transcriptome':      '.build',
    'discriminating_regions':   '.junctions',
    'exon_order':               '.junctions',
    'flatten_isoforms':         '.iso_flatten',
    'get_exon_seqs':            '.transcript_seq',
    'get_flattened_seqs':       '.transcript_seq',
    'get_intron_seqs':          '.transcript_seq',
    'get_spliced_seq':          '.transcript_seq',
    'junction_offsets':         '.junctions',
    'junction_probes':          '.junctions',
    'junction_windows':         '.junctions',
    'local_to_genomic':         '.transcript_seq',
    'mine_exons':               '.mine_transcripts',
    'mine_flattened_gene':      '.mine_transcripts',
    'mine_introns':             '.mine_transcripts',
    'mine_spliced_transcript':  '.mine_transcripts',
    'parse_interval_label':     '.transcript_seq',
    'spans_junction':           '.junctions',
    'transcript_ids':           '.build',
    'transcript_lengths':       '.build',
}

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS)
