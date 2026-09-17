from ..._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    'align_to_bed':        '.stream_align',
    'bam_to_bed':          '.bam_to_bed',
    'bowtie_align':        '.bowtie_align',
    'bowtie_build':        '.bowtie_build',
    'build_bowtie2_cmd':   '.bowtie_align',
    'check_index_exists':  '.bowtie_align',
    'get_fasta':           '.get_fasta',
    'process_alignments':  '.process_alignments',
    'trim_bed_coords':     '.trim_bed_coords',
    'validate_index':      '.bowtie_align',
    'write_fastq':         '.stream_align',
}

# submodules reachable as attributes of this package
_SUBMODULES = ('bowtie_presets', 'duplex')

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS, _SUBMODULES)
