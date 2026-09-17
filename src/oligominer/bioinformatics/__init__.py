"""
Bioinformatics sub-package.

Provides I/O utilities for common bioinformatics file formats (FASTA, FASTQ,
BED, SAM/BAM) and sequence classification tools. Transcriptome analysis
utilities are planned for future releases.
"""

from .._lazy import lazy_exports

# submodules reachable as attributes of this package
_SUBMODULES = ('file_io',)

__getattr__, __dir__, __all__ = lazy_exports(__name__, {}, _SUBMODULES)
