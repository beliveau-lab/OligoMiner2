"""
Specificity analysis sub-package.

Tools for assessing probe specificity through genome alignment (Bowtie2),
kmer frequency analysis (Jellyfish), and duplex stability prediction
(XGBoost and legacy LDA models).
"""

from .._lazy import lazy_exports

# submodules reachable as attributes of this package
_SUBMODULES = ("alignment", "duplex_stability", "kmers")

__getattr__, __dir__, __all__ = lazy_exports(__name__, {}, _SUBMODULES)
