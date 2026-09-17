"""
Specificity analysis sub-package.

Tools for assessing probe specificity through genome alignment (Bowtie2),
kmer frequency analysis (Jellyfish), and duplex stability prediction
(XGBoost and legacy LDA models).
"""

from . import alignment
from . import kmers
from . import duplex_stability
