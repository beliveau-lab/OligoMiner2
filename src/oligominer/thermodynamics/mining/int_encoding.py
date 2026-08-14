"""int_encoding.py — Integer encoding of DNA sequences for vectorized operations.

Converts DNA sequence strings into numpy uint8 arrays using a fixed encoding:
A=0, C=1, G=2, T=3, and all other characters (N, ambiguous IUPAC bases, etc.)
map to 4. This encoding enables fast vectorized lookups into the nearest-neighbor
thermodynamic parameter tables (nn_tables.py) and efficient boolean masking for
N-containing probe filtering.

The encoding is performed via a 256-element ASCII lookup table that maps each
possible byte value to its integer code, making the conversion a single numpy
advanced-indexing operation with no Python-level loops.

Two tables are provided. DNA_ASCII_LUT folds case, so soft-masked sequence
encodes identically to unmasked sequence. SOFTMASK_ASCII_LUT maps lowercase
bases to 4 instead, which the miner's ambiguous-base filter then excludes.
"""

import numpy as np

# ACGT → 0,1,2,3 in either case; all other characters (N, ambiguous bases) → 4
DNA_ASCII_LUT = np.full(256, 4, dtype=np.uint8)
DNA_ASCII_LUT[[ord(base) for base in 'ACGTacgt']] = [0, 1, 2, 3, 0, 1, 2, 3]

# upper-case ACGT → 0,1,2,3; lowercase (soft-masked) bases join N at 4
SOFTMASK_ASCII_LUT = np.full(256, 4, dtype=np.uint8)
SOFTMASK_ASCII_LUT[[ord(base) for base in 'ACGT']] = [0, 1, 2, 3]


def seq_to_8bit(seq, mask_soft=False):
    """
    Returns an 8-bit integer representation of the input sequence.

    Args:
        seq (str): the input DNA sequence.
        mask_soft (bool): if True, treat lowercase (soft-masked) bases as
            ambiguous so the miner excludes them. If False, case is folded and
            masked sequence is mined like any other.

    Returns:
        nuc_array (numpy.ndarray): the encoded DNA sequence.
    """
    lut = SOFTMASK_ASCII_LUT if mask_soft else DNA_ASCII_LUT

    # encode fasta sequence as 8bit integer array
    nuc_array = lut[np.frombuffer(bytes(str(seq), 'utf-8'), dtype=np.uint8)]

    # success
    return nuc_array


