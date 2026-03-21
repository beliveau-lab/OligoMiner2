"""int_encoding.py — Integer encoding of DNA sequences for vectorized operations.

Converts DNA sequence strings into numpy uint8 arrays using a fixed encoding:
A=0, C=1, G=2, T=3, and all other characters (N, ambiguous IUPAC bases, etc.)
map to 4. This encoding enables fast vectorized lookups into the nearest-neighbor
thermodynamic parameter tables (nn_tables.py) and efficient boolean masking for
N-containing probe filtering.

The encoding is performed via a 256-element ASCII lookup table (DNA_ASCII_LUT)
that maps each possible byte value to its integer code, making the conversion
a single numpy advanced-indexing operation with no Python-level loops.
"""

import numpy as np

# create ascii lookup table to encode fasta characters as ints
# ACGT → 0,1,2,3; all other characters (N, ambiguous bases, etc.) → 4
DNA_ASCII_LUT = np.full(256, 4, dtype=np.uint8)
DNA_ASCII_LUT[[ord(base) for base in 'ACGTacgt']] = [0, 1, 2, 3, 0, 1, 2, 3]

def seq_to_8bit(seq):
    """
    Returns an 8-bit integer representation of the input sequence.
    
    Args:
        seq (str): the input DNA sequence.
    
    Returns:
        nuc_array (numpy.ndarray): the encoded DNA sequence.    
    """    

    # encode fasta sequence as 8bit integer array
    nuc_array = DNA_ASCII_LUT[np.frombuffer(bytes(str(seq), 'utf-8'), dtype=np.uint8)]

    # success
    return nuc_array


