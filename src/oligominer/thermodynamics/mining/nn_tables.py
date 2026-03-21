"""nn_tables.py — Nearest-neighbor thermodynamic parameter lookup tables.

Stores the unified nearest-neighbor dH (enthalpy) and dS (entropy) parameters
from SantaLucia (1998) PNAS 95:1460-1465 as 4x4 numpy arrays, indexed by the
integer-encoded 5' and 3' bases of each dinucleotide (A=0, C=1, G=2, T=3).

Six lookup tables are provided:
  - DINUC_DH_LUT / DINUC_DS_LUT: base dH and dS for the 16 dinucleotides
  - TERMINAL_5_DH_LUT / TERMINAL_5_DS_LUT: initiation corrections keyed by
    the 5' terminal dinucleotide
  - TERMINAL_3_DH_LUT / TERMINAL_3_DS_LUT: initiation corrections keyed by
    the 3' terminal dinucleotide (transpose of the 5' tables)

The terminal corrections account for the asymmetric initiation parameters
for A/T vs G/C terminal base pairs. All values are sourced from the DNA_NN3
table in Biopython's Bio.SeqUtils.MeltingTemp module.
"""

import numpy as np

# configure datatype for nearest-neighbor Tm calculations
NN_DTYPE = np.float64

#
# all values below are from the DNA_NN3 table in biopython's Bio.SeqUtils.MeltingTemp module:
# https://github.com/biopython/biopython/blob/23aacea3e0ec1efc32ae01d4ccd729c70296fa0d/Bio/SeqUtils/MeltingTemp.py#L207
#

#
# base dH values for the 16 dinucleotides
#

# lookups are done via DINUC_DH_LUT[first_base][second_base] using int encoded sequences
DINUC_DH_LUT = np.array([
    [ -7.9,  -8.4,  -7.8,  -7.2], # AA, AC, AG, AT
    [ -8.5,  -8.0, -10.6,  -7.8], # CA, CC, CG, CT
    [ -8.2,  -9.8,  -8.0,  -8.4], # GA, GC, GG, GT
    [ -7.2,  -8.2,  -8.5,  -7.9], # TA, TC, TG, TT
], dtype=NN_DTYPE)


#
# base dS values for the 16 dinucleotides
#

# lookups are done via DINUC_DH_LUT[first_base][second_base] using int encoded sequences
DINUC_DS_LUT = np.array([
    [-22.2, -22.4, -21.0, -20.4], # AA, AC, AG, AT
    [-22.7, -19.9, -27.2, -21.0], # CA, CC, CG, CT
    [-22.2, -24.4, -19.9, -22.4], # GA, GC, GG, GT
    [-21.3, -22.2, -22.7, -22.2], # TA, TC, TG, TT
], dtype=NN_DTYPE)


#
# terminal nucleotide dH adjustments
#

# table with respect to 5' base in the dinucleotide
TERMINAL_5_DH_LUT = np.array([
    [2.3, 2.3, 2.3, 2.3], #  AA,  AC,  AG,  AT
    [0.1, 0.1, 0.1, 0.1], # *CA, *CC, *CG, *CT
    [0.1, 0.1, 0.1, 0.1], # *GA, *GC, *GG, *GT
    [2.3, 2.3, 2.3, 2.3], #  TA,  TC,  TG,  TT
], dtype=NN_DTYPE)

# table with respect to 3' base in the dinucleotide
#  AA,  AC*, AG*, AT
#  CA,  CC*, CG*, CT
#  GA,  GC*, GG*, GT
#  TA,  TC*, TG*, TT
TERMINAL_3_DH_LUT = TERMINAL_5_DH_LUT.T


#
# terminal nucleotide dS adjustments
#

# table with respect to 5' base in the dinucleotide
TERMINAL_5_DS_LUT = np.array([
    [ 4.1,  4.1,  4.1,  4.1], #  AA,  AC,  AG,  AT
    [-2.8, -2.8, -2.8, -2.8], # *CA, *CC, *CG, *CT
    [-2.8, -2.8, -2.8, -2.8], # *GA, *GC, *GG, *GT
    [ 4.1,  4.1,  4.1,  4.1], #  TA,  TC,  TG,  TT
], dtype=NN_DTYPE)

# table with respect to 3' base in the dinucleotide
#  AA,  AC*, AG*, AT
#  CA,  CC*, CG*, CT
#  GA,  GC*, GG*, GT
#  TA,  TC*, TG*, TT
TERMINAL_3_DS_LUT = TERMINAL_5_DS_LUT.T



