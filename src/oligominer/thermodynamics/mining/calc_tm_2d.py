"""calc_tm_2d.py — Vectorized nearest-neighbor melting temperature calculation.

Computes a 2D grid of Tm values for all candidate probes in a sequence,
where rows correspond to start positions and columns correspond to probe
lengths (from min_length to max_length). This enables simultaneous
evaluation of every position x length combination in a single pass.

The nearest-neighbor model uses the unified thermodynamic parameters from
SantaLucia (1998) PNAS 95:1460-1465 (the DNA_NN3 table in Biopython).
Salt corrections follow von Ahsen et al. (2001), using a Na-equivalent
concentration that accounts for K+, Tris, Mg2+, and dNTPs.

Tm = (dH * 1000) / (dS + R * ln(C_T)) - 273.15

where dH and dS are the sums of nearest-neighbor enthalpy and entropy
parameters (with terminal corrections), R is the gas constant, and C_T
is the effective total strand concentration.
"""

import numpy as np

from .nn_tables import (
    DINUC_DS_LUT,
    DINUC_DH_LUT,
    TERMINAL_5_DS_LUT,
    TERMINAL_5_DH_LUT,
    TERMINAL_3_DS_LUT,
    TERMINAL_3_DH_LUT
)

def get_dinuc_grid(nuc_array, max_length):
    """Build a 2D grid of dinucleotide pairs for all sliding windows.

    Uses numpy stride tricks (zero-copy views) to create an array where
    each row represents a candidate probe starting position and each column
    contains the dinucleotide pairs within that probe window.

    Args:
        nuc_array (numpy.ndarray): 1D uint8 array of integer-encoded bases
            (A=0, C=1, G=2, T=3, N=4).
        max_length (int): maximum probe length. Determines the number of
            dinucleotide columns (max_length - 1) per window.

    Returns:
        dinuc_grid (numpy.ndarray): array of shape (N_positions, 2, max_length - 1),
            where axis 1 holds the two bases of each dinucleotide pair (index 0 =
            5' base, index 1 = 3' base) and axis 2 spans the dinucleotides
            within the probe window. N_positions = len(nuc_array) - max_length + 1.
    """
    dinuc_array = np.lib.stride_tricks.sliding_window_view(nuc_array, 2)
    dinuc_grid =  np.lib.stride_tricks.sliding_window_view(dinuc_array, max_length - 1, axis=0)
    return dinuc_grid

def get_dS_grid(dinuc_grid, min_length, max_length, Na=50, K=0, Tris=0, Mg=0, dNTPs=0):
    """Compute the entropy (dS) grid for all candidate probes with salt correction.

    Looks up nearest-neighbor dS values from the SantaLucia (1998) table,
    applies 5' and 3' terminal corrections, accumulates via cumulative sum,
    and then applies the salt correction from von Ahsen et al. (2001):

        dS_salt = 0.368 * (L - 1) * ln([Mon])

    where [Mon] is the monovalent cation equivalent concentration (molar)
    and L is the probe length. If divalent ions are present, free Mg2+
    (after dNTP chelation) is converted to a Na-equivalent using:

        [Mon] = [Na] + [K] + [Tris]/2 + 120 * sqrt([Mg] - [dNTPs])

    Args:
        dinuc_grid (numpy.ndarray): dinucleotide grid from get_dinuc_grid(),
            shape (N_positions, 2, max_length - 1).
        min_length (int): minimum probe length (inclusive).
        max_length (int): maximum probe length (inclusive).
        Na (float): sodium concentration in mM (default 50).
        K (float): potassium concentration in mM (default 0).
        Tris (float): Tris buffer concentration in mM (default 0).
        Mg (float): magnesium concentration in mM (default 0).
        dNTPs (float): dNTP concentration in mM (default 0). dNTPs chelate
            Mg2+, reducing the effective free Mg2+ concentration.

    Returns:
        dS_grid (numpy.ndarray): array of shape (N_positions, N_lengths)
            where N_lengths = max_length - min_length + 1. Each value is
            the total dS (cal/mol/K) for the corresponding probe, including
            terminal and salt corrections.

    Raises:
        ValueError: if the total monovalent ion concentration is zero.
    """

    # look up dS values in lookup table
    dS_grid = DINUC_DS_LUT[dinuc_grid[:,0,:],dinuc_grid[:,1,:]]

    # adjust dS values based on 5' termini
    dS_grid[:,0] += TERMINAL_5_DS_LUT[dinuc_grid[:,0,0],dinuc_grid[:,1,0]]

    # compute cumulative sum across each row (in-place to avoid temporary memory allocation)
    np.cumsum(dS_grid, axis=1, out=dS_grid)

    # adjust dS values based on 3' termini, for every probe length to be designed
    for i in range(min_length - 2, (max_length + 1) - 2):
        dS_grid[:,i] += TERMINAL_3_DS_LUT[dinuc_grid[:,0,i],dinuc_grid[:,1,i]]

    # trim to needed columns and copy so the (N, max_length-1) backing array is freed
    dS_grid = dS_grid[:,min_length - 2:].copy()

    #
    # salt corrections to dS values
    #

    # adapted from https://github.com/biopython/biopython/blob/master/Bio/SeqUtils/MeltingTemp.py#L486
    Mon = Na + K + Tris / 2.0  # Note: all these values are millimolar
    mg = Mg * 1e-3  # Lowercase ions (mg, mon, dntps) are molar
    # Na equivalent according to von Ahsen et al. (2001):
    if sum((K, Mg, Tris, dNTPs)) > 0 and dNTPs < Mg:
        # dNTPs bind Mg2+ strongly. If [dNTPs] is larger or equal than
        # [Mg2+], free Mg2+ is considered not to be relevant.
        Mon += 120 * np.sqrt(Mg - dNTPs)
    mon = Mon * 1e-3
    if not mon:
        raise ValueError(
            "Total ion concentration of zero is not allowed in this method."
        )
    seq_lengths = np.arange(min_length, max_length + 1)
    salt_correction = 0.368 * (seq_lengths - 1) * np.log(mon)
    dS_grid += salt_correction

    # success
    return dS_grid

def get_dH_grid(dinuc_grid, min_length, max_length):
    """Compute the enthalpy (dH) grid for all candidate probes.

    Looks up nearest-neighbor dH values from the SantaLucia (1998) table,
    applies 5' and 3' terminal corrections, and accumulates via cumulative
    sum. No salt correction is applied to dH (salt effects are captured
    entirely in dS).

    Args:
        dinuc_grid (numpy.ndarray): dinucleotide grid from get_dinuc_grid(),
            shape (N_positions, 2, max_length - 1).
        min_length (int): minimum probe length (inclusive).
        max_length (int): maximum probe length (inclusive).

    Returns:
        dH_grid (numpy.ndarray): array of shape (N_positions, N_lengths)
            where N_lengths = max_length - min_length + 1. Each value is
            the total dH (kcal/mol) for the corresponding probe, including
            terminal corrections.
    """

    # look up dH values in lookup table
    dH_grid = DINUC_DH_LUT[dinuc_grid[:,0,:],dinuc_grid[:,1,:]]

    # adjust dH values based on 5' termini
    dH_grid[:,0] += TERMINAL_5_DH_LUT[dinuc_grid[:,0,0],dinuc_grid[:,1,0]]

    # compute cumulative sum across each row (in-place to avoid temporary memory allocation)
    np.cumsum(dH_grid, axis=1, out=dH_grid)

    # adjust dH values based on 3' termini, for every probe length to be designed
    for i in range(min_length - 2, (max_length + 1) - 2):
        dH_grid[:,i] += TERMINAL_3_DH_LUT[dinuc_grid[:,0,i],dinuc_grid[:,1,i]]

    # trim to needed columns and copy so the (N, max_length-1) backing array is freed
    dH_grid = dH_grid[:,min_length - 2:].copy()

    # success
    return dH_grid


def get_tm_grid(nuc_array, config):
    """Compute the full 2D Tm grid for all candidate probes in a sequence.

    This is the main entry point for Tm calculation. It orchestrates the
    full nearest-neighbor pipeline:

    1. Build dinucleotide grid from the encoded sequence (get_dinuc_grid)
    2. Compute dS grid with salt corrections (get_dS_grid)
    3. Add the concentration-dependent entropy term: R * ln(C_T)
    4. Compute dH grid with terminal corrections (get_dH_grid)
    5. Calculate Tm = (dH * 1000) / (dS + R*ln(C_T)) - 273.15
    6. Apply formamide correction if configured

    The concentration term uses C_T = [dnac1] - [dnac2]/2 (the effective
    total strand concentration for non-self-complementary sequences).

    Args:
        nuc_array (numpy.ndarray): 1D uint8 array of integer-encoded bases
            (A=0, C=1, G=2, T=3). N positions (encoded as 4) should be
            replaced with 0 before calling — Tm values at those positions
            will be masked out downstream by the N-filter in process_chunk.
        config (dict): mining configuration dict containing at minimum:
            - min_length, max_length: probe length range
            - Na, K, Tris, Mg, dNTPs: ion concentrations in mM
            - dnac1, dnac2: DNA strand concentrations in nM
            - pct_formamide: percent formamide (0 to disable)
            - formamide_factor: degrees C of Tm depression per percent
              formamide (typically 0.65)

    Returns:
        tm (numpy.ndarray): array of shape (N_positions, N_lengths) where
            N_positions = len(nuc_array) - max_length + 1 and
            N_lengths = max_length - min_length + 1. Each value is the
            predicted Tm in degrees Celsius for the corresponding probe.
    """

    # get a grid from the dinucleotide array using a sliding window
    dinuc_grid = get_dinuc_grid(nuc_array, config['max_length'])

    R = 1.987  # universal gas constant in Cal/degrees C*Mol
    k = (config['dnac1'] - (config['dnac2'] / 2.0)) * 1e-9

    # compute dS first, add concentration term in-place
    denom = get_dS_grid(dinuc_grid, config['min_length'], config['max_length'],
                        Na=config['Na'], K=config['K'], Tris=config['Tris'],
                        Mg=config['Mg'], dNTPs=config['dNTPs'])
    denom += R * np.log(k)

    # compute dH, scale and divide in-place to avoid temporary arrays
    tm = get_dH_grid(dinuc_grid, config['min_length'], config['max_length'])
    tm *= 1000.0
    tm /= denom
    del denom
    tm -= 273.15 # convert from Kelvin to Celsius

    if config['pct_formamide']:
        tm -= (config['pct_formamide'] * config['formamide_factor'])

    # success
    return tm
