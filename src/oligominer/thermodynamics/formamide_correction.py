"""
Formamide correction for DNA melting temperatures.

Formamide is a denaturant commonly used in FISH hybridization buffers. It
lowers the effective melting temperature of DNA duplexes in a concentration-
dependent manner, allowing hybridization at lower temperatures.
"""


def formamide_correction(temp, pct_fmd, fmd_factor=0.65):
    """
    Apply a formamide adjustment to a DNA melting temperature.

    Args:
        temp (float): initial temperature.
        pct_fmd (int): percent formamide as an integer (e.g. 60 for 60%).
        fmd_factor (float): how much Tm decreases per % formamide.

    Returns:
        temp_adj (float): the formamide-adjusted temperature.
    """
    temp_adj = temp + (pct_fmd * fmd_factor)
    return temp_adj
