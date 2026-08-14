"""
Formamide corrections for DNA melting temperatures.

Formamide is a denaturant commonly used in FISH hybridization buffers. It lowers
the melting temperature of DNA duplexes in a concentration-dependent manner,
allowing hybridization at lower temperatures.

Formamide enters a calculation in two opposite directions, and the two must not
be confused:

- ``formamide_correction`` lowers a melting temperature to say what a duplex's Tm
  becomes once formamide is present. Use it when reporting or filtering on Tm.
- ``effective_hyb_temperature`` raises a hybridization temperature to say which
  formamide-free temperature is thermodynamically equivalent to it. Use it when
  configuring a model that has no formamide term of its own, such as a NUPACK
  model. Hybridizing at 37 C in 50% formamide behaves like hybridizing at 69.5 C
  without it.

Passing a depressed Tm where an effective temperature belongs simulates a far
more permissive condition than the experiment, which inflates every duplex
probability rather than raising an error.
"""


def formamide_correction(temp, pct_fmd, fmd_factor=0.65):
    """
    Lower a melting temperature to account for formamide in the buffer.

    Args:
        temp (float): the melting temperature without formamide, in Celsius.
        pct_fmd (float): percent formamide as a float, e.g. 50.0 for 50%.
        fmd_factor (float): degrees Celsius of Tm depression per percent
            formamide.

    Returns:
        temp_adj (float): the formamide-depressed melting temperature.
    """
    temp_adj = temp - (pct_fmd * fmd_factor)

    # success
    return temp_adj


def effective_hyb_temperature(temp, pct_fmd, fmd_factor=0.65):
    """
    Return the formamide-free temperature equivalent to a hybridization condition.

    A model with no formamide term is run at this temperature to represent
    hybridization carried out in formamide at ``temp``.

    Args:
        temp (float): the hybridization temperature in Celsius.
        pct_fmd (float): percent formamide as a float, e.g. 50.0 for 50%.
        fmd_factor (float): degrees Celsius of Tm depression per percent
            formamide.

    Returns:
        temp_eff (float): the equivalent temperature without formamide.
    """
    temp_eff = temp + (pct_fmd * fmd_factor)

    # success
    return temp_eff
