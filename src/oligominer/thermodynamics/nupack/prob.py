
import nupack

from .config import DEFAULT_NUPACK_MODEL


def calc_prob(seq, model=None):
    """
    Calculate the probability that a sequence is entirely unpaired.

    Returns the NUPACK structure probability for the all-unpaired
    secondary structure (all dots).

    Args:
        seq (str): the input DNA sequence.
        model (nupack.Model, optional): a NUPACK thermodynamic model.
            Defaults to DEFAULT_NUPACK_MODEL.

    Returns:
        prob (float): probability of the all-unpaired structure.
    """
    if model is None:
        model = DEFAULT_NUPACK_MODEL

    prob = nupack.structure_probability(
        strands=[seq], structure='.' * len(seq), model=model
    )

    # success
    return prob
