
import nupack

from .config import DEFAULT_NUPACK_MODEL
from oligominer.utils import seq_utils


def calc_pdup(seq_a, seq_b=None, conc_a=1e-6, conc_b=1e-12, model=None):
    """
    Calculate the probability of duplex formation between two sequences.

    Uses a NUPACK test tube simulation to determine the fraction of the
    lower-concentration strand that is duplexed at equilibrium. When seq_b
    is omitted, the reverse complement of seq_a is used.

    Args:
        seq_a (str): the first DNA sequence.
        seq_b (str, optional): the second DNA sequence. Defaults to the
            reverse complement of seq_a.
        conc_a (float): molar concentration of seq_a.
        conc_b (float): molar concentration of seq_b.
        model (nupack.Model, optional): a NUPACK thermodynamic model.
            Defaults to DEFAULT_NUPACK_MODEL.

    Returns:
        pdup (float): fraction of the lower-concentration strand in duplex form.
    """
    if not seq_a or (seq_b is not None and not seq_b):
        raise ValueError("sequences must be non-empty strings")

    if model is None:
        model = DEFAULT_NUPACK_MODEL

    if seq_b is None:
        seq_b = seq_utils.rev_comp(seq_a)

    # define strands
    a = nupack.Strand(seq_a, name="a")
    b = nupack.Strand(seq_b, name="b")

    # define duplex
    ab = nupack.Complex([a, b], name="ab")

    # define test tube
    tube = nupack.Tube(
        {a: conc_a, b: conc_b},
        complexes=nupack.SetSpec(max_size=2, include=[ab]),
        name="tube"
    )

    # run nupack test tube simulation
    tube_result = nupack.tube_analysis(tubes=[tube], model=model)

    # get duplex concentration
    ab_conc = tube_result.tubes[tube].complex_concentrations[ab]
        
    # calculate pdup as the fraction of the lower-concentration strand in duplex form
    pdup = ab_conc / min(conc_a, conc_b)

    # success
    return pdup




def calc_competitive_pdup(seq_a, seq_b, target_seq, conc_a=1e-6, conc_b=1e-6, target_conc=1e-12, model=None):
    """
    Calculate the probability of duplex formation for two probes competing
    for the same target sequence.

    Extends the pDup concept to a competitive binding scenario: two probe
    sequences, seq_a and seq_b, both expected to hybridize to target_seq,
    are placed in the same NUPACK test tube simulation. The partition
    function is used to determine the equilibrium concentration of each
    probe-target duplex, yielding the fraction of target bound by each
    probe. This is useful for assessing probe specificity in the presence
    of a competing sequence, e.g. designing competitive blockers to
    suppress off-target binding.

    Args:
        seq_a (str): the first probe DNA sequence.
        seq_b (str): the second probe DNA sequence.
        target_seq (str): the target DNA sequence that both probes are
            expected to hybridize to.
        conc_a (float): molar concentration of seq_a.
        conc_b (float): molar concentration of seq_b.
        target_conc (float): molar concentration of target_seq.
        model (nupack.Model, optional): a NUPACK thermodynamic model.
            Defaults to DEFAULT_NUPACK_MODEL.

    Returns:
        at_pdup (float): fraction of the minority strand in the
            seq_a–target duplex at equilibrium.
        bt_pdup (float): fraction of the minority strand in the
            seq_b–target duplex at equilibrium.
    """
    if not seq_a or not seq_b or not target_seq:
        raise ValueError("sequences must be non-empty strings")

    if model is None:
        model = DEFAULT_NUPACK_MODEL

    # define strands
    a = nupack.Strand(seq_a, name="a")
    b = nupack.Strand(seq_b, name="b")
    target = nupack.Strand(target_seq, name="target")

    # define duplex
    # duplex_ab = nupack.Complex([a, b], name="ab")
    duplex_at = nupack.Complex([a, target], name="at")
    duplex_bt = nupack.Complex([b, target], name="bt")

    # define test tube
    tube = nupack.Tube(
        {a: conc_a, b: conc_b, target: target_conc},
        complexes=nupack.SetSpec(max_size=2, include=[duplex_at, duplex_bt]),
        name="tube"
    )

    # run nupack test tube simulation
    tube_result = nupack.tube_analysis(tubes=[tube], model=model)

    # get duplex concentrations
    at_conc = tube_result.tubes[tube].complex_concentrations[duplex_at]
    bt_conc = tube_result.tubes[tube].complex_concentrations[duplex_bt]

    # calculate pdup as the fraction of the lower-concentration strand in duplex form
    at_pdup = at_conc / min(conc_a, target_conc)
    bt_pdup = bt_conc / min(conc_b, target_conc)

    # success
    return at_pdup, bt_pdup
