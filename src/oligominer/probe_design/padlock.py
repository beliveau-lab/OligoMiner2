"""
Padlock (split-homology) probe mining, built on OM2's existing 2D Tm grid.

# The topology, and why one contiguous window is enough

A padlock probe is a single linear oligo whose two ends carry the target homology and
whose middle is a non-hybridizing backbone::

    5'--[arm_5p]--[backbone]--[arm_3p]--3'

Probe and target are antiparallel, so the probe's 5' arm pairs with the *downstream*
target segment and its 3' arm with the *upstream* one. Writing the target footprint as
``T1`` then ``T2`` (5'->3' along the target)::

    arm_5p = revcomp(T2)
    arm_3p = revcomp(T1)
    arm_5p + arm_3p == revcomp(T1 + T2)

**The two arms concatenated in probe order are exactly the reverse complement of one
contiguous genomic window.** That is what lets the whole homology go into specificity
analysis as a single probe sequence: the split is bookkeeping (a junction offset and an
optional gap), not a different molecule. Downstream code that screens k-mers, aligns to
a genome or scores a duplex never has to know a padlock is involved.

When the arms are separated by a gap on the target, the identity above holds for the
two arms alone -- ``probe_seq`` is what is *synthesized*, and the gap is filled by
polymerase before ligation. The gap's target sequence is recorded as ``gap_seq`` so the
ligated product can be reconstructed, but this module does not assert which end the
polymerase extends from; that is a chemistry choice, not a mining one. For ``gap == 0``
the arms are directly adjacent and ``probe_seq == revcomp(target window)`` exactly,
which ``check_identity()`` asserts.

# How the arms are mined

Each arm is mined by OM2's own ``process_chunk`` in exhaustive mode, once per arm with
that arm's own length / Tm / GC / homopolymer parameters. Reusing OM2's filter code
rather than reimplementing it means there is no parity risk between padlock mining and
probe mining -- they are the same filters, called twice.

The two candidate sets are then joined on adjacency::

    T2.start == T1.stop + gap

which is a vectorized lookup, not a nested loop.

Note on Tm: nearest-neighbour Tm is invariant under reverse-complement, so an arm's Tm
computed on the target segment is the Tm of the arm/target duplex. The arms are mined in
target coordinates and reverse-complemented once, at emission.
"""

import numpy as np

from oligominer.thermodynamics.mining.config import GET_DEFAULT_MINING_CONFIG
from oligominer.thermodynamics.mining.int_encoding import seq_to_8bit
from oligominer.thermodynamics.mining.mine_probes import process_chunk
from oligominer.utils.seq_utils import rev_comp

# one row per padlock; arm_5p/arm_3p are probe-strand sequences, coordinates are target
PADLOCK_COLUMNS = [
    "seq_id", "start", "stop", "probe_seq",
    "arm_5p", "arm_3p", "arm5_len", "arm3_len", "gap_len", "gap_seq",
    "arm5_tm", "arm3_tm", "junction",
]

DEFAULT_ARM = {
    "min_length": 15,
    "max_length": 25,
    "min_tm": 45,
    "max_tm": 65,
    "tm_target": None,
    "min_gc": 35,
    "max_gc": 65,
    "max_homopolymer": 4,
    "prohibited_seqs": None,
}


def arm_params(**overrides):
    """
    Build one arm's parameter dict, starting from the padlock defaults.

    Args:
        **overrides: any of min_length, max_length, min_tm, max_tm, tm_target,
            min_gc, max_gc, max_homopolymer, prohibited_seqs.

    Returns:
        params (dict): a complete arm parameter dict.
    """
    params = dict(DEFAULT_ARM)
    unknown = set(overrides) - set(DEFAULT_ARM)
    if unknown:
        raise ValueError(f"unknown arm parameter(s): {sorted(unknown)}")
    params.update(overrides)

    # success
    return params


def _arm_config(params, thermo):
    """
    Turn an arm parameter dict into an OM2 mining config in exhaustive mode.

    Args:
        params (dict): as returned by arm_params().
        thermo (dict): salt / concentration / formamide settings shared by both arms.

    Returns:
        config (dict): a mining config accepted by process_chunk().
    """
    config = GET_DEFAULT_MINING_CONFIG()
    config.update(params)
    config.update(thermo)

    # exhaustive so every valid (start, length) survives for the adjacency join;
    # selection between competing arms happens after the join, not before it
    config["exhaustive"] = True

    if params.get("prohibited_seqs"):
        config["_prohibited_encoded"] = [seq_to_8bit(p) for p in params["prohibited_seqs"]]

    # success
    return config


def _candidates(nuc_arr, config, offset):
    """
    Run OM2's chunk filters and return every valid (start, length, tm) for one arm.

    Args:
        nuc_arr (numpy.ndarray): 1D uint8 encoded sequence.
        config (dict): an exhaustive-mode mining config.
        offset (int): coordinate of nuc_arr[0] in the parent sequence.

    Returns:
        starts (numpy.ndarray): candidate start coordinates.
        lengths (numpy.ndarray): candidate lengths.
        tms (numpy.ndarray): candidate melting temperatures.
    """
    coords, tms = process_chunk("arm", nuc_arr, offset, offset + nuc_arr.size, config)
    if len(coords) == 0:
        empty = np.zeros(0, dtype=int)
        return empty, empty, np.zeros(0)

    starts = coords[:, 0]
    lengths = coords[:, 1] - coords[:, 0]

    # success
    return starts, lengths, tms


def _best_per_start(starts, lengths, tms, tm_target, span):
    """
    Reduce many candidates per start position to one, and index them by position.

    Args:
        starts (numpy.ndarray): candidate starts.
        lengths (numpy.ndarray): candidate lengths.
        tms (numpy.ndarray): candidate Tms.
        tm_target (float or None): if set, keep the candidate whose Tm is closest to it;
            otherwise keep the shortest.
        span (int): size of the position-indexed output arrays.

    Returns:
        has (numpy.ndarray): bool array, True where a candidate starts.
        length_at (numpy.ndarray): chosen length per position (0 where none).
        tm_at (numpy.ndarray): chosen Tm per position (nan where none).
    """
    has = np.zeros(span, dtype=bool)
    length_at = np.zeros(span, dtype=int)
    tm_at = np.full(span, np.nan)

    if len(starts) == 0:
        return has, length_at, tm_at

    # rank so the preferred candidate for each start sorts first, then take the first
    # occurrence of each start -- one pass, no grouping loop
    criterion = np.abs(tms - tm_target) if tm_target is not None else lengths.astype(float)
    order = np.lexsort((criterion, starts))
    s_sorted = starts[order]
    first = np.ones(len(s_sorted), dtype=bool)
    first[1:] = s_sorted[1:] != s_sorted[:-1]
    keep = order[first]

    in_range = (starts[keep] >= 0) & (starts[keep] < span)
    keep = keep[in_range]

    has[starts[keep]] = True
    length_at[starts[keep]] = lengths[keep]
    tm_at[starts[keep]] = tms[keep]

    # success
    return has, length_at, tm_at


def mine_padlock_sequence(seq, seq_id="seq", arm5=None, arm3=None, gap=0,
                          spacing=0, overlap=False, tm_target_join=None,
                          Na=390, K=0, Tris=0, Mg=0, dNTPs=0,
                          dnac1=25, dnac2=25, pct_formamide=0,
                          formamide_factor=0.65):
    """
    Mine padlock probes -- two homology arms, optionally separated by a fillable gap.

    Each arm is mined with its own length, Tm and GC constraints, then arms are paired
    by adjacency on the target. The emitted ``probe_seq`` is the full homology as one
    sequence, ready for k-mer screening, alignment and duplex scoring with no padlock
    awareness anywhere downstream.

    Args:
        seq (str): the target DNA sequence (ACGTN).
        seq_id (str): identifier for this sequence.
        arm5 (dict or None): parameters for the probe's 5' arm, from arm_params().
            This arm binds the DOWNSTREAM target segment. Defaults to DEFAULT_ARM.
        arm3 (dict or None): parameters for the probe's 3' arm, which binds the
            UPSTREAM target segment. Defaults to DEFAULT_ARM.
        gap (int): bases on the target between the two arms, filled by polymerase
            before ligation. 0 makes the arms directly adjacent.
        spacing (int): minimum bases between the footprints of adjacent padlocks.
        overlap (bool): if False, padlock footprints may not overlap.
        tm_target_join (float or None): when several arm lengths are valid at a
            position, prefer the one whose Tm is closest to this. None picks shortest.
        Na (float): sodium concentration in mM.
        K (float): potassium concentration in mM.
        Tris (float): Tris concentration in mM.
        Mg (float): magnesium concentration in mM.
        dNTPs (float): dNTP concentration in mM.
        dnac1 (float): probe strand concentration in nM.
        dnac2 (float): target strand concentration in nM.
        pct_formamide (int): percent formamide in the hybridization buffer. **Defaults
            to 0 here, unlike ``mine_sequence``, which defaults to 50.** Padlock
            hybridization and ligation is not run in a 50% formamide FISH buffer, and
            the correction is large enough to change the answer rather than shade it:
            measured on sacCer3 chrI, 16-24mers span Tm 7.4-37.2 C at 50% formamide and
            39.9-69.7 C at 0%. An arm Tm window written for one scale silently returns
            zero padlocks on the other.
        formamide_factor (float): degrees C of Tm depression per percent formamide.

    Returns:
        rows (list): list of tuples in PADLOCK_COLUMNS order.
    """
    if gap < 0:
        raise ValueError(f"gap must be >= 0, got {gap}")

    arm5 = dict(DEFAULT_ARM) if arm5 is None else arm5
    arm3 = dict(DEFAULT_ARM) if arm3 is None else arm3

    thermo = {
        "Na": Na, "K": K, "Tris": Tris, "Mg": Mg, "dNTPs": dNTPs,
        "dnac1": dnac1, "dnac2": dnac2,
        "pct_formamide": pct_formamide, "formamide_factor": formamide_factor,
    }

    seq_str = seq.upper()
    nuc_arr = seq_to_8bit(seq_str)
    span = nuc_arr.size

    # T1 (upstream on the target) becomes the probe's 3' arm; T2 (downstream) the 5' arm
    s1, l1, tm1 = _candidates(nuc_arr, _arm_config(arm3, thermo), 0)
    s2, l2, tm2 = _candidates(nuc_arr, _arm_config(arm5, thermo), 0)

    if len(s1) == 0 or len(s2) == 0:
        return []

    tj = tm_target_join if tm_target_join is not None else arm5.get("tm_target")
    has2, len2_at, tm2_at = _best_per_start(s2, l2, tm2, tj, span)

    # the adjacency join: a T1 candidate is usable iff a T2 candidate starts exactly
    # gap bases after it ends
    junction = s1 + l1 + gap
    in_bounds = junction < span
    s1, l1, tm1, junction = s1[in_bounds], l1[in_bounds], tm1[in_bounds], junction[in_bounds]

    paired = has2[junction]
    s1, l1, tm1, junction = s1[paired], l1[paired], tm1[paired], junction[paired]
    if len(s1) == 0:
        return []

    l2_sel = len2_at[junction]
    tm2_sel = tm2_at[junction]

    stop = junction + l2_sel                      # end of the full target footprint
    fits = stop <= span
    s1, l1, tm1, junction = s1[fits], l1[fits], tm1[fits], junction[fits]
    l2_sel, tm2_sel, stop = l2_sel[fits], tm2_sel[fits], stop[fits]
    if len(s1) == 0:
        return []

    # one padlock per start: prefer the shortest total footprint, then the one whose
    # 3' arm Tm best matches the request
    total = stop - s1
    tie = np.abs(tm1 - arm3["tm_target"]) if arm3.get("tm_target") is not None else total
    order = np.lexsort((tie, total, s1))
    s_sorted = s1[order]
    first = np.ones(len(s_sorted), dtype=bool)
    first[1:] = s_sorted[1:] != s_sorted[:-1]
    pick = order[first]

    # emit in coordinate order, applying the overlap / spacing policy greedily
    pick = pick[np.argsort(s1[pick], kind="stable")]

    rows = []
    last_stop = -1
    for i in pick:
        t_start, t_stop = int(s1[i]), int(stop[i])
        if not overlap and t_start < last_stop:
            continue
        if spacing > 0 and t_start < last_stop + spacing:
            continue

        t1_stop = t_start + int(l1[i])
        t1 = seq_str[t_start:t1_stop]
        t2 = seq_str[int(junction[i]):t_stop]
        gap_seq = seq_str[t1_stop:int(junction[i])]
        if "N" in t1 or "N" in t2 or "N" in gap_seq:
            continue

        arm_3p = rev_comp(t1)
        arm_5p = rev_comp(t2)

        rows.append((
            seq_id, t_start, t_stop,
            arm_5p + arm_3p,                    # the synthesized homology, probe 5'->3'
            arm_5p, arm_3p,
            int(l2_sel[i]), int(l1[i]), int(gap), gap_seq,
            round(float(tm2_sel[i]), 2), round(float(tm1[i]), 2),
            int(len(arm_5p)),                   # junction offset within probe_seq
        ))
        last_stop = t_stop

    # success
    return rows


def padlocks_to_df(rows):
    """
    Convert padlock tuples to a DataFrame.

    Args:
        rows (list): tuples in PADLOCK_COLUMNS order.

    Returns:
        df (pandas.DataFrame): one row per padlock.
    """
    import pandas as pd

    # success
    return pd.DataFrame(rows, columns=PADLOCK_COLUMNS)


def check_identity(row, target_seq):
    """
    Assert the invariant that makes downstream padlock-blindness safe.

    Verifies ``arm_5p + arm_3p == revcomp(T1 + T2)`` for a gapless padlock, and that
    the junction offset splits probe_seq back into the two arms.

    Args:
        row (dict or pandas.Series): one padlock record.
        target_seq (str): the sequence the padlock was mined from.

    Returns:
        ok (bool): True if every invariant holds.

    Raises:
        AssertionError: naming the invariant that failed.
    """
    probe = row["probe_seq"]
    j = int(row["junction"])

    assert probe[:j] == row["arm_5p"], "junction does not recover arm_5p"
    assert probe[j:] == row["arm_3p"], "junction does not recover arm_3p"
    assert len(row["arm_5p"]) == int(row["arm5_len"]), "arm5_len disagrees with arm_5p"
    assert len(row["arm_3p"]) == int(row["arm3_len"]), "arm3_len disagrees with arm_3p"

    if int(row["gap_len"]) == 0:
        window = target_seq[int(row["start"]):int(row["stop"])].upper()
        assert probe == rev_comp(window), (
            "arms concatenated are NOT revcomp of the contiguous target window -- "
            "the whole padlock-blind downstream path depends on this")

    # success
    return True
