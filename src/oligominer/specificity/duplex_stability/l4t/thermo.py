"""Nearest-neighbour duplex free energy as an ML feature, parameterised from dna04.2 itself.

WHY THIS EXISTS. Every model currently in the zoo is built on bowtie2's alignment score, and
that score is *representationally* incapable of the thing we are predicting. `--ma` is a flat
per-base bonus: a 20/20 all-AT match and a 20/20 all-GC match score IDENTICALLY. At 74.5 C
those two duplexes sit tens of kcal/mol apart -- one is not bound at all and the other is
essentially irreversible. No amount of model capacity recovers a distinction the input never
encoded. This module supplies the missing axis.

PARAMETER SOURCE -- THE POINT OF THIS MODULE. The numbers are read at import time out of the
NUPACK 4.0.2.0 `dna04.2` parameter file, which is this repo's canonical ground truth
(`knowledge/concepts/nupack-ground-truth.md`) and the very parameter set the pDup labels we
train on were computed under. They are NOT hand-transcribed out of a paper. A transcription
drifts; a load does not. Equally important, loading the file gets us the tables a paper
transcription cannot have: the 36-entry stack table including G.T-containing contexts, and the
96-entry `terminal_mismatch` table -- which is precisely one of the tables dna04.2 rewrote
wholesale relative to dna04.1.

Note `default_wobble_pairing` is False for DNA in this file: G.T is a MISMATCH here, not a pair.
It appears in the stack table only as the context of a mismatch column, never as a helix step.

ON THE NN BACKBONE BEING UNCHANGED. Sourcing from dna04.2 does not conflict with reproducing
older labels' stacking. `studies/20260718_nupack_nn_table_provenance/` finding 5 established
that the nearest-neighbour stacks are 36/36 IDENTICAL across dna1998 == dna04 == dna04.2, i.e.
SantaLucia (1998) Table 1 *is* dna04.2's stack backbone. Rather than take that on faith we keep
a 10-value SantaLucia cross-check table below and ASSERT the loaded values against it in the
self-check -- turning that study finding into a live, re-runnable tripwire. If some future
NUPACK release ever does move the stacks, this module fails loudly instead of silently changing
what a trained model was fit on.

dS IS RE-DERIVED, NOT SOURCED. NUPACK stores dG37 and dH; dS is implied. We recover it as
dS = (dH - dG37) / 310.15 and then extrapolate dG(T) = dH - T*dS. That is the standard
two-state assumption (dH, dS temperature-independent) and it is an inference we make, not a
number the file hands us. Everything the file does hand us is used verbatim.

WHAT THIS IS STILL NOT. Only Watson-Crick helices and their terminating columns are modelled.
NOT modelled: bulge loops, interior loops, hairpins, multiloops, coaxial stacking between
adjacent helices, dangling ends, and the whole partition function over alternative structures.
An 'I'/'D' bulge column and an 'S' soft-clip column are charged NOTHING. The returned number is
therefore NOT a NUPACK dG_assoc and must never be reported or thresholded as one -- for any
imperfect duplex it is biased toward "too stable", one-directionally, because every omitted
term is destabilising. It is an ML FEATURE: a monotone-ish measure of how much good helix is
present, computed under the same parameters as the label.

Input is the column form from `studies/20260720_a_realgenomic_dataset_construction/code/duplex.py`.
Its representation is SAMENESS, not complementarity -- see `_bottom` and `_run_energy`.

PROVENANCE
    studies/20260718_nupack_nn_table_provenance/data/nupack4020_dna04.2.json
    (NUPACK 4.0.2.0, material DNA, parameter set dna04.2; committed in-repo so this module
     never reaches into site-packages and never depends on an installed NUPACK.)
"""
from __future__ import annotations

import json
import math
import os
from functools import lru_cache
from pathlib import Path
from typing import Iterable, NamedTuple

import pandas as pd

GAP = "-"

# ---------------------------------------------------------------------------------------------
# Parameter loading
# ---------------------------------------------------------------------------------------------
# The parameter file is the FIRST-CLASS swappable knob for OM2's DNA params. Default: the corrected,
# literature-verified dna04.2 (studies/20260718 proved its numbers ARE the public NNDB/SantaLucia
# tables). Override with env var OM2_THERMO_PARAMS=/path/to/params.json to run under a different set
# (e.g. the classic dna1998/dna04.1 for legacy reproduction, or a lit-typed file ≡ dna04.2). The
# INVARIANT still holds: whatever the FEATURES use, the LABELS must be generated under the same set.
# VENDORED into this study 2026-07-27. The parameter file lives in THIS directory (`data/`), not in
# a sibling study, so the directory is self-contained: copy it out of the repo and the physics still
# loads. Byte-identical to `studies/20260718_nupack_nn_table_provenance/data/nupack4020_dna04.2.json`
# (sha256 asserted by `code/verify_vendored.py`), which is the NUPACK 4.0.2.0 dna04.2 export the
# pDup labels were themselves computed under. Same parameters for features and labels is the
# invariant; vendoring the file is what makes it checkable from inside this directory alone.
_DEFAULT_PARAM_PATH = Path(__file__).resolve().parent / "nupack4020_dna04.2.json"
PARAM_PATH = Path(os.environ.get("OM2_THERMO_PARAMS", _DEFAULT_PARAM_PATH))

# The reference temperature dG37 in the file is quoted at. dS derivation hangs off this exact
# value; do not "round" it to 310.
T37 = 310.15


def _load_params(path: Path = PARAM_PATH) -> dict:
    if not path.exists():
        raise FileNotFoundError(
            f"dna04.2 parameter file not found at {path}.\n"
            "model_zoo/thermo.py is parameterised directly from the committed NUPACK 4.0.2.0\n"
            "artifact in studies/20260718_nupack_nn_table_provenance/data/. Restore that file\n"
            "(it is tracked in-repo) -- do NOT substitute a copy from site-packages or a\n"
            "hand-typed table, because the whole point of this module is that the feature and\n"
            "the label come from the same parameter set."
        )
    with path.open() as fh:
        d = json.load(fh)
    for section in ("dG", "dH"):
        for table in ("stack", "terminal_penalty", "terminal_mismatch", "join_penalty"):
            if table not in d.get(section, {}):
                raise KeyError(f"{path.name} is missing d['{section}']['{table}'] -- "
                               "not a dna04.2-shaped parameter file.")
    return d


_PARAMS = _load_params()

# --- key conventions, established by inspecting the file (see module self-check) --------------
#
# A stack key is FOUR characters describing two adjacent duplex columns:
#
#       5' - k0  k1 - 3'      top strand
#       3' - k3  k2 - 5'      bottom strand
#
# so k0 pairs with k3, k1 pairs with k2, and the key reads "top 5'->3' then bottom 3'->5'
# backwards", i.e. `key = step + revcomp(step)` for a Watson-Crick step. That is why the table
# has 36 entries and not 16: the same layout also encodes columns containing G.T.
#
# `terminal_mismatch` uses the IDENTICAL layout, with the constraint that the INNER column
# (k1,k2) is the helix's closing pair (one of AT/TA/CG/GC/GT/TG -> 6 options) and the OUTER
# column (k0,k3) is the mismatched column sitting 5'-of-the-closing-pair on the top strand
# (all 16 combinations) -> 6 * 16 = 96 entries. A helix's other end is handled by rotating the
# two columns 180 degrees, which maps key k0k1k2k3 -> k2k3k0k1; the stack table is exactly
# invariant under that rotation, which is the check that the layout is read correctly.
#
# `terminal_penalty` is keyed on the TWO characters of a single terminal column (top, bottom):
# 0.05 kcal/mol dG37 (2.2 dH) for A.T and G.T ends, 0.0 for G.C.
STACK_dG = _PARAMS["dG"]["stack"]
STACK_dH = _PARAMS["dH"]["stack"]
TERM_MM_dG = _PARAMS["dG"]["terminal_mismatch"]
TERM_MM_dH = _PARAMS["dH"]["terminal_mismatch"]
TERM_PEN_dG = _PARAMS["dG"]["terminal_penalty"]
TERM_PEN_dH = _PARAMS["dH"]["terminal_penalty"]

# Charged once per helix, the cost of bringing two strands together (rotational/translational
# entropy). NUPACK charges it per join of two strands; we charge it per helix, which is the
# nearest available analogue in a model that has no loop terms to attach initiation to.
JOIN_dG = float(_PARAMS["dG"]["join_penalty"])
JOIN_dH = float(_PARAMS["dH"]["join_penalty"])

_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N", GAP: GAP}


def _revcomp(s: str) -> str:
    return "".join(_COMPLEMENT[c] for c in reversed(s))


def _dS(dH: float, dG37: float) -> float:
    """Implied entropy, kcal/(mol*K). See module docstring: this is a re-derivation."""
    return (dH - dG37) / T37


# Precompute (dH, dS) per stack key so the hot loop is one dict lookup and two adds.
_STACK_HS = {k: (STACK_dH[k], _dS(STACK_dH[k], STACK_dG[k])) for k in STACK_dG}
_TERM_MM_HS = {k: (TERM_MM_dH[k], _dS(TERM_MM_dH[k], TERM_MM_dG[k])) for k in TERM_MM_dG}
# See `_terminal_dG`: +10 in this table means "this context is not a terminal mismatch", not
# "this context costs 10 kcal/mol". Detect it once here rather than magic-numbering it inline.
_TERM_MM_PROHIBITED = {k for k, v in TERM_MM_dG.items() if v >= 9.9}
_TERM_PEN_HS = {k: (TERM_PEN_dH[k], _dS(TERM_PEN_dH[k], TERM_PEN_dG[k])) for k in TERM_PEN_dG}
_JOIN_HS = (JOIN_dH, _dS(JOIN_dH, JOIN_dG))

# ---------------------------------------------------------------------------------------------
# SantaLucia (1998) PNAS 95:1460-1465 Table 1, unified parameters, 1 M NaCl: (dH, dG37).
# THIS IS NOT USED TO COMPUTE ANYTHING. It exists solely so the self-check can assert that the
# loaded dna04.2 stacks still equal the 1998 unified set -- study 20260718 finding 5, kept as a
# live tripwire rather than a claim in a markdown file. If this assertion ever fires, a NUPACK
# release has moved the NN backbone and every feature/label comparison in the repo needs review.
# ---------------------------------------------------------------------------------------------
_SANTALUCIA_1998 = {
    "AA": (-7.9, -1.00),
    "AT": (-7.2, -0.88),
    "TA": (-7.2, -0.58),
    "CA": (-8.5, -1.45),
    "GT": (-8.4, -1.44),
    "CT": (-7.8, -1.28),
    "GA": (-8.2, -1.30),
    "CG": (-10.6, -2.17),
    "GC": (-9.8, -2.24),
    "GG": (-8.0, -1.84),
}


def santalucia_crosscheck(tol: float = 0.01) -> list[str]:
    """Compare the loaded dna04.2 stacks against SantaLucia 1998. Returns a list of mismatches."""
    bad = []
    for step, (dh, dg) in _SANTALUCIA_1998.items():
        key = step + _revcomp(step)
        if abs(STACK_dG[key] - dg) > tol:
            bad.append(f"{key} dG37 {STACK_dG[key]} != SL98 {dg}")
        if abs(STACK_dH[key] - dh) > tol:
            bad.append(f"{key} dH {STACK_dH[key]} != SL98 {dh}")
    return bad


# Salt correction coefficient, SantaLucia (1998) eq. 6: dS_salt = dS + 0.368 * N * ln[Na+].
# CONVENTION: N is counted PER PHOSPHATE, and we take N = (n_paired - 1) -- the number of
# phosphates per strand in an n_paired helix -- rather than the total over both strands. This is
# the convention the oligo-melting literature uses when quoting per-duplex values, and it is what
# reproduces published Tm's for the unified parameters. It matters that this is stated: the
# 2*(n-1) reading would roughly double the (destabilising, since ln[Na+] < 1 M is negative) salt
# term. NOTE this term is NOT from the parameter file -- NUPACK 4 handles salt separately -- so
# it is the one place a non-dna04.2 number enters. Units here are kcal/(mol*K), hence the 1e-3.
SALT_COEF = 0.368e-3

# Defaults are the assay conditions this repo's pDup ground truth is computed at: 74.5 C
# hybridisation, 0.390 M Na+. Keep them in sync with the NUPACK settings in the nupack-pdup
# skill or the feature will describe a different experiment than the label does.
DEFAULT_CELSIUS = 74.5
DEFAULT_SODIUM = 0.390


class _Helix(NamedTuple):
    """One contiguous run of '=' columns, already scored."""
    dG: float
    n_bp: int
    gc_stacks: int
    at_stacks: int
    term_mm_dG: float     # of dG, how much came from terminal-mismatch columns
    n_term_mm: int


def _runs(ops: str) -> Iterable[tuple[int, int]]:
    """Yield (start, stop) half-open spans of contiguous '=' columns.

    Anything that is not '=' ends a helix: 'X' is a mismatch, 'I'/'D' are bulges, 'S' is a base
    that is physically present but outside the seed. All three interrupt base stacking, so none
    of them may be walked through -- a helix is by construction uninterrupted.
    """
    i, n = 0, len(ops)
    while i < n:
        if ops[i] != "=":
            i += 1
            continue
        j = i
        while j < n and ops[j] == "=":
            j += 1
        yield i, j
        i = j


def _bottom(target_char: str) -> str:
    """The physical bottom-strand base of a column, given what `duplex.py` stored.

    SAMENESS: a matching column holds two IDENTICAL characters, because target_seq was oriented
    onto the probe. So the base actually facing the probe is the COMPLEMENT of the stored target
    character -- for an '=' column that recovers the Watson-Crick partner, and the same rule is
    what makes an 'X' column's real mismatched pair recoverable at all.
    """
    return _COMPLEMENT.get(target_char, "N")


def _terminal_dG(probe_aln: str, target_aln: str, start: int, stop: int,
                 ops: str, five_prime: bool) -> tuple[float, float, bool]:
    """(dH, dS, is_mismatch) for ONE end of the helix spanning [start, stop).

    THE BRANCH. A helix end is one of two physically distinct things and dna04.2 prices them
    differently:

      * the helix is terminated by a MISMATCH column ('X'). The unpaired bases stack on the
        closing pair and are worth real, context-dependent energy -- often several tenths of a
        kcal/mol of STABILISATION. dna04.2 tabulates all 96 such contexts, and rewrote every one
        of them relative to dna04.1, which is exactly why sourcing this table rather than
        transcribing a 1998 paper is the point of this module.
      * the helix ends at the end of the molecule, or abuts a gap ('I'/'D') or a soft clip
        ('S'). There is no defined stacking partner column: a bulge's flanking geometry is an
        interior/bulge loop we do not model, and a soft clip is unpaired flanking sequence. Here
        the correct minimal charge is the terminal_penalty on the closing pair itself -- the
        A.T / G.T fraying cost, zero for G.C.

    Only the first branch consults terminal_mismatch; conflating the two would hand a bulge the
    stabilisation of a stacked mismatch and make bulged duplexes look better than clean ones.
    """
    if five_prime:
        i = start - 1                       # column immediately 5' of the helix on the top strand
        close_top = probe_aln[start]
    else:
        i = stop                            # column immediately 3' of the helix on the top strand
        close_top = probe_aln[stop - 1]
    close_bot = _COMPLEMENT.get(close_top, "N")

    if 0 <= i < len(ops) and ops[i] == "X":
        mm_top = probe_aln[i]
        mm_bot = _bottom(target_aln[i])
        if five_prime:
            # outer column is 5' of the closing pair on the top strand -> native key layout
            key = mm_top + close_top + close_bot + mm_bot
        else:
            # outer column is 3' on the top strand; rotate the two columns 180 degrees, which
            # maps k0k1k2k3 -> k2k3k0k1 and puts the mismatch back on the 5' side
            key = mm_bot + close_bot + close_top + mm_top
        if key in _TERM_MM_PROHIBITED:
            # dna04.2 stores +10 kcal/mol as a PROHIBITION sentinel, not an energy, for the 20
            # keys whose outer column is itself Watson-Crick -- those are a stack, not a
            # terminal mismatch. Under SAMENESS an 'X' column can never produce one (differing
            # characters cannot complement), so this is a guard against a malformed ops string
            # rather than a real branch; charging +10 would be a silent 10 kcal/mol lie.
            pass
        else:
            hs = _TERM_MM_HS.get(key)
            if hs is not None:
                return hs[0], hs[1], True
        # An 'N' or other non-ACGT base leaves no defined context; fall through to the penalty
        # rather than invent an energy.

    hs = _TERM_PEN_HS.get(close_top + close_bot)
    if hs is None:
        return 0.0, 0.0, False
    return hs[0], hs[1], False


def _run_energy(probe_aln: str, target_aln: str, ops: str, start: int, stop: int,
                t_kelvin: float, ln_na: float) -> _Helix:
    """Score one helix.

    THE SAMENESS POINT. Within an '=' run the probe strand ALONE is the top strand, read 5'->3',
    and the parameter table supplies the partner: the base pair at column i is
    (probe[i], complement(probe[i])). We do NOT read the two aligned strings as the two strands
    of the helix -- doing so would treat a run of 'G' columns as a G/G stack and invert the whole
    energy model. target_aln carries no extra information inside a run, which is why it is only
    consulted at the ENDS, where an 'X' column's two bases genuinely differ.
    """
    seq = probe_aln[start:stop]
    n_bp = len(seq)

    dH = dS = 0.0
    gc_stacks = at_stacks = 0

    # A 1-bp "helix" has no adjacent pair and so no stacking term. It still gets both end terms
    # and the join penalty below, which for a lone A.T is net positive -- correctly, an isolated
    # base pair is not a stabilising element. No special-casing needed beyond the empty loop.
    for k in range(n_bp - 1):
        step = seq[k:k + 2]
        hs = _STACK_HS.get(step + _revcomp(step)) if "N" not in step else None
        if hs is None:
            # Non-ACGT character. Skip the step rather than guessing an energy; a silently
            # invented stack is worse than a slightly under-counted helix.
            continue
        dH += hs[0]
        dS += hs[1]
        # Classify the step by its own composition so a model can see the asymmetry directly
        # without having to reconstruct it from dG.
        if step[0] in "GC" and step[1] in "GC":
            gc_stacks += 1
        elif step[0] in "AT" and step[1] in "AT":
            at_stacks += 1

    # Helix initiation: one join penalty for the association event, plus one end term per end.
    dH += _JOIN_HS[0]
    dS += _JOIN_HS[1]

    term_mm_dG = 0.0
    n_term_mm = 0
    for five_prime in (True, False):
        e_dH, e_dS, is_mm = _terminal_dG(probe_aln, target_aln, start, stop, ops, five_prime)
        dH += e_dH
        dS += e_dS
        if is_mm:
            term_mm_dG += e_dH - t_kelvin * e_dS
            n_term_mm += 1

    # Salt. ln[Na+] is negative below 1 M, so this pushes dS down and (since dG subtracts T*dS)
    # makes the helix LESS stable -- the expected direction for sub-molar sodium.
    dS += SALT_COEF * max(n_bp - 1, 0) * ln_na

    return _Helix(dG=dH - t_kelvin * dS, n_bp=n_bp, gc_stacks=gc_stacks, at_stacks=at_stacks,
                  term_mm_dG=term_mm_dG, n_term_mm=n_term_mm)


def _helices(duplex, celsius: float, sodium: float) -> list[_Helix]:
    t_kelvin = celsius + 273.15
    ln_na = math.log(sodium)
    p, t, o = duplex.probe_aln, duplex.target_aln, duplex.ops
    return [_run_energy(p, t, o, a, b, t_kelvin, ln_na) for a, b in _runs(o)]


def duplex_dG(duplex, celsius: float = DEFAULT_CELSIUS, sodium: float = DEFAULT_SODIUM) -> float:
    """Summed nearest-neighbour dG (kcal/mol) over all Watson-Crick helices in the alignment.

    THIS IS AN UPPER BOUND ON STABILITY, NOT A dG_assoc. Only '=' columns and the columns
    immediately flanking them contribute. A mismatch column contributes only through the
    terminal_mismatch stacking of the helices it caps -- it is never charged an interior-loop
    penalty, and bulge ('I'/'D') and soft-clip ('S') columns are charged nothing at all. There is
    no coaxial stacking between adjacent helices, no dangling ends, no multiloop term, and no
    ensemble sum. A duplex broken into three short helices is scored as three independent
    helices that happen to be nearby, which physically it is not.

    The consequence is systematic and one-directional: for any imperfect duplex the returned
    value is MORE NEGATIVE than the truth, because every omitted term is destabilising. That is
    acceptable, and arguably desirable, for an ML feature -- a learned model can fit the missing
    loop penalties from the lesion-count features it already has. It must never be reported or
    thresholded as a physical association free energy.

    Returns 0.0 for an alignment with no '=' columns at all.
    """
    return sum(h.dG for h in _helices(duplex, celsius, sodium))


def duplex_features(duplex, celsius: float = DEFAULT_CELSIUS,
                    sodium: float = DEFAULT_SODIUM) -> dict:
    """A small thermodynamic feature block for one duplex.

    `max_helix_dG` is deliberately separate from `nn_dG`: nucleation is a single-helix event, so
    the best contiguous run predicts whether a duplex forms at all better than the sum does. Two
    duplexes with identical `nn_dG` -- one 20 bp clean, one four scattered 5 bp stubs -- behave
    completely differently, and only the max separates them.

    `term_mismatch_dG` is exposed separately from `nn_dG` (of which it is a component) because a
    model should be able to learn its weight independently: it is the only term here that is
    sensitive to WHICH base is mismatched rather than merely how many are, and it is the term
    that distinguishes a mismatch at a helix end from a mismatch in open flank.
    """
    helices = _helices(duplex, celsius, sodium)

    # Probe length excludes 'D' columns, where probe_aln holds a GAP: normalising by the number
    # of real probe bases keeps the per-base feature comparable across probe lengths.
    probe_len = sum(1 for c in duplex.probe_aln if c != GAP)

    nn_dG = sum(h.dG for h in helices)
    return {
        "nn_dG": nn_dG,
        "nn_dG_per_base": nn_dG / probe_len if probe_len else 0.0,
        "max_helix_dG": min((h.dG for h in helices), default=0.0),  # min = most negative
        "n_helices": len(helices),
        "longest_helix_bp": max((h.n_bp for h in helices), default=0),
        "gc_stacks": sum(h.gc_stacks for h in helices),
        "at_stacks": sum(h.at_stacks for h in helices),
        "term_mismatch_dG": sum(h.term_mm_dG for h in helices),
        "n_term_mismatches": sum(h.n_term_mm for h in helices),
    }


class _DuplexLite(NamedTuple):
    """Structural stand-in for duplex.Duplex, so this module imports nothing from the study.

    Deliberate: `thermo.py` lives in the zoo and must not depend on a dated study directory. It
    is duck-typed on the three fields, so a real `Duplex` works unchanged.
    """
    probe_aln: str
    target_aln: str
    ops: str


def stacking_profile(duplex, width: int, celsius: float = DEFAULT_CELSIUS,
                     sodium: float = DEFAULT_SODIUM):
    """Per-COLUMN stacking free energy: the vector whose sum `nn_dG` throws away.

    WHY THIS EXISTS. `duplex_features` reports the SUM of the stacking energies. The controlled-
    degradation sweep (`studies/20260720_e_controlled_degradation`) then measured that WHERE the
    energy is lost dominates HOW MUCH: at four substitutions, pDup is 0.875 with the lesions
    clustered at the 5' end and 0.046 with them spread evenly -- a 19x difference at identical
    lesion count, and (nearly) identical total stacking loss. A scalar sum is blind to the entire
    effect.

    Three independent literatures reach the same conclusion. The microarray PDNN model
    (Zhang et al.) multiplies each nearest-neighbour term by a positional weight omega(k) rather
    than summing them flat; the siRNA and antisense-oligo efficacy literatures both find
    position-dependent contribution profiles. This function supplies the raw material for a model
    to learn omega(k) itself instead of us fitting it.

    RETURNS a float vector of length `width`, entry k holding the dG of the stacking step BETWEEN
    aligned columns k and k+1, or 0.0 where that step does not exist -- because the run ends
    there, because a lesion interrupts it, or because the duplex is shorter than `width`. Zero is
    the honest encoding for "no stack here": it is exactly the contribution such a column makes
    to the sum.

    NOTE the entries are step energies only -- no join penalty, no terminal-mismatch term, no
    salt correction. Those are per-HELIX, not per-column, and cannot be attributed to a single
    position without inventing an attribution. They remain available as scalars from
    `duplex_features`, so a model gets both views and neither is double-counted.
    """
    import numpy as np

    t_kelvin = celsius + 273.15
    out = np.zeros(width, dtype=np.float32)
    p_aln, ops = duplex.probe_aln, duplex.ops

    for start, stop in _runs(ops):
        seq = p_aln[start:stop]
        for k in range(len(seq) - 1):
            col = start + k
            if col >= width:
                break
            step = seq[k:k + 2]
            hs = _STACK_HS.get(step + _revcomp(step)) if "N" not in step else None
            if hs is None:
                continue                       # non-ACGT: skip rather than invent an energy
            out[col] = hs[0] - t_kelvin * hs[1]
    return out


FEATURE_COLUMNS = [
    "nn_dG", "nn_dG_per_base", "max_helix_dG",
    "n_helices", "longest_helix_bp", "gc_stacks", "at_stacks",
    "term_mismatch_dG", "n_term_mismatches",
]


# Off-target alignments repeat heavily -- a handful of repeat families account for a large share
# of rows, and identical (probe, target, ops) triples recur thousands of times. Cache on the
# strings themselves rather than on a row index so the hit rate survives shuffling and chunking.
# 2**20 entries is a few hundred MB worst case; drop it if memory-bound.
@lru_cache(maxsize=1 << 20)
def _features_cached(probe_aln: str, target_aln: str, ops: str,
                     celsius: float, sodium: float) -> tuple:
    d = _DuplexLite(probe_aln, target_aln, ops)
    f = duplex_features(d, celsius=celsius, sodium=sodium)
    return tuple(f[k] for k in FEATURE_COLUMNS)


def dG_for_rows(probe_aln_list, target_aln_list, ops_list,
                celsius: float = DEFAULT_CELSIUS, sodium: float = DEFAULT_SODIUM,
                index=None) -> pd.DataFrame:
    """Apply `duplex_features` over parallel sequences of alignment columns -> DataFrame.

    Pure Python per unique row; the win at scale comes from the lru_cache, not from vectorising
    the inner loop, because the inner loop is a string walk that numpy cannot help with. Pass
    `index` to keep alignment with the source frame (e.g. `df.index`).
    """
    rows = [
        _features_cached(p, t, o, celsius, sodium)
        for p, t, o in zip(probe_aln_list, target_aln_list, ops_list)
    ]
    return pd.DataFrame(rows, columns=FEATURE_COLUMNS, index=index)


# ---------------------------------------------------------------------------------------------
# Self-check
# ---------------------------------------------------------------------------------------------
def _perfect(seq: str) -> _DuplexLite:
    """A perfect duplex in SAMENESS form: both strands identical, all columns '='."""
    return _DuplexLite(seq, seq, "=" * len(seq))


def _with_mismatch(seq: str, pos: int) -> _DuplexLite:
    """Same duplex with column `pos` turned into an 'X' by substituting the TARGET base."""
    sub = "A" if seq[pos] != "A" else "T"
    return _DuplexLite(seq, seq[:pos] + sub + seq[pos + 1:],
                       "=" * pos + "X" + "=" * (len(seq) - pos - 1))


if __name__ == "__main__":
    # 1. The tripwire. Study 20260718 finding 5, as a live assertion.
    bad = santalucia_crosscheck()
    assert not bad, "dna04.2 stacks have MOVED off SantaLucia 1998:\n  " + "\n  ".join(bad)
    print(f"SantaLucia-1998 cross-check: {len(_SANTALUCIA_1998)}/{len(_SANTALUCIA_1998)} unique "
          f"WC steps agree with dna04.2 to 0.01 kcal/mol (dG37 and dH)")
    print(f"loaded: {PARAM_PATH.name}  "
          f"stack={len(STACK_dG)}  terminal_mismatch={len(TERM_MM_dG)}  "
          f"terminal_penalty={len(TERM_PEN_dG)}  join_penalty={JOIN_dG}")

    # 2. Layout sanity: the stack table must be invariant under the 180-degree rotation that
    #    `_terminal_dG` relies on to score the 3' end of a helix. If this fails, the four-char
    #    key layout has been misread and every number below is wrong.
    assert all(abs(STACK_dG[k] - STACK_dG[k[2:] + k[:2]]) < 1e-9 for k in STACK_dG)

    gc = _perfect("GC" * 15)
    at = _perfect("AT" * 15)
    dg_gc, dg_at = duplex_dG(gc), duplex_dG(at)

    # 3. The headline: bowtie2 scores these two identically. Thermodynamics does not. And the
    #    SAMENESS handling is load-bearing here -- read the two strings as complementary strands
    #    and all-GC stops being the stable one.
    assert dg_gc < -20.0, dg_gc
    assert dg_at > 0.0, dg_at        # 30 bp of AT is not bound at 74.5 C
    assert dg_gc < dg_at - 15.0, (dg_gc, dg_at)

    # 4. A central mismatch splits one 30 bp helix into two, paying a second join penalty and
    #    losing a stack -- less stable, even though no interior-loop penalty is charged.
    mid = _with_mismatch("GC" * 15, 15)
    dg_mid = duplex_dG(mid)
    assert dg_mid > dg_gc, (dg_mid, dg_gc)

    # 5. THE POINT OF THE REWRITE. Same duplex, same mismatch count, mismatch moved from the
    #    middle to column 0. At the end it caps a single long helix with a tabulated terminal
    #    mismatch instead of severing the helix in two, so it must be markedly more stable --
    #    and it must actually register terminal-mismatch energy.
    end = _with_mismatch("GC" * 15, 0)
    dg_end = duplex_dG(end)
    f_end, f_mid = duplex_features(end), duplex_features(mid)
    assert dg_end < dg_mid - 2.0, (dg_end, dg_mid)
    assert f_end["n_term_mismatches"] == 1
    assert f_end["term_mismatch_dG"] != 0.0
    # A mismatch that abuts nothing but molecule-end (perfect duplex) gets terminal_penalty, not
    # terminal_mismatch -- the branch in `_terminal_dG`.
    assert duplex_features(gc)["n_term_mismatches"] == 0

    # 6. A 1 bp helix must not crash and must not be stabilising.
    lone = _DuplexLite("AATAA", "CGTCG", "XX=XX")
    assert duplex_dG(lone) > 0, duplex_dG(lone)

    labels = ["30bp all-GC", "30bp all-AT", "GC, mismatch at col 15 (middle)",
              "GC, mismatch at col 0 (end)", "1bp A.T island, mismatch-flanked"]
    cases = [gc, at, mid, end, lone]
    rows = dG_for_rows([c.probe_aln for c in cases], [c.target_aln for c in cases],
                       [c.ops for c in cases], index=labels)
    pd.set_option("display.width", 200)
    print(f"\nconditions: {DEFAULT_CELSIUS} C, {DEFAULT_SODIUM} M Na+ | params: dna04.2\n")
    print(rows.round(3).to_string())
    print("\nGC - AT dG gap:        {:+.2f} kcal/mol  <- invisible to bowtie2's --ma score".format(
        dg_gc - dg_at))
    print("end - middle mismatch: {:+.2f} kcal/mol  <- invisible to a mismatch COUNT".format(
        dg_end - dg_mid))
    print("\nall self-checks passed")
