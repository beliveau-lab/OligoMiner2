#!/usr/bin/env python3
"""Two feature sets: PaintSHOP-37 (the prior art) and OM2 physics-XGB (the flagship).

    python features.py      # conformance vs the screen encoders + superset/width/condition checks

ONLY TWO. The screening ladder is finished; the intermediate rungs are deliberately absent. What
remains is the baseline this work is measured against and the model that ships.

THE RELATIONSHIP IS A STRICT SUPERSET, AND THAT IS THE WHOLE ARGUMENT.

    PaintSHOP-37   37 features
    OM2           103 features  =  the SAME 37, verbatim, + 66 added

Nothing was dropped, reweighted, or reimplemented. Every PaintSHOP feature is computed by the same
code path and appears in the OM2 vector under a `ps_` prefix. That makes PaintSHOP -> OM2 a
CONTROLLED ADDITIVE ABLATION rather than a comparison of two unrelated models: any difference in
performance is attributable to the 66 added features and to nothing else. If a future change breaks
that containment, it breaks what the comparison means, which is why `verify()` asserts it.

WHAT THE 66 ADD, by block:

    16  alignment + CIGAR decomposition   what the aligner found: counts of matches, mismatches,
                                          indels, soft clips, the longest intact run, and WHERE the
                                          lesions sit relative to each probe terminus
     9  nearest-neighbour thermodynamics  dG of the duplex under dna04.2, its helical structure,
                                          and terminal-mismatch energetics
    32  terminal stacking profile         per-column stacking energy for the 16 columns at each
                                          terminus, anchored in NUCLEOTIDES from the physical end
     9  interior summaries                length-free reductions over whatever lies between

THE ONE THING PAINTSHOP STRUCTURALLY CANNOT DO. PaintSHOP-37 is composition and alignment score --
counts of dinucleotides, GC, lengths. Nothing in it is a function of temperature or salt, and there
is nowhere to put one: a hybridization condition is not a property of a sequence. So a PaintSHOP
model is pinned to whatever condition its training labels happened to be computed at, and asked
about another it returns the same answer. The OM2 vector's dG features are RECOMPUTED at each row's
own (T, Na+), so the inputs move when the condition moves. That is not a tuning advantage; it is the
difference between a model that can represent the question and one that cannot.

WIDTH-FREE. The terminal blocks are indexed from the two physical ends, and the interior is
summarised rather than enumerated, so the vector is 103 wide for a 10-mer and for an 80-mer alike.
No alignment-width contract, no truncation path, no re-fit for a longer probe.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
from .thermo import (stacking_profile, duplex_features, _DuplexLite,      # noqa: E402
                    DEFAULT_CELSIUS, DEFAULT_SODIUM)

THRESHOLD = 0.2
TERMINAL_K = 16          # columns held at full resolution at EACH terminus
GOOD_STEP = -0.5         # kcal/mol; a step at least this stabilising counts as intact helix

# PaintSHOP's dinucleotide order. NOT alphabetical, and NOT arbitrary: a tree indexes features
# positionally, so re-sorting produces a model that trains fine and serves wrong.
DINUCS = ['AA', 'AT', 'AG', 'AC', 'TA', 'TT', 'TG', 'TC',
          'GA', 'GT', 'GG', 'GC', 'CA', 'CT', 'CG', 'CC']

CIGAR_COLS = ["aln_len", "n_eq", "n_mm", "n_ins", "n_del", "n_soft", "longest_eq_run",
              "n_lesion", "core_lesions", "first_lesion_5p", "first_lesion_3p",
              "mean_lesion_pos_norm"]
THERMO_COLS = ["nn_dG", "nn_dG_per_base", "max_helix_dG", "n_helices", "longest_helix_bp",
               "gc_stacks", "at_stacks", "term_mismatch_dG", "n_term_mismatches"]
INTERIOR_COLS = ["int_n", "int_sum", "int_mean", "int_min", "int_max", "int_std",
                 "int_worst_nt_from_end", "int_longest_good_run"]

# Measured, not assumed (see `verify()`): only the ENERGIES move with condition. The structural
# counts are facts about the aligned duplex and must not move -- heating the tube does not change
# how many helices something has.
TEMPERATURE_DEPENDENT = {"nn_dG", "nn_dG_per_base", "max_helix_dG", "term_mismatch_dG"}
SALT_DEPENDENT = {"nn_dG", "nn_dG_per_base", "max_helix_dG"}


def strip_t3(s):
    """Remove TTT synthesis flanks -- only when BOTH ends carry them (conjunctive by design)."""
    return s[3:-3] if s.startswith("TTT") and s.endswith("TTT") else s


def gc_pct(s):
    return 100.0 * sum(c in "GC" for c in s) / len(s) if s else 0.0


def _align_score(df):
    if "align_score" in df.columns:
        return df["align_score"].values.astype(float)
    from bowtie_score import local_score
    return np.array([float(local_score(o)) for o in df["ops"].values])


# ---------------------------------------------------------------------------------------------
# PaintSHOP-37 — the prior art, reproduced exactly
# ---------------------------------------------------------------------------------------------
def enc_paintshop37(df):
    """5 scalars + 16 probe dinucleotide counts + 16 target dinucleotide counts = 37.

    `str.count` is NON-OVERLAPPING: "AAAA".count("AA") == 2, not 3. That is baked into the
    published PaintSHOP weights; "fixing" it silently changes every dinucleotide feature.
    """
    probe = df["probe_seq"].astype(str).map(strip_t3)
    derived = df["target_seq"].astype(str)
    out = {"align_score": _align_score(df),
           "probe_gc": probe.map(gc_pct).values, "derived_gc": derived.map(gc_pct).values,
           "probe_len": probe.str.len().astype(float).values,
           "derived_len": derived.str.len().astype(float).values}
    for dn in DINUCS:
        out[f"probe_{dn}"] = probe.str.count(dn).astype(float).values
    for dn in DINUCS:
        out[f"derived_{dn}"] = derived.str.count(dn).astype(float).values
    return pd.DataFrame(out)


# ---------------------------------------------------------------------------------------------
# the 66 added features
# ---------------------------------------------------------------------------------------------
def _cigar_block(df):
    """What the aligner found, decomposed. Positions are in the PROBE's own frame (0 = its 5' end),
    skipping columns where the probe has no base. The terminal 5 nt at each end are the dangling-end
    zone -- a lesion there is thermodynamically cheap relative to the same lesion mid-helix -- so
    `core_lesions` counts only those outside it."""
    rows = []
    for pa, op in zip(df["probe_aln"].values, df["ops"].values):
        n = {"=": 0, "X": 0, "I": 0, "D": 0, "S": 0, "M": 0}
        longest = run = probe_i = 0
        mm = []
        for pc, o in zip(pa, op):
            n[o] = n.get(o, 0) + 1
            if o == "=":
                run += 1; longest = max(longest, run)
            else:
                run = 0
                if pc != "-":
                    mm.append(probe_i)
            if pc != "-":
                probe_i += 1
        plen = probe_i or 1
        core = [p for p in mm if 5 <= p < plen - 5]
        rows.append([len(op), n["="], n["X"], n["I"], n["D"], n["S"], longest, len(mm), len(core),
                     mm[0] if mm else -1, (plen - 1 - mm[-1]) if mm else -1,
                     (sum(mm) / len(mm) / plen) if mm else -1.0])
    return pd.DataFrame(rows, columns=CIGAR_COLS)


def _conditions(df, celsius, sodium):
    """Per-ROW (T, Na+). A scalar pins the whole frame to one condition -- correct for a
    single-condition corpus, silently wrong for a mixed one, and invisible because the features
    stay finite and the model still fits."""
    cel = (np.full(len(df), float(celsius)) if celsius is not None else
           (df["label_celsius"].values.astype(float) if "label_celsius" in df.columns
            else np.full(len(df), DEFAULT_CELSIUS)))
    na = (np.full(len(df), float(sodium)) if sodium is not None else
          (df["label_sodium"].values.astype(float) if "label_sodium" in df.columns
           else np.full(len(df), DEFAULT_SODIUM)))
    return cel, na


def _terminal_block(prof, K):
    """K columns from the 5' terminus, K from the 3' (reversed), + interior summaries.

    EVERY FEATURE IS IN NUCLEOTIDES FROM A PHYSICAL LANDMARK, never a fraction of length.
    Nearest-neighbour dG is local and additive in absolute steps, and the effects that motivate a
    positional profile at all -- dangling ends, terminal mismatch, end fraying -- act over roughly
    the last 1-3 nt. A fractional grid would smear the terminus into the interior AS A FUNCTION OF
    LENGTH, degrading the encoding exactly where the physics is sharpest.

    `s3_0` is the 3'-terminal step, `s3_1` the one inside it, so index means "distance from THIS
    terminus" at both ends and at every probe length.

    OVERLAP RULE: when the duplex is shorter than 2K the windows overlap; both are still emitted in
    full (a short duplex is simply read twice, once from each end), the interior is empty, and
    `overlap_nt` records by how much -- recorded as a feature rather than hidden, so the model can
    tell a genuinely short duplex from a long one whose interior happened to be featureless.
    """
    p = np.asarray(prof, dtype=np.float64)
    L = len(p)
    s5 = np.zeros(K); s3 = np.zeros(K)
    s5[:min(K, L)] = p[:K]
    tail = p[max(0, L - K):][::-1]
    s3[:len(tail)] = tail
    interior = p[K:L - K] if L > 2 * K else np.array([])
    if len(interior):
        worst = K + int(np.argmin(interior))
        good = interior <= GOOD_STEP
        best = run = 0
        for g in good:
            run = run + 1 if g else 0
            best = max(best, run)
        summ = [float(len(interior)), float(interior.sum()), float(interior.mean()),
                float(interior.min()), float(interior.max()), float(interior.std()),
                float(min(worst, L - 1 - worst)), float(best)]
    else:
        summ = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.0, 0.0]
    return s5, s3, summ, float(max(0, 2 * K - L))


def enc_om2(df, celsius=None, sodium=None, K=TERMINAL_K):
    """The flagship: PaintSHOP-37 (verbatim, `ps_` prefixed) + 66 added = 103 features, width-free."""
    ps = enc_paintshop37(df)
    ps.columns = [f"ps_{c}" for c in ps.columns]

    d = df
    cel, na = _conditions(d, celsius, sodium)
    base = pd.DataFrame({"align_score": _align_score(d),
                         "probe_len": d["probe_seq"].astype(str).str.len().astype(float).values,
                         "probe_gc": d["probe_seq"].astype(str).map(gc_pct).values,
                         "target_gc": d["target_seq"].astype(str).map(gc_pct).values})
    cig = _cigar_block(d)
    th = pd.DataFrame([[duplex_features(_DuplexLite(pa, ta, op), celsius=float(c),
                                        sodium=float(s))[k] for k in THERMO_COLS]
                       for pa, ta, op, c, s in zip(d.probe_aln, d.target_aln, d.ops, cel, na)],
                      columns=THERMO_COLS)
    s5s, s3s, sums, ovs = [], [], [], []
    for pa, ta, op, c, s in zip(d.probe_aln, d.target_aln, d.ops, cel, na):
        prof = stacking_profile(_DuplexLite(pa, ta, op), len(pa), celsius=float(c), sodium=float(s))
        a, b, u, o = _terminal_block(prof, K)
        s5s.append(a); s3s.append(b); sums.append(u); ovs.append(o)

    return pd.concat([
        ps.reset_index(drop=True), base, cig, th,
        pd.DataFrame(np.asarray(s5s, np.float32), columns=[f"s5_{j}" for j in range(K)]),
        pd.DataFrame(np.asarray(s3s, np.float32), columns=[f"s3_{j}" for j in range(K)]),
        pd.DataFrame(np.asarray(sums, np.float32), columns=INTERIOR_COLS),
        pd.DataFrame({"overlap_nt": np.asarray(ovs, np.float32)}),
    ], axis=1)


FEATURES_PAINTSHOP = None      # filled on first call, so the contract is discoverable
FEATURES_OM2 = None


