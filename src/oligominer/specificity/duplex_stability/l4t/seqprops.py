#!/usr/bin/env python3
"""Sequence properties of a PROBE: length, GC%, and nearest-neighbour melting temperature.

    python seqprops.py        # self-check against an independent published value

PROBE-LEVEL, NOT DUPLEX-LEVEL. Everything here is a property of the probe oligo BEFORE it is
aligned to anything: what the molecule is, not what it found in a genome. That distinction is the
whole reason the corpus figures have two halves -- these describe the SUBSTRATE, and pDup describes
the OUTCOME of putting that substrate against a genome.

Tm IS COMPUTED AGAINST THE PROBE'S PERFECT COMPLEMENT, at 1 uM strand and 0.39 M Na+ (2x SSC, the
standard DNA-FISH condition, 2x SSC). It is therefore a property of the probe, NOT a prediction about any
genomic target -- a probe with Tm 80 C may still have pDup ~0 against every off-target site it
actually hits. Quoting it as if it predicted binding is the misreading to avoid.

PARAMETERS COME FROM THE SAME dna04.2 TABLES AS THE LABELS (`thermo.py`, vendored in this
directory). That is deliberate: a Tm computed from a hand-typed SantaLucia table and a pDup label
computed by NUPACK 4.0.2.0 would differ for reasons that have nothing to do with the biology.
"""
import math
import sys
from pathlib import Path

import numpy as np

from .thermo import _STACK_HS, TERM_PEN_dG, TERM_PEN_dH, JOIN_dG, JOIN_dH, _dS, _revcomp, T37

R = 1.9872                     # cal / (mol K)
CT = 1e-6                      # total strand concentration, M -- the standard 1 uM working concentration
NA = 0.390                     # M sodium, 2x SSC


def gc_pct(s):
    return 100.0 * sum(c in "GC" for c in s) / len(s) if s else float("nan")


def tm(seq, ct=CT, sodium=NA):
    """Nearest-neighbour Tm (C) of `seq` against its PERFECT COMPLEMENT.

    Tm = dH / (dS + R ln(CT/4))  -- the non-self-complementary form. CT/4 (not CT) because the
    duplex is a heterodimer: at Tm, [A] = [B] = CT/4 with half the strands duplexed.

    Salt correction is Owczarzy et al. 2004 eq. 22, which is a function of fractional GC as well as
    ln[Na+] -- a GC-independent correction misprices AT-rich oligos, and this corpus deliberately
    spans a wide GC range.
    """
    s = seq.upper()
    if len(s) < 2 or any(c not in "ACGT" for c in s):
        return float("nan")

    dH = JOIN_dH
    dG = JOIN_dG
    for i in range(len(s) - 1):
        step = s[i:i + 2]
        key = step + _revcomp(step)              # thermo.py's 4-char stack key layout
        h, _ = _STACK_HS[key]
        dH += h
        from thermo import STACK_dG
        dG += STACK_dG[key]
    for end in (s[0], s[-1]):                    # terminal penalty at BOTH helix ends
        k = end + _revcomp(end)
        dH += TERM_PEN_dH.get(k, 0.0)
        dG += TERM_PEN_dG.get(k, 0.0)

    dS = _dS(dH, dG)                             # kcal/(mol K), derived exactly as thermo.py does
    denom = dS * 1000.0 + R * math.log(ct / 4.0)
    if denom >= 0:
        return float("nan")
    tm_k = (dH * 1000.0) / denom

    # Owczarzy 2004 salt correction, applied from the 1 M reference the tables are quoted at.
    f_gc = gc_pct(s) / 100.0
    ln_na = math.log(max(sodium, 1e-9))
    inv = 1.0 / tm_k + (4.29 * f_gc - 3.95) * 1e-5 * ln_na + 9.40e-6 * ln_na ** 2
    return float(1.0 / inv - 273.15)


def describe(seqs):
    """(length, gc, tm) arrays for a list of probe sequences."""
    seqs = [str(s).upper() for s in seqs]
    return (np.array([len(s) for s in seqs], dtype=float),
            np.array([gc_pct(s) for s in seqs], dtype=float),
            np.array([tm(s) for s in seqs], dtype=float))


if __name__ == "__main__":
    # SELF-CHECK against an INDEPENDENT value: `20260726_f_om2_ssot_datasets/DATASETS.md` reports
    # chrX-30spot Tm median 79.6 C and Genomic-grid 78.87 C, computed by a different script
    # (its `seqprops.py`, itself validated against Biopython to within 0.8 C). Agreement to a
    # degree or two means this implementation is right; a large gap means it is not, and the
    # figures built on it would be wrong in a way no plot would reveal.
    import pandas as pd
    P = "/net/beliveau/vol1/project/conor/om2_ssot/datasets/probes"
    targets = {"chrx_30spot_probes.parquet": 79.6, "grid_probes.parquet": 78.87}
    ok = True
    for f, want in targets.items():
        df = pd.read_parquet(f"{P}/{f}")
        col = next(c for c in ("probe_seq", "sequence", "seq", "homology") if c in df.columns)
        # RANDOM SAMPLE, NEVER head(). `grid_probes.parquet` is SORTED BY LENGTH, so head(3000) is
        # 2,633 10-mers plus 367 20-mers -- the shortest probes in the set -- and returns a median
        # Tm of 36.9 C against a true 78.5. That looked exactly like a broken Tm implementation and
        # was a broken TEST. Any check on this file must sample randomly.
        samp = df[col] if len(df) <= 4000 else df[col].sample(4000, random_state=0)
        _, _, t = describe(samp)
        got = float(np.nanmedian(t))
        d = abs(got - want)
        ok &= d < 2.5
        print(f"  {f:34s} median Tm {got:6.2f} C  (independent: {want})  "
              f"delta {d:.2f}  {'OK' if d < 2.5 else 'MISMATCH'}")
    print("\nSELF-CHECK PASS" if ok else "\nSELF-CHECK FAILED — do not build figures on this")
    sys.exit(0 if ok else 1)
