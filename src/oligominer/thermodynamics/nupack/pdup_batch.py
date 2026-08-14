"""
# Batched pDup

Computes the same pDup values as ``calc_pdup`` for many duplexes at once.

``calc_pdup`` builds a Tube and calls ``tube_analysis`` per site. One tube of five
small complexes has nothing to parallelize, so that path is pinned at its
single-site cost no matter how many cores are available. Submitting the partition
functions as one job list instead releases the GIL and uses NUPACK's own executor,
which is what allows more than one core to be used at all.

Three levers, all exact:

1. skip the Tube and Result object graph and run the concentration solve directly
2. batch every partition function into one submit, which saturates around 64
   complexes so the work streams rather than needing to be materialized
3. when one probe is screened against many targets, compute the probe-only terms
   Q(A) and Q(A.A) once instead of once per target

Two properties of the low-level path produce wrong numbers rather than errors, and
both are handled here:

- a homodimer built from two Strand objects with different names is two
  distinguishable strands, so its 2-fold rotational symmetry goes uncorrected.
  One Strand object per distinct sequence is required, not merely faster.
- ``submit(...).get()`` returns a map keyed by sequence content rather than by
  submission order, so results are looked up by key. Zipping against the job list
  misaligns as soon as two complexes in a batch share sequences, which happens
  constantly in an off-target sweep.

The low-level NUPACK API is not part of its documented public surface, so
``low_level_available()`` reports whether it can be used and the batch entry points
fall back to ``tube_analysis`` when it cannot.
"""

import numpy as np

from .config import DEFAULT_NUPACK_MODEL

# concentrations of the probe and target strands, matching calc_pdup
CONC_A = 1e-6
CONC_B = 1e-12

# complexes whose partition functions determine pDup, in solve order:
# A.B, A, B, A.A, B.B
_COMPLEX_INDICES = [[0, 1], [0], [1], [0, 0], [1, 1]]

# batch throughput saturates near this many complexes per submit
BATCH_COMPLEXES = 320

_STRAND_CACHE = {}


def low_level_available():
    """
    Report whether NUPACK's low-level batch API can be used.

    Returns:
        available (bool): True when the batch submit path is importable.
    """
    try:
        from nupack import config                                      # noqa: F401
        from nupack.concentration import solve_complex_concentrations  # noqa: F401
        from nupack.core import SequenceList                           # noqa: F401
        from nupack.thermo import ComputeOptions, Job, PFJob, submit   # noqa: F401
    except ImportError:
        return False

    # success
    return True


def _strand(seq):
    """
    Return one Strand object per distinct sequence.

    Reusing a single object per sequence is required for correctness: a homodimer
    built from two differently named Strand objects is treated as two
    distinguishable strands and its rotational symmetry goes uncorrected.

    Args:
        seq (str): the sequence.

    Returns:
        strand (nupack.Strand): the cached strand object.
    """
    import nupack

    strand = _STRAND_CACHE.get(seq)
    if strand is None:
        strand = _STRAND_CACHE[seq] = nupack.Strand(seq, name=f's{len(_STRAND_CACHE)}')

    # success
    return strand


def _complex(*seqs):
    """
    Build a complex from sequences, reusing cached strand objects.

    Args:
        *seqs (str): the sequences in the complex.

    Returns:
        complex (nupack.Complex): the complex.
    """
    import nupack

    # success
    return nupack.Complex([_strand(s) for s in seqs])


def _options():
    """
    Build the compute options carrying NUPACK's executor and cache budget.

    Rebuilt per batch because the options hold the executor handle, and a stale
    handle serializes the batch.

    Returns:
        options (nupack.thermo.ComputeOptions): the options.
    """
    from nupack import config
    from nupack.thermo import ComputeOptions

    # success
    return ComputeOptions(max_bytes=int(float(config.cache) * 1e9))


def _logq(result, complex_):
    """
    Return one complex's symmetry-corrected log partition function.

    Args:
        result: the map returned by submit(...).get().
        complex_ (nupack.Complex): the complex to look up.

    Returns:
        logq (float): the corrected log partition function.
    """
    from nupack.core import SequenceList

    # success
    return (result[SequenceList(complex_)].pfunc.value().logq
            - np.log(int(complex_.symmetry())))


def _solve(logq, model):
    """
    Solve the two-strand equilibrium for pDup.

    Args:
        logq (list): log partition functions for A.B, A, B, A.A and B.B.
        model (nupack.Model): the thermodynamic model, whose temperature the
            equilibrium solve must run at.

    Returns:
        pdup (float): the fraction of the minority strand in duplex form.
    """
    from nupack.concentration import solve_complex_concentrations

    conc = solve_complex_concentrations(
        _COMPLEX_INDICES, logq, [CONC_A, CONC_B],
        kelvin=model.temperature, as_strands=True, rotational_correction=False)

    # success
    return float(conc[0] / CONC_B)


def calc_pdup_many(pairs, model=None, batch_complexes=BATCH_COMPLEXES):
    """
    Compute pDup for many probe and target pairs.

    Args:
        pairs (sequence): (probe_seq, target_seq) tuples.
        model (nupack.Model, optional): the thermodynamic model. Defaults to
            DEFAULT_NUPACK_MODEL.
        batch_complexes (int): complexes submitted per batch, which bounds peak
            memory. Throughput saturates well below the default.

    Returns:
        values (list): pDup per input pair, in input order.
    """
    from .pdup import calc_pdup

    model = model or DEFAULT_NUPACK_MODEL
    pairs = list(pairs)
    if not pairs:
        return []

    if not low_level_available():
        return [calc_pdup(a, b, conc_a=CONC_A, conc_b=CONC_B, model=model)
                for a, b in pairs]

    # five complexes per pair, so convert the complex budget into a pair budget
    per_batch = max(1, batch_complexes // 5)

    values = []
    for start in range(0, len(pairs), per_batch):
        values.extend(_pdup_batch(pairs[start:start + per_batch], model))

    # success
    return values


def _pdup_batch(pairs, model):
    """
    Compute pDup for one batch of pairs in a single submit.

    Args:
        pairs (list): (probe_seq, target_seq) tuples.
        model (nupack.Model): the thermodynamic model.

    Returns:
        values (list): pDup per pair.
    """
    from nupack.thermo import Job, PFJob, submit

    complexes, spans = [], []
    for probe, target in pairs:
        spans.append(len(complexes))
        complexes += [_complex(probe, target), _complex(probe), _complex(target),
                      _complex(probe, probe), _complex(target, target)]

    result = submit(model, [Job(c, PFJob()) for c in complexes], _options()).get()
    logq = [_logq(result, c) for c in complexes]

    # success
    return [_solve(logq[i:i + 5], model) for i in spans]


def calc_pdup_one_to_many(probe, targets, model=None,
                          batch_complexes=BATCH_COMPLEXES):
    """
    Compute pDup for one probe against many targets.

    Q(A) and Q(A.A) depend only on the probe, so they are computed once for the
    whole target list rather than once per target. A.A is a two-strand complex, so
    this removes one of the three expensive terms per site.

    Args:
        probe (str): the probe sequence.
        targets (sequence): the target sequences.
        model (nupack.Model, optional): the thermodynamic model. Defaults to
            DEFAULT_NUPACK_MODEL.
        batch_complexes (int): complexes submitted per batch.

    Returns:
        values (list): pDup per target, in input order.
    """
    model = model or DEFAULT_NUPACK_MODEL
    targets = list(targets)
    if not targets:
        return []

    if not low_level_available():
        return calc_pdup_many([(probe, t) for t in targets], model=model)

    # three complexes per target, plus two probe-only terms per batch
    per_batch = max(1, (batch_complexes - 2) // 3)

    values = []
    for start in range(0, len(targets), per_batch):
        values.extend(_pdup_one_to_many_batch(
            probe, targets[start:start + per_batch], model))

    # success
    return values


def _pdup_one_to_many_batch(probe, targets, model):
    """
    Compute pDup for one probe against one batch of targets.

    Args:
        probe (str): the probe sequence.
        targets (list): the target sequences.
        model (nupack.Model): the thermodynamic model.

    Returns:
        values (list): pDup per target.
    """
    from nupack.thermo import Job, PFJob, submit

    probe_complexes = [_complex(probe), _complex(probe, probe)]
    per_target = [(_complex(probe, t), _complex(t), _complex(t, t)) for t in targets]

    complexes = probe_complexes + [c for triple in per_target for c in triple]
    result = submit(model, [Job(c, PFJob()) for c in complexes], _options()).get()
    logq = [_logq(result, c) for c in complexes]

    logq_a, logq_aa = logq[0], logq[1]
    values = []
    for i in range(len(targets)):
        logq_ab, logq_b, logq_bb = logq[2 + 3 * i], logq[3 + 3 * i], logq[4 + 3 * i]
        values.append(_solve([logq_ab, logq_a, logq_b, logq_aa, logq_bb], model))

    # success
    return values


def add_pdup_batch(merged_df, model=None, probe_col='probe_seq',
                   target_col='derived_seq', out_col='pdup'):
    """
    Add a pDup column to a duplex table.

    Rows are grouped by probe so the probe-only partition functions are shared
    across every target that probe aligns to.

    Args:
        merged_df (pandas.DataFrame): duplex table carrying probe_col and
            target_col.
        model (nupack.Model, optional): the thermodynamic model.
        probe_col (str): column holding the probe sequence.
        target_col (str): column holding the target sequence.
        out_col (str): column to write pDup into.

    Returns:
        out (pandas.DataFrame): a copy with out_col added.
    """
    model = model or DEFAULT_NUPACK_MODEL
    out = merged_df.copy()
    values = np.zeros(len(out), dtype=np.float64)

    for probe, block in out.groupby(probe_col, sort=False):
        positions = out.index.get_indexer(block.index)
        values[positions] = calc_pdup_one_to_many(
            str(probe), [str(t) for t in block[target_col]], model=model)

    out[out_col] = values

    # success
    return out
