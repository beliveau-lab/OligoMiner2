"""
A vectorized L4t encoder, bit-identical to the vendored one.

# The problem

`features.enc_om2` (the vendored flagship implementation, and what physics-XGB IS) runs two
explicit Python loops over rows, one calling `thermo.stacking_profile` and one calling
`thermo.duplex_features`. Measured on 20,000 real rows:

    PaintSHOP-37 (enc_paintshop37)    37 feat     23.8 us/row
    OM2 L4t      (enc_om2)           103 feat    318.4 us/row

That is 13.4x slower than the incumbent's encoder, and it is 57% of the whole OM2 pipeline's
cost. Some of the gap is real -- computing nearest-neighbour free energies at each row's own
condition is more work than counting dinucleotides, and it is what buys condition-awareness --
but a per-row Python loop over a (rows x width) computation is not.

# What is vectorized, and what is not

**`stacking_profile` is vectorized across rows.** Its output depends on three things: the
dinucleotide step at each column, whether that column sits inside an uninterrupted `=` run and
is not its last column, and the temperature. All three become arrays:

- the alignment strings become a `(rows, width)` uint8 matrix;
- the step at column k is `probe[k], probe[k+1]`, so a pair code indexes a lookup table of
  (dH, dS) built once from `_STACK_HS`;
- run membership is `ops[k] == '=' and ops[k+1] == '='`, an elementwise AND on the ops matrix;
- `dG = dH - T*dS` with T broadcast per row.

**`duplex_features` is NOT vectorized.** Its inner work is a string walk that numpy cannot
help with -- the vendored module says so itself. Instead it is routed through the `lru_cache`
that `thermo` already provides and that `enc_om2` bypasses by constructing `_DuplexLite`
inline. Off-target alignments repeat heavily, so the cache is where that win lives.

# Bit-identity is the shipping gate

Float32 output means arithmetic order matters, so this cannot be "close enough". `verify()`
asserts every one of the 103 columns matches the vendored encoder exactly, on real rows, and
the fast path refuses to be used if it does not. A faster encoder that changes a feature value
is not an optimization, it is a different model.
"""

import numpy as np
import pandas as pd

from .l4t import features as _F                                                           # noqa: E402
from .l4t import thermo as _T                                                             # noqa: E402

# 128x128 lookup so a step is indexed by its two ASCII bytes directly. Entries absent from
# _STACK_HS (any step containing N, or a non-ACGT byte) stay NaN and are masked out, which
# reproduces the vendored `continue` rather than inventing an energy for them.
#
# THE TABLE IS BUILT BY THE SAME RULE THE VENDORED CODE LOOKS UP BY: key = step + revcomp(step).
# It must NOT be built by iterating _STACK_HS and slicing `key[:2]`. That table holds 36 keys,
# not 16 -- it carries mismatch and wobble steps too, so 2-char prefixes COLLIDE ('GG' prefixes
# four different keys). Slicing lets a mismatch entry overwrite the Watson-Crick one, which
# silently substitutes the wrong physics: 'GA' came out as dH -1.3 instead of -8.2. Inside a '='
# run the target is the complement by construction, which is why the probe step alone determines
# the key.
_DH = np.full((128, 128), np.nan, dtype=np.float64)
_DS = np.full((128, 128), np.nan, dtype=np.float64)
for _a in "ACGT":
    for _b in "ACGT":
        _step = _a + _b
        _hs = _T._STACK_HS.get(_step + _T._revcomp(_step))
        if _hs is not None:
            _DH[ord(_a), ord(_b)] = _hs[0]
            _DS[ord(_a), ord(_b)] = _hs[1]


def _as_matrix(strings, width):
    """
    Pack a sequence of strings into a (rows, width) uint8 matrix, zero-padded.

    Args:
        strings (iterable): the alignment strings.
        width (int): columns to keep.

    Returns:
        mat (numpy.ndarray): uint8, shape (rows, width).
    """
    rows = len(strings)
    mat = np.zeros((rows, width), dtype=np.uint8)
    for i, s in enumerate(strings):
        b = np.frombuffer(s.encode("ascii", "replace")[:width], dtype=np.uint8)
        mat[i, :len(b)] = b

    # success
    return mat


def stacking_profile_batch(probe_alns, ops_list, celsius, width=None):
    """
    Per-column stacking free energy for many duplexes at once.

    Reproduces ``thermo.stacking_profile`` exactly, including its treatment of run boundaries
    and non-ACGT steps, but computes every row in a handful of array operations.

    Args:
        probe_alns (sequence): probe alignment strings.
        ops_list (sequence): CIGAR-op strings, one character per alignment column.
        celsius (array-like): per-row temperature.
        width (int or None): output width. Defaults to the longest probe alignment, which is
            what ``enc_om2`` passes (``len(pa)`` per row) -- see the note below.

    Returns:
        prof (numpy.ndarray): float32, shape (rows, width).
    """
    lengths = np.fromiter((len(s) for s in probe_alns), dtype=np.int64,
                          count=len(probe_alns))
    if width is None:
        width = int(lengths.max()) if len(lengths) else 0

    probe = _as_matrix(probe_alns, width)
    ops = _as_matrix(ops_list, width)

    # a step exists at column k iff columns k and k+1 are both inside the same '=' run
    eq = ops == ord("=")
    step_ok = np.zeros((len(probe_alns), width), dtype=bool)
    step_ok[:, :-1] = eq[:, :-1] & eq[:, 1:]

    # the vendored loop breaks at `col >= width`, and a step needs column k+1 to exist within
    # the row's own length
    col_idx = np.arange(width)[None, :]
    step_ok &= (col_idx + 1) < lengths[:, None]

    left = probe[:, :-1]
    right = probe[:, 1:]
    dh = np.full((len(probe_alns), width), np.nan, dtype=np.float64)
    ds = np.full((len(probe_alns), width), np.nan, dtype=np.float64)
    dh[:, :-1] = _DH[left, right]
    ds[:, :-1] = _DS[left, right]

    # a step whose energy is unknown (contains N, or a non-ACGT byte) contributes nothing,
    # exactly as the vendored `continue` leaves out[col] at its initial 0.0
    valid = step_ok & np.isfinite(dh) & np.isfinite(ds)

    t_kelvin = np.asarray(celsius, dtype=np.float64)[:, None] + 273.15
    out = np.zeros((len(probe_alns), width), dtype=np.float32)
    # boolean indexing, NOT np.putmask: putmask indexes `values` by the FLAT POSITION in `out`
    # and cycles when values is shorter, so a mask-length value array gets scattered to the
    # wrong cells. It produces plausible numbers from the wrong rows, which is worse than an
    # error -- the bit-identity gate is what caught it.
    out[valid] = (dh - t_kelvin * ds).astype(np.float32)[valid]

    # success
    return out


def terminal_block_batch(prof, lengths, K):
    """
    The terminal/interior block for every row at once.

    Reproduces ``features._terminal_block`` exactly. The trick that makes it cheap is that the
    loops run over COLUMNS (width ~40) rather than over rows (millions): a column loop broadcast
    across all rows is ~1000x fewer Python iterations than a row loop.

    Everything is computed in float64, matching the vendored ``np.asarray(prof, np.float64)``,
    because the output is float32 and arithmetic order is observable at that precision.

    Args:
        prof (numpy.ndarray): (rows, width) stacking profile, zero-padded past each row's length.
        lengths (numpy.ndarray): per-row profile length.
        K (int): terminal window.

    Returns:
        s5 (numpy.ndarray): (rows, K) from the 5' terminus.
        s3 (numpy.ndarray): (rows, K) from the 3' terminus, index 0 = terminal step.
        summ (numpy.ndarray): (rows, 8) interior summaries in INTERIOR_COLS order.
        overlap (numpy.ndarray): (rows,) nt by which the two windows overlap.
    """
    n_rows = len(lengths)
    if n_rows == 0 or np.asarray(prof).shape[-1] == 0:
        # a reduction over an empty axis has no identity, so an empty batch is
        # returned in the declared shapes rather than raising
        return (np.zeros((n_rows, K), dtype=np.float64),
                np.zeros((n_rows, K), dtype=np.float64),
                np.zeros((n_rows, 8), dtype=np.float64),
                np.zeros(n_rows, dtype=np.float64))

    p = np.asarray(prof, dtype=np.float64)
    rows, width = p.shape
    L = np.asarray(lengths, dtype=np.int64)
    col = np.arange(width)[None, :]

    # s5: the first min(K, L) entries. The profile is zero past each row's length already, so a
    # straight slice reproduces the vendored zero-fill.
    s5 = np.zeros((rows, K), dtype=np.float64)
    take = min(K, width)
    s5[:, :take] = p[:, :take]
    s5[np.arange(K)[None, :] >= L[:, None]] = 0.0

    # s3: entry j is p[L-1-j], i.e. the profile read backwards from the 3' terminus
    j = np.arange(K)[None, :]
    idx = np.clip(L[:, None] - 1 - j, 0, width - 1)
    s3 = np.take_along_axis(p, idx, axis=1)
    s3[j >= L[:, None]] = 0.0

    # interior: columns [K, L-K), only where the duplex is longer than both windows
    interior_mask = (col >= K) & (col < (L[:, None] - K)) & (L[:, None] > 2 * K)
    has_interior = interior_mask.any(axis=1)

    masked = np.where(interior_mask, p, np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        n_int = interior_mask.sum(axis=1).astype(np.float64)
        i_sum = np.nansum(masked, axis=1)
        i_mean = np.where(has_interior, i_sum / np.where(n_int > 0, n_int, 1), 0.0)
        i_min = np.nanmin(np.where(interior_mask, p, np.inf), axis=1)
        i_max = np.nanmax(np.where(interior_mask, p, -np.inf), axis=1)
        # population std (ddof=0), matching numpy's default that the vendored code relies on
        dev = np.where(interior_mask, p - i_mean[:, None], 0.0)
        i_std = np.sqrt(np.where(has_interior,
                                 (dev ** 2).sum(axis=1) / np.where(n_int > 0, n_int, 1), 0.0))

    # worst = K + argmin(interior), reported as distance to the NEARER end
    worst_col = np.argmin(np.where(interior_mask, p, np.inf), axis=1)
    worst_nt = np.minimum(worst_col, L - 1 - worst_col).astype(np.float64)

    # longest run of "good" steps inside the interior. The loop is over WIDTH, not rows.
    good = interior_mask & (p <= _F.GOOD_STEP)
    run = np.zeros(rows, dtype=np.int64)
    best = np.zeros(rows, dtype=np.int64)
    for c in range(width):
        g = good[:, c]
        run = np.where(g, run + 1, 0)
        best = np.maximum(best, run)

    summ = np.zeros((rows, 8), dtype=np.float64)
    summ[:, 0] = np.where(has_interior, n_int, 0.0)
    summ[:, 1] = np.where(has_interior, i_sum, 0.0)
    summ[:, 2] = np.where(has_interior, i_mean, 0.0)
    summ[:, 3] = np.where(has_interior, i_min, 0.0)
    summ[:, 4] = np.where(has_interior, i_max, 0.0)
    summ[:, 5] = np.where(has_interior, i_std, 0.0)
    # the vendored code reports -1.0 for "no interior", which is a sentinel, not a distance
    summ[:, 6] = np.where(has_interior, worst_nt, -1.0)
    summ[:, 7] = np.where(has_interior, best.astype(np.float64), 0.0)

    overlap = np.maximum(0, 2 * K - L).astype(np.float64)

    # success
    return s5, s3, summ, overlap


def enc_om2_fast(df, celsius=None, sodium=None, K=None):
    """
    The L4t encoding, vectorized. Same 103 columns, same values, same order.

    Args:
        df (pandas.DataFrame): duplex frame with ``probe_aln``, ``target_aln``, ``ops``.
        celsius (float or None): pinned condition, or None for each row's own.
        sodium (float or None): pinned condition, or None for each row's own.
        K (int or None): terminal window. Defaults to the vendored ``TERMINAL_K``.

    Returns:
        X (pandas.DataFrame): 103 features.
    """
    K = _F.TERMINAL_K if K is None else K

    ps = _F.enc_paintshop37(df)
    ps.columns = [f"ps_{c}" for c in ps.columns]

    cel, na = _F._conditions(df, celsius, sodium)
    base = pd.DataFrame({
        "align_score": _F._align_score(df),
        "probe_len": df["probe_seq"].astype(str).str.len().astype(float).values,
        "probe_gc": df["probe_seq"].astype(str).map(_F.gc_pct).values,
        "target_gc": df["target_seq"].astype(str).map(_F.gc_pct).values,
    })
    cig = _F._cigar_block(df)

    # duplex_features stays per row -- its inner work is a string walk -- but goes through the
    # cache the vendored encoder bypasses. Repeat alignments are common, so this is the win.
    th = pd.DataFrame(
        [_T._features_cached(pa, ta, op, float(c), float(s))
         for pa, ta, op, c, s in zip(df.probe_aln, df.target_aln, df.ops, cel, na)],
        columns=_T.FEATURE_COLUMNS,
    )[_F.THERMO_COLS]

    # the stacking profile IS vectorizable, and so is the terminal block that reads it
    probe_alns = list(df.probe_aln)
    prof = stacking_profile_batch(probe_alns, list(df.ops), cel)
    lengths = np.fromiter((len(s) for s in probe_alns), dtype=np.int64, count=len(probe_alns))
    s5s, s3s, sums, ovs = terminal_block_batch(prof, lengths, K)

    # success
    return pd.concat([
        ps.reset_index(drop=True), base, cig, th,
        pd.DataFrame(np.asarray(s5s, np.float32), columns=[f"s5_{j}" for j in range(K)]),
        pd.DataFrame(np.asarray(s3s, np.float32), columns=[f"s3_{j}" for j in range(K)]),
        pd.DataFrame(np.asarray(sums, np.float32), columns=_F.INTERIOR_COLS),
        pd.DataFrame({"overlap_nt": np.asarray(ovs, np.float32)}),
    ], axis=1)


def verify(df, verbose=True):
    """
    Assert the fast encoder is bit-identical to the vendored one.

    This is a shipping gate, not a smoke test. Float32 output means arithmetic order matters,
    so "close" is not good enough: an encoder that changes a feature value is a different
    model wearing the same name.

    Args:
        df (pandas.DataFrame): real duplex rows to check on.
        verbose (bool): print the per-column verdict.

    Returns:
        report (dict): column counts, the worst absolute difference, and pass/fail.

    Raises:
        AssertionError: naming the first column that differs.
    """
    slow = _F.enc_om2(df)
    fast = enc_om2_fast(df)

    assert list(slow.columns) == list(fast.columns), (
        f"column names differ: {set(slow.columns) ^ set(fast.columns)}")

    a = slow.to_numpy(dtype=np.float64)
    b = fast.to_numpy(dtype=np.float64)
    both_nan = np.isnan(a) & np.isnan(b)
    diff = np.where(both_nan, 0.0, np.abs(a - b))
    worst = float(np.nanmax(diff)) if diff.size else 0.0

    bad = np.argwhere(diff > 0)
    if len(bad):
        r, c = bad[0]
        raise AssertionError(
            f"fast encoder differs from vendored at column {slow.columns[c]!r} row {r}: "
            f"{a[r, c]!r} vs {b[r, c]!r} ({len(bad)} cells differ, worst {worst:g}). "
            f"NOT shippable -- this is a different feature space.")

    report = {"n_rows": int(len(df)), "n_columns": int(slow.shape[1]),
              "max_abs_diff": worst, "bit_identical": True}
    if verbose:
        print(f"  bit-identical on {report['n_rows']:,} rows x "
              f"{report['n_columns']} columns (max |diff| = {worst:g})")

    # success
    return report
