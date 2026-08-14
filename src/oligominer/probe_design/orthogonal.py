"""
# Orthogonal sequence design

Finds which sequences in a set could cross-hybridize, so a barcode or readout set
can be made mutually orthogonal.

## Why nominating with k-mers is exact, not a heuristic

Comparing every pair at every offset is quadratic: 240,000 barcodes is 2.9e10
pairs. But two sequences can only cross-hybridize where they share a
complementary run, and a shared complementary run of length k or more is exactly
a shared k-mer between one sequence and the reverse complement of the other.

So indexing k-mers and joining returns the same pair set an all-pairs screen
would, for the stated minimum run length. `verify_against_screen` asserts that
equality on a small set rather than arguing for it.

## The three stages

**Nominate** with k-mers, which is exact for runs of length >= k.
**Rank** the nominated pairs with a registered model, which is cheap.
**Verify** the worst survivors with exact NUPACK, and act on that.

## Rank the model output, do not threshold it

A model fitted on genomic duplexes sees a very different population here. On
random barcode pairs its output concentrates well above the binder threshold
while the true binder rate is near zero, so an absolute cutoff evicts most of a
panel. Use the ordering to choose what to verify, and let the physics decide.
"""

import numpy as np
import pandas as pd

from oligominer.utils.seq_utils import rev_comp

# minimum complementary run treated as a cross-hybridization risk
DEFAULT_MIN_RUN = 7

PAIR_COLUMNS = ['i', 'j', 'shared_kmer', 'n_shared']


def _kmers(seq, k):
    """
    Return the distinct k-mers of a sequence.

    Args:
        seq (str): the sequence.
        k (int): k-mer length.

    Returns:
        kmers (set): the distinct k-mers.
    """
    seq = seq.upper()

    # success
    return {seq[i:i + k] for i in range(len(seq) - k + 1)}


def nominate(seqs, min_run=DEFAULT_MIN_RUN):
    """
    Return every pair of sequences sharing a complementary run.

    A pair is nominated when one sequence contains a k-mer that also occurs in
    the reverse complement of the other, with k = min_run. That is the exact
    condition for a complementary run of that length, so no qualifying pair is
    missed.

    Args:
        seqs (sequence): the sequences to screen.
        min_run (int): the shortest complementary run treated as a risk.

    Returns:
        pairs (pandas.DataFrame): one row per nominated pair with PAIR_COLUMNS,
            where i and j are indices into seqs and i < j.
    """
    seqs = [str(s).upper() for s in seqs]

    # every k-mer of every sequence, and of every reverse complement
    forward = {}
    for index, seq in enumerate(seqs):
        for kmer in _kmers(seq, min_run):
            forward.setdefault(kmer, set()).add(index)

    hits = {}
    for index, seq in enumerate(seqs):
        for kmer in _kmers(rev_comp(seq), min_run):
            for other in forward.get(kmer, ()):
                if other == index:
                    continue
                key = (min(index, other), max(index, other))
                entry = hits.setdefault(key, set())
                entry.add(kmer)

    rows = [{'i': i, 'j': j, 'shared_kmer': sorted(kmers)[0],
             'n_shared': len(kmers)}
            for (i, j), kmers in sorted(hits.items())]

    # success
    return pd.DataFrame(rows, columns=PAIR_COLUMNS)


def screen_bruteforce(seqs, min_run=DEFAULT_MIN_RUN):
    """
    Return the same pair set by direct comparison, for verification.

    Every pair is compared at every ungapped offset. This is the definition the
    k-mer join is claimed to reproduce, and is far too slow for a real panel.

    Args:
        seqs (sequence): the sequences to screen.
        min_run (int): the shortest complementary run treated as a risk.

    Returns:
        pairs (set): (i, j) tuples with i < j.
    """
    seqs = [str(s).upper() for s in seqs]
    complements = [rev_comp(s) for s in seqs]

    pairs = set()
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            if _longest_common_substring(seqs[i], complements[j]) >= min_run:
                pairs.add((i, j))

    # success
    return pairs


def _longest_common_substring(a, b):
    """
    Return the length of the longest substring shared by two sequences.

    Args:
        a (str): the first sequence.
        b (str): the second sequence.

    Returns:
        longest (int): the length of the longest shared substring.
    """
    previous = [0] * (len(b) + 1)
    longest = 0
    for i in range(1, len(a) + 1):
        current = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                current[j] = previous[j - 1] + 1
                if current[j] > longest:
                    longest = current[j]
        previous = current

    # success
    return longest


def verify_against_screen(seqs, min_run=DEFAULT_MIN_RUN):
    """
    Assert the k-mer join returns exactly the brute-force pair set.

    Args:
        seqs (sequence): the sequences to screen.
        min_run (int): the shortest complementary run treated as a risk.

    Returns:
        report (dict): the two set sizes, whether they are equal, and any pairs
            present in one and not the other.
    """
    nominated = {(int(row.i), int(row.j))
                 for row in nominate(seqs, min_run).itertuples()}
    screened = screen_bruteforce(seqs, min_run)

    # success
    return {
        'n_nominated': len(nominated),
        'n_screened': len(screened),
        'equal': nominated == screened,
        'missed_by_nominator': sorted(screened - nominated),
        'extra_in_nominator': sorted(nominated - screened),
    }


def worst_pairs(pairs, scores, n=None):
    """
    Order nominated pairs by a score, worst first.

    Args:
        pairs (pandas.DataFrame): nominated pairs.
        scores (array-like): one score per pair, higher meaning more risk.
        n (int, optional): keep only this many.

    Returns:
        ordered (pandas.DataFrame): the pairs with a score column, worst first.
    """
    ordered = pairs.copy()
    ordered['score'] = np.asarray(scores, dtype=np.float64)
    ordered = ordered.sort_values('score', ascending=False)

    # success
    return ordered.head(n) if n else ordered


def verify_pairs(seqs, pairs, model=None):
    """
    Compute exact pDup for a set of nominated pairs.

    Args:
        seqs (sequence): the sequences the pair indices refer to.
        pairs (pandas.DataFrame): nominated pairs carrying i and j.
        model (nupack.Model, optional): the thermodynamic model.

    Returns:
        out (pandas.DataFrame): the pairs with an exact pdup column.
    """
    from oligominer.thermodynamics.nupack import calc_pdup_many

    seqs = [str(s).upper() for s in seqs]
    duplexes = [(seqs[int(row.i)], seqs[int(row.j)]) for row in pairs.itertuples()]

    out = pairs.copy()
    out['pdup'] = calc_pdup_many(duplexes, model=model)

    # success
    return out


def evict(seqs, pairs, keep_first=True):
    """
    Choose sequences to drop so no nominated pair survives.

    Walks the pairs in the order given and drops one member of each pair that is
    still intact, which is a greedy vertex cover rather than a minimum one.

    Args:
        seqs (sequence): the sequences.
        pairs (pandas.DataFrame): the pairs to break, worst first.
        keep_first (bool): drop the higher-indexed member, keeping the earlier
            sequence. Set False to drop the lower-indexed one.

    Returns:
        kept (list): indices of the surviving sequences.
        dropped (list): indices of the evicted sequences.
    """
    dropped = set()
    for row in pairs.itertuples():
        i, j = int(row.i), int(row.j)
        if i in dropped or j in dropped:
            continue
        dropped.add(j if keep_first else i)

    kept = [index for index in range(len(seqs)) if index not in dropped]

    # success
    return kept, sorted(dropped)
