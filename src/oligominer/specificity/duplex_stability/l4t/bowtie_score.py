"""A bowtie2-consistent alignment score for CONSTRUCTED duplexes.

WHY THIS EXISTS. The controlled-degradation generator (study e) leaves `align_score = NaN` on
every generated row — an honest choice for its purpose, but it is exactly what makes OligoMiner's
LDA unscorable on constructed data (F43). This study's design decision (with CKC) is to give each
degraded duplex the score bowtie2 WOULD assign if that duplex were a genomic candidate, so the LDA
is scorable across the whole pDup range with no NaN — and to VALIDATE that computed score against
real bowtie2 rather than assert it.

THE SCHEME. bowtie2 `--local` (the args big_run.py uses to build the real corpus) scores a local
alignment with match bonus +2 (`--ma 2`) and, for FASTA input at max quality, mismatch penalty −6
(`--mp 6,2`). Local means it may soft-clip a low-identity end if clipping raises the score. For a
substitution-only, equal-length duplex the alignment columns are fixed (column i of probe vs
column i of target), so the local-optimal score is the maximum-sum contiguous window of per-column
scores — a Kadane scan, O(L). Indel families (bulges/truncations) are out of scope here exactly as
in make_augmented (a column map cannot represent them); they are dropped, not mis-scored.

`local_score(ops)` takes the ops string ('=' match, 'X' mismatch) and returns the bowtie2-style
AS. Validation (`validate.py`) aligns real probes to hg38 with the real args and checks this equals
bowtie2's reported AS on the alignments bowtie2 itself chose.
"""

MATCH_BONUS = 2     # bowtie2 --ma 2
MISMATCH_PEN = 6    # bowtie2 --mp 6,2 -> 6 at max (FASTA) quality


def local_score(ops: str, match_bonus: int = MATCH_BONUS, mismatch_pen: int = MISMATCH_PEN) -> int:
    """bowtie2 --local alignment score for a substitution-only column alignment.

    Maximum-sum contiguous window of per-column scores (+match_bonus on '=', -mismatch_pen on 'X'),
    floored at 0 — i.e. bowtie2 will soft-clip both ends to whatever maximizes the score, and an
    all-mismatch alignment scores 0 (fully clipped), never negative. This is Smith-Waterman's local
    rule specialized to a fixed column alignment.
    """
    best = 0
    cur = 0
    for o in ops:
        cur += match_bonus if o == "=" else -mismatch_pen
        if cur < 0:
            cur = 0
        elif cur > best:
            best = cur
    return best


if __name__ == "__main__":
    # tiny self-checks (only values computed by hand)
    assert local_score("=" * 30) == 60                      # perfect 30mer: 2*30
    assert local_score("X" * 30) == 0                       # all mismatch: fully clipped
    assert local_score("X" + "=" * 29) == 58                # leading mismatch clipped away: 2*29
    assert local_score("=" * 15 + "X" + "=" * 14) == 52     # 29 matches (58) minus one mismatch (6), no clip
    print("self-checks pass |",
          "perfect", local_score("=" * 30),
          "| 1 mid mm", local_score("=" * 15 + "X" + "=" * 14),
          "| end mm", local_score("X" + "=" * 29),
          "| 5 spread", local_score("=====X" * 5))
