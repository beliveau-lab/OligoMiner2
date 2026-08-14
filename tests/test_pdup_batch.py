"""Tests for batched pDup.

The batch path is an optimization, not an approximation, so the tests that matter
assert it reproduces calc_pdup exactly. Also covers the two low-level traps that
produce wrong numbers rather than errors: homodimer symmetry, and results being
keyed by sequence content rather than submission order.
"""

import random

import numpy as np
import pytest

nupack = pytest.importorskip('nupack')

from oligominer.thermodynamics.nupack import (          # noqa: E402
    add_pdup_batch,
    calc_pdup,
    calc_pdup_many,
    calc_pdup_one_to_many,
    low_level_available,
)
from oligominer.utils.seq_utils import rev_comp         # noqa: E402

# tolerance for an equilibrium solve reached by two different routes
EXACT = 1e-9


@pytest.fixture(scope='module')
def duplex_pairs():
    """Probe and target pairs spanning perfect matches through heavy mismatch."""
    random.seed(5)
    pairs = []
    for n_mismatch in (0, 1, 3, 6):
        for _ in range(3):
            probe = ''.join(random.choice('ACGT') for _ in range(30))
            target = list(rev_comp(probe))
            for _ in range(n_mismatch):
                j = random.randrange(len(target))
                target[j] = random.choice([b for b in 'ACGT' if b != target[j]])
            pairs.append((probe, ''.join(target)))

    # success
    return pairs


class TestAvailability:

    def test_low_level_api_is_reported(self):
        assert isinstance(low_level_available(), bool)


class TestExactness:
    """The batch path must reproduce tube_analysis, not approximate it."""

    def test_batched_matches_calc_pdup_pair_for_pair(self, duplex_pairs):
        batched = calc_pdup_many(duplex_pairs)
        reference = [calc_pdup(a, b, conc_a=1e-6, conc_b=1e-12) for a, b in duplex_pairs]

        assert len(batched) == len(reference)
        worst = max(abs(x - y) for x, y in zip(batched, reference))
        assert worst < EXACT, f'largest deviation {worst:.3e} over {len(batched)} pairs'

    def test_one_to_many_matches_calc_pdup(self, duplex_pairs):
        probe = duplex_pairs[0][0]
        targets = [t for _, t in duplex_pairs]

        batched = calc_pdup_one_to_many(probe, targets)
        reference = [calc_pdup(probe, t, conc_a=1e-6, conc_b=1e-12) for t in targets]

        worst = max(abs(x - y) for x, y in zip(batched, reference))
        assert worst < EXACT, f'largest deviation {worst:.3e}'

    def test_a_self_complementary_probe_is_handled(self):
        """Homodimer symmetry must be corrected or pDup comes out low."""
        palindrome = 'ACGTACGTACGTACGTACGTACGTACGTAC'
        target = rev_comp(palindrome)

        batched = calc_pdup_many([(palindrome, target)])[0]
        reference = calc_pdup(palindrome, target, conc_a=1e-6, conc_b=1e-12)

        assert abs(batched - reference) < EXACT

    def test_repeated_pairs_stay_aligned_to_their_positions(self, duplex_pairs):
        """Results are keyed by sequence content, so duplicates must not misalign."""
        perfect = duplex_pairs[0]
        mismatched = duplex_pairs[-1]
        repeated = [perfect, mismatched, perfect, mismatched]
        values = calc_pdup_many(repeated)

        # a duplicated pair must return the same value at both of its positions
        assert values[0] == pytest.approx(values[2], abs=EXACT)
        assert values[1] == pytest.approx(values[3], abs=EXACT)
        # and the two distinct pairs must not have been collapsed onto each other
        assert abs(values[0] - values[1]) > 1e-3

    def test_batching_boundary_does_not_change_the_answer(self, duplex_pairs):
        small = calc_pdup_many(duplex_pairs, batch_complexes=10)
        large = calc_pdup_many(duplex_pairs, batch_complexes=10_000)

        assert np.allclose(small, large, atol=EXACT)

    def test_one_to_many_batching_boundary_is_stable(self, duplex_pairs):
        probe = duplex_pairs[0][0]
        targets = [t for _, t in duplex_pairs]

        small = calc_pdup_one_to_many(probe, targets, batch_complexes=8)
        large = calc_pdup_one_to_many(probe, targets, batch_complexes=10_000)

        assert np.allclose(small, large, atol=EXACT)


class TestBehaviour:

    def test_a_perfect_match_beats_a_mismatched_one(self):
        probe = 'ACGTACGTACGTACGTACGTACGTACGTAC'
        perfect = rev_comp(probe)
        broken = list(perfect)
        for j in range(0, len(broken), 4):
            broken[j] = random.choice([b for b in 'ACGT' if b != broken[j]])

        values = calc_pdup_many([(probe, perfect), (probe, ''.join(broken))])

        assert values[0] > values[1]

    def test_values_are_probabilities(self, duplex_pairs):
        values = calc_pdup_many(duplex_pairs)
        # a converged concentration ratio can land a hair above 1
        assert all(0.0 <= v <= 1.0 + 1e-6 for v in values)

    def test_empty_input_returns_empty(self):
        assert calc_pdup_many([]) == []
        assert calc_pdup_one_to_many('ACGTACGTACGTACGTACGTACGT', []) == []

    def test_model_temperature_is_honored(self, duplex_pairs):
        """The equilibrium solve must run at the model's own temperature."""
        cold = nupack.Model(material='dna', ensemble='stacking',
                            celsius=25.0, sodium=0.39, magnesium=0.0)
        hot = nupack.Model(material='dna', ensemble='stacking',
                           celsius=85.0, sodium=0.39, magnesium=0.0)

        pair = [duplex_pairs[0]]
        assert calc_pdup_many(pair, model=cold)[0] > calc_pdup_many(pair, model=hot)[0]

    def test_temperature_matches_calc_pdup_at_that_temperature(self, duplex_pairs):
        cold = nupack.Model(material='dna', ensemble='stacking',
                            celsius=25.0, sodium=0.39, magnesium=0.0)
        batched = calc_pdup_many(duplex_pairs[:4], model=cold)
        reference = [calc_pdup(a, b, conc_a=1e-6, conc_b=1e-12, model=cold)
                     for a, b in duplex_pairs[:4]]

        assert np.allclose(batched, reference, atol=EXACT)


class TestDataFrameEntryPoint:

    def test_pdup_column_is_added(self, duplex_pairs):
        import pandas as pd

        df = pd.DataFrame([{'probe_seq': a, 'derived_seq': b} for a, b in duplex_pairs])
        out = add_pdup_batch(df)

        assert 'pdup' in out.columns
        assert len(out) == len(df)
        assert out['pdup'].between(0, 1 + 1e-6).all()

    def test_values_match_the_pairwise_reference(self, duplex_pairs):
        import pandas as pd

        df = pd.DataFrame([{'probe_seq': a, 'derived_seq': b} for a, b in duplex_pairs])
        out = add_pdup_batch(df)
        reference = [calc_pdup(a, b, conc_a=1e-6, conc_b=1e-12) for a, b in duplex_pairs]

        assert np.allclose(out['pdup'].to_numpy(), reference, atol=EXACT)

    def test_rows_stay_with_their_probe_when_grouped(self, duplex_pairs):
        """Grouping by probe must not permute results back onto the wrong rows."""
        import pandas as pd

        rows = []
        for probe, target in duplex_pairs:
            rows.append({'probe_seq': probe, 'derived_seq': target})
            rows.append({'probe_seq': duplex_pairs[0][0], 'derived_seq': target})
        df = pd.DataFrame(rows)

        out = add_pdup_batch(df)
        reference = [calc_pdup(r.probe_seq, r.derived_seq, conc_a=1e-6, conc_b=1e-12)
                     for r in df.itertuples()]

        assert np.allclose(out['pdup'].to_numpy(), reference, atol=EXACT)

    def test_a_non_default_index_is_preserved(self, duplex_pairs):
        import pandas as pd

        df = pd.DataFrame([{'probe_seq': a, 'derived_seq': b} for a, b in duplex_pairs],
                          index=[f'row{i}' for i in range(len(duplex_pairs))])
        out = add_pdup_batch(df)
        reference = [calc_pdup(a, b, conc_a=1e-6, conc_b=1e-12) for a, b in duplex_pairs]

        assert list(out.index) == list(df.index)
        assert np.allclose(out['pdup'].to_numpy(), reference, atol=EXACT)
