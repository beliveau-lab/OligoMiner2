"""Tests that the documented mining config IS the mining API.

``GET_DEFAULT_MINING_CONFIG()`` describes what the miner accepts and
``mine_sequence()`` is how those parameters are passed. These tests assert the two
agree on the set of parameters and on their defaults, and that each documented
parameter both reaches the miner and changes its output.
"""

import inspect

import pytest

from oligominer.thermodynamics.mining import mine_sequence, mine_fasta
from oligominer.thermodynamics.mining.config import GET_DEFAULT_MINING_CONFIG
from oligominer.utils.exceptions import ConfigurationError


# parameters of mine_sequence that describe the input rather than the mining regime
_NOT_CONFIG = {'seq', 'seq_id'}


def _signature_params():
    """Return the mining parameters reachable through mine_sequence()."""
    params = set(inspect.signature(mine_sequence).parameters) - _NOT_CONFIG

    # success
    return params


class TestConfigMatchesSignature:

    def test_every_documented_key_is_accepted(self):
        documented = set(GET_DEFAULT_MINING_CONFIG())
        rejected = documented - _signature_params()
        assert rejected == set(), (
            f"documented but unreachable through mine_sequence(): {sorted(rejected)}"
        )

    def test_every_accepted_param_is_documented(self):
        undocumented = _signature_params() - set(GET_DEFAULT_MINING_CONFIG())
        assert undocumented == set(), (
            f"accepted by mine_sequence() but absent from the config: {sorted(undocumented)}"
        )

    @pytest.mark.parametrize("key", sorted(GET_DEFAULT_MINING_CONFIG()))
    def test_each_documented_key_can_actually_be_passed(self, key):
        """Every documented key is accepted as a keyword argument."""
        value = GET_DEFAULT_MINING_CONFIG()[key]
        mine_sequence('ACGT' * 40, seq_id='t', **{key: value})

    def test_defaults_agree_in_value_not_just_in_name(self):
        """The documented default and the signature default are the same value."""
        documented = GET_DEFAULT_MINING_CONFIG()
        sig = inspect.signature(mine_sequence).parameters

        mismatched = {
            key: (documented[key], sig[key].default)
            for key in documented
            if key in sig and sig[key].default is not inspect.Parameter.empty
            and sig[key].default != documented[key]
        }
        assert mismatched == {}, f"config and signature defaults disagree: {mismatched}"


class TestGcIsReachable:
    """The GC bounds reach the filter and change which probes are returned."""

    # 30 nt of pure GC (100% GC) and pure AT (0% GC), both well outside 20-80
    GC_RICH = 'GCGCGCGCGCGCGCGCGCGCGCGCGCGCGC'
    AT_RICH = 'ATATATATATATATATATATATATATATAT'

    def test_widening_gc_bounds_changes_the_result(self):
        """A GC-rich sequence is filtered at the default bounds and kept when widened."""
        params = dict(seq_id='t', min_length=30, max_length=30,
                      min_tm=0, max_tm=200, max_homopolymer=None)

        default_bounds = mine_sequence(self.GC_RICH, **params)
        widened = mine_sequence(self.GC_RICH, min_gc=0, max_gc=100, **params)

        assert len(default_bounds) == 0
        assert len(widened) > 0

    def test_narrowing_gc_bounds_filters_probes_out(self):
        params = dict(seq_id='t', min_length=30, max_length=30,
                      min_tm=0, max_tm=200, max_homopolymer=None)

        kept = mine_sequence(self.AT_RICH, min_gc=0, max_gc=100, **params)
        filtered = mine_sequence(self.AT_RICH, min_gc=40, max_gc=60, **params)

        assert len(kept) > 0
        assert len(filtered) == 0

    def test_gc_bounds_can_be_disabled(self):
        params = dict(seq_id='t', min_length=30, max_length=30,
                      min_tm=0, max_tm=200, max_homopolymer=None)

        probes = mine_sequence(self.GC_RICH, min_gc=None, max_gc=None, **params)

        assert len(probes) > 0

    def test_gc_reaches_through_mine_fasta(self, example_fasta_path):
        """The kwargs path through mine_fasta must forward it too."""
        wide = mine_fasta(example_fasta_path, min_gc=0, max_gc=100)
        narrow = mine_fasta(example_fasta_path, min_gc=49, max_gc=51)

        assert len(wide) > len(narrow)


class TestAllowOverlap:

    SEQ = 'ACGT' * 60

    def _params(self):
        return dict(seq_id='t', min_length=30, max_length=30,
                    min_tm=0, max_tm=200, max_homopolymer=None,
                    min_gc=0, max_gc=100)

    def test_allow_overlap_actually_changes_behaviour(self):
        overlapping = mine_sequence(self.SEQ, allow_overlap=True, **self._params())
        non_overlapping = mine_sequence(self.SEQ, allow_overlap=False, **self._params())
        assert len(overlapping) > len(non_overlapping)

    def test_non_overlapping_probes_do_not_overlap(self):
        probes = mine_sequence(self.SEQ, allow_overlap=False, **self._params())
        spans = sorted((start, stop) for _, start, stop, _, _ in probes)
        assert all(a[1] <= b[0] for a, b in zip(spans, spans[1:]))


class TestExhaustiveGuard:

    SEQ = 'ACGT' * 60

    def test_exhaustive_with_no_overlap_raises(self):
        with pytest.raises(ConfigurationError):
            mine_sequence(self.SEQ, seq_id='t', exhaustive=True, allow_overlap=False)

    def test_exhaustive_with_spacing_raises(self):
        with pytest.raises(ConfigurationError):
            mine_sequence(self.SEQ, seq_id='t', exhaustive=True, spacing=5)
