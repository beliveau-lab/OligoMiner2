"""Tests for probe appending and scoring.

Covers: all appending schemes (same, unique, multiple, custom),
SABER appending, MERFISH barcodes, master table, scoring functions,
and the ProbeSet appending API.
"""

import numpy as np
import pandas as pd
import pytest

from oligominer.probe_design.appending.appending import append_sequences
from oligominer.utils.exceptions import InvalidInputError

from oligominer.data.appending import (
    load_bridges,
    load_outer_forward,
    load_inner_forward,
    load_saber_1x,
    load_merfish_bridges,
)
from oligominer.probe_design.appending import (
    append_same,
    append_unique,
    append_multiple,
    append_custom,
    append_saber,
    append_barcodes,
    build_appending_table,
)
from oligominer.probe_design.scoring import label_on_target, score_probes
from oligominer import ProbeSet


# ---------------------------------------------------------------------------
# appending schemes
# ---------------------------------------------------------------------------

class TestAppendSame:

    def test_prepends_to_all(self, synthetic_probe_df):
        primers = load_outer_forward()
        result, entries = append_same(synthetic_probe_df, primers, left=True, rc=False)
        primer_seq = primers['seq'].iloc[0]
        # every probe should start with the primer sequence
        for seq in result['sequence']:
            assert primer_seq in seq

    def test_entries_are_uniform(self, synthetic_probe_df):
        primers = load_outer_forward()
        _, entries = append_same(synthetic_probe_df, primers, left=True)
        assert entries.nunique() == 1

    def test_append_right(self, synthetic_probe_df):
        primers = load_outer_forward()
        result, _ = append_same(synthetic_probe_df, primers, left=False, rc=False)
        primer_seq = primers['seq'].iloc[0]
        for seq in result['sequence']:
            assert seq.endswith(primer_seq)


class TestAppendUnique:

    def test_one_per_target(self, synthetic_probe_df):
        primers = load_inner_forward()
        result, entries = append_unique(
            synthetic_probe_df, primers,
            target_column='target', left=True, rc=False,
        )
        # each target should get a different primer
        target_entries = {}
        for target in synthetic_probe_df['target'].unique():
            mask = result['target'] == target
            target_entries[target] = entries[mask].iloc[0]
        assert len(set(target_entries.values())) == len(target_entries)


class TestAppendMultiple:

    def test_n_per_target(self, synthetic_probe_df):
        bridges = load_bridges()
        n_per_target = 3
        result, entries = append_multiple(
            synthetic_probe_df, bridges, n_per_target,
            target_column='target', left=False, rc=False,
        )
        for target in synthetic_probe_df['target'].unique():
            mask = result['target'] == target
            unique = entries[mask].unique()
            assert len(unique) <= n_per_target


class TestAppendCustom:

    def test_custom_ranges(self, synthetic_probe_df):
        primers = load_outer_forward()
        ranges = ["1-5", "6-10", "11-15", "16-20"]
        result, entries = append_custom(
            synthetic_probe_df, primers, ranges, left=False, rc=False,
        )
        assert len(result) == len(synthetic_probe_df)
        # different ranges should get different primers
        assert entries.iloc[0] != entries.iloc[5]


# ---------------------------------------------------------------------------
# SABER appending
# ---------------------------------------------------------------------------

class TestSaberAppending:

    def test_saber_unique(self, synthetic_probe_df):
        saber = load_saber_1x()
        result, entries = append_saber(
            synthetic_probe_df, saber,
            scheme='unique', target_column='target',
        )
        assert len(result) == len(synthetic_probe_df)
        # sequences should be longer after appending
        for orig, appended in zip(
            synthetic_probe_df['sequence'], result['sequence']
        ):
            assert len(appended) > len(orig)


# ---------------------------------------------------------------------------
# MERFISH barcodes
# ---------------------------------------------------------------------------

class TestMerfishBarcodes:

    def test_barcode_appending(self, synthetic_probe_df):
        bridges = load_merfish_bridges()
        barcodes = pd.DataFrame({
            'barcode': [
                "1001000100010000",
                "0100100010000100",
                "0010010001000010",
                "0001001000100001",
            ]
        })
        result = append_barcodes(
            synthetic_probe_df, bridges, barcodes,
            target_column='refseq',
        )
        assert len(result) == len(synthetic_probe_df)
        for orig, appended in zip(
            synthetic_probe_df['sequence'], result['sequence']
        ):
            assert len(appended) > len(orig)


# ---------------------------------------------------------------------------
# master table
# ---------------------------------------------------------------------------

class TestAppendingTable:

    def test_build_appending_table(self, synthetic_probe_df):
        primers = load_outer_forward()
        _, entries_a = append_same(synthetic_probe_df, primers, left=True)
        _, entries_b = append_same(synthetic_probe_df, primers, left=False)

        table = build_appending_table(
            synthetic_probe_df,
            {'left_primer': entries_a, 'right_primer': entries_b},
        )
        assert 'left_primer' in table.columns
        assert 'right_primer' in table.columns
        assert len(table) == len(synthetic_probe_df)

    def test_build_appending_table(self, synthetic_probe_df):
        """build_appending_table assembles a tracking table from entries."""
        primers = load_outer_forward()
        _, entries = append_same(synthetic_probe_df, primers, left=True)

        table = build_appending_table(
            synthetic_probe_df, {'primer': entries},
        )
        assert 'primer' in table.columns


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------

class TestScoring:

    def test_label_on_target(self, synthetic_merged_df):
        labeled = label_on_target(synthetic_merged_df)
        assert 'on_target' in labeled.columns
        # 5 probes, each with 1 on-target alignment
        assert labeled['on_target'].sum() == 5

    def test_score_probes(self, synthetic_merged_df):
        labeled = label_on_target(synthetic_merged_df)
        scores = score_probes(labeled, pred_column='duplex_pred')
        assert 'on_target_score' in scores.columns
        assert 'off_target_score' in scores.columns
        assert len(scores) == 5  # 5 unique probes

    def test_on_target_higher_than_off_target(self, synthetic_merged_df):
        labeled = label_on_target(synthetic_merged_df)
        scores = score_probes(labeled, pred_column='duplex_pred')
        # each probe has on_target=0.85, off_target=2*0.05=0.10
        for _, row in scores.iterrows():
            assert row['on_target_score'] > row['off_target_score']


# ---------------------------------------------------------------------------
# ProbeSet appending API
# ---------------------------------------------------------------------------

class TestProbeSetAppending:

    def test_append_method(self, synthetic_probe_df):
        ps = ProbeSet(synthetic_probe_df)
        primers = load_outer_forward()
        result = ps.append(
            primers, scheme='same', label='outer_fwd', left=True, rc=False,
        )
        assert result is ps  # fluent interface
        assert 'sequence' in ps.df.columns
        assert ps.master_table is not None

    def test_chained_appending(self, synthetic_probe_df):
        ps = ProbeSet(synthetic_probe_df)
        of_primers = load_outer_forward()
        if_primers = load_inner_forward()

        ps.append(of_primers, scheme='same', label='outer_fwd',
                  left=True, rc=True)
        ps.append(if_primers, scheme='unique', label='inner_fwd',
                  target_column='target', left=True, rc=True)

        master = ps.master_table
        assert 'outer_fwd' in master.columns
        assert 'inner_fwd' in master.columns


class TestSchemeArgumentValidation:
    """
    Each scheme reaches for its required argument deep inside the assignment,
    where a missing one surfaced as a KeyError on None or an arithmetic error
    rather than as a statement of what the caller left out.
    """

    @pytest.fixture
    def probes(self):
        return pd.DataFrame({'sequence': ['ACGT', 'TTTT'],
                             'target': ['a', 'b']})

    @pytest.fixture
    def sequences(self):
        return pd.DataFrame({'id': ['s1'], 'seq': ['GGG']})

    def test_unique_without_a_target_column(self, probes, sequences):
        with pytest.raises(InvalidInputError, match='target_column'):
            append_sequences(probes, sequences, 'unique')

    def test_multiple_without_a_count(self, probes, sequences):
        with pytest.raises(InvalidInputError, match='n_per_target'):
            append_sequences(probes, sequences, 'multiple',
                             target_column='target')

    def test_multiple_without_a_target_column(self, probes, sequences):
        with pytest.raises(InvalidInputError, match='target_column'):
            append_sequences(probes, sequences, 'multiple', n_per_target=1)

    def test_custom_without_ranges(self, probes, sequences):
        with pytest.raises(InvalidInputError, match='ranges'):
            append_sequences(probes, sequences, 'custom')

    def test_same_needs_none_of_them(self, probes, sequences):
        result, _ = append_sequences(probes, sequences, 'same')
        assert result['sequence'].tolist() == ['GGGTTTACGT', 'GGGTTTTTTT']

    def test_an_unknown_scheme_is_still_rejected(self, probes, sequences):
        with pytest.raises(InvalidInputError, match='Unknown appending'):
            append_sequences(probes, sequences, 'nonsense')


class TestBarcodeShapeValidation:
    """
    A malformed MHD4 barcode surfaced as an IndexError from deep inside bridge
    selection, which says nothing about the barcode being wrong.
    """

    @pytest.fixture
    def probes(self):
        return pd.DataFrame({'sequence': ['ACGT', 'TTTT'],
                             'refseq': ['a', 'b']})

    @pytest.fixture
    def bridges(self):
        return pd.DataFrame({'id': [f'b{i}' for i in range(16)],
                             'seq': ['ACGT'] * 16})

    def _barcodes(self, code):
        return pd.DataFrame({'barcode': [code, code]})

    def test_a_valid_codeword_is_accepted(self, probes, bridges):
        from oligominer.probe_design.appending.paintshop_appending import (
            append_barcodes,
        )

        result = append_barcodes(probes, bridges,
                                 self._barcodes('1111000000000000'))
        assert len(result) == len(probes)

    def test_a_barcode_narrower_than_the_bridge_set_is_rejected(self, probes,
                                                                bridges):
        from oligominer.probe_design.appending.paintshop_appending import (
            append_barcodes,
        )

        with pytest.raises(InvalidInputError, match='16 bridges'):
            append_barcodes(probes, bridges, self._barcodes('1000'))

    def test_a_non_binary_barcode_is_rejected(self, probes, bridges):
        from oligominer.probe_design.appending.paintshop_appending import (
            append_barcodes,
        )

        with pytest.raises(InvalidInputError, match='not binary'):
            append_barcodes(probes, bridges,
                            self._barcodes('abcd000000000000'))

    def test_a_codeword_of_the_wrong_weight_is_rejected(self, probes, bridges):
        from oligominer.probe_design.appending.paintshop_appending import (
            append_barcodes,
        )

        with pytest.raises(InvalidInputError, match='exactly 4'):
            append_barcodes(probes, bridges,
                            self._barcodes('1100000000000000'))

    def test_collect_indices_returns_the_set_positions(self):
        from oligominer.probe_design.appending.paintshop_appending import (
            collect_indices,
        )

        assert collect_indices('1000000010001001') == [1, 9, 13, 16]
