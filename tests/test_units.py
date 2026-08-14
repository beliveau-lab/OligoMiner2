"""Tests for targeting units.

A split-FISH or HCR pair only produces signal when both members bind, so
dropping one member removes the site while leaving the other in the order,
costing synthesis and contributing off-target binding for no signal. Filtering
therefore has to act on units rather than rows.
"""

import pandas as pd
import pytest

from oligominer.probe_design import (
    assign_units,
    drop_incomplete_units,
    filter_units,
    unit_sizes,
    units_to_orders,
)


@pytest.fixture
def pairs():
    """Three two-member units, one row per oligo."""
    rows = []
    for site, (left, right) in {
        's1': ('AAA', 'TTT'),
        's2': ('CCC', 'GGG'),
        's3': ('ACA', 'TGT'),
    }.items():
        rows.append({'site': site, 'sequence': left})
        rows.append({'site': site, 'sequence': right})

    # success
    return assign_units(pd.DataFrame(rows), by='site', roles=['left', 'right'])


class TestAssignUnits:

    def test_rows_sharing_a_site_share_a_unit(self, pairs):
        assert pairs.groupby('site')['unit_id'].nunique().eq(1).all()

    def test_different_sites_get_different_units(self, pairs):
        assert pairs['unit_id'].nunique() == 3

    def test_roles_are_assigned_in_row_order(self, pairs):
        assert list(pairs['unit_role']) == ['left', 'right'] * 3

    def test_units_can_be_keyed_on_several_columns(self):
        df = pd.DataFrame([{'chrom': 'chr1', 'start': 1, 'sequence': 'A'},
                           {'chrom': 'chr1', 'start': 1, 'sequence': 'T'},
                           {'chrom': 'chr1', 'start': 9, 'sequence': 'C'}])
        out = assign_units(df, by=['chrom', 'start'])
        assert out['unit_id'].nunique() == 2

    def test_roles_are_optional(self):
        df = pd.DataFrame([{'site': 's1', 'sequence': 'A'}])
        assert 'unit_role' not in assign_units(df, by='site').columns


class TestUnitSizes:

    def test_sizes_are_reported_per_unit(self, pairs):
        assert set(unit_sizes(pairs)) == {2}

    def test_a_partial_unit_reports_its_smaller_size(self, pairs):
        assert set(unit_sizes(pairs.iloc[:3])) == {2, 1}


class TestFilterUnits:
    """A row-wise mask must be promoted to a unit-wise decision."""

    def test_a_unit_survives_when_all_members_pass(self, pairs):
        kept, dropped = filter_units(pairs, keep=[True] * 6)
        assert len(kept) == 6
        assert dropped == 0

    def test_one_failing_member_removes_its_whole_unit(self, pairs):
        kept, dropped = filter_units(pairs, keep=[True, False, True, True, True, True])
        assert dropped == 1
        assert len(kept) == 4
        assert 0 not in set(kept['unit_id'])

    def test_no_orphan_is_left_behind(self, pairs):
        """The surviving member of a broken pair must not remain in the output."""
        kept, _ = filter_units(pairs, keep=[True, False, True, True, True, True])
        assert 'AAA' not in set(kept['sequence'])

    def test_surviving_units_keep_every_member(self, pairs):
        kept, _ = filter_units(pairs, keep=[False, False, True, True, True, True])
        assert set(unit_sizes(kept)) == {2}

    def test_failing_every_unit_returns_nothing(self, pairs):
        kept, dropped = filter_units(pairs, keep=[False] * 6)
        assert len(kept) == 0
        assert dropped == 3

    def test_a_row_wise_filter_would_have_left_an_orphan(self, pairs):
        """This is the failure the unit-wise filter exists to prevent."""
        mask = [True, False, True, True, True, True]
        row_wise = pairs[pd.Series(mask, index=pairs.index)]
        assert (unit_sizes(row_wise) == 1).any()

        unit_wise, _ = filter_units(pairs, keep=mask)
        assert not (unit_sizes(unit_wise) == 1).any()


class TestDropIncompleteUnits:

    def test_an_orphaned_member_is_removed(self, pairs):
        orphaned = pairs.iloc[:5]           # unit 2 lost its second member
        kept, dropped = drop_incomplete_units(orphaned)
        assert dropped == 1
        assert set(unit_sizes(kept)) == {2}

    def test_complete_units_are_untouched(self, pairs):
        kept, dropped = drop_incomplete_units(pairs)
        assert len(kept) == len(pairs)
        assert dropped == 0

    def test_the_expected_size_can_be_stated(self, pairs):
        kept, dropped = drop_incomplete_units(pairs, expected_size=3)
        assert len(kept) == 0
        assert dropped == 3

    def test_an_empty_table_is_handled(self):
        empty = pd.DataFrame(columns=['unit_id', 'sequence'])
        kept, dropped = drop_incomplete_units(empty)
        assert len(kept) == 0
        assert dropped == 0


class TestUnitsToOrders:

    def test_one_row_per_unit(self, pairs):
        assert len(units_to_orders(pairs)) == 3

    def test_each_role_becomes_a_column(self, pairs):
        orders = units_to_orders(pairs)
        assert 'left' in orders.columns
        assert 'right' in orders.columns

    def test_members_stay_with_their_unit(self, pairs):
        orders = units_to_orders(pairs).set_index('unit_id')
        first = pairs[pairs['unit_id'] == 0]
        assert orders.loc[0, 'left'] == first.iloc[0]['sequence']
        assert orders.loc[0, 'right'] == first.iloc[1]['sequence']
