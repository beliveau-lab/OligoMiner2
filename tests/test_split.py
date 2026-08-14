"""Tests for split-architecture probe design."""

import pandas as pd
import pytest

from oligominer.probe_design.split import (
    DOWNSTREAM, UPSTREAM, SplitArchitecture, assemble_split, design_split,
    hcr3, pair_probes, pair_summary, split_fish,
)
from oligominer.probe_design.units import (
    ROLE_COLUMN, UNIT_COLUMN, drop_incomplete_units, filter_units,
)
from oligominer.utils.exceptions import InvalidInputError


def probes(spans, seq_id='chr1', seq=None):
    """Build a probe table from (start, stop) spans."""
    return pd.DataFrame({
        'seq_id': seq_id,
        'start': [s for s, _ in spans],
        'stop': [e for _, e in spans],
        'probe_seq': [seq or 'A' * (e - s) for s, e in spans],
    })


class TestSplitArchitecture:

    def test_rejects_inverted_gap_window(self):
        with pytest.raises(InvalidInputError):
            SplitArchitecture('x', min_gap=5, max_gap=2)

    def test_rejects_negative_gap(self):
        with pytest.raises(InvalidInputError):
            SplitArchitecture('x', min_gap=-1)

    def test_append_for_returns_the_declared_sequence(self):
        arch = SplitArchitecture('x', upstream_3p='GGG', downstream_5p='CCC')
        assert arch.append_for(UPSTREAM, '3p') == 'GGG'
        assert arch.append_for(DOWNSTREAM, '5p') == 'CCC'
        assert arch.append_for(UPSTREAM, '5p') == ''

    def test_rejects_unknown_role_and_end(self):
        arch = SplitArchitecture('x')
        with pytest.raises(InvalidInputError):
            arch.append_for('middle', '5p')
        with pytest.raises(InvalidInputError):
            arch.append_for(UPSTREAM, 'left')


class TestPairProbes:

    def test_pairs_probes_inside_the_gap_window(self):
        arch = SplitArchitecture('x', min_gap=0, max_gap=2)
        pairs = pair_probes(probes([(0, 25), (27, 52)]), arch)

        assert len(pairs) == 2
        assert pairs['gap'].tolist() == [2, 2]
        assert pairs[UNIT_COLUMN].tolist() == [0, 0]

    def test_leaves_probes_outside_the_window_unpaired(self):
        arch = SplitArchitecture('x', min_gap=0, max_gap=2)
        assert pair_probes(probes([(0, 25), (40, 65)]), arch).empty

    def test_leaves_overlapping_probes_unpaired(self):
        arch = SplitArchitecture('x', min_gap=0, max_gap=2)
        assert pair_probes(probes([(0, 25), (20, 45)]), arch).empty

    def test_roles_follow_target_coordinate(self):
        arch = SplitArchitecture('x', max_gap=2)
        pairs = pair_probes(probes([(30, 55), (0, 25), (27, 52)]), arch)

        upstream = pairs[pairs[ROLE_COLUMN] == UPSTREAM].iloc[0]
        downstream = pairs[pairs[ROLE_COLUMN] == DOWNSTREAM].iloc[0]
        assert upstream['start'] < downstream['start']

    def test_does_not_pair_across_targets(self):
        arch = SplitArchitecture('x', max_gap=2)
        df = pd.concat([probes([(0, 25)], seq_id='chr1'),
                        probes([(26, 51)], seq_id='chr2')],
                       ignore_index=True)

        assert pair_probes(df, arch).empty

    def test_consumes_both_members_before_pairing_again(self):
        # three adjacent probes yield one pair, not two overlapping ones
        arch = SplitArchitecture('x', min_gap=0, max_gap=2)
        pairs = pair_probes(probes([(0, 25), (26, 51), (52, 77)]), arch)

        assert len(pairs) == 2
        assert pairs['start'].tolist() == [0, 26]

    def test_pairs_every_member_of_a_dense_run(self):
        arch = SplitArchitecture('x', min_gap=0, max_gap=2)
        spans = [(i * 26, i * 26 + 25) for i in range(6)]
        pairs = pair_probes(probes(spans), arch)

        assert len(pairs) == 6
        assert pairs[UNIT_COLUMN].tolist() == [0, 0, 1, 1, 2, 2]

    def test_skips_a_probe_that_blocks_a_valid_pair(self):
        # the middle probe is too far from the first but adjacent to the last
        arch = SplitArchitecture('x', min_gap=0, max_gap=2)
        pairs = pair_probes(probes([(0, 25), (100, 125), (126, 151)]), arch)

        assert pairs['start'].tolist() == [100, 126]

    def test_min_gap_excludes_abutting_probes(self):
        arch = SplitArchitecture('x', min_gap=2, max_gap=4)
        assert pair_probes(probes([(0, 25), (25, 50)]), arch).empty

    def test_empty_input_returns_the_labelled_columns(self):
        arch = SplitArchitecture('x')
        out = pair_probes(probes([]), arch)

        assert out.empty
        assert UNIT_COLUMN in out.columns
        assert ROLE_COLUMN in out.columns

    def test_requires_the_coordinate_columns(self):
        arch = SplitArchitecture('x')
        with pytest.raises(InvalidInputError):
            pair_probes(pd.DataFrame({'seq_id': ['chr1']}), arch)


class TestAssembleSplit:

    def test_each_member_carries_its_own_appends(self):
        arch = SplitArchitecture('x', upstream_3p='GGG', downstream_5p='CCC')
        out = design_split(probes([(0, 25), (26, 51)]), arch)

        upstream = out[out[ROLE_COLUMN] == UPSTREAM].iloc[0]
        downstream = out[out[ROLE_COLUMN] == DOWNSTREAM].iloc[0]
        assert upstream['sequence'] == 'A' * 25 + 'GGG'
        assert downstream['sequence'] == 'CCC' + 'A' * 25

    def test_linker_sits_between_append_and_homology(self):
        arch = SplitArchitecture('x', upstream_3p='GGG', linker='AA')
        out = design_split(probes([(0, 25), (26, 51)]), arch)

        upstream = out[out[ROLE_COLUMN] == UPSTREAM].iloc[0]
        assert upstream['sequence'] == 'A' * 25 + 'AA' + 'GGG'

    def test_no_linker_where_nothing_is_appended(self):
        # the downstream member appends nothing, so it must be bare homology
        arch = SplitArchitecture('x', upstream_3p='GGG', linker='AA')
        out = design_split(probes([(0, 25), (26, 51)]), arch)

        downstream = out[out[ROLE_COLUMN] == DOWNSTREAM].iloc[0]
        assert downstream['sequence'] == 'A' * 25

    def test_requires_pairing_first(self):
        arch = SplitArchitecture('x')
        with pytest.raises(InvalidInputError):
            assemble_split(probes([(0, 25)]), arch)

    def test_requires_the_homology_column(self):
        arch = SplitArchitecture('x')
        pairs = pair_probes(probes([(0, 25), (26, 51)]), arch)
        with pytest.raises(InvalidInputError):
            assemble_split(pairs.drop(columns=['probe_seq']), arch)

    def test_empty_input_returns_the_assembled_columns(self):
        arch = SplitArchitecture('x')
        out = design_split(probes([(0, 25), (500, 525)]), arch)

        assert out.empty
        assert 'sequence' in out.columns


class TestChemistries:

    def test_hcr3_halves_meet_at_the_gap(self):
        # the upstream member's 3' end and the downstream member's 5' end are
        # the ends that face each other once both are bound
        arch = hcr3(initiator_5p='CCCCC', initiator_3p='GGGGG', spacer='')
        out = design_split(probes([(0, 25), (26, 51)]), arch)

        upstream = out[out[ROLE_COLUMN] == UPSTREAM].iloc[0]
        downstream = out[out[ROLE_COLUMN] == DOWNSTREAM].iloc[0]
        assert upstream['sequence'].endswith('GGGGG')
        assert downstream['sequence'].startswith('CCCCC')

    def test_hcr3_spacer_defaults_to_two_bases(self):
        arch = hcr3('CCCCC', 'GGGGG')
        out = design_split(probes([(0, 25), (26, 51)]), arch)

        upstream = out[out[ROLE_COLUMN] == UPSTREAM].iloc[0]
        assert upstream['sequence'] == 'A' * 25 + 'AA' + 'GGGGG'

    def test_split_fish_halves_sit_on_the_outer_ends(self):
        arch = split_fish(bridge_5p='CCCCC', bridge_3p='GGGGG')
        out = design_split(probes([(0, 25), (26, 51)]), arch)

        upstream = out[out[ROLE_COLUMN] == UPSTREAM].iloc[0]
        downstream = out[out[ROLE_COLUMN] == DOWNSTREAM].iloc[0]
        assert upstream['sequence'].startswith('CCCCC')
        assert downstream['sequence'].endswith('GGGGG')

    def test_architecture_name_is_recorded(self):
        arch = hcr3('C', 'G', name='B1')
        out = design_split(probes([(0, 25), (26, 51)]), arch)

        assert set(out['architecture']) == {'B1'}


class TestUnitsIntegration:

    def test_dropping_one_member_removes_the_whole_pair(self):
        arch = SplitArchitecture('x', max_gap=2)
        out = design_split(probes([(0, 25), (26, 51), (60, 85), (86, 111)]),
                           arch)
        out['keep'] = [True, False, True, True]

        survivors, n_dropped = filter_units(out, out['keep'])
        assert len(survivors) == 2
        assert survivors[UNIT_COLUMN].nunique() == 1
        assert n_dropped == 1

    def test_a_row_wise_filter_leaves_a_widow_that_is_then_dropped(self):
        arch = SplitArchitecture('x', max_gap=2)
        out = design_split(probes([(0, 25), (26, 51), (60, 85), (86, 111)]),
                           arch)

        row_wise = out.drop(index=1)
        kept, n_dropped = drop_incomplete_units(row_wise)
        assert len(kept) == 2
        assert n_dropped == 1


class TestPairSummary:

    def test_counts_paired_and_unpaired_probes(self):
        arch = SplitArchitecture('x', max_gap=2)
        df = probes([(0, 25), (26, 51), (500, 525)])
        summary = pair_summary(df, pair_probes(df, arch))

        assert summary['n_probes'] == 3
        assert summary['n_pairs'] == 1
        assert summary['n_unpaired_probes'] == 1
        assert summary['fraction_paired'] == pytest.approx(2 / 3)

    def test_empty_table_reports_zero_rather_than_dividing(self):
        arch = SplitArchitecture('x')
        df = probes([])
        assert pair_summary(df, pair_probes(df, arch))['fraction_paired'] == 0.0
