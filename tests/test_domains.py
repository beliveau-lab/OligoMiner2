"""Tests for domain-based oligo assembly.

The edge-anchored appending interface puts a sequence at the 5' or 3' end of one
homology region, which covers the PaintSHOP anatomy. A padlock's backbone sits
between two homology arms, so the insertion point is interior and no value of a
left/right flag produces it. These tests cover the ordered-domain generalization
that does, including per-join linkers, which matter because the arms must abut on
the target for ligase to seal the circle.
"""

import random

import pandas as pd
import pytest

from oligominer.probe_design import (
    DomainAssembly,
    assemble_padlock,
    check_backbone_placement,
    mine_padlock_sequence,
    padlocks_to_df,
)
from oligominer.probe_design.appending.config import LINKER
from oligominer.probe_design.domains import LINEAR_LAYOUT, PADLOCK_LAYOUT
from oligominer.utils.exceptions import InvalidInputError


@pytest.fixture
def index():
    return pd.RangeIndex(3)


@pytest.fixture(scope='module')
def padlock_df():
    random.seed(3)
    target = ''.join(random.choice('ACGT') for _ in range(6000))

    # success
    return padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1'))


class TestDomainAssembly:

    def _seamless(self, asm):
        """Clear the default linker from every join in a layout."""
        for left, right in zip(asm.layout, asm.layout[1:]):
            asm.set_linker(left, right, '')
        return asm

    def test_domains_join_in_layout_order(self, index):
        asm = DomainAssembly(index, ['a', 'b', 'c'])
        asm.set_domain('a', 'AAA')
        asm.set_domain('b', 'CCC')
        asm.set_domain('c', 'GGG')
        assert self._seamless(asm).assemble().iloc[0] == 'AAACCCGGG'

    def test_an_empty_domain_makes_its_neighbours_adjacent(self):
        """The join uses the boundary after the last non-empty domain, which is
        configurable; the a-to-c boundary is not, because they are not adjacent."""
        index = pd.RangeIndex(3)
        asm = DomainAssembly(index, ['a', 'b', 'c'])
        asm.set_domain('a', 'AAA')
        asm.set_domain('c', 'GGG')
        asm.set_linker('a', 'b', '')
        assert asm.assemble().iloc[0] == 'AAAGGG'

    def test_an_empty_domain_contributes_no_sequence_of_its_own(self):
        index = pd.RangeIndex(3)
        asm = DomainAssembly(index, ['a', 'b', 'c'])
        asm.set_domain('a', 'AAA')
        asm.set_domain('c', 'GGG')
        asm.set_linker('a', 'b', '')
        assembled = asm.assemble().iloc[0]
        assert len(assembled) == len('AAA') + len('GGG')

    def test_per_row_sequences_are_honored(self, index):
        asm = DomainAssembly(index, ['a', 'b'])
        asm.set_domain('a', pd.Series(['A', 'AA', 'AAA'], index=index))
        asm.set_domain('b', 'T')
        assert list(self._seamless(asm).assemble()) == ['AT', 'AAT', 'AAAT']

    def test_the_default_linker_is_applied_at_every_join(self, index):
        asm = DomainAssembly(index, ['a', 'b', 'c'])
        for name, seq in (('a', 'AAA'), ('b', 'CCC'), ('c', 'GGG')):
            asm.set_domain(name, seq)
        assert asm.assemble().iloc[0] == f'AAA{LINKER}CCC{LINKER}GGG'

    def test_a_linker_is_inserted_at_its_join_only(self, index):
        asm = DomainAssembly(index, ['a', 'b', 'c'])
        for name, seq in (('a', 'AAA'), ('b', 'CCC'), ('c', 'GGG')):
            asm.set_domain(name, seq)
        asm.set_linker('a', 'b', 'GG')
        asm.set_linker('b', 'c', '')
        assert asm.assemble().iloc[0] == 'AAAGGCCCGGG'

    def test_joins_can_carry_different_linkers(self, index):
        asm = DomainAssembly(index, ['a', 'b', 'c'])
        for name, seq in (('a', 'AAA'), ('b', 'CCC'), ('c', 'GGG')):
            asm.set_domain(name, seq)
        asm.set_linker('a', 'b', 'TT')
        asm.set_linker('b', 'c', '')
        assert asm.assemble().iloc[0] == 'AAATTCCCGGG'

    def test_a_duplicate_domain_name_raises(self, index):
        with pytest.raises(InvalidInputError, match='duplicate domain names'):
            DomainAssembly(index, ['a', 'a'])

    def test_an_unknown_domain_raises(self, index):
        asm = DomainAssembly(index, ['a'])
        with pytest.raises(Exception):
            asm.set_domain('nope', 'AAA')

    def test_the_linear_layout_reproduces_the_paintshop_anatomy(self, index):
        asm = DomainAssembly(index, LINEAR_LAYOUT)
        asm.set_domain('homology', 'ACGTACGT')
        asm.set_domain('inner_5p', 'CC')
        asm.set_domain('outer_5p', 'AA')
        asm.set_domain('inner_3p', 'GG')
        asm.set_domain('outer_3p', 'TT')
        assert (self._seamless(asm).assemble().iloc[0]
                == 'AA' + 'CC' + 'ACGTACGT' + 'GG' + 'TT')


class TestPadlockAssembly:

    def test_the_backbone_lands_between_the_arms(self, padlock_df):
        out, _ = assemble_padlock(padlock_df, backbone='TTTGCTAGCTAGCTAGCAAA')
        assert all(check_backbone_placement(row) for _, row in out.iterrows())

    def test_the_layout_is_arm_backbone_arm(self):
        assert PADLOCK_LAYOUT == ['arm_5p', 'backbone', 'arm_3p']

    def test_oligo_length_is_the_arms_plus_the_backbone(self, padlock_df):
        backbone = 'TTTGCTAGCTAGCTAGCAAA'
        out, _ = assemble_padlock(padlock_df, backbone=backbone)
        expected = out['arm5_len'] + out['arm3_len'] + len(backbone)
        assert (out['full_oligo'].str.len() == expected).all()

    def test_ligation_junctions_are_seamless_by_default(self, padlock_df):
        """A linker across the ligation junction would stop ligase sealing the circle."""
        backbone = 'TTTGCTAGCTAGCTAGCAAA'
        out, _ = assemble_padlock(padlock_df, backbone=backbone)
        row = out.iloc[0]
        assert row['full_oligo'] == row['arm_5p'] + backbone + row['arm_3p']

    def test_linkers_can_be_added_when_wanted(self, padlock_df):
        out, _ = assemble_padlock(padlock_df, backbone='GGGG',
                                  linker_5p='TT', linker_3p='AA')
        row = out.iloc[0]
        assert row['full_oligo'] == row['arm_5p'] + 'TT' + 'GGGG' + 'AA' + row['arm_3p']

    def test_the_arms_are_not_disturbed(self, padlock_df):
        out, _ = assemble_padlock(padlock_df, backbone='GGGG')
        assert out['full_oligo'].str.startswith(out['arm_5p'].iloc[0]).iloc[0]

    def test_an_appending_table_is_returned(self, padlock_df):
        _, table = assemble_padlock(padlock_df, backbone='GGGG', backbone_id='bb1')
        assert 'backbone' in table.columns
        assert len(table) == len(padlock_df)

    def test_a_per_row_backbone_is_supported(self, padlock_df):
        backbones = pd.Series(['AAAA'] * len(padlock_df), index=padlock_df.index)
        out, _ = assemble_padlock(padlock_df, backbone=backbones)
        assert out['full_oligo'].str.contains('AAAA').all()


class TestSplitArchitectures:
    """A split-FISH or HCR pair is two oligos, each its own domain layout."""

    def test_a_split_pair_can_be_built_from_one_padlock_record(self, padlock_df):
        left = DomainAssembly(padlock_df.index, ['homology', 'initiator_half'])
        left.set_domain('homology', padlock_df['arm_5p'])
        left.set_domain('initiator_half', 'GAGGAGGGCAGCAAACGG')

        right = DomainAssembly(padlock_df.index, ['initiator_half', 'homology'])
        right.set_domain('homology', padlock_df['arm_3p'])
        right.set_domain('initiator_half', 'AAGAGTCTTCCTTTACG')

        left_oligos, right_oligos = left.assemble(), right.assemble()

        assert len(left_oligos) == len(right_oligos) == len(padlock_df)
        # each member carries its own half of the readout on the correct side
        assert left_oligos.iloc[0].startswith(padlock_df['arm_5p'].iloc[0])
        assert right_oligos.iloc[0].endswith(padlock_df['arm_3p'].iloc[0])
