"""
Split-architecture probe design.

A split architecture places one functional element across two oligos that bind
adjacent sites on the same target. The element is only complete when both
members are bound, so background from a single stray binding event does not
produce signal. HCR 3.0 split-initiator pairs and split-FISH bridge pairs are
both this shape, differing only in which sequences the members carry and how
close the two footprints must be.

The design has two steps. `pair_probes` finds the adjacent pairs a probe table
supports, subject to a gap window measured on the target. `assemble_split`
attaches the architecture's domains to each member and concatenates them.
`design_split` runs both.

An architecture declares, per member, what is appended at each end:

    upstream oligo:    5'--[upstream_5p]--[homology]--[upstream_3p]--3'
    downstream oligo:  5'--[downstream_5p]--[homology]--[downstream_3p]--3'

`upstream` is the member whose footprint starts earlier on the target, so the
roles are fixed by coordinate and do not depend on which strand the probes were
mined from. A chemistry that needs the two halves of its element adjacent in the
assembled construct puts them on the ends that face each other once both are
bound; `hcr3` and `split_fish` below build the architecture for each.

Pairs are targeting units: both members are needed for one site to give signal.
`pair_probes` labels them with `units.UNIT_COLUMN` and `units.ROLE_COLUMN`, so
`units.filter_units` and `units.drop_incomplete_units` keep pairs intact through
any later row-wise filtering.
"""

import pandas as pd

from oligominer.probe_design.domains import DomainAssembly
from oligominer.probe_design.units import ROLE_COLUMN, UNIT_COLUMN
from oligominer.utils.exceptions import InvalidInputError

# the two members of a pair, named by their position on the target
UPSTREAM = 'upstream'
DOWNSTREAM = 'downstream'
ROLES = [UPSTREAM, DOWNSTREAM]

# domain layout each member is assembled from
SPLIT_LAYOUT = ['append_5p', 'homology', 'append_3p']


class SplitArchitecture:
    """
    What each member of a split pair carries, and how close the pair must be.

    Attributes:
        name (str): the chemistry this describes, recorded on the output.
        appends (dict): (role, end) -> sequence appended there, where role is
            'upstream' or 'downstream' and end is "5p" or "3p".
        linker (str): inserted between an appended domain and the homology.
        min_gap (int): fewest bases allowed between the two footprints.
        max_gap (int): most bases allowed between the two footprints.
    """

    def __init__(self, name, upstream_5p='', upstream_3p='',
                 downstream_5p='', downstream_3p='', linker='',
                 min_gap=0, max_gap=2):
        if min_gap > max_gap:
            raise InvalidInputError(
                f'min_gap {min_gap} exceeds max_gap {max_gap}')
        if min_gap < 0:
            raise InvalidInputError(f'min_gap must be >= 0, got {min_gap}')

        self.name = name
        self.appends = {
            (UPSTREAM, '5p'): upstream_5p,
            (UPSTREAM, '3p'): upstream_3p,
            (DOWNSTREAM, '5p'): downstream_5p,
            (DOWNSTREAM, '3p'): downstream_3p,
        }
        self.linker = linker
        self.min_gap = min_gap
        self.max_gap = max_gap

    def append_for(self, role, end):
        """
        Return the sequence this architecture appends to one end of one member.

        Args:
            role (str): 'upstream' or 'downstream'.
            end (str): "5p" or "3p".

        Returns:
            sequence (str): the appended sequence, "" when nothing is appended.
        """
        if role not in ROLES:
            raise InvalidInputError(f'unknown role {role!r}; expected {ROLES}')
        if end not in ('5p', '3p'):
            raise InvalidInputError(f"unknown end {end!r}; expected '5p' or '3p'")

        # success
        return self.appends[(role, end)]

    def __repr__(self):
        return (f'SplitArchitecture({self.name!r}, '
                f'gap {self.min_gap}-{self.max_gap})')


def hcr3(initiator_5p, initiator_3p, name='hcr3', spacer='AA',
         min_gap=0, max_gap=2):
    """
    Build the architecture for an HCR 3.0 split-initiator pair.

    The amplifier's initiator is split in two. The upstream member carries the
    3' half on its 3' end and the downstream member carries the 5' half on its
    5' end, so once both are bound the two halves meet at the gap between the
    footprints and present a complete initiator to the hairpins.

    Initiator sequences are amplifier-specific and are supplied by the caller,
    so a set can be used without waiting for it to be added to this package.

    Args:
        initiator_5p (str): the initiator half carried by the downstream
            member, written 5'->3'.
        initiator_3p (str): the initiator half carried by the upstream member,
            written 5'->3'.
        name (str): recorded on the output, e.g. the amplifier's name.
        spacer (str): inserted between an initiator half and the homology.
        min_gap (int): fewest bases between the two footprints.
        max_gap (int): most bases between the two footprints.

    Returns:
        architecture (SplitArchitecture): ready for design_split.
    """
    # success
    return SplitArchitecture(
        name, upstream_3p=initiator_3p, downstream_5p=initiator_5p,
        linker=spacer, min_gap=min_gap, max_gap=max_gap)


def split_fish(bridge_5p, bridge_3p, name='split-fish', spacer='',
               min_gap=0, max_gap=1):
    """
    Build the architecture for a split-FISH bridge pair.

    Each member carries half of the readout bridge on the end facing away from
    its partner, so the two halves sit at the outer ends of the bound pair and
    the bridge is only complete at a site where both members bound.

    Args:
        bridge_5p (str): the bridge half carried by the upstream member, placed
            on its 5' end.
        bridge_3p (str): the bridge half carried by the downstream member,
            placed on its 3' end.
        name (str): recorded on the output.
        spacer (str): inserted between a bridge half and the homology.
        min_gap (int): fewest bases between the two footprints.
        max_gap (int): most bases between the two footprints.

    Returns:
        architecture (SplitArchitecture): ready for design_split.
    """
    # success
    return SplitArchitecture(
        name, upstream_5p=bridge_5p, downstream_3p=bridge_3p,
        linker=spacer, min_gap=min_gap, max_gap=max_gap)


def pair_probes(df, architecture, by='seq_id', start_column='start',
                end_column='stop'):
    """
    Pair probes whose footprints are adjacent on the target.

    Within each target the probes are taken in coordinate order and paired
    greedily: the first probe pairs with the nearest following probe whose gap
    falls in the architecture's window, both are consumed, and the search
    resumes at the next unpaired probe. Greedy pairing in coordinate order
    yields the most pairs when the gap window is a single interval, which is
    what an architecture specifies.

    Probes that overlap, or that sit further apart than the window allows, are
    left unpaired and do not appear in the result.

    Args:
        df (pandas.DataFrame): probe table carrying the target, start and end
            columns.
        architecture (SplitArchitecture): supplies the gap window.
        by (str): column identifying the target a probe was mined from.
        start_column (str): column holding each probe's start coordinate.
        end_column (str): column holding each probe's end coordinate, treated
            as exclusive.

    Returns:
        pairs (pandas.DataFrame): two rows per pair, carrying unit and role
            columns and a gap column, ordered by target and coordinate.
    """
    for column in (by, start_column, end_column):
        if column not in df.columns:
            raise InvalidInputError(f'probe table has no {column!r} column')

    if df.empty:
        out = df.copy()
        for column, dtype in ((UNIT_COLUMN, 'int64'), (ROLE_COLUMN, object),
                              ('gap', 'int64'), ('architecture', object)):
            out[column] = pd.Series(dtype=dtype)
        return out

    keep = []
    gaps = []

    for _, group in df.groupby(by, sort=True):
        ordered = group.sort_values([start_column, end_column])
        starts = ordered[start_column].to_numpy()
        ends = ordered[end_column].to_numpy()
        labels = ordered.index.to_numpy()

        i = 0
        while i < len(ordered) - 1:
            gap = int(starts[i + 1] - ends[i])
            if architecture.min_gap <= gap <= architecture.max_gap:
                keep.extend([labels[i], labels[i + 1]])
                gaps.extend([gap, gap])
                i += 2
            else:
                i += 1

    out = df.loc[keep].copy()
    out['gap'] = gaps

    # each consecutive kept row belongs to one pair, in the order they were kept
    out[UNIT_COLUMN] = [i // 2 for i in range(len(out))]
    out[ROLE_COLUMN] = [ROLES[i % 2] for i in range(len(out))]
    out['architecture'] = architecture.name

    # success
    return out.reset_index(drop=True)


def assemble_split(pairs, architecture, homology_column='probe_seq'):
    """
    Attach the architecture's domains to each member of each pair.

    Args:
        pairs (pandas.DataFrame): the output of pair_probes.
        architecture (SplitArchitecture): what each member carries.
        homology_column (str): column holding each probe's homology sequence.

    Returns:
        out (pandas.DataFrame): a copy with a sequence column holding the
            assembled oligo, and append_5p/append_3p recording what was added.
    """
    if homology_column not in pairs.columns:
        raise InvalidInputError(
            f'pair table has no {homology_column!r} column')
    if ROLE_COLUMN not in pairs.columns:
        raise InvalidInputError(
            f'pair table has no {ROLE_COLUMN!r} column; run pair_probes first')

    roles = pairs[ROLE_COLUMN]
    five = roles.map(lambda role: architecture.append_for(role, '5p'))
    three = roles.map(lambda role: architecture.append_for(role, '3p'))

    assembly = DomainAssembly(pairs.index, SPLIT_LAYOUT)
    assembly.set_domain('append_5p', five)
    assembly.set_domain('homology', pairs[homology_column])
    assembly.set_domain('append_3p', three)
    assembly.set_linker('append_5p', 'homology', architecture.linker)
    assembly.set_linker('homology', 'append_3p', architecture.linker)

    out = pairs.copy()
    out['append_5p'] = five
    out['append_3p'] = three
    out['sequence'] = assembly.assemble()

    # success
    return out


def design_split(df, architecture, by='seq_id', start_column='start',
                 end_column='stop', homology_column='probe_seq'):
    """
    Pair adjacent probes and assemble both members of every pair.

    Args:
        df (pandas.DataFrame): probe table to design from.
        architecture (SplitArchitecture): the chemistry to build.
        by (str): column identifying the target a probe was mined from.
        start_column (str): column holding each probe's start coordinate.
        end_column (str): column holding each probe's end coordinate.
        homology_column (str): column holding each probe's homology sequence.

    Returns:
        out (pandas.DataFrame): two rows per pair, assembled and unit-labelled.
    """
    pairs = pair_probes(df, architecture, by=by, start_column=start_column,
                        end_column=end_column)
    if pairs.empty:
        out = pairs.copy()
        for column in ('append_5p', 'append_3p', 'sequence'):
            out[column] = pd.Series(dtype=object)
        return out

    # success
    return assemble_split(pairs, architecture,
                          homology_column=homology_column)


def pair_summary(df, pairs):
    """
    Report how much of a probe table the pairing consumed.

    Args:
        df (pandas.DataFrame): the probe table pairing ran on.
        pairs (pandas.DataFrame): the output of pair_probes.

    Returns:
        summary (dict): probe and pair counts, and the fraction of probes that
            found a partner.
    """
    n_probes = len(df)
    n_pairs = len(pairs) // 2

    # success
    return {
        'n_probes': n_probes,
        'n_pairs': n_pairs,
        'n_paired_probes': 2 * n_pairs,
        'n_unpaired_probes': n_probes - 2 * n_pairs,
        'fraction_paired': (2 * n_pairs / n_probes) if n_probes else 0.0,
    }
