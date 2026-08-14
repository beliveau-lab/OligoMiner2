"""
# Targeting units

A targeting unit is the set of oligos that must all be present for one target
site to produce signal.

For a conventional probe the unit is one oligo, and every filter can act on rows
independently. Split architectures break that: a split-FISH pair or an HCR 3.0
split-initiator pair only produces signal when both members bind, so dropping one
member does not weaken the site, it silently removes it while leaving the other
member in the order. The surviving oligo then costs synthesis and contributes
off-target binding while contributing no signal.

Filtering therefore has to act on units rather than rows. `filter_units` keeps a
unit only when every member passes, and `drop_incomplete_units` removes units
that lost members to an earlier row-wise stage.
"""

import pandas as pd

# column naming the unit a row belongs to
UNIT_COLUMN = 'unit_id'

# column naming a row's role within its unit
ROLE_COLUMN = 'unit_role'


def assign_units(df, by, roles=None, unit_column=UNIT_COLUMN,
                 role_column=ROLE_COLUMN):
    """
    Label each row with the unit it belongs to.

    Args:
        df (pandas.DataFrame): the oligo table.
        by (str or list): column or columns identifying a target site. Rows
            sharing these values form one unit.
        roles (list, optional): role names assigned in row order within each
            unit, e.g. ['left', 'right']. Recorded in role_column when given.
        unit_column (str): column to write the unit identifier into.
        role_column (str): column to write the role into.

    Returns:
        out (pandas.DataFrame): a copy carrying the unit and role columns.
    """
    out = df.copy()
    keys = [by] if isinstance(by, str) else list(by)

    out[unit_column] = out.groupby(keys, sort=False).ngroup()

    if roles is not None:
        position = out.groupby(unit_column, sort=False).cumcount()
        out[role_column] = [roles[i % len(roles)] for i in position]

    # success
    return out


def unit_sizes(df, unit_column=UNIT_COLUMN):
    """
    Return how many oligos each unit currently has.

    Args:
        df (pandas.DataFrame): the oligo table carrying unit_column.
        unit_column (str): the unit identifier column.

    Returns:
        sizes (pandas.Series): unit identifier mapped to its member count.
    """
    # success
    return df.groupby(unit_column, sort=False).size()


def filter_units(df, keep, unit_column=UNIT_COLUMN):
    """
    Keep only the units whose every member passes.

    A row-wise mask is promoted to a unit-wise decision: a unit survives when all
    of its members are marked keep, and is removed entirely otherwise.

    Args:
        df (pandas.DataFrame): the oligo table carrying unit_column.
        keep (array-like): boolean per row, True for rows that pass.
        unit_column (str): the unit identifier column.

    Returns:
        kept (pandas.DataFrame): the members of surviving units.
        n_units_dropped (int): how many units were removed.
    """
    keep = pd.Series(list(keep), index=df.index)
    all_pass = keep.groupby(df[unit_column]).transform('all')

    kept = df[all_pass].copy()
    n_units_dropped = (df[unit_column].nunique() - kept[unit_column].nunique()
                       if len(kept) else df[unit_column].nunique())

    # success
    return kept, int(n_units_dropped)


def drop_incomplete_units(df, expected_size=None, unit_column=UNIT_COLUMN):
    """
    Remove units that do not have all of their members.

    Use after a stage that filtered rows without unit awareness, to remove the
    orphaned members it left behind.

    Args:
        df (pandas.DataFrame): the oligo table carrying unit_column.
        expected_size (int, optional): the number of oligos a complete unit has.
            The largest size present is used when None, since a unit can only
            lose members to an earlier stage, never gain them.
        unit_column (str): the unit identifier column.

    Returns:
        kept (pandas.DataFrame): the members of complete units.
        n_units_dropped (int): how many units were removed.
    """
    if len(df) == 0:
        return df.copy(), 0

    sizes = unit_sizes(df, unit_column=unit_column)
    if expected_size is None:
        expected_size = int(sizes.max())

    complete = sizes[sizes == expected_size].index
    kept = df[df[unit_column].isin(complete)].copy()

    # success
    return kept, int(sizes.size - len(complete))


def units_to_orders(df, sequence_column='sequence', unit_column=UNIT_COLUMN,
                    role_column=ROLE_COLUMN):
    """
    Return one row per unit with each member's sequence in its own column.

    This is the shape an oligo order takes, where a unit's members are ordered
    together and must be tracked together.

    Args:
        df (pandas.DataFrame): the oligo table carrying unit and role columns.
        sequence_column (str): column holding each oligo's sequence.
        unit_column (str): the unit identifier column.
        role_column (str): the role column.

    Returns:
        orders (pandas.DataFrame): one row per unit, one column per role.
    """
    orders = df.pivot(index=unit_column, columns=role_column,
                      values=sequence_column)
    orders.columns.name = None

    # success
    return orders.reset_index()
