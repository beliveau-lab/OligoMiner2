"""
Domain-aware appending -- the generalization the padlock backbone forces.

# Why the existing interface cannot express a padlock

OM2's appending prototype (``probe_design/appending/appending.py``) is edge-anchored:
every function takes ``left=True/False`` and puts the appended sequence at the 5' or 3'
end of one ``sequence`` column. That covers the PaintSHOP anatomy::

    5'--[Outer]--[Inner]--[Homology]--[Inner]--[Outer]--3'

because every layer wraps symmetrically around a single homology core.

A padlock is not that shape. Its backbone sits **between** the two homology arms::

    5'--[arm_5p]--[backbone]--[arm_3p]--3'

There is no value of ``left`` that produces this, because the insertion point is
interior. So the boolean generalizes to an **ordered list of named domains**, and
appending targets a slot rather than an edge. The old behaviour is the two-domain
special case: ``left=True`` is "insert before the homology domain".

# The linker is per-join, not global

``appending.config.LINKER = "TTT"`` is inserted at every join by ``_join()``. For a
padlock that is wrong in a specific and consequential way: the arms must **abut on the
target** for ligase to seal the circle, so a linker must never appear at the ligation
junction. Linkers therefore attach to a join, not to the assembly, and a join can carry
the empty linker.

# What is preserved

The four assignment schemes (``same`` / ``unique`` / ``multiple`` / ``custom``) and the
``entries`` tracking contract are unchanged -- this module delegates to the vendored
implementations for assignment and only owns *where* the result lands. That keeps the
PaintSHOP appending path bit-identical while making padlocks expressible.
"""

import pandas as pd

from oligominer.probe_design.appending.appending import (
    append_sequences, build_appending_table,
)
from oligominer.probe_design.appending.config import LINKER
from oligominer.utils.exceptions import InvalidInputError

# the domain layout a padlock record produces; homology is split, backbone is interior
PADLOCK_LAYOUT = ["arm_5p", "backbone", "arm_3p"]

# the classic single-homology layout, for probes that are not padlocks
LINEAR_LAYOUT = ["outer_5p", "inner_5p", "homology", "inner_3p", "outer_3p"]


class DomainAssembly:
    """
    An ordered set of named sequence domains, assembled into one oligo.

    Each domain holds a per-row sequence (a pandas Series) or is empty. Joins between
    adjacent domains carry their own linker, so a ligation junction can be linker-free
    while other joins are not.

    Attributes:
        layout (list): domain names in 5'->3' order.
        domains (dict): domain name -> pandas.Series of sequences.
        linkers (dict): (left_name, right_name) -> linker string.
        entries (dict): label -> pandas.Series of appending tracking strings.
    """

    def __init__(self, index, layout):
        """
        Args:
            index (pandas.Index): the row index every domain is aligned to.
            layout (list): domain names in 5'->3' order.
        """
        if len(set(layout)) != len(layout):
            raise InvalidInputError(f"duplicate domain names in layout: {layout}")

        self.index = index
        self.layout = list(layout)
        # object dtype throughout: pandas 3 str-dtype and object do not concatenate,
        # and domain values arrive from callers under both
        self.domains = {name: pd.Series("", index=index, dtype=object)
                        for name in self.layout}
        self.linkers = {}
        self.entries = {}

    def set_domain(self, name, values):
        """
        Fill a domain with per-row sequences.

        Args:
            name (str): a domain name from the layout.
            values (pandas.Series or str): sequences, or one sequence for every row.

        Returns:
            self (DomainAssembly): for chaining.
        """
        self._require(name)
        if isinstance(values, str):
            values = pd.Series(values, index=self.index, dtype=object)
        self.domains[name] = (pd.Series(values).astype(object)
                              .reindex(self.index).fillna(""))

        # success
        return self

    def set_linker(self, left, right, linker):
        """
        Set the linker inserted between two adjacent domains.

        Args:
            left (str): the 5' domain name.
            right (str): the 3' domain name.
            linker (str): the linker sequence; "" for a seamless join.

        Returns:
            self (DomainAssembly): for chaining.
        """
        self._require(left)
        self._require(right)
        if self.layout.index(right) != self.layout.index(left) + 1:
            raise InvalidInputError(
                f"{left!r} and {right!r} are not adjacent in {self.layout}")
        self.linkers[(left, right)] = linker

        # success
        return self

    def append_into(self, name, sequences, scheme, label=None, rc=False,
                    target_column=None, n_per_target=None, ranges=None,
                    probes=None):
        """
        Assign sequences into a domain using one of the four appending schemes.

        Delegates assignment to the vendored ``append_sequences`` so scheme behaviour
        and the ``entries`` tracking strings stay bit-identical to the package; only the
        destination differs -- a named domain instead of an edge.

        Args:
            name (str): the domain to fill.
            sequences (pandas.DataFrame): appending sequences with ``id`` and ``seq``.
            scheme (str): ``"same"``, ``"unique"``, ``"multiple"`` or ``"custom"``.
            label (str or None): key for this step in the appending table. Defaults to
                the domain name.
            rc (bool): reverse-complement the sequences before assigning.
            target_column (str, optional): grouping column for unique/multiple.
            n_per_target (int, optional): sequences per target for multiple.
            ranges (list, optional): range strings for custom.
            probes (pandas.DataFrame, optional): the frame carrying target_column.
                Required for the unique and multiple schemes.

        Returns:
            self (DomainAssembly): for chaining.
        """
        self._require(name)

        # the vendored functions assign onto a 'sequence' column; give them an empty one
        # so what comes back IS the assigned sequence rather than a concatenation
        frame = pd.DataFrame(index=self.index)
        frame["sequence"] = ""
        if probes is not None and target_column is not None:
            frame[target_column] = probes[target_column].reindex(self.index)

        result, entries = append_sequences(
            frame, sequences, scheme,
            target_column=target_column, n_per_target=n_per_target,
            ranges=ranges, left=True, rc=rc, linker="",
        )

        self.domains[name] = result["sequence"].astype(object)
        self.entries[label or name] = entries

        # success
        return self

    def assemble(self):
        """
        Concatenate the domains in layout order, inserting each join's linker.

        A join's linker is skipped when either neighbour is empty on that row, so an
        unused domain does not leave a dangling linker in the product.

        When a domain is empty its neighbours become adjacent, and the join uses the
        linker of the first boundary crossed -- the one immediately after the last
        non-empty domain. Only adjacent boundaries can be configured, so this keeps
        every linker actually used one the caller is able to set.

        Returns:
            seqs (pandas.Series): the assembled oligo per row.
        """
        seqs = pd.Series("", index=self.index, dtype=object)
        previous_index = None

        for position, name in enumerate(self.layout):
            current = self.domains[name].astype(object)
            if previous_index is None:
                seqs = current.copy()
                if (current.str.len() > 0).any():
                    previous_index = position
                continue

            # the boundary immediately after the last non-empty domain
            boundary = (self.layout[previous_index], self.layout[previous_index + 1])
            linker = self.linkers.get(boundary, LINKER)

            both = (seqs.str.len() > 0) & (current.str.len() > 0)
            joiner = pd.Series("", index=self.index, dtype=object)
            joiner[both] = linker
            seqs = seqs + joiner + current

            # an empty domain must not become the left neighbour of the next join,
            # or its linker would be attributed to the wrong boundary
            if (current.str.len() > 0).any():
                previous_index = position

        # success
        return seqs

    def table(self):
        """
        Build the appending table recording what went into each domain.

        Returns:
            table (pandas.DataFrame): one row per oligo, one column per appending step.
        """
        # success
        return build_appending_table(pd.DataFrame(index=self.index), self.entries)

    def _require(self, name):
        """Raise if name is not in the layout."""
        if name not in self.domains:
            raise InvalidInputError(
                f"unknown domain {name!r}; layout is {self.layout}")

        # success
        return True


def assemble_padlock(padlock_df, backbone, linker_5p="", linker_3p="",
                     backbone_id="backbone"):
    """
    Insert a backbone between the two homology arms of every padlock.

    This is the operation the edge-anchored interface cannot express. Both joins default
    to a **seamless** linker: the arm ends must abut their target for ligase to seal the
    circle, and an inserted TTT would sit across the ligation junction.

    Args:
        padlock_df (pandas.DataFrame): padlock records with ``arm_5p`` and ``arm_3p``,
            as produced by ``padlock.mine_padlock_sequence``.
        backbone (str or pandas.Series): the backbone sequence, one for all rows or
            one per row.
        linker_5p (str): linker between arm_5p and the backbone. Default seamless.
        linker_3p (str): linker between the backbone and arm_3p. Default seamless.
        backbone_id (str): identifier recorded in the appending table.

    Returns:
        result (pandas.DataFrame): a copy of padlock_df with ``full_oligo`` added.
        table (pandas.DataFrame): the appending table.
    """
    asm = DomainAssembly(padlock_df.index, PADLOCK_LAYOUT)
    asm.set_domain("arm_5p", padlock_df["arm_5p"])
    asm.set_domain("arm_3p", padlock_df["arm_3p"])
    asm.set_domain("backbone", backbone)
    asm.set_linker("arm_5p", "backbone", linker_5p)
    asm.set_linker("backbone", "arm_3p", linker_3p)

    bb_repr = backbone if isinstance(backbone, str) else "per-row"
    asm.entries["backbone"] = pd.Series(f"{backbone_id}_{bb_repr}",
                                        index=padlock_df.index, dtype=object)

    result = padlock_df.copy()
    result["full_oligo"] = asm.assemble()

    # success
    return result, asm.table()


def check_backbone_placement(row):
    """
    Assert the backbone landed between the arms and did not disturb them.

    Args:
        row (dict or pandas.Series): a row from ``assemble_padlock``'s result, plus the
            backbone sequence under key ``_backbone``.

    Returns:
        ok (bool): True if the assembled oligo starts with arm_5p and ends with arm_3p.

    Raises:
        AssertionError: naming the invariant that failed.
    """
    oligo = row["full_oligo"]
    assert oligo.startswith(row["arm_5p"]), "5' arm is not at the 5' end of the oligo"
    assert oligo.endswith(row["arm_3p"]), "3' arm is not at the 3' end of the oligo"
    assert len(oligo) > len(row["arm_5p"]) + len(row["arm_3p"]), \
        "nothing was inserted between the arms"

    # success
    return True
