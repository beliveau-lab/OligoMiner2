"""
Bundled PaintSHOP appending sequence data and convenience loaders.

Provides helper functions to load appending datasets (primers, bridges,
SABER concatemers, MERFISH sequences) from package data via
importlib.resources, so that scripts and notebooks can access them
without constructing file paths manually.
"""

from importlib.resources import files

import pandas as pd


_PKG = "oligominer.data.appending.paintshop"


def load_appending_data(name):
    """
    Load a bundled PaintSHOP appending TSV file by name.

    Args:
        name (str): the base name of the data file (without extension).
            For example, ``"ps_bridges"`` loads ``ps_bridges.tsv.gz``.

    Returns:
        df (pandas.DataFrame): two-column DataFrame with ``id`` and
            ``seq`` columns.
    """
    filename = f"{name}.tsv.gz"
    ref = files(_PKG).joinpath(filename)
    df = pd.read_csv(str(ref), sep="\t", compression="gzip")

    # success
    return df


# ------------------------------------------------------------------
# named convenience loaders
# ------------------------------------------------------------------

def load_bridges():
    """
    Load the PaintSHOP bridge set (800 orthogonal bridges).

    Returns:
        df (pandas.DataFrame): bridge sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("ps_bridges")


def load_outer_forward():
    """
    Load the PaintSHOP outer forward primer set (10 primers).

    Returns:
        df (pandas.DataFrame): primer sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("ps_of")


def load_outer_reverse():
    """
    Load the PaintSHOP outer reverse primer set (10 primers).

    Returns:
        df (pandas.DataFrame): primer sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("ps_or")


def load_inner_forward():
    """
    Load the PaintSHOP inner forward primer set (74 primers).

    Returns:
        df (pandas.DataFrame): primer sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("ps_if")


def load_inner_reverse():
    """
    Load the PaintSHOP inner reverse primer set (74 primers).

    Returns:
        df (pandas.DataFrame): primer sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("ps_ir")


def load_saber_1x():
    """
    Load the SABER 1x concatemer set (50 sequences).

    Returns:
        df (pandas.DataFrame): SABER sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("saber_1x")


def load_saber_2x():
    """
    Load the SABER 2x concatemer set (50 sequences).

    Returns:
        df (pandas.DataFrame): SABER sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("saber_2x")


def load_merfish_bridges():
    """
    Load the MERFISH bridge set (16 bridges).

    Returns:
        df (pandas.DataFrame): bridge sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("merfish_bridges")


def load_merfish_primers():
    """
    Load the MERFISH primer set (318 primers).

    Returns:
        df (pandas.DataFrame): primer sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("merfish_primers")


def load_kishi_bridges():
    """
    Load the Kishi et al. 2019 bridge set (84 bridges).

    Returns:
        df (pandas.DataFrame): bridge sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("Kishi2019_bridges")


def load_mateo_bridges():
    """
    Load the Mateo et al. 2019 bridge set (199 bridges).

    Returns:
        df (pandas.DataFrame): bridge sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("Mateo2019_bridges")


def load_xia_bridges():
    """
    Load the Xia et al. 2019 bridge set (70 bridges).

    Returns:
        df (pandas.DataFrame): bridge sequences with ``id`` and ``seq``
            columns.
    """
    return load_appending_data("Xia2019_bridges")
