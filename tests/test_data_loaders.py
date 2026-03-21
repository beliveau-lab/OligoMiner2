"""Tests for bundled data loaders.

Covers: test_data loaders (example FASTA, test probes) and all
PaintSHOP appending data loaders.
"""

import os

import pandas as pd
import pytest

from oligominer import get_example_fasta_path, load_example_fasta, ProbeSet
from oligominer.data.test_data import load_test_probes
from oligominer.data.appending import (
    load_appending_data,
    load_bridges,
    load_outer_forward,
    load_outer_reverse,
    load_inner_forward,
    load_inner_reverse,
    load_saber_1x,
    load_saber_2x,
    load_merfish_bridges,
    load_merfish_primers,
    load_kishi_bridges,
    load_mateo_bridges,
    load_xia_bridges,
)


# ---------------------------------------------------------------------------
# test data loaders
# ---------------------------------------------------------------------------

class TestExampleFasta:

    def test_path_exists(self):
        path = get_example_fasta_path()
        assert os.path.exists(path)
        assert path.endswith('.fa')

    def test_load_example_fasta(self):
        fasta = load_example_fasta()
        assert len(fasta.keys()) > 0


class TestLoadTestProbes:

    def test_returns_probe_set(self):
        ps = load_test_probes()
        assert isinstance(ps, ProbeSet)
        assert len(ps) > 0

    def test_has_expected_columns(self):
        ps = load_test_probes()
        assert 'probe_seq' in ps.df.columns
        assert 'tm' in ps.df.columns
        assert 'seqid' in ps.df.columns


# ---------------------------------------------------------------------------
# appending data loaders
# ---------------------------------------------------------------------------

APPENDING_LOADERS = [
    ("bridges", load_bridges),
    ("outer_forward", load_outer_forward),
    ("outer_reverse", load_outer_reverse),
    ("inner_forward", load_inner_forward),
    ("inner_reverse", load_inner_reverse),
    ("saber_1x", load_saber_1x),
    ("saber_2x", load_saber_2x),
    ("merfish_bridges", load_merfish_bridges),
    ("merfish_primers", load_merfish_primers),
    ("kishi_bridges", load_kishi_bridges),
    ("mateo_bridges", load_mateo_bridges),
    ("xia_bridges", load_xia_bridges),
]


@pytest.mark.parametrize("name,loader", APPENDING_LOADERS, ids=[n for n, _ in APPENDING_LOADERS])
def test_appending_loader(name, loader):
    """Each appending data loader should return a non-empty DataFrame with id and seq columns."""
    df = loader()
    assert isinstance(df, pd.DataFrame)
    assert len(df) > 0
    assert 'id' in df.columns
    assert 'seq' in df.columns
    # sequences should be non-empty strings
    assert df['seq'].str.len().min() > 0
