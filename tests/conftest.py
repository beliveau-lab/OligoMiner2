"""Shared fixtures for OligoMiner2 tests."""

import numpy as np
import pandas as pd
import pytest

from oligominer import get_example_fasta_path, mine_fasta, mine_sequence, ProbeSet


# ---------------------------------------------------------------------------
# paths
# ---------------------------------------------------------------------------

@pytest.fixture
def example_fasta_path():
    """Path to the bundled example FASTA file."""
    return get_example_fasta_path()


# ---------------------------------------------------------------------------
# sequences
# ---------------------------------------------------------------------------

SHORT_SEQ = "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"

@pytest.fixture
def short_seq():
    """A short synthetic DNA sequence (~50 bp) for fast unit tests."""
    return SHORT_SEQ


# ---------------------------------------------------------------------------
# mining results
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def example_probes(tmp_path_factory):
    """Mine probes from the bundled example FASTA (session-scoped for speed)."""
    fasta_path = get_example_fasta_path()
    probes = mine_fasta(fasta_path)
    return probes


@pytest.fixture
def example_probe_set(example_probes):
    """A ProbeSet built from the bundled example FASTA mining results."""
    return ProbeSet(example_probes)


# ---------------------------------------------------------------------------
# synthetic data for appending / scoring tests
# ---------------------------------------------------------------------------

@pytest.fixture
def synthetic_probe_df():
    """A small synthetic probe DataFrame with target/refseq columns."""
    rng = np.random.default_rng(42)
    bases = list("ATCG")
    n_probes = 20
    n_targets = 4

    seqs = [''.join(rng.choice(bases, size=30)) for _ in range(n_probes)]

    target_names = [f"chr1:region_{i}" for i in range(n_targets)]
    refseq_names = [f"NM_00{i+1}" for i in range(n_targets)]

    targets, refseqs, starts = [], [], []
    for i in range(n_targets):
        for j in range(n_probes // n_targets):
            targets.append(target_names[i])
            refseqs.append(refseq_names[i])
            starts.append(i * 10000 + j * 100)

    df = pd.DataFrame({
        'seq_id': ['chr1'] * n_probes,
        'start': starts,
        'stop': [s + 30 for s in starts],
        'probe_seq': seqs,
        'tm': [44.5 + rng.normal(0, 1) for _ in range(n_probes)],
        'sequence': seqs,
        'target': targets,
        'refseq': refseqs,
    })

    return df


@pytest.fixture
def synthetic_merged_df():
    """A small synthetic merged DataFrame with on/off-target alignments."""
    rows = []
    for i in range(5):
        seqid = f"chr1:{i*100}-{i*100+30}"
        rows.append({
            'seqid': seqid,
            'seq_id': 'chr1',
            'start': i * 100,
            'align_seqid': 'chr1',
            'align_start': i * 100,
            'align_stop': i * 100 + 30,
            'duplex_pred': 0.85,
        })
        for j in range(2):
            rows.append({
                'seqid': seqid,
                'seq_id': 'chr1',
                'start': i * 100,
                'align_seqid': f'chr{j+2}',
                'align_start': 50000 + j * 200,
                'align_stop': 50000 + j * 200 + 30,
                'duplex_pred': 0.05,
            })

    return pd.DataFrame(rows)
