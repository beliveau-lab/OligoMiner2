"""Tests that the public surface stays importable and composable.

The primitives are meant to be called individually and composed, including by
callers outside this repo. A capability that exists but can only be reached
through a private module path is not usable, and a rename that breaks an import
should fail here rather than in a downstream caller.
"""

import importlib

import pytest

# every verb a caller is expected to reach for, and where it lives
PUBLIC_SURFACE = {
    'oligominer': [
        'mine_sequence', 'mine_fasta', 'write_probes', 'ProbeSet', '__version__',
    ],
    'oligominer.utils.cores': ['resolve_cores', 'describe'],
    'oligominer.thermodynamics': [
        'formamide_correction', 'effective_hyb_temperature',
    ],
    'oligominer.thermodynamics.nupack': [
        'calc_pdup', 'calc_pdup_many', 'calc_pdup_one_to_many', 'add_pdup_batch',
    ],
    'oligominer.specificity.kmers': [
        'build_index', 'max_kmer', 'resolve_backend', 'KmerIndex',
    ],
    'oligominer.specificity.alignment': [
        'bowtie_build', 'build_bowtie2_cmd', 'align_to_bed', 'duplex',
    ],
    'oligominer.specificity.duplex_stability.frames': [
        'build_duplex_frame', 'build_aln', 'expand_cigar',
    ],
    'oligominer.specificity.triage': ['triage', 'triage_summary'],
    'oligominer.models': [
        'load', 'load_all', 'spec', 'available', 'REGISTRY',
        'build_flat_corpus', 'retrain',
    ],
    'oligominer.probe_design': [
        'ProbeSet', 'exclude_intervals', 'mine_padlock_sequence', 'padlocks_to_df',
        'DomainAssembly', 'assemble_padlock', 'assign_units', 'filter_units',
        'nominate', 'verify_against_screen',
    ],
    'oligominer.probe_design.schema': [
        'new_manifest', 'record_stage', 'attrition_summary', 'SCHEMA_VERSION',
    ],
    'oligominer.bioinformatics.file_io': [
        'load_fasta', 'write_fasta', 'load_gtf', 'load_bam_file',
    ],
    'oligominer.bioinformatics.transcriptome': [
        'mine_exons', 'flatten_isoforms', 'get_spliced_seq',
        'junction_windows', 'discriminating_regions', 'build_transcriptome',
    ],
}


@pytest.mark.parametrize('module_name', sorted(PUBLIC_SURFACE))
def test_module_imports(module_name):
    importlib.import_module(module_name)


@pytest.mark.parametrize(
    'module_name,attribute',
    [(m, a) for m, names in sorted(PUBLIC_SURFACE.items()) for a in names])
def test_attribute_is_exported(module_name, attribute):
    module = importlib.import_module(module_name)
    assert hasattr(module, attribute), (
        f'{module_name}.{attribute} is not exported; a capability reachable only '
        f'through a private path is not usable')


class TestOptionalDependencies:
    """Heavyweight dependencies must not be required to import the package."""

    def test_torch_is_not_imported_at_module_scope(self):
        import oligominer.specificity.duplex_stability.bilstm_arch as arch
        source = open(arch.__file__).read()
        for line in source.splitlines():
            if line.startswith('import torch') or line.startswith('from torch'):
                pytest.fail('torch is imported at module scope; it is an optional extra')

    def test_the_tree_models_load_without_torch(self):
        from oligominer.models import load
        assert load('physics-xgb').outputs_pdup


class TestComposition:
    """The primitives compose in the order a probe design actually runs."""

    def test_mine_then_wrap_then_export(self, example_fasta_path, tmp_path):
        from oligominer import ProbeSet, mine_fasta

        probes = mine_fasta(example_fasta_path, min_tm=42, max_tm=47)
        probe_set = ProbeSet(probes)

        csv_path = tmp_path / 'probes.csv'
        json_path = tmp_path / 'probes.json'
        probe_set.to_csv(csv_path)
        probe_set.to_json(json_path)

        assert csv_path.exists()
        assert len(ProbeSet.from_json(json_path)) == len(probe_set)

    def test_mine_then_screen_kmers(self, example_fasta_path, tmp_path):
        from oligominer import mine_fasta
        from oligominer.specificity.kmers import build_index, max_kmer

        index = tmp_path / 'ref.npz'
        build_index(str(example_fasta_path), str(index), k=12, backend='numpy')

        probes = mine_fasta(example_fasta_path, min_tm=42, max_tm=47)[:50]
        counts = max_kmer(index, [p[3] for p in probes], k=12, backend='numpy')
        assert len(counts) == len(probes)

    def test_padlock_then_assemble_then_unit(self):
        import pandas as pd

        from oligominer.probe_design import (
            assemble_padlock, assign_units, mine_padlock_sequence, padlocks_to_df,
        )

        import random
        random.seed(12)
        target = ''.join(random.choice('ACGT') for _ in range(4000))

        padlocks = padlocks_to_df(mine_padlock_sequence(target, seq_id='chr1'))
        assembled, _ = assemble_padlock(padlocks, backbone='TTTGCTAGCTAGCTAGCAAA')
        united = assign_units(assembled, by=['seq_id', 'start'])

        assert 'full_oligo' in assembled.columns
        assert united['unit_id'].nunique() == len(assembled)
