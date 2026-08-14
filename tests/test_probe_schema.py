"""Tests for the probe set schema, provenance and JSON round-trip.

Covers manifest construction and validation, the attrition record that explains
a short probe set, version handling including refusal of a future version, and
that a probe set survives a write and read unchanged.
"""

import json

import pandas as pd
import pytest

from oligominer import ProbeSet
from oligominer.probe_design import schema


class TestManifest:

    def test_new_manifest_is_at_the_current_version(self):
        assert schema.new_manifest()['schema_version'] == schema.SCHEMA_VERSION

    def test_new_manifest_validates(self):
        schema.validate_manifest(schema.new_manifest())

    def test_manifest_records_the_package_version(self):
        assert schema.new_manifest()['oligominer_version']

    def test_name_is_carried(self):
        assert schema.new_manifest('my probes')['name'] == 'my probes'


class TestStages:

    def test_a_generative_stage_records_output_only(self):
        manifest = schema.new_manifest()
        schema.record_stage(manifest, 'mine', n_out=100)
        stage = manifest['stages'][0]
        assert stage['n_out'] == 100
        assert 'n_dropped' not in stage

    def test_a_filtering_stage_records_what_it_dropped(self):
        manifest = schema.new_manifest()
        schema.record_stage(manifest, 'screen', n_in=100, n_out=70)
        assert manifest['stages'][0]['n_dropped'] == 30

    def test_parameters_are_recorded(self):
        manifest = schema.new_manifest()
        schema.record_stage(manifest, 'screen', params={'max_kmer': 5})
        assert manifest['stages'][0]['params'] == {'max_kmer': 5}

    def test_a_stage_cannot_emit_more_than_it_received(self):
        manifest = schema.new_manifest()
        with pytest.raises(schema.SchemaError, match='cannot emit more'):
            schema.record_stage(manifest, 'bogus', n_in=10, n_out=99)

    def test_stages_accumulate_in_order(self):
        manifest = schema.new_manifest()
        schema.record_stage(manifest, 'mine', n_out=100)
        schema.record_stage(manifest, 'screen', n_in=100, n_out=70)
        assert [s['stage'] for s in manifest['stages']] == ['mine', 'screen']


class TestAttrition:

    def _manifest(self):
        manifest = schema.new_manifest()
        schema.record_stage(manifest, 'mine', n_out=1000)
        schema.record_stage(manifest, 'kmer', n_in=1000, n_out=780)
        schema.record_stage(manifest, 'specificity', n_in=780, n_out=402)
        return manifest

    def test_start_and_end_counts(self):
        summary = schema.attrition_summary(self._manifest())
        assert summary['n_in'] == 1000
        assert summary['n_out'] == 402

    def test_losses_are_attributed_per_stage(self):
        summary = schema.attrition_summary(self._manifest())
        assert summary['by_stage'] == {'kmer': 220, 'specificity': 378}

    def test_the_generative_stage_is_not_counted_as_a_loss(self):
        assert 'mine' not in schema.attrition_summary(self._manifest())['by_stage']

    def test_losses_sum_to_the_total_drop(self):
        summary = schema.attrition_summary(self._manifest())
        assert sum(summary['by_stage'].values()) == summary['n_in'] - summary['n_out']

    def test_an_empty_manifest_has_no_attrition(self):
        assert schema.attrition_summary(schema.new_manifest()) == {}


class TestVersioning:

    def test_a_manifest_without_a_version_is_refused(self):
        with pytest.raises(schema.SchemaError, match='no schema_version'):
            schema.validate_manifest({'created': 'x', 'n_probes': 0,
                                      'columns': [], 'stages': []})

    def test_a_future_version_is_refused_rather_than_guessed_at(self):
        manifest = schema.new_manifest()
        manifest['schema_version'] = schema.SCHEMA_VERSION + 99
        with pytest.raises(schema.SchemaError, match='cannot be read'):
            schema.validate_manifest(manifest)

    def test_a_missing_required_field_is_refused(self):
        manifest = schema.new_manifest()
        del manifest['stages']
        with pytest.raises(schema.SchemaError, match='stages'):
            schema.validate_manifest(manifest)

    def test_json_round_trip_preserves_the_manifest(self):
        manifest = schema.new_manifest('demo')
        schema.record_stage(manifest, 'mine', n_out=10)
        assert schema.loads(schema.dumps(manifest)) == manifest

    def test_upgrade_is_a_no_op_at_the_current_version(self):
        manifest = schema.new_manifest()
        assert schema.upgrade_manifest(manifest) is manifest


class TestProbeSetProvenance:

    @pytest.fixture
    def probe_set(self, example_probes):
        return ProbeSet(example_probes)

    def test_construction_is_recorded(self, probe_set):
        assert probe_set.manifest['stages'][0]['stage'] == 'create'

    def test_manifest_tracks_the_probe_count(self, probe_set):
        assert probe_set.manifest['n_probes'] == len(probe_set)

    def test_manifest_tracks_the_columns(self, probe_set):
        assert probe_set.manifest['columns'] == list(probe_set.df.columns)

    def test_describe_mentions_every_stage(self, probe_set):
        probe_set._record('screen', n_in=len(probe_set), n_out=1)
        text = probe_set.describe()
        assert 'create' in text
        assert 'screen' in text

    def test_mining_records_its_parameters(self, example_fasta_path):
        probe_set = ProbeSet.from_fasta(example_fasta_path, min_tm=42, max_tm=47)
        mine = [s for s in probe_set.manifest['stages'] if s['stage'] == 'mine'][0]
        assert mine['params']['min_tm'] == 42
        assert mine['params']['max_tm'] == 47

    def test_mining_records_the_source_fasta(self, example_fasta_path):
        probe_set = ProbeSet.from_fasta(example_fasta_path)
        assert probe_set.manifest['target']['source_fasta'] == str(example_fasta_path)


class TestJsonIo:

    @pytest.fixture
    def probe_set(self, example_probes):
        probe_set = ProbeSet(example_probes)
        probe_set._record('screen', params={'max_kmer': 5},
                          n_in=len(probe_set), n_out=len(probe_set))
        return probe_set

    def test_round_trip_preserves_the_probes(self, probe_set, tmp_path):
        path = tmp_path / 'ps.json'
        probe_set.to_json(path)
        back = ProbeSet.from_json(path)
        pd.testing.assert_frame_equal(
            back.df.reset_index(drop=True), probe_set.df.reset_index(drop=True),
            check_dtype=False)

    def test_round_trip_preserves_the_stages(self, probe_set, tmp_path):
        path = tmp_path / 'ps.json'
        probe_set.to_json(path)
        back = ProbeSet.from_json(path)
        assert ([s['stage'] for s in back.manifest['stages']]
                == [s['stage'] for s in probe_set.manifest['stages']])

    def test_the_document_is_valid_json_with_a_version(self, probe_set, tmp_path):
        path = tmp_path / 'ps.json'
        probe_set.to_json(path)
        document = json.loads(path.read_text())
        assert document['schema_version'] == schema.SCHEMA_VERSION
        assert len(document['probes']) == len(probe_set)

    def test_a_future_version_on_disk_is_refused(self, probe_set, tmp_path):
        path = tmp_path / 'ps.json'
        probe_set.to_json(path)
        document = json.loads(path.read_text())
        document['schema_version'] = schema.SCHEMA_VERSION + 99
        path.write_text(json.dumps(document))
        with pytest.raises(schema.SchemaError):
            ProbeSet.from_json(path)
