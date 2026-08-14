"""Tests for the re-training API.

The rules that matter are not about fit quality: a registered artifact must never
be overwritten, every output must carry a card identifying its corpus, and fitting
a model away from the regime it ships at must say so.
"""

import json
import random
import warnings

import numpy as np
import pandas as pd
import pytest

from oligominer.models import artifact_path, build_flat_corpus, retrain, write_card
from oligominer.specificity.duplex_stability.frames import build_duplex_frame


def cigar_from(probe, target):
    """Build an xeq CIGAR from two equal-length sequences."""
    ops = ['=' if a == b else 'X' for a, b in zip(probe, target)]
    out, run, current = [], 0, ops[0]
    for op in ops:
        if op == current:
            run += 1
        else:
            out.append(f'{run}{current}')
            current, run = op, 1
    out.append(f'{run}{current}')
    return ''.join(out)


@pytest.fixture
def corpus():
    """A duplex frame with a pdup label spread across every decile."""
    random.seed(8)
    rows = []
    for i in range(400):
        probe = ''.join(random.choice('ACGT') for _ in range(32))
        target = list(probe)
        for j in random.sample(range(32), random.randint(0, 10)):
            target[j] = random.choice([b for b in 'ACGT' if b != probe[j]])
        target = ''.join(target)
        rows.append({'probe_seq': probe, 'derived_seq': target,
                     'align_cigar': cigar_from(probe, target),
                     'align_score': -5, 'pdup': (i % 100) / 100.0})

    # success
    return build_duplex_frame(pd.DataFrame(rows), celsius=69.5, sodium=0.39)


class TestFlatCorpus:

    def test_every_decile_gets_the_same_number_of_rows(self, corpus, tmp_path):
        path = build_flat_corpus(corpus, tmp_path / 'c.parquet')
        cut = pd.read_parquet(path)
        deciles = np.clip((cut['pdup'] * 10).astype(int), 0, 9)
        assert deciles.value_counts().nunique() == 1

    def test_a_sidecar_records_the_corpus_checksum(self, corpus, tmp_path):
        path = build_flat_corpus(corpus, tmp_path / 'c.parquet')
        sidecar = json.loads(path.with_suffix('.corpus.json').read_text())
        assert len(sidecar['sha256']) == 64
        assert sidecar['flat_in_pdup'] is True

    def test_the_seed_is_recorded(self, corpus, tmp_path):
        path = build_flat_corpus(corpus, tmp_path / 'c.parquet', seed=7)
        assert json.loads(path.with_suffix('.corpus.json').read_text())['seed'] == 7

    def test_a_missing_label_column_raises(self, corpus, tmp_path):
        with pytest.raises(KeyError, match='pdup'):
            build_flat_corpus(corpus.drop(columns='pdup'), tmp_path / 'c.parquet')

    def test_sampling_is_reproducible(self, corpus, tmp_path):
        a = pd.read_parquet(build_flat_corpus(corpus, tmp_path / 'a.parquet', seed=3))
        b = pd.read_parquet(build_flat_corpus(corpus, tmp_path / 'b.parquet', seed=3))
        assert list(a['probe_seq']) == list(b['probe_seq'])


class TestArtifactProtection:

    def test_writing_over_a_registered_artifact_is_refused(self, corpus):
        registered = str(artifact_path('physics-xgb'))
        with pytest.raises(ValueError, match='registered artifact'):
            retrain('physics-xgb', corpus, registered)

    def test_a_baseline_model_cannot_be_refitted_here(self, corpus, tmp_path):
        with pytest.raises(ValueError, match='kept for comparison'):
            retrain('ps-xgb', corpus, tmp_path / 'x.pkl')


class TestRegimeWarnings:

    def test_pooling_temperatures_for_the_bilstm_warns(self, corpus, tmp_path):
        with pytest.warns(RuntimeWarning, match='no condition channel'):
            with pytest.raises(ValueError, match='torch'):
                retrain('duplex-BiLSTM', corpus, tmp_path / 'b.pt', celsius=None)

    def test_fitting_the_bilstm_off_its_shipping_condition_warns(self, corpus, tmp_path):
        with pytest.warns(RuntimeWarning, match='ships fitted at'):
            with pytest.raises(ValueError, match='torch'):
                retrain('duplex-BiLSTM', corpus, tmp_path / 'b.pt', celsius=47.0)

    def test_the_bilstm_is_not_fitted_by_this_api(self, corpus, tmp_path):
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            with pytest.raises(ValueError, match='torch'):
                retrain('duplex-BiLSTM', corpus, tmp_path / 'b.pt', celsius=69.5)


class TestFitting:

    def test_a_tree_model_fits_and_writes_a_card(self, corpus, tmp_path):
        out = tmp_path / 'refit.pkl'
        info = retrain('physics-xgb', corpus, out, celsius=69.5)

        assert out.exists()
        card = json.loads(tmp_path.joinpath('refit.card.json').read_text())
        assert card['based_on'] == 'physics-xgb'
        assert card['role'] == 'RETRAINED'
        assert card['celsius'] == 69.5
        assert info['n_features'] == 103

    def test_the_card_records_the_corpus(self, corpus, tmp_path):
        built = build_flat_corpus(corpus, tmp_path / 'c.parquet')
        retrain('physics-xgb', str(built), tmp_path / 'refit.pkl', celsius=69.5)
        card = json.loads(tmp_path.joinpath('refit.card.json').read_text())
        assert card['corpus'] == str(built)
        assert len(card['corpus_sha256']) == 64

    def test_the_refitted_model_loads_and_predicts(self, corpus, tmp_path):
        import pickle

        out = tmp_path / 'refit.pkl'
        retrain('physics-xgb', corpus, out, celsius=69.5)

        with out.open('rb') as handle:
            bundle = pickle.load(handle)
        assert 'model' in bundle
        assert len(bundle['features']) == 103

    def test_an_in_memory_corpus_records_no_checksum(self, corpus, tmp_path):
        """A fit from an in-memory frame cannot be reproduced from its card alone."""
        retrain('physics-xgb', corpus, tmp_path / 'refit.pkl', celsius=69.5)
        card = json.loads(tmp_path.joinpath('refit.card.json').read_text())
        assert card['corpus_sha256'] is None
        assert card['corpus'] == 'in-memory frame'

    def test_write_card_can_be_called_directly(self, tmp_path):
        path = write_card('physics-xgb', tmp_path / 'm.pkl', tmp_path / 'c.parquet',
                          celsius=69.5, seed=0)
        assert json.loads(path.read_text())['based_on'] == 'physics-xgb'
