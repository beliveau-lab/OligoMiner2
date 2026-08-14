"""
# Model loading and scoring

Loads a registered model and gives every model the same ``predict`` interface, so
a caller can iterate over the zoo without knowing what kind each entry is.

Four properties of these artifacts change the numbers if they are handled wrongly,
and none of them raise when they are:

- the XGBoost pickles are dicts keyed ``'model'``, not bare estimators
- those models are fitted with a logit link, so the inverse sigmoid is required
  before any value is read as a probability; skipping it leaves the ranking intact
  and every calibration and threshold number wrong
- the boosters carry no feature names, so names come from the encoder
- ``LinearDiscriminantAnalysis.predict`` returns hard class labels, so ``om1-lda``
  is scored through ``decision_function``, which is the ordering it provides

``duplex-BiLSTM`` is checked to have no condition channel on load. An artifact
with one is a different model and would shift every input.
"""

import json
import pickle

import numpy as np

from oligominer.utils.exceptions import MissingDependency
from .registry import artifact_path, card_path, spec

# alignment width assumed when a model card does not record one
DEFAULT_ALN_WIDTH = 42


def _read_card(name):
    """
    Read a model's card.

    Args:
        name (str): a registry key.

    Returns:
        card (dict): the card contents, or an empty dict when absent.
    """
    path = card_path(name)
    card = json.loads(path.read_text()) if path.is_file() else {}

    # success
    return card


def _encode(encoding, df, width):
    """
    Encode a duplex frame into the feature matrix a tabular model expects.

    Args:
        encoding (str): 'L4t', 'paintshop-37feat' or 'alignment-3feat'.
        df (pandas.DataFrame): a duplex frame carrying the aligned columns.
        width (int): alignment width, for width-bound encodings.

    Returns:
        X (pandas.DataFrame): the feature matrix.

    Raises:
        KeyError: on an unknown encoding.
    """
    from oligominer.specificity.duplex_stability import features_fast
    from oligominer.specificity.duplex_stability.l4t import features as l4t_features

    if encoding == 'L4t':
        X = features_fast.enc_om2_fast(df)
    elif encoding == 'paintshop-37feat':
        X = l4t_features.enc_paintshop37(df)
    elif encoding == 'alignment-3feat':
        X = _enc_alignment_3feat(df)
    else:
        raise KeyError(f'unknown encoding {encoding!r}')

    # success
    return X


def _enc_alignment_3feat(df):
    """
    Encode the three alignment features the OligoMiner1 model reads.

    Args:
        df (pandas.DataFrame): a duplex frame with probe_seq, derived_seq and
            align_score columns.

    Returns:
        X (pandas.DataFrame): probe length, probe GC percent and alignment score.
    """
    import pandas as pd

    from oligominer.utils.seq_utils import calc_gc

    probe = df['probe_seq'].astype(str)

    # success
    return pd.DataFrame({
        'probe_len': probe.str.len().to_numpy(dtype=np.float64),
        'probe_gc': np.array([calc_gc(s) for s in probe], dtype=np.float64),
        'align_score': df['align_score'].to_numpy(dtype=np.float64),
    })


class LoadedModel:
    """
    A loaded model exposing a uniform predict.

    Attributes:
        name (str): the registry key.
        family (str): 'om2' or 'baseline'.
        kind (str): 'xgb', 'bilstm' or 'lda'.
        outputs_pdup (bool): False for om1-lda, whose output is a decision score
            on an arbitrary scale rather than a probability.
        card (dict): the model card, empty when the artifact has none.
        encoding (str): the feature encoding key.
        width (int): the alignment width the model was fitted at.
        link (str or None): 'logit' when the model emits logits.
    """

    def __init__(self, name, entry, model, card=None):
        """
        Args:
            name (str): the registry key.
            entry (dict): the registry declaration.
            model: the loaded estimator.
            card (dict, optional): the model card.
        """
        self.name = name
        self.entry = entry
        self.model = model
        self.card = card or {}
        self.family = entry['family']
        self.kind = entry['kind']
        self.outputs_pdup = entry['outputs_pdup']

        # the card is authoritative where it exists, the registry is the fallback
        self.encoding = self.card.get('encoding') or entry['encoding']
        self.width = int(self.card.get('aln_width', DEFAULT_ALN_WIDTH))
        self.link = self.card.get('link', entry.get('link'))

    def __repr__(self):
        return (f'<LoadedModel {self.name} family={self.family} '
                f'encoding={self.encoding} outputs_pdup={self.outputs_pdup}>')

    def predict(self, df):
        """
        Score a duplex frame.

        Args:
            df (pandas.DataFrame): must carry the aligned columns the encoding
                reads: probe_aln, target_aln and ops, as built by
                oligominer.specificity.duplex_stability.frames.build_duplex_frame.
                Raw sequences alone are not sufficient.

        Returns:
            values (numpy.ndarray): pDup in [0, 1], or a decision score for
                om1-lda.
        """
        if self.kind == 'bilstm':
            values = self._predict_bilstm(df)
        else:
            values = self._predict_tabular(df)

        # success
        return values

    def _predict_tabular(self, df):
        """
        Encode a frame, predict, and invert the link when the model has one.

        Args:
            df (pandas.DataFrame): the duplex frame.

        Returns:
            values (numpy.ndarray): predictions.
        """
        X = _encode(self.encoding, df, self.width)
        features = X.values.astype(np.float32)

        if self.kind == 'lda':
            # decision_function is the ordering this model provides; predict would
            # collapse every row to a hard class label and destroy the ranking
            values = np.asarray(self.model.decision_function(features),
                                dtype=np.float64)
            return values

        values = self.model.predict(features).astype(np.float64)

        if self.link == 'logit':
            values = 1.0 / (1.0 + np.exp(-np.clip(values, -50, 50)))

        # success
        return np.clip(values, 0, 1) if self.outputs_pdup else values

    def _predict_bilstm(self, df):
        """
        Tokenize the aligned duplex and score it.

        Args:
            df (pandas.DataFrame): the duplex frame.

        Returns:
            values (numpy.ndarray): pDup in [0, 1].
        """
        from oligominer.specificity.duplex_stability.bilstm_arch import (
            encode_conditions,
            tokenize,
        )

        ncond = int(getattr(self.model, 'ncond', 0))
        tokens = tokenize(df, width=self.width)
        conditions = encode_conditions(df, ncond=ncond)
        values = np.asarray(self.model.predict(tokens, conditions), dtype=np.float64)

        # success
        return np.clip(values, 0, 1)

    def feature_names(self, df):
        """
        Return the encoder's column names.

        Args:
            df (pandas.DataFrame): a small duplex frame to encode.

        Returns:
            names (list): the feature names, empty for the BiLSTM which consumes
                tokens rather than a feature matrix.
        """
        if self.kind == 'bilstm':
            return []

        # success
        return list(_encode(self.encoding, df.head(2), self.width).columns)


def load(name, device=None):
    """
    Load a registered model.

    Args:
        name (str): a registry key. See registry.available().
        device (str, optional): torch device for the BiLSTM, e.g. 'cuda'.
            Defaults to CPU.

    Returns:
        model (LoadedModel): the loaded model with a uniform predict.

    Raises:
        ValueError: if the BiLSTM artifact carries a condition channel.
    """
    entry = spec(name)
    card = _read_card(name)
    path = artifact_path(name)

    if entry['kind'] == 'bilstm':
        model = _load_bilstm(path, device=device)
        if int(getattr(model, 'ncond', 0)) != 0:
            raise ValueError(
                f'{name} artifact reports ncond={model.ncond}; the registered model '
                f'has no condition channel and an artifact with one is a different '
                f'model')
    else:
        model = _load_pickle(path)

    # success
    return LoadedModel(name, entry, model, card=card)


def _load_pickle(path):
    """
    Load a pickled estimator.

    Args:
        path: the artifact path.

    Returns:
        model: the estimator, unwrapped when the pickle is a dict keyed 'model'.
    """
    with path.open('rb') as handle:
        obj = pickle.load(handle)

    if isinstance(obj, dict) and 'model' in obj:
        obj = obj['model']

    # success
    return obj


def _load_bilstm(path, device=None):
    """
    Load the BiLSTM checkpoint.

    Args:
        path: the artifact path.
        device (str, optional): torch device.

    Returns:
        model: the loaded network wrapper.

    Raises:
        MissingDependency: if torch is not installed.
    """
    from oligominer.specificity.duplex_stability.bilstm_arch import OM2BiLSTM

    try:
        import torch                                                   # noqa: F401
    except ImportError:
        raise MissingDependency(
            'torch (required by duplex-BiLSTM; install with pip install "oligominer[bilstm]", '
            'or use physics-xgb, which needs no extra dependency)')

    # success
    return OM2BiLSTM.load(str(path), device=device)


def load_all(include_baselines=True):
    """
    Load every registered model.

    Args:
        include_baselines (bool): include the incumbent baselines as well as the
            OligoMiner2 models.

    Returns:
        models (dict): registry key mapped to a LoadedModel.
    """
    from .registry import BASELINE_MODELS, OM2_MODELS

    names = list(OM2_MODELS)
    if include_baselines:
        names += list(BASELINE_MODELS)

    # success
    return {name: load(name) for name in names}
