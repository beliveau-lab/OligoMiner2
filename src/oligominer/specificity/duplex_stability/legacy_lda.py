"""
Legacy LDA-based duplex stability prediction from OligoMiner v1.

Predicts the probability that a probe has thermodynamically relevant
off-target binding using pre-fit Linear Discriminant Analysis models at
various temperatures. The model coefficients were extracted from the
original OligoMiner outputClean.py script.

See: https://github.com/beliveau-lab/OligoMiner/blob/master/outputClean.py

Available model temperatures: 32, 37, 42, 47, 52, 57 °C.
"""

import numpy as np
import pandas as pd
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

from oligominer.utils.seq_utils import calc_gc, clamp
from oligominer.utils.exceptions import ConfigurationError
from .config import LDA_TEMPERATURES as AVAILABLE_TEMPERATURES

# ordered feature columns expected by the model
FEATURE_COLUMNS = ('probe_len', 'align_score', 'probe_gc')

# pre-fit LDA coefficients from OligoMiner v1, indexed by temperature
# features are [probe_len, align_score, gc_percentage]
_MODEL_COEFS = {
    32: [[-0.14494789, 0.18791679, 0.02588474]],
    37: [[-0.13364364, 0.22510179, 0.05494031]],
    42: [[-0.09006122, 0.25660706, 0.1078303]],
    47: [[-0.01593182, 0.24498485, 0.15753649]],
    52: [[0.01860365, 0.1750174, 0.17003374]],
    57: [[0.03236755, 0.11624593, 0.24306498]],
}

_MODEL_INTERCEPTS = {
    32: -1.17545204,
    37: -5.40436344,
    42: -12.45549846,
    47: -19.32670233,
    52: -20.11992898,
    57: -23.98652919,
}

_MODEL_CLASSES = np.array([-1, 1])

# cache loaded models to avoid repeated construction
_model_cache = {}


def load_model(temperature):
    """
    Load a pre-fit OligoMiner v1 LDA model for a given temperature.

    Models are constructed from hardcoded coefficients extracted from the
    original OligoMiner outputClean.py and cached after first construction.

    Args:
        temperature (int): hybridization temperature in °C. Must be one of
            32, 37, 42, 47, 52, or 57.

    Returns:
        model (sklearn.discriminant_analysis.LinearDiscriminantAnalysis):
            the pre-fit LDA classifier.
    """
    if temperature not in AVAILABLE_TEMPERATURES:
        raise ConfigurationError(
            f"No model available for {temperature}°C. "
            f"Available temperatures: {AVAILABLE_TEMPERATURES}"
        )

    if temperature in _model_cache:
        return _model_cache[temperature]

    # build model from pre-fit coefficients
    clf = LinearDiscriminantAnalysis()
    clf.coef_ = np.array(_MODEL_COEFS[temperature])
    clf.intercept_ = np.array([_MODEL_INTERCEPTS[temperature]])
    clf.classes_ = _MODEL_CLASSES

    _model_cache[temperature] = clf

    # success
    return clf


def compute_features(df):
    """
    Compute sequence-derived features for LDA duplex prediction.

    Takes a DataFrame with probe_seq and align_score columns and returns
    a feature matrix matching the model's expected input.

    Args:
        df (pandas.DataFrame): must contain columns probe_seq and
            align_score.

    Returns:
        features (pandas.DataFrame): feature matrix with columns
            probe_len, align_score, and probe_gc.
    """
    features = pd.DataFrame()

    features['probe_len'] = df['probe_seq'].str.len().astype(float)
    features['align_score'] = df['align_score'].astype(float)

    # gc as percentage (0-100) to match original BioPython GC() used in training
    features['probe_gc'] = df['probe_seq'].apply(calc_gc, as_percent=True)

    # success
    return features


def predict_duplex_batch(df, temperature=42, normalize=True):
    """
    Predict off-target duplex probability for a batch of probe-target pairs.

    Computes features from the input DataFrame and runs the legacy LDA
    model to produce predictions. The model outputs the probability that
    a probe has thermodynamically relevant off-target binding.

    Args:
        df (pandas.DataFrame): must contain columns probe_seq and
            align_score.
        temperature (int): hybridization temperature in °C. Must be one of
            32, 37, 42, 47, 52, or 57.
        normalize (bool): if True, predictions are returned in 0.0-1.0
            range. If False, predictions are scaled to 0-100 for
            consistency with the XGBoost model interface.

    Returns:
        predictions (numpy.ndarray): predicted off-target duplex
            probabilities for each row.
    """
    model = load_model(temperature)
    features = compute_features(df)

    # probability of class 1 (having off-target binding)
    predictions = model.predict_proba(features.values)[:, 1]

    if normalize:
        predictions = clamp(predictions, 0.0, 1.0)
    else:
        predictions = clamp(predictions * 100.0, 0.0, 100.0)

    # success
    return predictions


def predict_duplex(probe_seq, align_score, temperature=42, normalize=True):
    """
    Predict off-target duplex probability for a single probe.

    Convenience wrapper around predict_duplex_batch for single-probe usage.

    Args:
        probe_seq (str): the probe DNA sequence.
        align_score (int): the Bowtie2 alignment score.
        temperature (int): hybridization temperature in °C. Must be one of
            32, 37, 42, 47, 52, or 57.
        normalize (bool): if True, prediction is returned in 0.0-1.0
            range. If False, scaled to 0-100.

    Returns:
        prediction (float): predicted off-target duplex probability.
    """
    df = pd.DataFrame({
        'probe_seq': [probe_seq],
        'align_score': [align_score],
    })

    predictions = predict_duplex_batch(
        df, temperature=temperature, normalize=normalize
    )

    # success
    return float(predictions[0])
