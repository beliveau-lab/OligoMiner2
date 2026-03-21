"""
XGBoost-based duplex stability prediction, inspired by the PaintSHOP pipeline.

Approximates NUPACK pDup values using pre-trained XGBoost models at various
temperatures. Each model was trained on NUPACK-computed duplex probabilities
and predicts from sequence-derived features (GC content, lengths, and
dinucleotide counts).

Available model temperatures: 37, 42, 47, 52, 60 °C.
"""

from importlib.resources import files

import pandas as pd
import xgboost as xgb

from oligominer.utils.seq_utils import rev_comp, calc_gc, clamp
from oligominer.utils.exceptions import ConfigurationError
from .config import XGBOOST_TEMPERATURES as AVAILABLE_TEMPERATURES

DINUCLEOTIDES = [
    'AA', 'AT', 'AG', 'AC',
    'TA', 'TT', 'TG', 'TC',
    'GA', 'GT', 'GG', 'GC',
    'CA', 'CT', 'CG', 'CC',
]

# ordered feature columns expected by the model
FEATURE_COLUMNS = (
    ['align_score', 'probe_gc', 'derived_gc', 'probe_len', 'derived_len']
    + [f'probe_{dn}' for dn in DINUCLEOTIDES]
    + [f'derived_{dn}' for dn in DINUCLEOTIDES]
)

# cache loaded models to avoid repeated disk reads
_model_cache = {}


def load_model(temperature):
    """
    Load a pre-trained PaintSHOP XGBoost model for a given temperature.

    Models are bundled as package data and cached after first load.

    Args:
        temperature (int): hybridization temperature in °C. Must be one of
            37, 42, 47, 52, or 60.

    Returns:
        model: the unpickled XGBoost model object.
    """
    if temperature not in AVAILABLE_TEMPERATURES:
        raise ConfigurationError(
            f"No model available for {temperature}°C. "
            f"Available temperatures: {AVAILABLE_TEMPERATURES}"
        )

    if temperature in _model_cache:
        return _model_cache[temperature]

    filename = f"{temperature}_all_fixed_xgb.ubj"

    # load xgboost model from ubj files
    model_ref = files("oligominer.data.models.paintshop").joinpath(filename)

    # Initialize an empty Booster object
    xgb_model = xgb.Booster()

    # Load the model from the .ubj file
    xgb_model.load_model(model_ref)

    _model_cache[temperature] = xgb_model

    # success
    return xgb_model



def compute_features(df):
    """
    Compute sequence-derived features for XGBoost duplex prediction.

    Takes a DataFrame with probe_seq, derived_seq, and align_score columns
    and returns a feature matrix matching the model's expected input.

    The derived sequence is reverse-complemented to match the orientation
    used during model training.

    Args:
        df (pandas.DataFrame): must contain columns probe_seq, derived_seq,
            and align_score.

    Returns:
        features (pandas.DataFrame): feature matrix with columns ordered
            for model input.
    """
    features = pd.DataFrame()

    features['align_score'] = df['align_score'].values

    # reverse complement derived to match training data orientation
    derived_rc = df['derived_seq'].apply(rev_comp)

    # gc content as percentage to match model training data
    features['probe_gc'] = df['probe_seq'].apply(calc_gc, as_percent=True)
    features['derived_gc'] = derived_rc.apply(calc_gc, as_percent=True)

    # sequence lengths
    features['probe_len'] = df['probe_seq'].str.len()
    features['derived_len'] = derived_rc.str.len()

    # dinucleotide counts for both sequences
    for dn in DINUCLEOTIDES:
        features[f'probe_{dn}'] = df['probe_seq'].str.count(dn)
        features[f'derived_{dn}'] = derived_rc.str.count(dn)

    # reorder to match model expectation
    features = features[FEATURE_COLUMNS]

    # success
    return features


def predict_duplex_batch(df, temperature=37, normalize=True):
    """
    Predict duplex formation probability for a batch of probe-target pairs.

    Computes features from the input DataFrame and runs the XGBoost model
    to produce predictions analogous to NUPACK pDup.

    Args:
        df (pandas.DataFrame): must contain columns probe_seq, derived_seq,
            and align_score.
        temperature (int): hybridization temperature in °C. Must be one of
            37, 42, 47, 52, or 60.
        normalize (bool): if True, scale predictions from the native 0-100
            range to 0.0-1.0 to match calc_pdup output.

    Returns:
        predictions (numpy.ndarray): predicted duplex probabilities for
            each row, in 0.0-1.0 range if normalize is True, else 0-100.
    """
    model = load_model(temperature)
    features = compute_features(df)

    predictions = model.predict(xgb.DMatrix(features.values))

    if normalize:
        predictions = clamp(predictions / 100.0, 0.0, 1.0)
    else:
        predictions = clamp(predictions, 0.0, 100.0)

    # success
    return predictions


def predict_duplex(probe_seq, derived_seq, align_score,
                   temperature=37, normalize=True):
    """
    Predict duplex formation probability for a single probe-target pair.

    Convenience wrapper around predict_duplex_batch for single-pair usage,
    providing a similar interface to calc_pdup.

    Args:
        probe_seq (str): the probe DNA sequence.
        derived_seq (str): the off-target (derived) DNA sequence.
        align_score (int): the Bowtie2 alignment score.
        temperature (int): hybridization temperature in °C. Must be one of
            37, 42, 47, 52, or 60.
        normalize (bool): if True, scale prediction from the native 0-100
            range to 0.0-1.0 to match calc_pdup output.

    Returns:
        prediction (float): predicted duplex formation probability.
    """
    df = pd.DataFrame({
        'probe_seq': [probe_seq],
        'derived_seq': [derived_seq],
        'align_score': [align_score],
    })

    predictions = predict_duplex_batch(
        df, temperature=temperature, normalize=normalize
    )

    # success
    return float(predictions[0])
