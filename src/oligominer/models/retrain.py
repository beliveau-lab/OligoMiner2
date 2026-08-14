"""
# Re-training the OligoMiner2 models

Fits `physics-xgb` or `duplex-BiLSTM` on a new corpus and writes a card recording
what was fitted and on what.

## Three rules this enforces

**A registered artifact is never overwritten.** Re-training produces a new file;
promoting it to the shipping model is a separate, deliberate act.

**Every output gets a card carrying the corpus checksum.** A model whose training
data cannot be identified is not reproducible, and identifying it from prose has
failed repeatedly.

**The training regime defaults to the model's own and warns when it is left.**
`duplex-BiLSTM` has no condition channel, so a corpus spanning several
temperatures shows it identical tokens at different labels and it correctly
learns to hedge. The API will do it if asked, and says so.

## Why a flat corpus

`build_flat_corpus` takes the same number of rows from each pDup decile. A real
genomic duplex population is bimodal, with most rows in the two extreme deciles,
so an unweighted draw teaches the base rate rather than the sequence signal.

A metric measured on such a cut is arithmetic about the cut, not a finding.
Evaluate on a random draw from aligner output instead.
"""

import hashlib
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from .registry import REGISTRY, artifact_path, spec

# the condition each model's shipping artifact was fitted at
SHIPPING_CELSIUS = {'duplex-BiLSTM': 69.5, 'physics-xgb': None}

# pDup deciles the flat corpus balances across
N_DECILES = 10


def _sha256(path):
    """
    Return a file's SHA-256 digest.

    Args:
        path (str or pathlib.Path): the file.

    Returns:
        digest (str): the hex digest.
    """
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)

    # success
    return digest.hexdigest()


def build_flat_corpus(duplexes, out, celsius=None, n_per_decile=None, seed=0,
                      label_col='pdup'):
    """
    Cut a corpus with equal representation in every pDup decile.

    Args:
        duplexes (pandas.DataFrame or str): duplex rows, or a parquet path. Must
            carry the label column and the aligned columns the encoders read.
        out (str): destination parquet path.
        celsius (float, optional): keep only rows at this temperature. None pools
            every temperature, which suits physics-xgb and not duplex-BiLSTM.
        n_per_decile (int, optional): rows per decile. The smallest decile's size
            is used when None, which is what makes the cut exactly flat.
        seed (int): sampling seed, recorded in the sidecar.
        label_col (str): the column holding pDup.

    Returns:
        path (pathlib.Path): the corpus written.

    Raises:
        KeyError: if the label column is absent.
    """
    df = pd.read_parquet(duplexes) if isinstance(duplexes, str) else duplexes.copy()
    if label_col not in df.columns:
        raise KeyError(f'no {label_col!r} column in the duplex table')

    if celsius is not None:
        for column in ('label_celsius', 'celsius', 't_eff'):
            if column in df.columns:
                df = df[np.isclose(df[column], celsius)].copy()
                break

    decile = np.clip((df[label_col] * N_DECILES).astype(int), 0, N_DECILES - 1)
    counts = decile.value_counts().sort_index()
    if n_per_decile is None:
        n_per_decile = int(counts.min())

    rng = np.random.default_rng(seed)
    picks = []
    for value, group in df.groupby(decile):
        take = min(n_per_decile, len(group))
        picks.append(group.iloc[rng.choice(len(group), take, replace=False)])
    cut = pd.concat(picks, ignore_index=True)

    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cut.to_parquet(out)

    sidecar = {
        'rows': int(len(cut)),
        'n_per_decile': int(n_per_decile),
        'celsius': celsius,
        'seed': seed,
        'label_col': label_col,
        'flat_in_pdup': True,
        'decile_counts': {int(k): int(v) for k, v in counts.items()},
        'sha256': _sha256(out),
        'built_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    }
    out.with_suffix('.corpus.json').write_text(json.dumps(sidecar, indent=1))

    # success
    return out


def _check_not_registered(out):
    """
    Refuse to write over a registered artifact.

    Args:
        out (pathlib.Path): the destination path.

    Returns:
        out (pathlib.Path): the same path.

    Raises:
        ValueError: if out is a registered model's artifact.
    """
    resolved = Path(out).resolve()
    for name in REGISTRY:
        if Path(str(artifact_path(name))).resolve() == resolved:
            raise ValueError(
                f'{out} is the registered artifact for {name!r}. Re-training writes a '
                f'new file; promoting it is a separate, deliberate act.')

    # success
    return Path(out)


def _warn_off_regime(name, celsius):
    """
    Warn when a model is being fitted away from its shipping regime.

    Args:
        name (str): the registry key.
        celsius (float or None): the corpus temperature, or None if pooled.

    Returns:
        ok (bool): True when the regime matches, False when a warning was issued.
    """
    expected = SHIPPING_CELSIUS.get(name, None)

    if name == 'duplex-BiLSTM' and celsius is None:
        warnings.warn(
            'duplex-BiLSTM has no condition channel, so a corpus pooling several '
            'temperatures shows it identical tokens at different labels and it will '
            'correctly learn to hedge. Pass celsius to fit at one condition.',
            RuntimeWarning, stacklevel=3)
        return False

    if expected is not None and celsius is not None and celsius != expected:
        warnings.warn(
            f'{name} ships fitted at {expected} C and is being fitted at {celsius} C. '
            f'The result is a different model and must be labeled as such.',
            RuntimeWarning, stacklevel=3)
        return False

    # success
    return True


def _corpus_digest(corpus):
    """
    Return the corpus checksum when the corpus is a file on disk.

    An in-memory frame has no stable identity to record, which is itself worth
    knowing: such a fit cannot be reproduced from the card alone.

    Args:
        corpus (str or pathlib.Path or pandas.DataFrame): the corpus.

    Returns:
        digest (str or None): the hex digest, or None for an in-memory frame.
    """
    if isinstance(corpus, (str, Path)) and Path(corpus).is_file():
        return _sha256(corpus)

    # success
    return None


def write_card(name, out, corpus, celsius, seed, extra=None):
    """
    Write the card recording what was fitted and on what.

    Args:
        name (str): the registry key the fit is based on.
        out (str or pathlib.Path): the artifact written.
        corpus (str or pathlib.Path): the corpus used.
        celsius (float or None): the condition fitted at.
        seed (int): the training seed.
        extra (dict, optional): further fields to record.

    Returns:
        path (pathlib.Path): the card written.
    """
    entry = spec(name)
    card = {
        'name': f'{name}-retrained',
        'based_on': name,
        'role': 'RETRAINED',
        'family': entry['family'],
        'kind': entry['kind'],
        'encoding': entry['encoding'],
        'outputs_pdup': entry['outputs_pdup'],
        'link': entry.get('link'),
        'celsius': celsius,
        'seed': seed,
        'corpus': str(corpus) if isinstance(corpus, (str, Path)) else 'in-memory frame',
        'corpus_sha256': _corpus_digest(corpus),
        'fitted_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    }
    card.update(extra or {})

    path = Path(out).with_suffix('.card.json')
    path.write_text(json.dumps(card, indent=1))

    # success
    return path


def retrain(name, corpus, out, celsius=None, seed=0, label_col='pdup',
            hyperparams=None):
    """
    Fit a registered model on a new corpus and write its card.

    Args:
        name (str): the registry key to re-fit, 'physics-xgb' or 'duplex-BiLSTM'.
        corpus (str or pandas.DataFrame): the training corpus.
        out (str): destination artifact path. Must not be a registered artifact.
        celsius (float, optional): the condition the corpus represents.
        seed (int): training seed.
        label_col (str): the column holding pDup.
        hyperparams (dict, optional): overrides for the fit.

    Returns:
        info (dict): the artifact and card paths, row count and metrics.

    Raises:
        ValueError: if out is a registered artifact, or the model cannot be
            re-fitted through this API.
    """
    entry = spec(name)
    if entry['family'] != 'om2':
        raise ValueError(
            f'{name!r} is a {entry["family"]} model kept for comparison; only '
            f'OligoMiner2 models are re-fitted here')

    out = _check_not_registered(out)
    _warn_off_regime(name, celsius)

    df = pd.read_parquet(corpus) if isinstance(corpus, str) else corpus
    if label_col not in df.columns:
        raise KeyError(f'no {label_col!r} column in the corpus')

    if entry['kind'] == 'xgb':
        info = _fit_xgb(name, df, out, label_col, seed, hyperparams)
    else:
        raise ValueError(
            f'{name!r} is fitted with torch on a GPU, which this API does not '
            f'orchestrate. Use the corpus builder here and fit it in a GPU job.')

    info['card'] = str(write_card(name, out, corpus, celsius, seed,
                                  extra={'n_train': int(len(df))}))

    # success
    return info


def _fit_xgb(name, df, out, label_col, seed, hyperparams):
    """
    Fit the tree model on a corpus and save it.

    Args:
        name (str): the registry key.
        df (pandas.DataFrame): the corpus.
        out (pathlib.Path): destination artifact.
        label_col (str): the label column.
        seed (int): training seed.
        hyperparams (dict or None): overrides.

    Returns:
        info (dict): artifact path and row count.
    """
    import pickle

    import xgboost as xgb

    from .loaders import _encode

    entry = spec(name)
    features = _encode(entry['encoding'], df, 42)
    labels = df[label_col].to_numpy(dtype=np.float64)

    # the model is fitted on the logit of pDup, matching its link
    if entry.get('link') == 'logit':
        clipped = np.clip(labels, 1e-6, 1 - 1e-6)
        labels = np.log(clipped / (1 - clipped))

    params = {'max_depth': 8, 'n_estimators': 200, 'learning_rate': 0.1,
              'subsample': 0.9, 'random_state': seed}
    params.update(hyperparams or {})

    model = xgb.XGBRegressor(**params)
    model.fit(features.values.astype(np.float32), labels)

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('wb') as handle:
        pickle.dump({'model': model, 'features': list(features.columns)}, handle)

    # success
    return {'artifact': str(out), 'n_train': int(len(df)),
            'n_features': int(features.shape[1]), 'params': params}
