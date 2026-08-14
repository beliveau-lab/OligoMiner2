"""
# Model registry

The single declaration of which duplex-stability models ship and what each one is.
Every consumer reads this rather than naming an artifact path, so a model's
identity has one source.

Four models are registered. Two are OligoMiner2's own and two are the incumbents
they are measured against; ``family`` distinguishes them so a comparison cannot
present a baseline as an OM2 result.

| name           | family   | what it is                                          |
|----------------|----------|-----------------------------------------------------|
| physics-xgb    | om2      | 103 L4t features, width-free, condition-aware        |
| duplex-BiLSTM  | om2      | reads the aligned duplex, no condition input         |
| ps-xgb         | baseline | the PaintSHOP model, 37 features, condition-blind    |
| om1-lda        | baseline | the OligoMiner1 model, emits a decision score        |

``physics-xgb``'s features are PaintSHOP's 37 verbatim, prefixed ``ps_``, plus 66
added physics features, so comparing it against ``ps-xgb`` is a controlled
ablation rather than a comparison of unrelated models.

``duplex-BiLSTM`` has no condition channel and was fitted at one condition
(69.5 C effective, 37 C with 50% formamide). It cannot be asked about another
temperature.

``om1-lda`` returns a decision score on an arbitrary scale, not a probability.
"""

from importlib.resources import files

MODEL_DATA = files('oligominer.data.models.zoo')

REGISTRY = {
    'physics-xgb': {
        'display': 'physics-XGB',
        'family': 'om2',
        'kind': 'xgb',
        'artifact': 'om2_xgb_flagship.pkl',
        'card': 'om2_xgb_flagship.card.json',
        'encoding': 'L4t',
        'n_features': 103,
        'width_free': True,
        'outputs_pdup': True,
        'link': 'logit',
        'condition_aware': True,
    },
    'duplex-BiLSTM': {
        'display': 'duplex-BiLSTM',
        'family': 'om2',
        'kind': 'bilstm',
        'artifact': 'om2_bilstm_flagship.pt',
        'card': 'om2_bilstm_flagship.card.json',
        'encoding': 'duplex-paircol',
        'ncond': 0,
        'outputs_pdup': True,
        'link': None,
        'condition_aware': False,
    },
    'ps-xgb': {
        'display': 'PaintSHOP-XGB',
        'family': 'baseline',
        'kind': 'xgb',
        'artifact': 'ps_xgb_baseline.pkl',
        'card': 'ps_xgb_baseline.card.json',
        'encoding': 'paintshop-37feat',
        'n_features': 37,
        'width_free': True,
        'outputs_pdup': True,
        'link': 'logit',
        'condition_aware': False,
    },
    'om1-lda': {
        'display': 'OligoMiner1-LDA',
        'family': 'baseline',
        'kind': 'lda',
        'artifact': 'om1_lda_baseline.pkl',
        'card': 'om1_lda_baseline.card.json',
        'encoding': 'alignment-3feat',
        'n_features': 3,
        'outputs_pdup': False,
        'link': None,
        'condition_aware': False,
    },
}

# the models this package ships as its own, in the order a comparison presents them
OM2_MODELS = [name for name, entry in REGISTRY.items() if entry['family'] == 'om2']

# the incumbents OM2 is measured against
BASELINE_MODELS = [name for name, entry in REGISTRY.items()
                   if entry['family'] == 'baseline']

# the model used when a caller does not name one
DEFAULT_MODEL = 'physics-xgb'


def spec(name):
    """
    Return one model's registry entry.

    Args:
        name (str): a registry key.

    Returns:
        entry (dict): the model's declaration.

    Raises:
        KeyError: naming the valid keys, so a typo cannot become a silent fallback.
    """
    if name not in REGISTRY:
        raise KeyError(f'unknown model {name!r}; the registered models are {sorted(REGISTRY)}')

    # success
    return REGISTRY[name]


def artifact_path(name):
    """
    Return the path to a model's artifact file.

    Args:
        name (str): a registry key.

    Returns:
        path (pathlib.Path): the artifact's location in the package data.
    """
    # success
    return MODEL_DATA / spec(name)['artifact']


def card_path(name):
    """
    Return the path to a model's card.

    Args:
        name (str): a registry key.

    Returns:
        path (pathlib.Path): the card's location in the package data.
    """
    # success
    return MODEL_DATA / spec(name)['card']


def available():
    """
    Return the registered model names.

    Returns:
        names (list): every key of the registry, sorted.
    """
    # success
    return sorted(REGISTRY)


def missing_artifacts():
    """
    Return any registered artifact or card that is not present on disk.

    Returns:
        missing (list): 'name:key -> filename' for each absent file.
    """
    missing = []
    for name, entry in REGISTRY.items():
        for key in ('artifact', 'card'):
            if not (MODEL_DATA / entry[key]).is_file():
                missing.append(f"{name}:{key} -> {entry[key]}")

    # success
    return missing
