from .._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    "BASELINE_MODELS": ".registry",
    "DEFAULT_MODEL": ".registry",
    "LoadedModel": ".loaders",
    "OM2_MODELS": ".registry",
    "REGISTRY": ".registry",
    "artifact_path": ".registry",
    "available": ".registry",
    "build_flat_corpus": ".retrain",
    "card_path": ".registry",
    "load": ".loaders",
    "load_all": ".loaders",
    "missing_artifacts": ".registry",
    "retrain": ".retrain",
    "spec": ".registry",
    "write_card": ".retrain",
}

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS)
