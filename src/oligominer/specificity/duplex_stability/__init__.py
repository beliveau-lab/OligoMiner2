from ..._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    'AVAILABLE_TEMPERATURES':  '.paintshop_xgboost',
    'compute_features':        '.paintshop_xgboost',
    'load_model':              '.paintshop_xgboost',
    'pdup_summary':            '.pdup',
    'predict_duplex':          '.paintshop_xgboost',
    'predict_duplex_batch':    '.paintshop_xgboost',
    'predict_pdup':            '.pdup',
}

# submodules reachable as attributes of this package
_SUBMODULES = ('legacy_lda',)

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS, _SUBMODULES)
