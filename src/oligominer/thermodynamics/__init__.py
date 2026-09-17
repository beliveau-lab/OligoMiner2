from .._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    "effective_hyb_temperature": ".formamide_correction",
    "formamide_correction": ".formamide_correction",
}

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS)
