from ..._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    "KmerIndex": ".numpy_index",
    "build_index": ".backends",
    "calc_max_kmer": ".jellyfish_query",
    "calc_max_kmer_multi": ".jellyfish_query",
    "have_jellyfish": ".backends",
    "jellyfish_build": ".jellyfish_build",
    "jellyfish_query": ".jellyfish_query",
    "max_kmer": ".backends",
    "read_metadata": ".backends",
    "resolve_backend": ".backends",
    "sidecar_path": ".backends",
    "validate_index": ".jellyfish_query",
    "write_metadata": ".backends",
}

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS)
