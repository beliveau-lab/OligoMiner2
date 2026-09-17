from ._version import __version__

from ._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    "ProbeSet": ".probe_design.probe_set",
    "get_example_fasta_path": ".data.test_data",
    "load_example_fasta": ".data.test_data",
    "mine_fasta": ".thermodynamics.mining",
    "mine_sequence": ".thermodynamics.mining",
    "write_probes": ".thermodynamics.mining",
}

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS)
