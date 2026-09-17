from ..._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    'read_align_csv':   '.csv_io',
    'read_probe_csv':   '.csv_io',
    'write_align_csv':  '.csv_io',
    'write_probe_csv':  '.csv_io',
}

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS)
