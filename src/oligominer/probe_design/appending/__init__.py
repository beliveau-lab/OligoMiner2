from ..._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    'LINKER':                 '.config',
    'add_bridges':            '.paintshop_appending',
    'append_barcodes':        '.paintshop_appending',
    'append_custom':          '.appending',
    'append_multiple':        '.appending',
    'append_saber':           '.paintshop_appending',
    'append_same':            '.appending',
    'append_sequences':       '.appending',
    'append_unique':          '.appending',
    'build_appending_table':  '.appending',
    'collect_indices':        '.paintshop_appending',
}

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS)
