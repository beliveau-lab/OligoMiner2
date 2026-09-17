"""# Lazy package exports

Builds the PEP 562 hooks that let a package re-export a name without importing
the module that defines it.

A package `__init__` that re-exports eagerly makes every import pay for the whole
subtree: `from oligominer.thermodynamics.mining import mine_sequence` runs the
`oligominer` `__init__`, which imports `probe_design`, which imports the model
loaders, which import xgboost. Mining candidate probes by nearest-neighbor
thermodynamics needs none of that, but it waits for all of it, once per process.
A pipeline that shards across a genome forks that cost rather than amortizing it.

Dropping the re-exports would fix it and would break every documented import
path. PEP 562 fixes it without: a module-level `__getattr__` runs on the first
access to a name the module does not already have, so `oligominer.ProbeSet`
imports `probe_design.probe_set` at that moment and not before. The resolved
value is cached in the module globals, so the hook runs once per name.

## What stays eager

A name whose module has an import-time side effect the caller must not miss, and
a name cheap enough that deferring it only adds indirection. `__version__` is the
second kind. The nupack subpackage is the first: its `__init__` checks the nupack
installation, and that check has to run when the subpackage is imported rather
than when one of its functions is first touched.
"""

import importlib
import sys
from types import ModuleType


class _LazyPackage(ModuleType):
    """A package whose re-exported names win against same-named submodules.

    Most lazy names need nothing but `__getattr__`, which runs only when normal
    lookup fails. A name that collides with a submodule defeats that: importing
    `oligominer.thermodynamics.formamide_correction` binds the *module* as an
    attribute of `oligominer.thermodynamics`, so the later
    `from oligominer.thermodynamics import formamide_correction` finds an
    attribute, never calls `__getattr__`, and hands the caller a module where the
    eager re-export used to hand back a function.

    Whether that happens depends on whether something else imported the submodule
    first, which makes it an ordering bug: the import works alone and breaks in a
    process that touched the submodule earlier.

    So the colliding names, and only those, resolve through `__getattribute__`,
    which runs unconditionally. Every other attribute takes the normal path.
    """

    def __getattribute__(self, name):
        namespace = object.__getattribute__(self, "__dict__")
        shadowed = namespace.get("_SHADOWED_EXPORTS")

        if shadowed and name in shadowed:
            value = namespace.get(name)
            if value is not None and not isinstance(value, ModuleType):
                # already resolved to the function
                return value
            module = importlib.import_module(shadowed[name], self.__name__)
            value = getattr(module, name)
            namespace[name] = value

            # success
            return value

        # success
        return ModuleType.__getattribute__(self, name)


def lazy_exports(package, exports, submodules=()):
    """Build the `__getattr__`, `__dir__` and `__all__` a lazy package needs.

    Args:
        package (str): the importing package's `__name__`.
        exports (dict): public name -> the relative module defining it, as
            '.module' or '.sub.module'.
        submodules (tuple): submodule names the package exposes as attributes.

    Returns:
        hooks (tuple): (`__getattr__`, `__dir__`, `__all__`), to be unpacked into
            the package's module globals.
    """
    names = sorted({*exports, *submodules})

    # names a same-named submodule would otherwise shadow
    shadowed = {
        name: module for name, module in exports.items() if module.rsplit(".", 1)[-1] == name
    }
    if shadowed:
        module = sys.modules[package]
        module._SHADOWED_EXPORTS = shadowed
        module.__class__ = _LazyPackage

    def __getattr__(name):
        if name in submodules:
            value = importlib.import_module(f"{package}.{name}")
        elif name in exports:
            module = importlib.import_module(exports[name], package)
            value = getattr(module, name)
        else:
            raise AttributeError(f"module {package!r} has no attribute {name!r}")

        # cache on the package so the hook runs once per name
        importlib.import_module(package).__dict__[name] = value

        # success
        return value

    def __dir__():
        return sorted({*importlib.import_module(package).__dict__, *names})

    # success
    return __getattr__, __dir__, names
