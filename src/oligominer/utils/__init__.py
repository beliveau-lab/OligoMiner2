from .._lazy import lazy_exports

# public name -> the module that defines it
_EXPORTS = {
    "ConfigurationError": ".exceptions",
    "DirectoryPathError": ".exceptions",
    "ExternalCommandFailed": ".exceptions",
    "InvalidInputError": ".exceptions",
    "MissingDependency": ".exceptions",
    "MissingDirectory": ".exceptions",
    "MissingInputFile": ".exceptions",
    "MissingOutputFile": ".exceptions",
    "OligominerError": ".exceptions",
    "PipelineStateError": ".exceptions",
    "ShellPipeline": ".shell_pipeline",
    "check_dir_exists": ".required_files",
    "check_input_exists": ".required_files",
    "check_output_exists": ".required_files",
    "ensure_executable": ".dependencies",
    "ensure_python_package": ".dependencies",
    "get_abs_path": ".file_paths",
    "get_dir_name": ".file_paths",
    "require_one_of": ".input_dispatch",
    "run_cmd": ".shell_pipeline",
}

# submodules reachable as attributes of this package
_SUBMODULES = ("seq_utils",)

__getattr__, __dir__, __all__ = lazy_exports(__name__, _EXPORTS, _SUBMODULES)
