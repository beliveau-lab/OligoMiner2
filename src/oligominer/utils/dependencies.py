"""
# Required dependency helper functions

Utilities for checking the presence of required dependencies for OligoMiner2.
"""

from shutil import which

from .exceptions import MissingDependency, DEFAULT_HELP_URL


def ensure_executable(executable, raise_on_missing=True, help_url=DEFAULT_HELP_URL):
    """
    Check if an executable is available in the system PATH.

    Args:
        executable (str): name of the executable to check.
        raise_on_missing (bool): whether to raise if the executable is not
            found. If False, the function returns False when the executable
            is missing.
        help_url (str): URL for installation instructions.

    Returns:
        exe_path (str or bool): path to the executable if found, or False
            if missing and raise_on_missing is False.

    Raises:
        MissingDependency: if the executable is not found and
            raise_on_missing is True.
    """
    # check whether executable is currently found in PATH
    exe_path = which(executable)
    if exe_path is None:
        if raise_on_missing:
            raise MissingDependency(executable, help_url=help_url)
        else:
            return False

    # success
    return exe_path


def ensure_python_package(package_name, raise_on_missing=True, help_url=DEFAULT_HELP_URL):
    """
    Check if a Python package can be imported.

    Args:
        package_name (str): name of the package to check.
        raise_on_missing (bool): whether to raise if the package is not
            found. If False, the function returns False when the package
            is missing.
        help_url (str): URL for installation instructions.

    Returns:
        pkg_exists (bool): True if the package can be imported, False
            if missing and raise_on_missing is False.

    Raises:
        MissingDependency: if the package cannot be imported and
            raise_on_missing is True.
    """
    try:
        __import__(package_name)
        pkg_exists = True
    except ImportError:
        pkg_exists = False
        if raise_on_missing:
            raise MissingDependency(package_name, help_url=help_url)

    # success
    return pkg_exists
