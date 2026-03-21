"""
# File path helpers

Utilities for working with file paths. These helpers normalize user-provided
paths to absolute paths, expanding user home directories as needed.
"""

import os
from pathlib import Path

from .exceptions import InvalidInputError


def get_abs_path(path):
    """Convert a path-like input to an absolute, normalized Path.

    Args:
        path (path-like): A path-like object (str, bytes, os.PathLike) that may be
            relative or absolute, may contain "~" or environment variables,
            and may or may not exist.

    Returns:
        abs_path (str): An absolute Path with "~" and environment variables
            expanded, and path components normalized. Symlinks are resolved when
            possible, but the path does not need to exist.
    """
    if path is None:
        raise InvalidInputError("path must be a path-like object, not None")

    # Support os.PathLike and bytes
    fs_path = os.fspath(path)

    # Expand env vars like $HOME / %USERPROFILE%
    fs_path = os.path.expandvars(fs_path)

    # Make Path and expand "~"
    p = Path(fs_path).expanduser()

    try:
        # Absolute, normalized, resolve symlinks where possible
        abs_path = p.resolve(strict=False)
    except RuntimeError:
        # Fallback for weird cases (e.g. symlink loops)
        abs_path = p.absolute()

    # convert final abs_path to string for consistency
    abs_path = str(abs_path)

    # success
    return abs_path


def get_dir_name(path):
    """Get the directory name component of a path-like input.

    Args:
        path (path-like): A path-like object (str, bytes, os.PathLike).

    Returns:
        dir_name (str): The directory name component of the input path.
    """
    # get absolute path
    abs_path = get_abs_path(path)

    # get directory name using os.path.dirname
    dir_name = os.path.dirname(abs_path)

    # success
    return dir_name
