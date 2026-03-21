"""Utility functions for validating the existence of required files and directories."""

import os

from .exceptions import MissingInputFile, MissingOutputFile, MissingDirectory, DirectoryPathError
from .file_paths import get_abs_path


def check_input_exists(file_path):
    """
    Checks that a required input file exists on disk.

    Args:
        file_path (str): path to the input file.

    Returns:
        result (bool): True if the file exists.

    Raises:
        MissingInputFile: if the file does not exist.
    """
    if not os.path.exists(file_path):
        raise MissingInputFile(file_path)

    # success
    return True


def check_output_exists(file_path):
    """
    Checks that an expected output file exists on disk.

    Args:
        file_path (str): path to the output file.

    Returns:
        result (bool): True if the file exists.

    Raises:
        MissingOutputFile: if the file does not exist.
    """
    if not os.path.exists(file_path):
        raise MissingOutputFile(file_path)

    # success
    return True


def check_dir_exists(dir_path, create=False, parent_dir=False):
    """
    Checks that a required directory exists on disk, optionally creating it.

    Args:
        dir_path (str): path to the target directory, or to a file whose
            parent directory should be checked (when parent_dir=True).
        create (bool): if True, create the directory if it does not exist.
        parent_dir (bool): if True, treat dir_path as a file path and check
            its parent directory instead.

    Returns:
        result (bool): True if the directory exists or was successfully created.

    Raises:
        DirectoryPathError: if a non-directory file already exists at the path.
        MissingDirectory: if the directory does not exist and create is False.
    """

    # optionally resolve to parent directory of a file path
    if parent_dir:
        dir_path = os.path.dirname(get_abs_path(dir_path))
    else:
        dir_path = get_abs_path(dir_path)

    # check whether this path exists on disk
    if os.path.exists(dir_path) and not os.path.isdir(dir_path):
        # an existing file (not dir) was found at the target path
        raise DirectoryPathError(dir_path)
    else:
        # this path does not exist
        if create:
            # create the target directory
            os.makedirs(dir_path, exist_ok=True)
        if not os.path.exists(dir_path):
            # the target directory does not exist
            raise MissingDirectory(dir_path)

    # success
    return True
