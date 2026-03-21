"""
# Alignment Exceptions

Custom exceptions for Bowtie2 index validation and alignment errors.
"""

from oligominer.utils import get_abs_path
from oligominer.utils.exceptions import OligominerError

class BowtieError(OligominerError):
    """Base class for exceptions in this module."""
    __module__ = OligominerError.__module__

class BowtieIndexError(BowtieError):
    """Exception raised for invalid index used to query a given Bowtie2 index file."""

    def __init__(self, file_path, error_msg):
        self.file_path = get_abs_path(file_path)
        self.error_msg = error_msg
        super().__init__(f"Invalid Bowtie2 index: {self.file_path}")

    def __str__(self):
        error_text = (
            f'Invalid index for use with Bowtie2 file:\n\n'
            f'{self.file_path}\n\nError message: {self.error_msg}.\n\nExiting...'
        )
        return error_text

class MissingBowtieIndexError(BowtieError):
    """Exception raised when a Bowtie2 index file is not found during an attempted query."""

    def __init__(self, file_path):
        self.file_path = get_abs_path(file_path)
        super().__init__(f"Bowtie2 index not found: {self.file_path}")

    def __str__(self):
        error_text = (
            f'A Bowtie2 file was not found at the specified path:\n\n{self.file_path}\n\n'
            f'For info on creating this file, run:\n\n  $ oligominer build_bowtie2 --help\n\nExiting...'
        )
        return error_text
