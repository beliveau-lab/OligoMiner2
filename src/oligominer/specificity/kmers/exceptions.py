"""
# Kmer Exceptions

Custom exceptions for Jellyfish index validation and query errors.
"""

from oligominer.utils import get_abs_path
from oligominer.utils.exceptions import OligominerError

class JellyfishError(OligominerError):
    """Base class for exceptions in this module."""
    __module__ = OligominerError.__module__


class JellyfishIndexError(JellyfishError):
    """Exception raised for invalid k value used to query a given Jellyfish index file."""

    def __init__(self, file_path, jf_k, error_k):
        self.file_path = get_abs_path(file_path)
        self.jf_k = jf_k
        self.error_k = error_k
        super().__init__(f"Invalid k value ({error_k}) for Jellyfish index: {self.file_path}")

    def __str__(self):
        error_text = (
            f'Invalid k value ({self.error_k}) for use with Jellyfish file:\n\n'
            f'{self.file_path}\n\nk value for this index is {self.jf_k}.\n\nExiting...'
        )
        return error_text


class MissingJellyfishIndexError(JellyfishError):
    """Exception raised when a Jellyfish index file is not found during an attempted query."""

    def __init__(self, file_path):
        self.file_path = get_abs_path(file_path)
        super().__init__(f"Jellyfish index not found: {self.file_path}")

    def __str__(self):
        error_text = (
            f'A Jellyfish file was not found at the specified path:\n\n{self.file_path}\n\n'
            f'For info on creating this file, run:\n\n  $ oligominer build_jellyfish --help\n\nExiting...'
        )
        return error_text
