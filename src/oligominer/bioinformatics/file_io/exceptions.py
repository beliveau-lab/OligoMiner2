"""
# File I/O Exceptions

Custom exceptions for file I/O operations.
"""

from oligominer.utils import get_abs_path, get_dir_name
from oligominer.utils.exceptions import OligominerError


class FileIOError(OligominerError):
    """Base class for exceptions in this module."""
    __module__ = OligominerError.__module__


class FastaPermissionError(FileIOError):
    """Raised when the directory containing a FASTA file is not writeable."""

    # don't print module full path in traceback
    __module__ = OligominerError.__module__

    def __init__(self, file_path):
        self.fasta_path = get_abs_path(file_path)
        self.fasta_dir = get_dir_name(file_path)
        super().__init__(f"Directory not writeable for FASTA index: {self.fasta_dir}")

    def __str__(self):
        error_text = f'Error creating .fai file.\n\nThe directory containing the '
        error_text += f'input fasta is not writeable:\n\n'
        error_text += f'  directory: {self.fasta_dir}\n\n'
        error_text += f'  fasta file: {self.fasta_path}\n\n'
        error_text += f'Please copy this fasta file to a writeable location or\n'
        error_text += f'modify permissions as needed and try again.\n\nExiting...'
        return error_text


class EmptyExportError(FileIOError):
    """Raised when an export is attempted with no records.

    Covers both FASTA and GTF exports. The format_name is included in
    the error message to help the user identify which pipeline step
    produced an empty result.
    """

    # don't print module full path in traceback
    __module__ = OligominerError.__module__

    def __init__(self, file_path, format_name='file'):
        self.file_path = get_abs_path(file_path)
        self.format_name = format_name
        super().__init__(
            f"No records for {format_name} export: {self.file_path}"
        )

    def __str__(self):
        error_text = f'No records present during attempted {self.format_name} export. This\n'
        error_text += f'can happen if you apply a filter such that all records are removed.\n\n'
        error_text += f'The target export path provided was:\n\n{self.file_path}\n\nExiting...'
        return error_text
