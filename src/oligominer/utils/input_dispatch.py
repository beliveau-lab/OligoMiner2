"""
# Input Dispatch

Validation helper for functions that accept exactly one of two mutually
exclusive input sources (e.g. a file path or in-memory data).
"""

from .exceptions import InvalidInputError


def require_one_of(a, b, name_a, name_b):
    """
    Validate that exactly one of two values is provided (not None).

    Args:
        a: first value.
        b: second value.
        name_a (str): parameter name for first value (used in error message).
        name_b (str): parameter name for second value (used in error message).

    Raises:
        InvalidInputError: if both or neither value is provided.
    """
    if (a is None) == (b is None):
        raise InvalidInputError(f"Exactly one of '{name_a}' or '{name_b}' must be provided.")
