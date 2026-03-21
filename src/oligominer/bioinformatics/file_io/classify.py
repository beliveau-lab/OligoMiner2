"""
Sequence record classification utilities.

Provides functions for classifying sequence IDs by user-defined regex
rules (e.g. canonical, alt, hap, fix, unlocalized, unplaced).
"""

import re

import pandas as pd

from .config import DEFAULT_CLASSIFICATION_RULES
from oligominer.utils import get_dir_name, check_dir_exists

def classify_seq_ids(seq_source, rules=None, default_category='canonical'):
    """
    Classify sequence IDs from any dict-like source by regex rules.

    Rules are evaluated in order; the first matching pattern wins.
    Sequence IDs that match no rule receive the default_category label.

    Args:
        seq_source (dict-like): any object with .keys() returning sequence IDs
            (e.g. pyfaidx.Fasta or dict).
        rules (dict or None): {category_name: regex_pattern} evaluated in order.
            None uses DEFAULT_CLASSIFICATION_RULES.
        default_category (str): label for IDs that match no rule.

    Returns:
        classifications (pandas.DataFrame): columns ['seq_id', 'category'],
            one row per sequence record in original key order.
    """
    if rules is None:
        rules = DEFAULT_CLASSIFICATION_RULES

    # pre-compile patterns
    compiled_rules = [(category, re.compile(pattern)) for category, pattern in rules.items()]

    # classify each sequence ID
    rows = []
    for seq_id in seq_source.keys():
        category = default_category
        for rule_category, pattern in compiled_rules:
            if pattern.search(seq_id):
                category = rule_category
                break
        rows.append((seq_id, category))

    classifications = pd.DataFrame(rows, columns=['seq_id', 'category'])

    # success
    return classifications


def classify_and_write(seq_source, output_path, rules=None,
                       default_category='canonical'):
    """
    Classify sequence IDs and write the results to a TSV file.

    Convenience wrapper around classify_seq_ids that also writes the
    classification table to disk. Returns the DataFrame for further use.

    Args:
        seq_source (dict-like): any object with .keys() returning sequence IDs.
        output_path (str): destination file path for TSV output.
        rules (dict or None): classification rules (see classify_seq_ids).
        default_category (str): label for IDs that match no rule.

    Returns:
        classifications (pandas.DataFrame): columns ['seq_id', 'category'].
    """
    classifications = classify_seq_ids(
        seq_source, rules=rules, default_category=default_category
    )

    # ensure output directory exists
    check_dir_exists(get_dir_name(output_path), create=True)

    # write TSV
    classifications.to_csv(output_path, sep='\t', index=False)

    # success
    return classifications
