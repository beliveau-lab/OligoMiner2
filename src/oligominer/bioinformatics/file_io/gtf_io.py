"""
GTF/GFF annotation file I/O utilities.

Loading parses standard GTF/GFF format into a pandas DataFrame. Filtering,
writing, splitting, and merging mirror the structure of fasta_io for a
consistent preprocessing interface.
"""

import os

import pandas as pd

from oligominer.utils import (
    get_abs_path, get_dir_name,
    check_dir_exists, check_input_exists, check_output_exists,
)
from .exceptions import EmptyExportError
from .config import GTF_EXTENSIONS, GTF_COLUMNS, merge_files_by_extension


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

def load_gtf(gtf_path):
    """
    Load a GTF/GFF annotation file into a pandas DataFrame.

    Skips comment lines (starting with '#') and parses the standard
    9-column GTF format. The attributes column is left unparsed; use
    parse_attributes to expand it.

    Args:
        gtf_path (str): path to the input GTF/GFF file.

    Returns:
        df (pandas.DataFrame): one row per annotation record with columns
            seqid, source, type, start, end, score, strand, phase, attributes.
    """
    gtf_path = get_abs_path(gtf_path)
    check_input_exists(gtf_path)

    df = pd.read_csv(
        gtf_path,
        sep='\t',
        comment='#',
        header=None,
        names=GTF_COLUMNS,
        dtype={'seqid': str, 'start': int, 'end': int},
    )

    # success
    return df


# ---------------------------------------------------------------------------
# Attribute parsing
# ---------------------------------------------------------------------------

def parse_attributes(df):
    """
    Expand the GTF attributes column into individual DataFrame columns.

    Parses the semicolon-delimited key-value pairs in the attributes field
    and adds gene_id, gene_name, transcript_id, and transcript_id_full
    columns.  The raw attributes column is dropped.

    Args:
        df (pandas.DataFrame): a GTF DataFrame with an 'attributes' column.

    Returns:
        df (pandas.DataFrame): the input DataFrame with attributes replaced
            by gene_id, gene_name, transcript_id (version stripped), and
            transcript_id_full (original versioned ID) columns.
    """
    df = df.copy()

    # parse each attribute string into a dict
    attr_dicts = df['attributes'].apply(_parse_attr_string)

    # extract commonly used fields
    df['gene_id'] = attr_dicts.apply(lambda x: x.get('gene_id', ''))
    df['gene_name'] = attr_dicts.apply(lambda x: x.get('gene_name', ''))
    df['transcript_id_full'] = attr_dicts.apply(
        lambda x: x.get('transcript_id', '')
    )
    df['transcript_id'] = df['transcript_id_full'].apply(
        lambda x: x.split('.')[0]
    )

    # drop the raw attributes column
    df = df.drop(columns=['attributes'])

    # success
    return df


# ---------------------------------------------------------------------------
# Filtering
# ---------------------------------------------------------------------------

def filter_gtf(df, chrom_names=None, feature_type='exon',
               incl_str=None, excl_str=None):
    """
    Filter a GTF DataFrame by chromosome names, feature type, and regex.

    When multiple filters are provided they are applied in order:
    chrom_names, then feature_type, then incl_str, then excl_str.

    Args:
        df (pandas.DataFrame): a GTF DataFrame (see load_gtf).
        chrom_names (list or None): keep only records whose seqid is in
            this list. None skips this filter.
        feature_type (str or None): keep only records whose type matches
            this value (e.g. 'exon'). None skips this filter.
        incl_str (str or None): regex — keep seqids that match.
        excl_str (str or None): regex — drop seqids that match.

    Returns:
        filtered (pandas.DataFrame): the filtered annotation records.
    """
    filtered = df

    if chrom_names is not None:
        filtered = filtered[filtered['seqid'].isin(chrom_names)]

    if feature_type is not None:
        filtered = filtered[filtered['type'] == feature_type]

    if incl_str is not None:
        filtered = filtered[filtered['seqid'].str.contains(incl_str, regex=True)]

    if excl_str is not None:
        filtered = filtered[~filtered['seqid'].str.contains(excl_str, regex=True)]

    # success
    return filtered.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Writing
# ---------------------------------------------------------------------------

def write_gtf(df, filepath):
    """
    Write a GTF DataFrame to a tab-separated file.

    If the DataFrame contains parsed attribute columns (gene_id,
    transcript_id, etc.) but no 'attributes' column, those columns are
    written as-is.  Otherwise the standard GTF 9-column format is used.

    Args:
        df (pandas.DataFrame): annotation records to write.
        filepath (str): destination file path.
    """
    if df.empty:
        raise EmptyExportError(filepath, 'GTF')

    check_dir_exists(filepath, parent_dir=True, create=True)
    df.to_csv(filepath, sep='\t', index=False, header=False)


def write_bed(df, filepath):
    """
    Write a GTF DataFrame as a BED file.

    Selects BED-compatible columns (seqid, start, end, transcript_id,
    score, strand, transcript_id_full, gene_id) and writes them
    tab-separated with no header.

    Args:
        df (pandas.DataFrame): annotation records with parsed attributes.
        filepath (str): destination file path.
    """
    if df.empty:
        raise EmptyExportError(filepath, 'GTF')

    bed_columns = [
        'seqid', 'start', 'end', 'transcript_id',
        'score', 'strand', 'transcript_id_full', 'gene_id',
    ]

    # only include columns that are present
    available = [c for c in bed_columns if c in df.columns]
    bed_df = df[available]

    check_dir_exists(filepath, parent_dir=True, create=True)
    bed_df.to_csv(filepath, sep='\t', index=False, header=False)


# ---------------------------------------------------------------------------
# Splitting
# ---------------------------------------------------------------------------

def split_gtf(df, target_dir, suffix='_filtered_gtf.tsv'):
    """
    Split a GTF DataFrame into per-chromosome files.

    Each chromosome's records are written to a separate TSV file named
    <seqid><suffix> inside target_dir. This is analogous to split_fasta
    for sequence files.

    Args:
        df (pandas.DataFrame): annotation records with a 'seqid' column.
        target_dir (str): output directory (created if needed).
        suffix (str): filename suffix appended to each chromosome name.

    Returns:
        written_paths (list): file paths written.
    """
    if df.empty:
        raise EmptyExportError(target_dir, 'GTF')

    check_dir_exists(target_dir, create=True)

    written_paths = []
    for chrom, chrom_df in df.groupby('seqid'):
        out_path = os.path.join(target_dir, f'{chrom}{suffix}')
        chrom_df.to_csv(out_path, sep='\t', index=False)
        written_paths.append(out_path)

    # success
    return written_paths


# ---------------------------------------------------------------------------
# Merging
# ---------------------------------------------------------------------------

def merge_annotation_beds(input_dir, output_path, extension='.bed'):
    """
    Merge per-chromosome annotation BED files into a single file.

    Concatenates all files matching the given extension inside input_dir,
    sorted alphabetically by filename. This is the annotation-side
    counterpart of merge_fastas for sequence files, intended for
    reassembling per-chromosome BED outputs (e.g. from split_gtf or
    isoform flattening) into a genome-wide annotation file.

    Args:
        input_dir (str): directory containing per-chromosome BED files.
        output_path (str): path for the merged output file.
        extension (str): file extension to match.

    Returns:
        output_path (str): path to the created file.
    """
    input_dir = get_abs_path(input_dir)
    check_dir_exists(input_dir)
    check_dir_exists(get_dir_name(output_path), create=True)

    merged = merge_files_by_extension(input_dir, output_path, {extension})
    if not merged:
        raise EmptyExportError(output_path, 'GTF')

    check_output_exists(output_path)

    # success
    return output_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_attr_string(data):
    """
    Parse a single GTF attributes string into a dict.

    Handles both GTF (space-separated key "value") and GFF3 (key=value)
    attribute formats.

    Args:
        data (str): the raw attributes string.

    Returns:
        attr_data (dict): {key: value} pairs.
    """
    attr_data = {}
    cleaned = data.replace('"', '').strip(' ;')

    for field in cleaned.split(';'):
        field = field.strip()
        if not field:
            continue

        # handle GFF3 key=value format
        if '=' in field:
            key, _, val = field.partition('=')
        else:
            # GTF key "value" format (space-separated)
            parts = field.split(' ', 1)
            key = parts[0]
            val = parts[1] if len(parts) > 1 else ''

        attr_data[key.strip()] = val.strip()

    # success
    return attr_data
