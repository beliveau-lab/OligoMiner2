"""
Bundled test data and convenience loaders.

Provides helper functions to load test datasets from package data via
importlib.resources, so that tests and demo scripts can access them
without constructing file paths manually. Similar to sklearn.datasets,
these loaders return ready-to-use objects or paths for the bundled
example genome and annotation files.
"""

from importlib.resources import files


def load_test_probes():
    """
    Load the bundled test probe CSV as a ProbeSet.

    Returns:
        probes (ProbeSet): the test probe set.
    """
    from oligominer.probe_design import ProbeSet

    csv_ref = files("oligominer.data.test_data.probes.csv").joinpath("test_probes.csv")
    probes = ProbeSet.from_csv(str(csv_ref))

    # success
    return probes


def load_example_fasta():
    """
    Load the bundled example FASTA file using pyfaidx.

    Returns:
        fasta (pyfaidx.Fasta): dict-like object (.keys(), [seq_id] -> sequence).
    """
    from oligominer.bioinformatics.file_io import load_fasta

    path = get_example_fasta_path()
    fasta = load_fasta(path)

    # success
    return fasta


def get_example_fasta_path():
    """
    Return the path to the bundled example FASTA file.

    Returns:
        path (str): absolute path to example.fa.
    """
    fasta_ref = files("oligominer.data.test_data.fasta").joinpath("example.fa")
    path = str(fasta_ref)

    # success
    return path


def load_example_gtf():
    """
    Load the bundled example GTF annotation file as a DataFrame.

    Returns:
        df (pandas.DataFrame): one row per annotation record with columns
            seqid, source, type, start, end, score, strand, phase,
            attributes.
    """
    from oligominer.bioinformatics.file_io import load_gtf

    path = get_example_gtf_path()
    df = load_gtf(path)

    # success
    return df


def get_example_gtf_path():
    """
    Return the path to the bundled example GTF annotation file.

    Returns:
        path (str): absolute path to example.gtf.
    """
    gtf_ref = files("oligominer.data.test_data.gtf").joinpath("example.gtf")
    path = str(gtf_ref)

    # success
    return path
