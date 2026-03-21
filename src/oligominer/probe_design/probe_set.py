"""
# ProbeSet

A convenience class that wraps a probe DataFrame and provides methods for
the common probe design workflow: alignment, kmer filtering, merging with
alignment results, pDup computation, sequence appending, and scoring.

A ProbeSet can be created from mining result tuples, loaded from a CSV on
disk, constructed directly from a DataFrame, or mined from a FASTA file.
"""

import pandas as pd

from oligominer.thermodynamics.mining import mine_fasta, probes_to_df, write_probes, PROBE_COLUMNS
from oligominer.utils.exceptions import PipelineStateError

from .pipeline import (
    _make_seqid,
    align_probes,
    add_max_kmer,
    merge_probes_alignments,
    add_pdup,
    add_duplex_pred,
)
from .probe_io import (
    read_probe_csv,
    write_probe_csv,
    read_align_csv,
    write_align_csv,
)
from .appending import (
    append_sequences,
    append_saber,
    append_barcodes as _append_barcodes,
    build_appending_table,
)
from .scoring import label_on_target, score_probes


class ProbeSet:
    """
    A set of candidate probes backed by a pandas DataFrame.

    Provides a fluent interface for the probe design pipeline steps:
    mining, alignment, kmer counting, merging, and pDup computation.

    Can be created from:
      - A list of probe tuples from mine_sequence() or mine_fasta()
      - A pandas DataFrame with probe columns
      - A FASTA file via ProbeSet.from_fasta()
      - A CSV file via ProbeSet.from_csv()

    Attributes:
        df (pandas.DataFrame): probe candidates with columns
            seq_id, start, stop, probe_seq, tm, seqid, and optionally
            max_kmer.
        align_df (pandas.DataFrame or None): alignment results, populated
            after calling align().
        merged_df (pandas.DataFrame or None): merged probe-alignment
            duplex table, populated after calling merge().
        score_df (pandas.DataFrame or None): per-probe on-target and
            off-target scores, populated after calling score().
    """

    def __init__(self, probes):
        """
        Create a ProbeSet from probe tuples or a probe DataFrame.

        When given a list of tuples (as returned by mine_sequence or
        mine_fasta), converts them to a DataFrame via probes_to_df().
        When given a DataFrame, uses it directly.

        If the DataFrame does not already contain a seqid column, one is
        built from the seq_id, start, and stop columns.

        Args:
            probes (list or pandas.DataFrame): either a list of
                (seq_id, start, stop, probe_seq, tm) tuples, or a
                DataFrame with at least those columns.
        """
        if isinstance(probes, pd.DataFrame):
            self.df = probes.copy()
        else:
            self.df = probes_to_df(probes)

        # ensure seqid column exists
        if 'seqid' not in self.df.columns:
            self.df['seqid'] = self.df.apply(_make_seqid, axis=1)

        self.align_df = None
        self.merged_df = None
        self.score_df = None
        self._master_entries = {}

    # ------------------------------------------------------------------
    # constructors
    # ------------------------------------------------------------------

    @classmethod
    def from_fasta(cls, input_fasta, cores=1, **mining_params):
        """
        Mine probes from a FASTA file and return a ProbeSet.

        Args:
            input_fasta (str): path to the input FASTA file.
            cores (int): number of CPU cores for parallel mining.
            **mining_params: forwarded to mine_fasta() / mine_sequence()
                (min_length, max_length, min_tm, max_tm, Na, Mg, etc.).

        Returns:
            probe_set (ProbeSet): the mined probe set.
        """
        probes = mine_fasta(input_fasta, cores=cores, **mining_params)

        # success
        return cls(probes)

    @classmethod
    def from_csv(cls, path):
        """
        Load a ProbeSet from a probe CSV file.

        Args:
            path (str): path to the probe CSV file.

        Returns:
            probe_set (ProbeSet): the loaded probe set.
        """
        probe_df = read_probe_csv(path)

        # success
        return cls(probe_df)

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------

    def to_csv(self, path):
        """
        Write the probe DataFrame to a CSV file.

        Args:
            path (str): output file path.
        """
        write_probe_csv(self.df, path)

    def to_bed(self, path):
        """
        Write probes to a BED file.

        Args:
            path (str): output file path.
        """
        tuples = self._to_tuples()
        write_probes(tuples, path, fmt='bed')

    def to_fastq(self, path):
        """
        Write probes to a FASTQ file.

        Args:
            path (str): output file path.
        """
        tuples = self._to_tuples()
        write_probes(tuples, path, fmt='fastq')

    def to_fasta(self, path):
        """
        Write probe sequences to a FASTA file.

        Each probe is written as a separate record with the seqid as the
        header and the probe sequence as the body.

        Args:
            path (str): output file path.
        """
        from oligominer.bioinformatics.file_io import write_fasta

        seqs = dict(zip(self.df['seqid'], self.df['probe_seq']))
        write_fasta(seqs, path)

    def _to_tuples(self):
        """Convert the probe DataFrame back to a list of tuples.

        Returns:
            tuples (list): list of (seq_id, start, stop, probe_seq, tm)
                tuples.
        """
        tuples = list(self.df[PROBE_COLUMNS[:5]].itertuples(index=False, name=None))

        # success
        return tuples

    def save_align_csv(self, path):
        """
        Write the alignment DataFrame to a CSV file.

        Args:
            path (str): output file path.

        Raises:
            PipelineStateError: if align() has not been called yet.
        """
        if self.align_df is None:
            raise PipelineStateError("No alignment data. Call align() first.")
        write_align_csv(self.align_df, path)

    def load_align_csv(self, path):
        """
        Load alignment results from a CSV file.

        Args:
            path (str): path to the alignment CSV file.

        Returns:
            self (ProbeSet): for method chaining.
        """
        self.align_df = read_align_csv(path)

        # success
        return self

    def save_merged_csv(self, path):
        """
        Write the merged duplex DataFrame to a CSV file.

        Args:
            path (str): output file path.

        Raises:
            PipelineStateError: if merge() has not been called yet.
        """
        if self.merged_df is None:
            raise PipelineStateError("No merged data. Call merge() first.")
        write_align_csv(self.merged_df, path)

    # ------------------------------------------------------------------
    # pipeline steps
    # ------------------------------------------------------------------

    def align(self, bt2_index, ref_fasta, preset=None, k=100,
              threads=1, verbose=False, **bt2_params):
        """
        Align probes to a reference genome.

        Populates self.align_df with alignment results.

        Args:
            bt2_index (str): path to the Bowtie2 index.
            ref_fasta (str): path to the reference FASTA file.
            preset (dict, optional): Bowtie2 preset parameters.
            k (int): number of alignments to report per read.
            threads (int): number of Bowtie2 threads.
            verbose (bool): if True, print alignment output.
            **bt2_params: additional parameters forwarded to bowtie_align().

        Returns:
            self (ProbeSet): for method chaining.
        """
        self.align_df = align_probes(
            self.df, bt2_index, ref_fasta,
            preset=preset, k=k, threads=threads, verbose=verbose,
            **bt2_params
        )

        # success
        return self

    def compute_max_kmer(self, jf_index, k=18, verbose=False):
        """
        Compute max kmer counts for each probe.

        Adds a max_kmer column to self.df.

        Args:
            jf_index (str): path to the Jellyfish index file.
            k (int): kmer length.
            verbose (bool): if True, print query output.

        Returns:
            self (ProbeSet): for method chaining.
        """
        self.df = add_max_kmer(
            self.df, jf_index, k=k, verbose=verbose
        )

        # success
        return self

    def merge(self):
        """
        Merge probe and alignment tables into a duplex table.

        Requires that align() (or load_align_csv()) has been called.
        Populates self.merged_df.

        Returns:
            self (ProbeSet): for method chaining.

        Raises:
            PipelineStateError: if no alignment data is available.
        """
        if self.align_df is None:
            raise PipelineStateError("No alignment data. Call align() first.")

        self.merged_df = merge_probes_alignments(self.df, self.align_df)

        # success
        return self

    def compute_pdup(self, model=None, conc_a=1e-6, conc_b=1e-12):
        """
        Compute pDup for each duplex in the merged table.

        Requires that merge() has been called. Adds a pdup column to
        self.merged_df.

        Args:
            model (nupack.Model, optional): a NUPACK thermodynamic model.
            conc_a (float): molar concentration of the probe strand.
            conc_b (float): molar concentration of the target strand.

        Returns:
            self (ProbeSet): for method chaining.

        Raises:
            PipelineStateError: if no merged data is available.
        """
        if self.merged_df is None:
            raise PipelineStateError("No merged data. Call merge() first.")

        self.merged_df = add_pdup(
            self.merged_df, model=model, conc_a=conc_a, conc_b=conc_b
        )

        # success
        return self

    def compute_duplex_pred(self, temperature=37, normalize=True):
        """
        Predict duplex stability for each duplex using the PaintSHOP XGBoost
        model.

        Requires that merge() has been called. Adds a duplex_pred column to
        self.merged_df. This is a fast alternative to compute_pdup() that
        does not require NUPACK.

        Args:
            temperature (int): hybridization temperature in °C. Must be one
                of 37, 42, 47, 52, or 60.
            normalize (bool): if True, scale predictions to 0.0-1.0 to match
                pdup output. If False, predictions are on the native 0-100
                scale.

        Returns:
            self (ProbeSet): for method chaining.

        Raises:
            PipelineStateError: if no merged data is available.
        """
        if self.merged_df is None:
            raise PipelineStateError("No merged data. Call merge() first.")

        self.merged_df = add_duplex_pred(
            self.merged_df, temperature=temperature, normalize=normalize
        )

        # success
        return self

    # ------------------------------------------------------------------
    # appending
    # ------------------------------------------------------------------

    def _ensure_sequence_column(self):
        """Initialize the sequence column from probe_seq if not present."""
        if 'sequence' not in self.df.columns:
            self.df['sequence'] = self.df['probe_seq']

    def append(self, sequences, scheme, label, target_column=None,
               n_per_target=None, ranges=None, left=True, rc=False):
        """
        Append sequences to probes using the specified scheme.

        Modifies the ``sequence`` column in df (creating it from
        ``probe_seq`` if it does not exist). Tracks appended sequences
        in the internal master table.

        Args:
            sequences (pandas.DataFrame): appending sequences with
                ``id`` and ``seq`` columns.
            scheme (str): one of ``"same"``, ``"unique"``,
                ``"multiple"``, or ``"custom"``.
            label (str): name for this appending step in the master
                table (e.g. ``"outer_forward"``).
            target_column (str, optional): column for grouping by
                target.
            n_per_target (int, optional): sequences per target for
                ``"multiple"`` scheme.
            ranges (list of str, optional): ranges for ``"custom"``.
            left (bool): if True, prepend to the 5' end.
            rc (bool): if True, reverse-complement sequences.

        Returns:
            self (ProbeSet): for method chaining.
        """
        self._ensure_sequence_column()

        self.df, entries = append_sequences(
            self.df, sequences, scheme,
            target_column=target_column,
            n_per_target=n_per_target,
            ranges=ranges,
            left=left, rc=rc,
        )
        self._master_entries[label] = entries

        # success
        return self

    def append_saber_seqs(self, sequences, scheme, label="saber",
                          target_column=None, n_per_target=None,
                          ranges=None):
        """
        Append SABER concatemer sequences to the 3' end.

        Args:
            sequences (pandas.DataFrame): SABER sequences with ``id``
                and ``seq`` columns.
            scheme (str): appending scheme.
            label (str): master table column name.
            target_column (str, optional): column for grouping by
                target.
            n_per_target (int, optional): sequences per target for
                ``"multiple"`` scheme.
            ranges (list of str, optional): ranges for ``"custom"``.

        Returns:
            self (ProbeSet): for method chaining.
        """
        self._ensure_sequence_column()

        self.df, entries = append_saber(
            self.df, sequences, scheme,
            target_column=target_column,
            n_per_target=n_per_target,
            ranges=ranges,
        )
        self._master_entries[label] = entries

        # success
        return self

    def append_merfish_barcodes(self, bridges, barcodes,
                                target_column="refseq"):
        """
        Append MERFISH barcode-encoded bridges.

        Args:
            bridges (pandas.DataFrame): the 16 MERFISH bridges with
                ``id`` and ``seq`` columns.
            barcodes (pandas.DataFrame): barcodes with a ``barcode``
                column.
            target_column (str): column for grouping by target.

        Returns:
            self (ProbeSet): for method chaining.
        """
        self._ensure_sequence_column()

        self.df = _append_barcodes(
            self.df, bridges, barcodes,
            target_column=target_column,
        )

        # success
        return self

    @property
    def master_table(self):
        """
        Return the master table tracking all appending steps.

        Returns:
            master (pandas.DataFrame or None): one row per probe, one
                column per appending step, or None if no appending has
                been performed.
        """
        if not self._master_entries:
            return None
        return build_appending_table(self.df, self._master_entries)

    # ------------------------------------------------------------------
    # scoring
    # ------------------------------------------------------------------

    def score(self, pred_column="duplex_pred"):
        """
        Compute on-target and off-target scores from duplex predictions.

        Requires that merge() and either compute_pdup() or
        compute_duplex_pred() have been called. Populates self.score_df.

        Args:
            pred_column (str): name of the prediction column in
                merged_df.

        Returns:
            self (ProbeSet): for method chaining.

        Raises:
            PipelineStateError: if no merged data is available.
        """
        if self.merged_df is None:
            raise PipelineStateError("No merged data. Call merge() first.")

        labeled = label_on_target(self.merged_df)
        self.score_df = score_probes(labeled, pred_column=pred_column)

        # success
        return self

    # ------------------------------------------------------------------
    # display
    # ------------------------------------------------------------------

    def __len__(self):
        return len(self.df)

    def __repr__(self):
        n_probes = len(self.df)
        n_align = len(self.align_df) if self.align_df is not None else 0
        n_merged = len(self.merged_df) if self.merged_df is not None else 0
        return (
            f"ProbeSet(probes={n_probes}, "
            f"alignments={n_align}, "
            f"duplexes={n_merged})"
        )
