"""
# Probe Design Pipeline

High-level functions that compose the lower-level OligoMiner2 primitives into
a probe design workflow:

  1. Mine candidate probes from a FASTA file
  2. Align probes to a reference genome
  3. Optionally compute per-probe max kmer counts
  4. Merge probe and alignment tables into a duplex table
  5. Optionally compute pDup for each duplex
  6. Optionally predict duplex stability via XGBoost (PaintSHOP model)
"""

import pandas as pd

from oligominer.thermodynamics.mining import mine_fasta, probes_to_df
from oligominer.bioinformatics.file_io import seqs_to_fastq
from oligominer.specificity.alignment import (
    bowtie_align, bowtie_presets, process_alignments
)
from oligominer.specificity.kmers import calc_max_kmer_multi


def _make_seqid(row):
    """Build a seqid string from a probe row."""
    return f"{row['seq_id']}:{row['start']}-{row['stop']}"


def mine_probe_candidates(input_fasta, cores=1, **mining_params):
    """
    Mine candidate probes from a FASTA file.

    Thin wrapper around mine_fasta that returns a DataFrame with an added
    seqid column for downstream merging.

    Args:
        input_fasta (str): path to the input FASTA file.
        cores (int): number of CPU cores for parallel mining.
        **mining_params: forwarded to mine_sequence() (min_length,
            max_length, min_tm, max_tm, etc.).

    Returns:
        probe_df (pandas.DataFrame): probe candidates with columns
            seq_id, start, stop, probe_seq, tm, seqid.
    """
    probe_tuples = mine_fasta(input_fasta, cores=cores, **mining_params)
    probe_df = probes_to_df(probe_tuples)

    # build seqid for linking to alignment results
    probe_df['seqid'] = probe_df.apply(_make_seqid, axis=1)

    # success
    return probe_df


def align_probes(probe_df, bt2_index, ref_fasta, preset=None, k=100,
                 threads=1, verbose=False, **bt2_params):
    """
    Align probe candidates to a reference genome.

    Converts probe sequences to FASTQ, aligns with Bowtie2, and processes
    the results into an alignment DataFrame with derived sequences.

    Args:
        probe_df (pandas.DataFrame): probe candidates as returned by
            mine_probe_candidates(). Must have seqid and probe_seq columns.
        bt2_index (str): path to the Bowtie2 index.
        ref_fasta (str): path to the reference FASTA file.
        preset (dict, optional): Bowtie2 preset parameters. Defaults to
            VERY_SENSITIVE_LOCAL.
        k (int): number of alignments to report per read.
        threads (int): number of Bowtie2 threads.
        verbose (bool): if True, print alignment output.
        **bt2_params: additional parameters forwarded to bowtie_align().

    Returns:
        align_df (pandas.DataFrame): alignment results with columns
            align_seqid, align_start, align_stop, seqid, align_score,
            align_strand, align_cigar, derived_seq.
    """
    if preset is None:
        preset = bowtie_presets.VERY_SENSITIVE_LOCAL

    # convert probes to fastq for alignment
    fastq_data = seqs_to_fastq(
        seq_list=probe_df['probe_seq'],
        seq_id_list=probe_df['seqid']
    )

    # align to reference genome
    sam_data = bowtie_align(
        bt2_index,
        input_data=fastq_data,
        preset=preset,
        k=k,
        threads=threads,
        verbose=verbose,
        **bt2_params
    )

    # process alignments into a dataframe with derived sequences
    align_df = process_alignments(sam_data=sam_data, ref_fasta=ref_fasta)

    # success
    return align_df


def add_max_kmer(probe_df, jf_index, k=18, verbose=False):
    """
    Add a max_kmer column to the probe DataFrame.

    Queries a Jellyfish kmer index to find the maximum kmer count for
    each probe sequence. Lower values indicate higher specificity.

    Args:
        probe_df (pandas.DataFrame): probe candidates with a probe_seq
            column.
        jf_index (str): path to the Jellyfish index file.
        k (int): kmer length. Must match the Jellyfish index.
        verbose (bool): if True, print query output.

    Returns:
        probe_df (pandas.DataFrame): the input DataFrame with an added
            max_kmer column.
    """
    seqs = probe_df['probe_seq'].tolist()
    max_kmer_values = calc_max_kmer_multi(jf_index, seqs, k=k, verbose=verbose)
    probe_df['max_kmer'] = max_kmer_values

    # success
    return probe_df


def merge_probes_alignments(probe_df, align_df):
    """
    Merge probe and alignment tables into a duplex table.

    Each row in the result represents one probe-to-genome alignment (a
    potential duplex). Probes with no alignments are dropped.

    Args:
        probe_df (pandas.DataFrame): probe candidates with a seqid column.
        align_df (pandas.DataFrame): alignment results with a seqid column.

    Returns:
        merged_df (pandas.DataFrame): one row per probe-alignment pair,
            containing all columns from both tables.
    """
    merged_df = probe_df.merge(align_df, on='seqid', how='inner')

    # success
    return merged_df


def add_pdup(merged_df, model=None, conc_a=1e-6, conc_b=1e-12):
    """
    Add a pdup column to the merged duplex DataFrame.

    Computes the duplex formation probability between each probe sequence
    and its derived off-target sequence using NUPACK.

    Args:
        merged_df (pandas.DataFrame): merged probe-alignment table with
            probe_seq and derived_seq columns.
        model (nupack.Model, optional): a NUPACK thermodynamic model.
            Defaults to DEFAULT_NUPACK_MODEL.
        conc_a (float): molar concentration of the probe strand.
        conc_b (float): molar concentration of the target strand.

    Returns:
        merged_df (pandas.DataFrame): the input DataFrame with an added
            pdup column.
    """
    from oligominer.thermodynamics.nupack import calc_pdup

    pdup_values = []
    for _, row in merged_df.iterrows():
        pdup = calc_pdup(
            row['probe_seq'], row['derived_seq'],
            conc_a=conc_a, conc_b=conc_b, model=model
        )
        pdup_values.append(pdup)

    merged_df['pdup'] = pdup_values

    # success
    return merged_df


def add_duplex_pred(merged_df, temperature=37, normalize=True):
    """
    Add a duplex_pred column to the merged duplex DataFrame.

    Uses a pre-trained PaintSHOP XGBoost model to predict duplex formation
    probability, providing a fast approximation of NUPACK pDup without
    requiring the NUPACK dependency.

    Args:
        merged_df (pandas.DataFrame): merged probe-alignment table with
            probe_seq, derived_seq, and align_score columns.
        temperature (int): hybridization temperature in °C. Must be one of
            37, 42, 47, 52, or 60.
        normalize (bool): if True, scale predictions to 0.0-1.0 to match
            calc_pdup output. If False, predictions are on the native
            0-100 scale.

    Returns:
        merged_df (pandas.DataFrame): the input DataFrame with an added
            duplex_pred column.
    """
    from oligominer.specificity.duplex_stability import predict_duplex_batch

    merged_df['duplex_pred'] = predict_duplex_batch(
        merged_df, temperature=temperature, normalize=normalize
    )

    # success
    return merged_df


def design_probes(input_fasta, bt2_index, ref_fasta,
                  jf_index=None, compute_pdup=False,
                  compute_duplex_pred=False,
                  preset=None, k_align=100, threads=1,
                  jf_k=18, cores=1, verbose=False,
                  nupack_model=None, conc_a=1e-6, conc_b=1e-12,
                  duplex_pred_temperature=37,
                  **mining_params):
    """
    End-to-end probe design pipeline.

    Mines candidate probes, aligns them to a reference genome, and
    optionally computes max kmer counts and duplex formation probabilities.

    Args:
        input_fasta (str): path to the input FASTA file for probe mining.
        bt2_index (str): path to the Bowtie2 index for alignment.
        ref_fasta (str): path to the reference FASTA file.
        jf_index (str, optional): path to a Jellyfish index. If provided,
            max_kmer values are computed for each probe.
        compute_pdup (bool): if True, compute pDup for each duplex using
            NUPACK.
        compute_duplex_pred (bool): if True, predict duplex stability using
            the PaintSHOP XGBoost model. Does not require NUPACK.
        preset (dict, optional): Bowtie2 preset parameters.
        k_align (int): number of alignments to report per probe.
        threads (int): number of Bowtie2 threads.
        jf_k (int): kmer length for Jellyfish queries.
        cores (int): number of CPU cores for probe mining.
        verbose (bool): if True, print progress output.
        nupack_model (nupack.Model, optional): NUPACK model for pDup.
        conc_a (float): probe strand concentration for pDup.
        conc_b (float): target strand concentration for pDup.
        duplex_pred_temperature (int): temperature in °C for the XGBoost
            model. Must be one of 37, 42, 47, 52, or 60.
        **mining_params: forwarded to mine_sequence().

    Returns:
        result (dict): a dictionary with keys:
            probe_df (pandas.DataFrame): probe candidates.
            align_df (pandas.DataFrame): alignment results.
            merged_df (pandas.DataFrame): merged duplex table.
    """
    # step 1: mine probe candidates
    probe_df = mine_probe_candidates(
        input_fasta, cores=cores, **mining_params
    )

    # step 2: align probes to reference genome
    align_df = align_probes(
        probe_df, bt2_index, ref_fasta,
        preset=preset, k=k_align, threads=threads, verbose=verbose
    )

    # step 3: optionally compute max kmer counts
    if jf_index is not None:
        probe_df = add_max_kmer(
            probe_df, jf_index, k=jf_k, verbose=verbose
        )

    # step 4: merge probe and alignment tables
    merged_df = merge_probes_alignments(probe_df, align_df)

    # step 5: optionally compute pDup
    if compute_pdup:
        merged_df = add_pdup(
            merged_df, model=nupack_model,
            conc_a=conc_a, conc_b=conc_b
        )

    # step 6: optionally predict duplex stability via xgboost
    if compute_duplex_pred:
        merged_df = add_duplex_pred(
            merged_df, temperature=duplex_pred_temperature
        )

    result = {
        'probe_df': probe_df,
        'align_df': align_df,
        'merged_df': merged_df,
    }

    # success
    return result
