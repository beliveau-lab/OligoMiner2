"""mine_probes.py — Oligonucleotide probe mining engine.

Provides the core probe mining pipeline for OligoMiner2. The central function,
mine_sequence(), takes a DNA sequence string and returns candidate probes as a
list of tuples. It is a pure compute function with no file I/O dependency.

The pipeline for each sequence:
  1. Integer-encode the DNA string (int_encoding.py)
  2. Split into overlapping chunks for memory-efficient processing
  3. For each chunk, compute a 2D Tm grid (calc_tm_2d.py) covering all
     positions x probe lengths simultaneously
  4. Apply vectorized filters: Tm range, N-masking, homopolymer runs,
     prohibited subsequences, and GC% — all using cumulative-sum tricks
     for O(1) range queries
  5. Select the best probe length per position (greedy shortest or
     closest-to-target Tm)
  6. Apply overlap and spacing constraints

Higher-level functions compose mine_sequence() with file I/O:
  - mine_fasta(): load a FASTA file, mine all sequences, return probe tuples
  - write_probes(): write probe tuples to BED, FASTQ, or CSV
  - probes_to_df(): convert probe tuples to a pandas DataFrame
"""

import multiprocessing

import numpy as np

from .int_encoding import seq_to_8bit
from .calc_tm_2d import get_tm_grid
from .config import GET_DEFAULT_MINING_CONFIG, WRITE_BUFFER_SIZE
from oligominer.utils import check_dir_exists
from oligominer.utils.exceptions import ConfigurationError

# Column order for probe result tuples: (seq_id, start, stop, probe_seq, tm)
PROBE_COLUMNS = ['seq_id', 'start', 'stop', 'probe_seq', 'tm']


def mine_sequence(
    seq,
    seq_id='seq',
    min_length=30,
    max_length=37,
    min_tm=42,
    max_tm=47,
    tm_target=None,
    chunk_size=100000,
    overlap=True,
    spacing=0,
    exhaustive=False,
    max_homopolymer=4,
    prohibited_seqs=None,
    cores=1,
    Na=390,
    K=0,
    Tris=0,
    Mg=0,
    dNTPs=0,
    dnac1=25,
    dnac2=25,
    pct_formamide=50,
    formamide_factor=0.65,
):
    """Mine candidate oligonucleotide probes from a DNA sequence.

    Pure compute — takes a sequence string, returns probe results as a list
    of tuples. No file I/O.

    Args:
        seq (str): DNA sequence string (ACGTN).
        seq_id (str): identifier for this sequence (used in output).
        min_length (int): minimum probe length.
        max_length (int): maximum probe length.
        min_tm (float): minimum melting temperature (°C).
        max_tm (float): maximum melting temperature (°C).
        tm_target (float or None): target Tm for closest-to-target selection.
            None = greedy shortest valid probe. Ignored when exhaustive=True.
        chunk_size (int): bases per processing chunk.
        overlap (bool): allow overlapping probes. Mutually exclusive with
            exhaustive.
        spacing (int): minimum spacing between adjacent probes. Mutually
            exclusive with exhaustive.
        exhaustive (bool): return all valid (position, length) combinations
            instead of selecting one probe per position. Mutually exclusive
            with overlap and spacing.
        max_homopolymer (int or None): max homopolymer run. None to disable.
        prohibited_seqs (list or None): list of subsequence strings to exclude.
        cores (int): number of CPU cores for parallel chunk processing.
        Na (float): sodium concentration in mM.
        K (float): potassium concentration in mM.
        Tris (float): Tris buffer concentration in mM.
        Mg (float): magnesium concentration in mM.
        dNTPs (float): dNTP concentration in mM. dNTPs chelate Mg2+,
            reducing the effective free Mg2+ concentration.
        dnac1 (float): concentration of the probe strand in nM.
        dnac2 (float): concentration of the target strand in nM.
        pct_formamide (int): percent formamide in the hybridization buffer
            (e.g. 50 for 50%). Set to 0 to disable formamide correction.
        formamide_factor (float): degrees C of Tm depression per percent
            formamide.

    Returns:
        probes (list): list of tuples (seq_id, start, stop, probe_seq, tm).
            Use probes_to_df() to convert to a pandas DataFrame.
    """

    if exhaustive and (not overlap or spacing > 0):
        raise ConfigurationError(
            "exhaustive mode is incompatible with overlap=False and "
            "spacing > 0. In exhaustive mode all valid probes are returned.")

    config = GET_DEFAULT_MINING_CONFIG()
    config['min_length'] = min_length
    config['max_length'] = max_length
    config['min_tm'] = min_tm
    config['max_tm'] = max_tm
    config['tm_target'] = tm_target
    config['chunk_size'] = chunk_size
    config['exhaustive'] = exhaustive
    config['max_homopolymer'] = max_homopolymer
    config['prohibited_seqs'] = prohibited_seqs
    config['Na'] = Na
    config['K'] = K
    config['Tris'] = Tris
    config['Mg'] = Mg
    config['dNTPs'] = dNTPs
    config['dnac1'] = dnac1
    config['dnac2'] = dnac2
    config['pct_formamide'] = pct_formamide
    config['formamide_factor'] = formamide_factor
    if prohibited_seqs:
        config['_prohibited_encoded'] = [seq_to_8bit(p) for p in prohibited_seqs]

    seq_str = seq.upper()
    nuc_array = seq_to_8bit(seq_str)
    chunks = chunk_generator(seq_id, nuc_array, config)

    if cores > 1:
        with multiprocessing.Pool(processes=cores, maxtasksperchild=64) as pool:
            chunk_results = pool.imap(_process_chunk_packed, chunks, chunksize=4)
            rows = _collect_probes(seq_id, seq_str, chunk_results, overlap, spacing)
    else:
        chunk_results = (process_chunk(*c) for c in chunks)
        rows = _collect_probes(seq_id, seq_str, chunk_results, overlap, spacing)

    return rows


def mine_fasta(input_file, cores=1, **mining_params):
    """Mine probes from all sequences in a FASTA file.

    Loads the FASTA via pyfaidx and calls mine_sequence() per sequence.
    Returns raw probe tuples — use ProbeSet() to convert to a DataFrame
    or to export results to disk.

    Args:
        input_file (str): path to input FASTA file.
        cores (int): number of CPU cores (passed to mine_sequence).
        **mining_params: forwarded to mine_sequence(). See mine_sequence()
            for the full list of parameters (min_length, max_length,
            min_tm, max_tm, Na, Mg, etc.).

    Returns:
        probes (list): list of tuples (seq_id, start, stop, probe_seq, tm).
            Pass to ProbeSet() to get a DataFrame or use export methods.
    """
    from oligominer.bioinformatics.file_io import load_fasta

    fasta = load_fasta(input_file)

    all_probes = []
    for seq_id in fasta.keys():
        seq_str = str(fasta[seq_id])
        all_probes.extend(
            mine_sequence(seq_str, seq_id=seq_id, cores=cores, **mining_params)
        )

    # success
    return all_probes


def write_probes(probes, output_file, fmt='bed'):
    """Write probe results to a file.

    Args:
        probes (list): list of (seq_id, start, stop, probe_seq, tm) tuples.
        output_file (str): destination file path.
        fmt (str): output format — 'bed', 'fastq', or 'csv'.
    """
    check_dir_exists(output_file, parent_dir=True, create=True)

    write_fn = _WRITERS[fmt]
    with open(output_file, 'w', buffering=WRITE_BUFFER_SIZE) as outfile:
        write_fn(probes, outfile)


def probes_to_df(probes):
    """Convert probe results (list of tuples) to a pandas DataFrame.

    Args:
        probes (list): list of (seq_id, start, stop, probe_seq, tm) tuples,
            as returned by mine_sequence().

    Returns:
        df (pandas.DataFrame): DataFrame with columns: seq_id, start, stop,
            probe_seq, tm.
    """
    import pandas as pd
    return pd.DataFrame(probes, columns=PROBE_COLUMNS)


# ---------------------------------------------------------------------------
# Probe file writers (internal)
# ---------------------------------------------------------------------------

def _write_probes_bed(probes, outfile):
    """Write probe tuples as BED records."""
    lines = []
    for seq_id, start, stop, probe_seq, tm in probes:
        lines.append(f'{seq_id}\t{start}\t{stop}\t{probe_seq}\t{tm:.2f}\n')
    outfile.writelines(lines)


def _write_probes_fastq(probes, outfile):
    """Write probe tuples as FASTQ records."""
    lines = []
    for seq_id, start, stop, probe_seq, tm in probes:
        lines.append(
            f"@{seq_id}:{start}-{stop}\n"
            f"{probe_seq}\n+\n{'~' * len(probe_seq)}\n"
        )
    outfile.writelines(lines)


def _write_probes_csv(probes, outfile):
    """Write probe tuples as CSV records."""
    outfile.write(','.join(PROBE_COLUMNS) + '\n')
    lines = []
    for seq_id, start, stop, probe_seq, tm in probes:
        lines.append(f'{seq_id},{start},{stop},{probe_seq},{tm:.2f}\n')
    outfile.writelines(lines)


_WRITERS = {
    'bed': _write_probes_bed,
    'fastq': _write_probes_fastq,
    'csv': _write_probes_csv,
}


def _collect_probes(seq_id, seq_str, chunk_results, allow_overlap, spacing):
    """Collect probe results from chunk processing into a list of row tuples.

    Args:
        seq_id (str): sequence identifier.
        seq_str (str): the full sequence string for subsequence extraction.
        chunk_results (iterable): iterable of (coord_result, tm_result) tuples
            from process_chunk().
        allow_overlap (bool): whether overlapping probes are permitted.
        spacing (int): minimum bases between adjacent probes.

    Returns:
        rows (list): list of (seq_id, start, stop, probe_seq, tm) tuples.
    """
    rows = []
    current_probe_stop = -1
    for coord_result, tm_result in chunk_results:
        starts = coord_result[:, 0].tolist()
        stops  = coord_result[:, 1].tolist()
        tms    = tm_result.tolist()
        for probe_start, probe_stop, tm in zip(starts, stops, tms):
            if not allow_overlap and probe_start < current_probe_stop:
                continue
            if spacing > 0 and probe_start < current_probe_stop + spacing:
                continue
            current_probe_stop = probe_stop
            probe_seq = seq_str[probe_start:probe_stop]
            rows.append((seq_id, probe_start, probe_stop, probe_seq, round(tm, 2)))
    return rows


def _process_chunk_packed(chunk):
    """Top-level wrapper so pool.imap can pickle it."""
    return process_chunk(*chunk)


def chunk_generator(seq_id, nuc_array, config):
    """Yield overlapping chunks of an encoded sequence for parallel processing.

    Chunks overlap by max_length - 1 bases so that probes spanning chunk
    boundaries are not missed.

    Args:
        seq_id (str): sequence identifier (passed through to chunk tuples).
        nuc_array (numpy.ndarray): 1D uint8 array of integer-encoded bases.
        config (dict): mining configuration dict (needs 'chunk_size' and
            'max_length' keys).

    Yields:
        chunk (tuple): (seq_id, chunk_nuc_arr, chunk_start, chunk_stop, config)
            for each chunk.
    """

    num_chunks = np.ceil(nuc_array.size / config['chunk_size']).astype(int)

    chunk_start = 0
    chunk_stop = chunk_start + config['chunk_size']
    for _ in range(num_chunks):
        # .copy() so the full nuc_array can be freed after chunking
        chunk = (seq_id, nuc_array[chunk_start:chunk_stop].copy(), chunk_start, chunk_stop, config)
        chunk_start = chunk_stop - config['max_length'] + 1
        chunk_stop = min(chunk_stop + config['chunk_size'], nuc_array.size)
        yield chunk


def process_chunk(_seq_id, chunk_nuc_arr, chunk_start, _chunk_stop, config):
    """Process a single chunk: compute Tm grid, apply filters, select probes.

    Args:
        _seq_id (str): sequence identifier (unused, passed for pickling compat).
        chunk_nuc_arr (numpy.ndarray): 1D uint8 array of integer-encoded bases
            for this chunk.
        chunk_start (int): start position of this chunk in the full sequence
            (used to convert chunk-relative coords to genome coords).
        _chunk_stop (int): stop position of this chunk (unused).
        config (dict): mining configuration dict.

    Returns:
        coord_result (numpy.ndarray): array of shape (N_probes, 2) with
            columns [genome_start, genome_stop] for each selected probe.
        tm_result (numpy.ndarray): 1D array of Tm values (°C) for each
            selected probe.
    """

    min_length = config['min_length']
    max_length = config['max_length']
    n_lengths  = max_length - min_length + 1

    result_length = chunk_nuc_arr.size - max_length + 1
    if result_length <= 0:
        return np.zeros((0, 2), dtype=int), np.zeros(0)

    # --- N-mask: identify positions with non-ACGT bases (encoded as 4) ---
    # Clip N→0 for Tm LUT lookups (Tm values for N-containing probes will be masked out anyway)
    has_n = (chunk_nuc_arr == 4)
    if has_n.any():
        n_cumsum = np.concatenate([[0], np.cumsum(has_n.astype(np.int32))])
        tm_nuc_arr = chunk_nuc_arr.copy()
        tm_nuc_arr[has_n] = 0
    else:
        n_cumsum = None
        tm_nuc_arr = chunk_nuc_arr

    # compute the nearest-neighbor Tm of every candidate seq in this chunk
    tm_grid = get_tm_grid(tm_nuc_arr, config)

    # valid_mask: inclusive Tm range
    valid_mask = (tm_grid >= config['min_tm']) & (tm_grid <= config['max_tm'])

    # apply N-mask: exclude any probe window containing an N
    if n_cumsum is not None:
        for j in range(n_lengths):
            L = min_length + j
            contains_n = (n_cumsum[L:L + result_length] - n_cumsum[:result_length]) > 0
            valid_mask[:, j] &= ~contains_n

    # apply homopolymer filter: reject probes with any run longer than max_homopolymer
    max_hpoly = config.get('max_homopolymer')
    if max_hpoly is not None:
        r = max_hpoly + 1  # minimum prohibited run length
        n_run_pos = chunk_nuc_arr.size - r + 1
        if n_run_pos > 0:
            # has_run[k] = True if chunk_nuc_arr[k:k+r] are all the same base
            has_run = np.ones(n_run_pos, dtype=bool)
            for i in range(1, r):
                has_run &= (chunk_nuc_arr[:n_run_pos] == chunk_nuc_arr[i:i + n_run_pos])
            run_cumsum = np.concatenate([[0], np.cumsum(has_run.astype(np.int32))])
            for j in range(n_lengths):
                L = min_length + j
                n_run_in_probe = L - r + 1  # positions in probe where a run of r can start
                if n_run_in_probe > 0:
                    valid_mask[:, j] &= (
                        run_cumsum[n_run_in_probe:n_run_in_probe + result_length]
                        - run_cumsum[:result_length]
                    ) == 0

    # apply prohibited_seqs filter: reject probes containing any prohibited substring
    prohibited_encoded = config.get('_prohibited_encoded')
    if prohibited_encoded:
        for p_enc in prohibited_encoded:
            pl = len(p_enc)
            if pl == 0 or chunk_nuc_arr.size < pl:
                continue
            n_match_pos = chunk_nuc_arr.size - pl + 1
            # match[k] = True if chunk_nuc_arr[k:k+pl] == p_enc
            match = chunk_nuc_arr[:n_match_pos] == p_enc[0]
            for i in range(1, pl):
                match &= chunk_nuc_arr[i:i + n_match_pos] == p_enc[i]
            match_cumsum = np.concatenate([[0], np.cumsum(match.astype(np.int32))])
            for j in range(n_lengths):
                L = min_length + j
                if pl > L:
                    continue
                n_start = L - pl + 1  # positions in probe where p can start
                valid_mask[:, j] &= (
                    match_cumsum[n_start:n_start + result_length]
                    - match_cumsum[:result_length]
                ) == 0

    if config.get('exhaustive'):
        # exhaustive mode: return all valid (position, length) combinations
        starts, length_idx = np.where(valid_mask)
    else:
        # select best probe length per start position
        if config['tm_target'] is None:
            # greedy minimum length: pick the first (shortest) valid probe per position
            length_idx = np.argmax(valid_mask, axis=1)
        else:
            # closest-to-target: among valid lengths, pick the one with Tm nearest tm_target
            masked_deltas = np.where(valid_mask, np.abs(tm_grid - config['tm_target']), np.inf)
            length_idx = np.argmin(masked_deltas, axis=1)

        # filter out positions where no valid probe exists
        any_valid = valid_mask.any(axis=1)
        starts     = np.where(any_valid)[0]
        length_idx = length_idx[any_valid]

    # apply GC% filter post-selection (only compute for chosen probes, not full grid)
    min_gc = config.get('min_gc')
    max_gc = config.get('max_gc')
    if (min_gc is not None or max_gc is not None) and len(starts):
        lengths = min_length + length_idx
        is_gc = ((tm_nuc_arr == 1) | (tm_nuc_arr == 2)).astype(np.int32)
        gc_cumsum = np.concatenate([[0], np.cumsum(is_gc)])
        gc_fracs = (gc_cumsum[starts + lengths] - gc_cumsum[starts]) / lengths
        gc_ok = np.ones(len(starts), dtype=bool)
        if min_gc is not None:
            gc_ok &= gc_fracs >= min_gc / 100.0
        if max_gc is not None:
            gc_ok &= gc_fracs <= max_gc / 100.0
        starts     = starts[gc_ok]
        length_idx = length_idx[gc_ok]

    if len(starts) == 0:
        return np.zeros((0, 2), dtype=int), np.zeros(0)

    # lookup probe Tm values
    tm_result = tm_grid[starts, length_idx]

    # convert chunk-relative start + length index → genome start + stop
    coord_result = np.empty((len(starts), 2), dtype=int)
    coord_result[:, 0] = starts + chunk_start                           # genome start
    coord_result[:, 1] = coord_result[:, 0] + min_length + length_idx  # genome stop

    return coord_result, tm_result
