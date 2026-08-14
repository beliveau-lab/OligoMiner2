"""config.py — Oligonucleotide probe mining configuration management.

The keys of this dict are exactly the keyword arguments ``mine_sequence()``
accepts, and their values are exactly its defaults. ``tests/test_mining_api.py``
asserts both, so the config stays a complete and accurate description of the
mining API.
"""


GET_DEFAULT_MINING_CONFIG = lambda: {

    # probe mining params
    'min_length': 30,
    'max_length': 37,
    'min_tm': 42,
    'max_tm': 47,
    'tm_target': None,  # None = greedy minimum length; set to float for closest-to-target behavior
    'allow_overlap': True,     # permit probes that overlap on the target
    'spacing': 0,              # minimum bases between adjacent probes (0 = no gap required)
    'exhaustive': False,       # return every valid (position, length) pair
    'chunk_size': 100000,
    'cores': None,             # None = resolve from the scheduler grant, see utils.cores

    # additional probe quality filters
    'min_gc': 20,              # minimum GC% (0–100)
    'max_gc': 80,              # maximum GC% (0–100)
    'mask_soft': False,        # exclude probes overlapping lowercase (repeat-masked) sequence
    'min_entropy': None,       # minimum normalized dinucleotide entropy (0-1), None to disable
    'max_homopolymer': 4,      # reject probes with any homopolymer run longer than this (None to disable)
    'prohibited_seqs': None,   # optional list of exact substring sequences to prohibit (e.g. ['AAAAA','TTTTT'])

    # tm calculation params
    'dnac1': 25,
    'dnac2': 25,
    'Na': 390,
    'K': 0,
    'Tris': 0,
    'Mg': 0,
    'dNTPs': 0,
    'pct_formamide': 50,
    'formamide_factor': 0.65,
    
}

WRITE_BUFFER_SIZE = 1000 * 1000  # 1 MiB buffer size for writing probe output to disk


