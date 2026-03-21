"""config.py — Oligonucleotide probe mining configuration management."""


GET_DEFAULT_MINING_CONFIG = lambda: {

    # probe mining params
    'min_length': 30,
    'max_length': 37,
    'min_tm': 42,
    'max_tm': 47,
    'tm_target': None,  # None = greedy minimum length; set to float for closest-to-target behavior
    'bed': False,
    'allow_overlap': True,
    'spacing': 0,              # minimum bases between adjacent probes (0 = no gap required)
    'chunk_size': 100000,
    'cores': 1,

    # additional probe quality filters
    'min_gc': 20,              # minimum GC% (0–100)
    'max_gc': 80,              # maximum GC% (0–100)
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


