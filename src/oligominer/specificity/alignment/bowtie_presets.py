"""
# Bowtie2 Presets

Named parameter presets and shared constants for Bowtie2. Presets are dicts of
seed and alignment parameters that can be passed directly to bowtie_align().
See bowtie_align.py for the alignment function.
"""

# expected bowtie2 index file extensions
BT2_INDEX_EXTENSIONS = ['.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2']

VERY_FAST = {'D': 5, 'R': 1, 'N': 0, 'L': 22, 'i': 'S,0,2.50', 'local': False}
FAST = {'D': 10, 'R': 2, 'N': 0, 'L': 22, 'i': 'S,0,2.50', 'local': False}
SENSITIVE = {'D': 15, 'R': 2, 'N': 0, 'L': 22, 'i': 'S,1,1.15', 'local': False}
VERY_SENSITIVE = {'D': 20, 'R': 3, 'N': 0, 'L': 20, 'i': 'S,1,0.50', 'local': False}
VERY_FAST_LOCAL = {'D': 5, 'R': 1, 'N': 0, 'L': 25, 'i': 'S,1,2.00', 'local': True}
FAST_LOCAL = {'D': 10, 'R': 2, 'N': 0, 'L': 22, 'i': 'S,1,1.75', 'local': True}
SENSITIVE_LOCAL = {'D': 15, 'R': 2, 'N': 0, 'L': 20, 'i': 'S,1,0.75', 'local': True}
VERY_SENSITIVE_LOCAL = {'D': 20, 'R': 3, 'N': 0, 'L': 20, 'i': 'S,1,0.50', 'local': True}
