"""
Default NUPACK model.

The default represents the standard FISH condition: hybridization at 37 C in 50%
formamide, at 390 mM sodium. A NUPACK model has no formamide term, so it is run at
the equivalent formamide-free temperature of 69.5 C.
"""

import nupack

from oligominer.thermodynamics import effective_hyb_temperature

# the standard FISH hybridization condition
HYB_CELSIUS = 37.0
PCT_FORMAMIDE = 50.0
FORMAMIDE_FACTOR = 0.65
SODIUM_M = 0.390

GET_DEFAULT_NUPACK_MODEL = lambda: nupack.Model(
    material='dna',
    ensemble='stacking',
    celsius=effective_hyb_temperature(HYB_CELSIUS, PCT_FORMAMIDE, FORMAMIDE_FACTOR),
    sodium=SODIUM_M,
    magnesium=0.0
)

# pre-instantiated default model for repeated use
DEFAULT_NUPACK_MODEL = GET_DEFAULT_NUPACK_MODEL()
