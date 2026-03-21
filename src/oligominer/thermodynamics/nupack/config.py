
import nupack

from oligominer.thermodynamics import formamide_correction

GET_DEFAULT_NUPACK_MODEL = lambda: nupack.Model(
    material='dna',
    ensemble='stacking',
    celsius=formamide_correction(37.0, 50, 0.65),
    sodium=0.390,
    magnesium=0.0
)

# pre-instantiated default model for repeated use
DEFAULT_NUPACK_MODEL = GET_DEFAULT_NUPACK_MODEL()
