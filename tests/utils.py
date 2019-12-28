from hypothesis.strategies import floats

from transit_chem.config import LARGE_NUMBER, SMALL_NUMBER

reasonable_floats = floats(
    allow_infinity=False,
    allow_nan=False,
    min_value=-LARGE_NUMBER,
    max_value=LARGE_NUMBER,
)

reasonable_pos_floats = floats(
    allow_infinity=False,
    allow_nan=False,
    min_value=SMALL_NUMBER,
    max_value=LARGE_NUMBER,
)
