import numpy as np
from hypothesis.strategies import floats

from transit_chem.config import LARGE_NUMBER, SMALL_NUMBER

reasonable_floats = floats(
    allow_infinity=False, allow_nan=False, min_value=-LARGE_NUMBER, max_value=LARGE_NUMBER
)

reasonable_pos_floats = floats(
    allow_infinity=False, allow_nan=False, min_value=SMALL_NUMBER, max_value=LARGE_NUMBER
)


def is_diagonal(x: np.ndarray) -> bool:
    n, m = x.shape
    off_diags = x[~np.eye(n, m, dtype=bool)]
    return np.allclose(off_diags, 0)


def is_identity(x: np.ndarray) -> bool:
    n, m = x.shape
    return np.allclose(np.eye(n, m), x)


def is_hermitian(x: np.ndarray) -> bool:
    return np.allclose(x, np.conj(x).T)
