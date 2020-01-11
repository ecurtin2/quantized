import numpy as np
from hypothesis.strategies import floats

from typing import Iterable

from quantized.config import conf
from quantized.utils import isclose

reasonable_floats = floats(
    allow_infinity=False, allow_nan=False, min_value=-conf.large_number, max_value=conf.large_number
)

reasonable_pos_floats = floats(
    allow_infinity=False, allow_nan=False, min_value=conf.small_number, max_value=conf.large_number
)


def allclose(x: Iterable[float], y: Iterable[float]) -> bool:
    return all(isclose(a, b) for a, b in zip(x, y))


def is_diagonal(x: np.ndarray) -> bool:
    n, m = x.shape
    off_diags = x[~np.eye(n, m, dtype=bool)]
    return np.allclose(off_diags, 0)


def is_identity(x: np.ndarray) -> bool:
    n, m = x.shape
    return np.allclose(np.eye(n, m), x)


def is_hermitian(x: np.ndarray) -> bool:
    return np.allclose(x, np.conj(x).T)
