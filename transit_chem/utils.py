from __future__ import annotations

from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from itertools import combinations_with_replacement, product
from math import acos
from math import isclose as math_isclose
from math import sqrt
from typing import Callable, Dict, Iterable, Sequence, Tuple, TypeVar

import attr
import numpy as np
from joblib import Memory
from tqdm import tqdm

from transit_chem.config import conf
from transit_chem.validation import Range, not_inf, not_nan

__all__ = ["pairwise_array_from_func", "Parabola", "LinearComb", "cache"]


memory = Memory(conf.cache_dir, verbose=conf.joblib_verbosity)
cache = memory.cache

T = TypeVar("T")


def isclose(x: float, y: float) -> bool:
    return math_isclose(x, y, abs_tol=conf.small_number)


class LinearComb:
    def __init__(self, c, f):
        self.c = c
        self.f = f

    def __call__(self, x):
        return sum(ci * fi(x) for ci, fi in zip(self.c, self.f))

    def __repr__(self):
        return " + ".join(f"{c} * {f}" for c, f in zip(self.c, self.f))

    def __kinetic__(self):
        return LinearComb(c=self.c, f=[f.__kinetic__() for f in self.f])


@cache
def pairwise_array_from_func(
    items: Sequence[T], func: Callable[[T, T], float], symmetric=False
) -> np.ndarray:
    """Create a pairwise array by applying a function to all pairs of items.


    Parameters
    -----------
    items
        A container from which pairs will be generated. Must support len() and integer indexing over range(len(items))
    func
        A function f(first, second, *args, **kwargs) which takes 2 items and returns a float.
    symmetric
        Whether the resulting matrix should be symmetric. If true,
        will only compute each (i, j) pair once and set both [i, j] and [j, i] to that value.

    Returns
    --------
    np.array
        The resulting matrix

    Examples
    ---------

    >>> from transit_chem.utils import pairwise_array_from_func
    >>> def distance(i, j):
    ...     return abs(i - j)
    ...
    >>> pairwise_array_from_func([1, 2, 4], distance)
    array([[0., 1., 3.],
           [1., 0., 2.],
           [3., 2., 0.]])
    >>> pairwise_array_from_func([1, 2, 4], distance, symmetric=True)
    array([[0., 1., 3.],
           [1., 0., 2.],
           [3., 2., 0.]])


    """
    n = len(items)
    result = np.zeros((n, n))

    fut_res: Dict[Future, Tuple] = {}

    if symmetric:
        combs = list(combinations_with_replacement(range(n), 2))
    else:
        combs = list(product(range(n), repeat=2))

    with ProcessPoolExecutor() as exc:
        for i, j in combs:
            fut_res[exc.submit(func, items[i], items[j])] = (i, j)

        with tqdm(total=len(combs), disable=not conf.enable_progressbar, ncols=100) as pbar:
            for future in as_completed(fut_res.keys()):
                res = future.result()
                i, j = fut_res[future]
                if symmetric:
                    result[i, j] = result[j, i] = res
                else:
                    result[i, j] = res
                pbar.update(1)

    return result


@attr.s(frozen=True)
class Parabola:
    a: float = attr.ib(validator=[not_nan, not_inf, Range(-conf.large_number, conf.large_number)])
    b: float = attr.ib(validator=[not_nan, not_inf, Range(-conf.large_number, conf.large_number)])
    c: float = attr.ib(validator=[not_nan, not_inf, Range(-conf.large_number, conf.large_number)])

    def __call__(self, x):
        return self.a * x ** 2 + self.b * x + self.c

    @property
    def has_vertex(self) -> bool:
        return not isclose(self.a, 0.0)

    @property
    def vertex(self) -> float:
        if not self.has_vertex:
            raise ValueError(f"{self} is a line and has no vertex.")
        return -self.b / (2.0 * self.a)

    @staticmethod
    def from_points(
        point1: Tuple[float, float], point2: Tuple[float, float], point3: Tuple[float, float]
    ) -> Parabola:
        """Create a Parabola passing through 3 points.

        Parameters
        -----------
        point1
            x, y point
        point2
            x, y point
        point3
            x, y point

        Returns
        --------
        Parabola
            Parabola with coefficients fit to the points.

        """
        x1, y1 = point1
        x2, y2 = point2
        x3, y3 = point3

        # The polynomial coefficients are the solution to the
        # linear equation Ax = b where
        # x = [a, b, c]
        A = np.array([[x1 ** 2, x1, 1], [x2 ** 2, x2, 1], [x3 ** 2, x3, 1]])
        b = np.array([y1, y2, y3])
        coeffs = np.linalg.solve(A, b)

        return Parabola(a=coeffs[0], b=coeffs[1], c=coeffs[2])


def dot(v1: Iterable[float], v2: Iterable[float]) -> float:
    return sum((a * b) for a, b in zip(v1, v2))


def length(v: Iterable[float]) -> float:
    return sqrt(dot(v, v))


def angle(v1: Iterable[float], v2: Iterable[float]) -> float:
    return acos(dot(v1, v2) / (length(v1) * length(v2)))
