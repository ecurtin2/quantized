from __future__ import annotations
from math import isclose
from typing import Callable, Sequence, Tuple, TypeVar

import attr
import numpy as np

from transit_chem.validation import not_inf, not_nan, Range
from transit_chem.config import SMALL_NUMBER, LARGE_NUMBER


__all__ = ["pairwise_array_from_func", "Parabola"]


T = TypeVar("T")


def pairwise_array_from_func(
    items: Sequence[T], func: Callable[[T, T], float], *args, symmetric=False, **kwargs
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

    >>> from transit_chem import pairwise_array_from_func
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

    if symmetric:
        for i in range(n):
            for j in range(i + 1):
                result[i, j] = result[j, i] = func(items[i], items[j])
    else:
        for i in range(n):
            for j in range(n):
                result[i, j] = func(items[i], items[j])
    return result


@attr.s
class Parabola:
    a: float = attr.ib(validator=[not_nan, not_inf, Range(-LARGE_NUMBER, LARGE_NUMBER)])
    b: float = attr.ib(validator=[not_nan, not_inf, Range(-LARGE_NUMBER, LARGE_NUMBER)])
    c: float = attr.ib(validator=[not_nan, not_inf, Range(-LARGE_NUMBER, LARGE_NUMBER)])

    def __call__(self, x):
        return self.a * x ** 2 + self.b * x + self.c

    @property
    def has_vertex(self) -> bool:
        return not isclose(self.a, 0.0, abs_tol=SMALL_NUMBER)

    @property
    def vertex(self) -> float:
        if not self.has_vertex:
            raise ValueError(f"{self} is a line and has no vertex.")
        return -self.b / (2.0 * self.a)

    @staticmethod
    def from_points(
        point1: Tuple[float, float],
        point2: Tuple[float, float],
        point3: Tuple[float, float],
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
