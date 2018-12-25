from typing import Tuple

import numpy as np


__all__ = [
    'pairwise_array_from_func',
    'parabola_from_points',

]


def pairwise_array_from_func(items, func, *args, symmetric=False, **kwargs):
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
                result[i, j] = result[j, i] = func(items[i], items[j], *args, **kwargs)
    else:
        for i in range(n):
            for j in range(n):
                result[i, j] = func(items[i], items[j], *args, **kwargs)
    return result


def parabola_from_points(
        point1: Tuple[float, float],
        point2: Tuple[float, float],
        point3: Tuple[float, float]
):
    """Return a coefficients for a parabola passing through 3 points.

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
    np.array
        The coefficients [a, b, c] for the parabola f(x) = ax^2 + bx + c.

    """
    x1, y1 = point1
    x2, y2 = point2
    x3, y3 = point3

    # The polynomial coefficients are the solution to the
    # linear equation Ax = b where
    # x = [a, b, c]
    A = np.array([
        [x1 ** 2, x1, 1],
        [x2 ** 2, x2, 1],
        [x3 ** 2, x3, 1]
    ])
    b = np.array([y1, y2, y3])
    coeffs = np.linalg.solve(A, b)
    return coeffs
