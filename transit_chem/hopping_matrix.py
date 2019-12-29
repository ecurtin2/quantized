from itertools import count
from typing import Iterable, Union

import numpy as np


def hopping_matrix(s1, s2, s3):
    occ_probs = [s1, s2, s3]

    def f(t: float, delta_t: float) -> np.array:
        d1 = s1(t) - s1(t - delta_t)
        d2 = s2(t) - s2(t - delta_t)
        d3 = s3(t) - s3(t - delta_t)
        increasing = [d > 0 for d in (d1, d2, d3)]

        d_occ_probs = [d1, d2, d3]

        two_states_increasing = d1 * d2 * d3 < 0

        A = np.zeros((3, 3))
        if two_states_increasing:
            # If the state is decreasing in probability, the diagonal is
            # determined
            for j in range(3):
                if not increasing[j]:
                    A[j, j] = occ_probs[j](t) / occ_probs[j](t - delta_t)
                    for k in range(3):
                        if k != j:
                            A[k, j] = d_occ_probs[k] / occ_probs[j](t)

                # If the state is increasing, diagonal is 1. Off diagonal
                # columns = 0
                elif increasing[j] > 0:
                    A[j, j] = 1
                    for k in range(3):
                        if k != j:
                            A[k, j] = 0

        else:
            for j in range(3):
                # If the state is decreasing in probability, the diagonal
                # is determined
                if not increasing[j]:
                    A[j, j] = occ_probs[j](t) / occ_probs[j](t - delta_t)
                    for k in range(3):
                        if k != j:
                            A[j, k] = 0

            # If the state is increasing, diagonal is 1. Off diagonal
            # column element is 1 - diagonal
            for j in range(3):
                if increasing[j]:
                    A[j, j] = 1
                    for k in range(3):
                        if k != j:
                            A[j, k] = 1 - A[k, k]

        return A

    # Wrap f such that it can be called with iterable of times
    def maybe_array_wrapper(t: Union[float, np.ndarray], delta_t: float):
        if isinstance(t, Iterable):
            result = [f(t_, delta_t) for t_ in t]
            if len(result) == 0:
                raise ValueError("Call to hopping matrix had no elements in the iterable.")
            return np.stack(result)
        else:
            return f(t, delta_t)

    return maybe_array_wrapper


def accumulate(iterable, func, *, initial=None):
    """https://docs.python.org/3/library/itertools.html#itertools.accumulate

    Stolen but not used since I ddidn have ython 3.8 handy and wanted to use the initial
    keyword argument

    """
    it = iter(iterable)
    total = initial
    if initial is None:
        try:
            total = next(it)
        except StopIteration:
            return
    yield total
    for element in it:
        total = func(total, element)
        yield total


def p_not_generator(acceptor: int, hopping_matrix: callable, p0: np.ndarray, delta_t):
    """Determine P_not, the probability than acceptor state has never been occupied.

    Parameters
    -----------
    acceptor
        The index of the acceptor well
    hopping_matrix
        A (N, 3, 3) array containing the hopping matrix over time
    p0
        An (3,) array containing probabilities of occupying each region initially

    """

    if not 0 <= acceptor <= 2:
        raise ValueError(f"Acceptor value out of range, must be 0, 1 or 2, got {acceptor}")
    if not p0.shape == (3,):
        raise ValueError(f"p0 shape must be (3,), fot {p0.shape}")

    #  Remove the acceptor rows and columns
    p0_tilde = np.delete(p0, acceptor)

    times = (delta_t * i for i in count())
    matrices = (hopping_matrix(t, delta_t) for t in times)
    A_tilde = (np.delete(np.delete(a, acceptor, axis=0), acceptor, axis=1) for a in matrices)

    #  Taking [a, b, c, d, ...] -> [(b @ a), (c @ b @ a), (d @ c @ b @ a), ...]
    a_prod = accumulate(A_tilde, lambda a, b: b @ a, initial=p0_tilde)
    pnot = map(np.sum, a_prod)
    return pnot
