import attr
from itertools import count, tee, takewhile
from math import isclose
from typing import Iterable, Union, Generator, Tuple, Callable, List, Any

import numpy as np

from transit_chem.time_evolution import TimeEvolvingState
from transit_chem.operators import Overlap


from loguru import logger


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


@attr.s(frozen=True)
class OccupancyProbabilites:
    s: Tuple = attr.ib(converter=tuple)

    @property
    def initial(self) -> np.array:
        return self.__call__(0)

    def __call__(self, t: float) -> np.array:
        return np.array([s(t) for s in self.s])

    def __getitem__(self, item) -> Callable:
        return self.s[item]

    def __iter__(self):
        return iter(self.s)

    @staticmethod
    def from_1d_state(
        state: TimeEvolvingState, borders: Tuple[float, ...]
    ) -> "OccupancyProbabilites":

        bounds = pairwise(borders)
        overlaps = [Overlap(a, b) for a, b in bounds]
        return OccupancyProbabilites(
            tuple(state.observable(ovlp, hermitian=True) for ovlp in overlaps)
        )


@attr.s(frozen=True)
class HoppingMatrix:
    occ_probs: OccupancyProbabilites = attr.ib()
    N: int = 3

    def at_time(self, t: float, delta_t: float) -> np.array:
        d_occ_probs = self.occ_probs(t) - self.occ_probs(t - delta_t)
        increasing = d_occ_probs > 0
        two_states_increasing = np.count_nonzero(increasing) == 2

        A = np.zeros((self.N, self.N))
        if two_states_increasing:
            # If the state is decreasing in probability, the diagonal is
            # determined
            for j in range(3):
                if not increasing[j]:
                    A[j, j] = self.occ_probs[j](t) / self.occ_probs[j](t - delta_t)
                    for k in range(3):
                        if k != j:
                            A[k, j] = d_occ_probs[k] / self.occ_probs[j](t)

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
                    A[j, j] = self.occ_probs[j](t) / self.occ_probs[j](t - delta_t)
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

    def __call__(self, t: Union[float, np.ndarray], delta_t: float):
        if isinstance(t, Iterable):
            result = [self.at_time(t_, delta_t) for t_ in t]
            if len(result) == 0:
                raise ValueError("Call to hopping matrix had no elements in the iterable.")
            return np.stack(result)
        else:
            return self.at_time(t, delta_t)


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


@attr.s(frozen=True)
class Pnot:
    acceptor: int = attr.ib()
    hopping_matrix: HoppingMatrix = attr.ib()
    delta_t: float = attr.ib()
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

    def __iter__(self):
        times = (self.delta_t * i for i in count())
        matrices = (self.hopping_matrix(t, self.delta_t) for t in times)
        A_tilde = (
            np.delete(np.delete(a, self.acceptor, axis=0), self.acceptor, axis=1) for a in matrices
        )
        p0_tilde = np.delete(self.hopping_matrix.occ_probs.initial, self.acceptor)
        #  Taking [a, nb, c, d, ...] -> [(b @ a), (c @ b @ a), (d @ c @ b @ a), ...]
        a_prod = accumulate(A_tilde, lambda a, b: b @ a, initial=p0_tilde)
        pnot = map(np.sum, a_prod)
        return ((self.delta_t * i, p) for i, p in enumerate(pnot))

    @property
    def acceptor_occ_prob(self) -> Callable:
        return self.hopping_matrix.occ_probs[self.acceptor]

    def until_time(self, t: float):
        return takewhile(lambda x: x[0] < t, self)

    def until_prob(self, p: float):
        return takewhile(lambda x: x[1] > p, self)

    @staticmethod
    def converged_with_timestep(
        a: HoppingMatrix, tau: float, tolerance: float, max_dt: float = 1.0, min_dt: float = 1e-4
    ) -> Tuple[List[float], List[float], "Pnot"]:
        if not max_dt > min_dt:
            raise ValueError(f"Min dt must be below max dt min_dt={min_dt}, max_dt={max_dt}")

        dts = [max_dt]
        last_dt = dts[-1]
        while last_dt > min_dt:
            dts.append(last_dt / 1.61)
            last_dt = dts[-1]

        logger.info(f"Delta_t set = {dts}")

        last_t = -(10.0 ** 6)

        for dt in dts:
            logger.info(f"Getting tau for delta_t={dt}")
            p_not = Pnot(acceptor=2, hopping_matrix=a, delta_t=dt)
            times, p_not_list = zip(*p_not.until_prob(tau))
            t = times[-1]
            logger.info(f"Tau = {t}")
            if isclose(times[-1], last_t, abs_tol=tolerance):
                break
            last_t = t
        else:
            raise ValueError(f"Tau={tau} did not converge with timestep to tolerance={tolerance}")
        logger.info(
            f"Yes!! tau ({tau}) successfully converged to {t} with difference of {t-last_t}"
            f", which is within the given tolerance, {tolerance}"
        )
        return times, p_not_list, p_not
