from functools import lru_cache
from itertools import chain, count, islice, takewhile, tee
from math import isclose
from typing import Callable, Generator, Iterable, Iterator, List, Tuple, Union

import attr
import numpy as np
from loguru import logger
from scipy.interpolate import interp1d
from scipy.optimize import root_scalar

from transit_chem.operators import Overlap
from transit_chem.time_evolution import TimeEvolvingObservable, TimeEvolvingState


def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


@attr.s(frozen=True)
class OccupancyProbabilites:
    s: Tuple[TimeEvolvingObservable] = attr.ib(converter=tuple)

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
    hopping_matrix: HoppingMatrix = attr.ib(hash=False)
    times: List[float] = attr.ib()
    values: List[float] = attr.ib()
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

    @property
    def min_time(self) -> float:
        return self.times[0]

    @property
    def max_time(self) -> float:
        return self.times[-1]

    @property
    def tau90(self) -> float:
        return self.time_when_equal_to(0.1)

    @property
    def initial_value(self) -> float:
        return self(self.min_time)

    @property
    def final_value(self) -> float:
        return self(self.max_time)

    @lru_cache()
    def time_when_equal_to(self, x: float) -> float:
        def f(t):
            return self(t) - x

        if not self.final_value < x < self.initial_value:
            raise ValueError(
                f"Can't solve for Pnot = {x}. It is out the calculated range: {self.initial_value}, {self.final_value}"
            )

        result = root_scalar(f, bracket=[self.min_time, self.max_time])
        return result.root

    @lru_cache()
    def _interpolate(self):
        return interp1d(self.times, self.values, bounds_error=False)

    def __call__(self, t: float) -> float:
        return self._interpolate()(t)

    @property
    def acceptor_occ_prob(self) -> Callable:
        return self.hopping_matrix.occ_probs[self.acceptor]

    @staticmethod
    def gen(
        delta_t: float, hopping_matrix: HoppingMatrix, acceptor: int
    ) -> Generator[Tuple, None, None]:
        times = (delta_t * i for i in count())
        matrices = (hopping_matrix(t, delta_t) for t in times)
        A_tilde = (np.delete(np.delete(a, acceptor, axis=0), acceptor, axis=1) for a in matrices)
        p0_tilde = np.delete(hopping_matrix.occ_probs.initial, acceptor)
        #  Taking [a, nb, c, d, ...] -> [(b @ a), (c @ b @ a), (d @ c @ b @ a), ...]
        a_prod = accumulate(A_tilde, lambda a, b: b @ a, initial=p0_tilde)
        pnot = map(np.sum, a_prod)
        return ((delta_t * i, p) for i, p in enumerate(pnot))

    @staticmethod
    def gen_until_time(
        t: float, delta_t: float, hopping_matrix: HoppingMatrix, acceptor: int
    ) -> Iterator[Tuple[float, float]]:

        gen = Pnot.gen(delta_t=delta_t, hopping_matrix=hopping_matrix, acceptor=acceptor)
        while_lt = takewhile(lambda x: x[0] < t, gen)
        extra_terms = islice(gen, 100)
        return chain(while_lt, extra_terms)

    @staticmethod
    def gen_until_prob(
        p: float, delta_t: float, hopping_matrix: HoppingMatrix, acceptor: int
    ) -> Iterator[Tuple[float, float]]:
        gen = Pnot.gen(delta_t=delta_t, hopping_matrix=hopping_matrix, acceptor=acceptor)
        while_gt = takewhile(lambda x: x[1] > p, gen)
        extra_terms = islice(gen, 100)
        return chain(while_gt, extra_terms)

    @staticmethod
    def converged_with_timestep(
        a: HoppingMatrix,
        acceptor: int,
        until_equals: float,
        tolerance: float,
        max_dt: float = 1.0,
        min_dt: float = 1e-4,
    ) -> "Pnot":
        if not max_dt > min_dt:
            raise ValueError(f"Min dt must be below max dt min_dt={min_dt}, max_dt={max_dt}")

        dts = [max_dt]
        last_dt = dts[-1]
        while last_dt > min_dt:
            dts.append(last_dt / 1.61)
            last_dt = dts[-1]

        last_t = -(10.0 ** 6)

        for dt in dts:

            times, p_not_list = zip(
                *Pnot.gen_until_prob(until_equals, acceptor=acceptor, delta_t=dt, hopping_matrix=a)
            )

            p_not = Pnot(
                acceptor=acceptor, times=times, values=p_not_list, delta_t=dt, hopping_matrix=a
            )

            t = p_not.time_when_equal_to(until_equals)
            logger.info(f"For delta_t = {dt:.5f}, Tau = {t:.5f}")
            if isclose(t, last_t, abs_tol=tolerance):
                break
            last_t = t
        else:
            raise ValueError(
                f"Time for Pnot to reach {until_equals} did not converge with timestep to tolerance={tolerance}"
            )
        logger.info(
            f"Yes!! Time for pnot to reach ({until_equals:6f}) successfully converged to {t:6f}"
            f" with difference of {t - last_t:3f}"
            f", which is within the given tolerance, {tolerance}"
        )
        return p_not
