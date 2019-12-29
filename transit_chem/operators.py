from __future__ import annotations

from abc import abstractmethod
from typing import Callable

import attr
import numpy as np
from scipy import integrate


# TODO: get the type hints to work for Operator throughout codebase
class Operator:
    @abstractmethod
    def __call__(self, first: Callable[[float], float], second: Callable[[float], float]) -> float:
        pass


def overlap(
    first: Callable,
    second: Callable,
    *args,
    lower_limit: float = -np.inf,
    upper_limit: float = np.inf,
) -> float:
    """Compute the overlap integral in 1 dimension over the specified range

    Parameters
    -----------
    first
        The first function
    second
        The second function
    *args
        Extra args to pass to **both** first and second function
    lower_limit
        The lower limit of integration
    upper_limit
        The upper limit of integration

    Returns
    --------
    overlap
        The value of the overlap integral.

    Examples
    ---------

    >>> from transit_chem import overlap
    >>> from numpy import sin, pi
    >>>
    >>> def f1(x):
    ...     return sin(x)
    ...
    >>> def f2(x):
    ...     return sin(2*x)
    ...
    >>> overlap(f1, f2, lower_limit=0, upper_limit=2*pi)
    6.869054119163646e-17
    >>> overlap(f1, f1, lower_limit=0, upper_limit=2*pi)
    3.141592653589793

    """

    def integrand(x, *args_):
        return first(x, *args_) * second(x, *args_)

    return integrate.quad(integrand, a=lower_limit, b=upper_limit, args=args)[0]


def kinetic(first, second) -> float:
    return overlap(first, second.__kinetic__())


@attr.s(frozen=True)
class Hamiltonian:
    potential: Callable[[float], float] = attr.ib()

    def __call__(self, first, second) -> float:
        return Potential(self.potential)(first, second) + kinetic(first, second)


@attr.s(frozen=True)
class Potential:
    potential: Callable[[float], float] = attr.ib()

    def __call__(self, first, second) -> float:
        def v(x):
            return self.potential(x) * second(x)

        return overlap(first, v)
