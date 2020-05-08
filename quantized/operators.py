from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Callable

import numpy as np
from scipy import integrate


from quantized.attr_wrapped import attrs, attrib, document_me
from quantized.molecule import Molecule
from quantized.utils import pairwise_array_from_func


__all__ = [
    "Operator",
    "Overlap",
    "Kinetic",
    "Hamiltonian",
    "Potential",
    "ExtendedHuckelHamiltonian",
]

# TODO: get the type hints to work for Operator throughout codebase
class Operator(ABC):
    @document_me
    def __call__(self, first: Callable[[float], float], second: Callable[[float], float]) -> float:
        pass

    @property
    @abstractmethod
    def hermitian(self) -> bool:
        return True

    def matrix(self, basis) -> np.array:
        """Return a matrix of the operator projected onto a basis."""
        return pairwise_array_from_func(basis, self, symmetric=self.hermitian)


@attrs(frozen=True)
class Overlap(Operator):
    lower_limit: float = attrib(default=-np.inf)
    upper_limit: float = attrib(default=np.inf)
    hermitian: bool = True
    """Compute the overlap integral in 1 dimension over the specified range

    Parameters
    -----------
    first
        The first function
    second
        The second function
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

    >>> from quantized import overlap
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

    @document_me
    def __call__(self, first: Callable, second: Callable) -> float:

        try:
            val = first.__overlap__(second, self.lower_limit, self.upper_limit)
            return val
        except (AttributeError, NotImplementedError):
            pass

        def integrand(x):
            return first(x) * second(x)

        return integrate.quad(integrand, a=self.lower_limit, b=self.upper_limit)[0]


@attrs()
class Kinetic(Operator):
    overlap: Overlap = attrib(default=Overlap())
    hermitian = True

    @document_me
    def __call__(self, first, second) -> float:
        return Overlap()(first, second.__kinetic__())


@attrs()
class Hamiltonian(Operator):
    potential: Callable[[float], float] = attrib()
    kinetic: Kinetic = attrib(default=Kinetic())
    hermitian = True

    def __attrs_post_init__(self):
        self.Potential = Potential(self.potential)

    @document_me
    def __call__(self, first, second) -> float:
        return self.Potential(first, second) + self.kinetic(first, second)


@attrs(frozen=True)
class Potential(Operator):
    potential: Callable[[float], float] = attrib()
    overlap: Overlap = attrib(default=Overlap())
    hermitian = True

    @document_me
    def __call__(self, first, second) -> float:
        def v(x):
            return self.potential(x) * second(x)

        return self.overlap(first, v)


@attrs(frozen=True)
class ExtendedHuckelHamiltonian(Operator):
    S: np.array = attrib()
    molecule: Molecule = attrib()
    hermitian = True

    def __attrs_post_init__(self):
        invalid = [a.element for a in self.molecule if a.element.voie is None]
        if invalid:
            raise ValueError(
                f"Could not make EHT Hamiltonian from elements {invalid}"
                f". They are not configured with VOIE."
            )

    def matrix(self, basis=None) -> np.array:
        """Create the Hamiltonian under the Extended Hueckel Approximation."""

        h = np.zeros(self.S.shape)
        n = self.S.shape[0]

        for i in range(n):
            h[i, i] = -self.molecule.atoms[i].element.voie

        for i in range(n):
            for j in range(n):
                if i != j:
                    h[i, j] = 1.75 * (h[i, i] + h[j, j]) * self.S[i, j] / 2.0
        return h
