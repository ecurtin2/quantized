from typing import Callable, Union

import attr
import numpy as np

from transit_chem.basis import EigenBasis, get_expansion_coeffs
from transit_chem.operators import Operator


@attr.s
class TimeEvolvingState:
    initial_state: Callable = attr.ib()
    eigen_basis: EigenBasis = attr.ib(repr=False)

    def __attrs_post_init__(self):
        self.expansion_coeffs = get_expansion_coeffs(self.initial_state, self.eigen_basis.states)

    def __call__(self, x: Union[float, np.ndarray], t: float) -> float:
        return sum(
            [
                c * np.exp(-1j * e * t) * state(x)
                for state, e, c in zip(
                    self.eigen_basis.states, self.eigen_basis.energies, self.expansion_coeffs
                )
            ]
        )

    def observable(self, operator: Callable, hermitian: bool) -> "TimeEvolvingObservable":
        return TimeEvolvingObservable(self, operator, hermitian=hermitian)


@attr.s()
class TimeEvolvingObservable:
    time_evolving_state: TimeEvolvingState = attr.ib()
    operator: Operator = attr.ib()
    hermitian: bool = attr.ib()

    def __attrs_post_init__(self):
        e = self.time_evolving_state.eigen_basis.energies
        N = len(self.time_evolving_state.eigen_basis)
        c = self.time_evolving_state.expansion_coeffs

        # TODO: this is awkward
        ao_matrix = self.operator.matrix(self.time_evolving_state.eigen_basis.ao_basis)
        eigen_basis_matrix = self.time_evolving_state.eigen_basis.transformed(ao_matrix)

        # The following code is equivalent to the definition of f below.
        # f is written that way to optimize for speed.
        # On my machine, it was roughly 7x faster
        #
        # def f(t: float):
        #     return sum(
        #         c[i] * c[j] * np.exp(1j * (e[j] - e[i]) * t) * eigen_basis_matrix[i, j]
        #         for i in range(N)
        #         for j in range(N)
        #     )

        # Precalculate reused terms
        # P is a real matrix, but we cast it to complex here so that it doesn't
        # get casted inside the function during the multiplication.
        self.P = np.zeros_like(eigen_basis_matrix, dtype=np.complex128)
        self.W = np.zeros_like(eigen_basis_matrix, dtype=np.complex128)
        for i in range(N):
            for j in range(N):
                self.P[i, j] = c[i] * c[j] * eigen_basis_matrix[i, j]
                self.W[i, j] = np.exp(1j * (e[j] - e[i]))

    def __call__(self, t: float) -> float:
        return np.abs(np.sum(self.P * (self.W ** t)))
