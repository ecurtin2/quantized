from typing import Callable, Union
from numbers import Real

import numpy as np

from quantized.attr_wrapped import attrib, attrs, document_me
from quantized.basis import EigenBasis, get_expansion_coeffs
from quantized.operators import Operator

__all__ = ["TimeEvolvingState", "TimeEvolvingObservable"]


@attrs()
class TimeEvolvingState:
    initial_state: Callable = attrib()
    eigen_basis: EigenBasis = attrib(repr=False)

    def __attrs_post_init__(self):
        self.expansion_coeffs = get_expansion_coeffs(self.initial_state, self.eigen_basis.states)

    @document_me
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


@attrs()
class TimeEvolvingObservable:
    time_evolving_state: TimeEvolvingState = attrib()
    operator: Operator = attrib()
    hermitian: bool = attrib()

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

    @document_me
    def __call__(self, t: Union[float, np.array]) -> Union[float, np.array]:
        if isinstance(t, Real):
            return np.abs(np.sum(self.P * (self.W ** t)))
        elif isinstance(t, np.ndarray):
            return np.array([np.abs(np.sum(self.P * (self.W ** t_))) for t_ in t])
