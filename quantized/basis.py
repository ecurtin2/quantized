from __future__ import annotations

import math
from itertools import count, takewhile
from typing import Callable, List, Tuple, Union

from quantized.attr_wrapped import attrib, attrs, document_me

# import numba
import attr
import numpy as np
from numpy.linalg import inv
from scipy import special
from scipy.linalg import eigh

from quantized.config import conf
from quantized.operators import Overlap
from quantized.potentials import Harmonic
from quantized.utils import LinearComb, Parabola, cache, isclose
from quantized.validation import Range

__all__ = [
    "HarmonicOscillator",
    "EigenBasis",
    "harmonic_basis_from_parabola",
    "get_expansion_coeffs",
]


# class ThreeDBasis:
#     def __init__(self, factory, x0, y0, z0, alpha):
#         """Create a basis object from some parameters and a function-factory function.
#
#         params
#         ======
#             factory: function g(x0, y0, z0, alpha) that returns a function f(x, y, z)
#             x0, y0, z0: centers of basis function
#             alpha: parameter for basis (decay coefficient or some other such nonsense).
#         """
#         self.f = factory(x0, y0, z0, alpha)
#         self.x0 = x0
#         self.y0 = y0
#         self.z0 = z0
#         self.alpha = alpha
#
#     def __call__(self, x, y, z):
#         return self.f(x, y, z)
#
#     def overlap(
#         self,
#         other,
#         distance_cutoff=5.0,
#         interval=5.0,
#         xmin=None,
#         xmax=None,
#         ymin=None,
#         ymax=None,
#         zmin=None,
#         zmax=None,
#     ):
#         """Return the overlap of this with another
#
#         This is to reduce the calculations performed by naively integrating over
#         all space for all combinations of basis functions.
#         Return 0.0 if the basis functions are separated by more than distance_cutoff.
#         Integrate from interval away from the centers. If bounds are not given.
#
#         If the midpoint of the centers of the functions is more than distance_cutoff from the boundaries,
#         0.0 is returned.
#
#         """
#         if not isinstance(other, ThreeDBasis):
#             raise TypeError(f"ThreeDBasis can't be overlapped with {type(other)}")
#         r = np.sqrt(
#             (self.x0 - other.x0) ** 2 + (self.y0 - other.y0) ** 2 + (self.z0 - other.z0) ** 2
#         )
#         if r > distance_cutoff:
#             return 0.0
#         else:
#             if xmin is None:
#                 xmin = min(self.x0, other.x0) - interval
#             if ymin is None:
#                 ymin = min(self.y0, other.y0) - interval
#             if zmin is None:
#                 zmin = min(self.z0, other.z0) - interval
#             if xmax is None:
#                 xmax = max(self.x0, other.x0) + interval
#             if ymax is None:
#                 ymax = max(self.y0, other.y0) + interval
#             if zmax is None:
#                 zmax = max(self.z0, other.z0) + interval
#
#             #  Check to see if the midpoint of the overlapped functions
#             # falls outside integration range.
#             #  The maximum value between the centers,
#             # and if it's sufficiently far from the region of integration
#             #  we can return 0.0 and not lose anything.
#
#             midpoint = (
#                 0.5 * (self.x0 + other.x0),
#                 0.5 * (self.y0 + other.y0),
#                 0.5 * (self.z0 + other.z0),
#             )
#             if midpoint[0] > (xmax + interval) or midpoint[0] < (xmin - interval):
#                 return 0.0
#             if midpoint[1] > (ymax + interval) or midpoint[1] < (ymin - interval):
#                 return 0.0
#             if midpoint[2] > (zmax + interval) or midpoint[2] < (zmin - interval):
#                 return 0.0
#
#             # Finally actually do the overlap if none of the zero-conditions apply
#             ovlp = overlap_factory(self.f, other.f)
#             integrator = integrate_region_gen_3D([[xmin, xmax], [ymin, ymax], [zmin, zmax]])
#             return integrator(ovlp)
#
#
# def overlap_factory(f1, f2):
#     """Generate a jitted function of f1(x, y, z) * f2(x, y, z). Only works on f1, f2 that are numba jitted."""
#
#     @numba.jit("float64(float64, float64, float64)", nopython=True)
#     def overlap_function(x, y, z):
#         return f1(x, y, z) * f2(x, y, z)
#
#     return overlap_function
#
#
# @attr.s()
# class PzOrbital:
#     x: float = attr.ib()
#     y: float = attr.ib()
#     z: float = attr.ib()
#     alpha: float = attr.ib()
#
#     @property
#     def function(self) -> Callable[[float, float, float], float]:
#         """Return a function for pz orbital function centered at (x0, y0, z0) and decay coeff alpha.
#
#         The pz orbital is compiled at runtime using numba's jit LLVM compiler. As such, it is fast!
#         """
#         N = self.alpha ** 2.5 / math.sqrt(math.pi)
#
#         @numba.jit("float64(float64, float64, float64)", nopython=True)
#         def pz_orbital(x, y, z):
#             r = math.sqrt((x - self.x) ** 2 + (y - self.y) ** 2 + (z - self.z) ** 2)
#             return N * (z - self.z) * math.exp(-self.alpha * r)
#
#         return pz_orbital


@attrs(frozen=True)
class HarmonicOscillator:
    """A 1D quantum harmonic oscillator wave function.

    ```
    >>> ho = HarmonicOscillator(n=2, center=0.5)
    >>> ho
    HarmonicOscillator(n=2, center=0.5, omega=1, mass=1.0)
    >>> ho(0.5)
    -0.5311259660135984
    >>> ho(1000)
    0.0
    >>> ho(-1000)
    0.0
    ```
    """

    n: int = attrib(validator=[Range(0, conf.harmonic_oscillator_max_n)], desc="The quantum number")
    center: float = attrib(
        validator=[Range(-conf.large_number, conf.large_number)], desc="The center of the function"
    )
    mass: float = attrib(
        default=1.0,
        validator=Range(conf.small_number, conf.large_number),
        desc="Mass of the particle",
    )
    omega: float = attrib(
        default=1.0,
        validator=Range(conf.small_number, conf.large_number),
        desc="Angular frequency of the oscillator",
    )

    @staticmethod
    def from_parabola(p: Parabola, n: int, mass: float = 1.0) -> HarmonicOscillator:
        """Create a harmonic oscillator, who's potential is defined by the given parabola"""
        # a = m/2 * w**2
        # 2a / m = w**2
        # sqrt(2a/m) = w
        w = np.sqrt(2 * p.a / mass)
        return HarmonicOscillator(n=n, center=p.vertex, omega=w, mass=mass)

    @staticmethod
    def from_potential_points(
        point1: Tuple[float, float],
        point2: Tuple[float, float],
        point3: Tuple[float, float],
        n: int,
        mass: float = 1.0,
    ) -> HarmonicOscillator:
        """Create a harmonic oscillator wave function from 3 samples of the potential.

        The three points are fit to a parabola (harmonic potential), then the parameters
        for the harmonic oscillator are determined, and the corresponding wave function
        generated and returned.

        **point1**
        A sample point of the potential

        **point2**
        A sample point of the potential

        **point3**
        A sample point of the potential

        **n**
        The quantum number of the resulting wave function

        **mass**
        The mass of the particle

        **Examples**

        ```python
        ho = HarmonicOscillator.from_potential_points(
        ...     point1=(0.5, 1),
        ...     point2=(2.0, 0.5),
        ...     point3=(3.0, 1.5),
        ...     n=0
        ... )
        ho
        HarmonicOscillator(n=0, center=1.5624999999999998, omega=1.0327955589886444, mass=1.0)
        ```

        """
        p = Parabola.from_points(point1, point2, point3)
        return HarmonicOscillator.from_parabola(p, n=n, mass=mass)

    @document_me
    def __kinetic__(self) -> Callable:
        """Return kinetic energy operator applied on this."""

        # K = \frac{p^2}{2m}
        # p = i * \sqrt{m w hbar/2}(a+ - a)
        #
        # k = -1/2 * (m w hbar / 2)[(a+ - a)^2]
        # [] = (a+^2 + a+a + aa+  - a^2)
        #
        # [] = sqrt(n+1)sqrt(n+2)| n + 2 >  always
        #         sqrt(n+1)sqrt(n+1)| n >      always
        #         sqrt(n)  sqrt(n)  | n >      if n == 0  0
        #         sqrt(n)  sqrt(n-1)| n - 2 >  if n <= 1  0
        #
        # k = - (m w hbar) / 4 * []

        def k(x):
            first = np.sqrt(self.n + 1) * np.sqrt(self.n + 2) * attr.evolve(self, n=self.n + 2)(x)

            if self.n == 0:
                # Lowering operator on n=0
                # means the second term is zero
                inner = (self.n + 1) * self(x)
            else:
                inner = (2 * self.n + 1) * self(x)

            if self.n <= 1:
                # Lowering operator twice on n=0 or n=1 will be 0
                # if x is numpy array, we want to return a numpy
                # array of zeros, so multiply x by zeros will work on
                # numpy array or float.
                last = 0.0 * x
            else:
                last = np.sqrt(self.n) * np.sqrt(self.n - 1) * attr.evolve(self, n=self.n - 2)(x)

            return -self.mass * self.omega * (first - inner + last) / 4.0

        return k

    @document_me
    def __overlap__(self, other, lower_limit: float, upper_limit: float) -> float:
        """Determine the overlap with some other function

        This specializes a generic overlap integral, and short circuits integral
        calculations if the integral is analytically known.
        """
        if lower_limit != -np.inf or upper_limit != np.inf:
            raise NotImplementedError
        if not isinstance(other, HarmonicOscillator):
            raise NotImplementedError

        if (
            isclose(self.center, other.center)
            and isclose(self.mass, other.mass)
            and isclose(self.omega, other.omega)
        ):
            return 1.0 if self.n == other.n else 0.0
        raise NotImplementedError

    @property
    def energy(self):
        """The energy of harmonic oscillator"""
        return (self.n + 0.5) * self.omega

    @property
    def potential(self) -> Harmonic:
        """The potential for this oscillator"""
        return Harmonic(center=self.center, mass=self.mass, omega=self.omega)

    @property
    def N(self):
        """The normalization constant"""
        return (
            1.0
            / math.sqrt(2 ** self.n * math.factorial(self.n))
            * ((self.mass * self.omega) / math.pi) ** 0.25
        )

    @property
    def _hermite(self):
        return special.hermite(self.n)

    @document_me
    def __call__(self, x: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
        """Return"""
        y = (np.sqrt(self.mass * self.omega)) * (x - self.center)
        return self.N * np.exp(-(y ** 2) / 2.0) * self._hermite(y)


def harmonic_basis_from_parabola(p: Parabola, cutoff_energy: float) -> List[HarmonicOscillator]:
    """Create a set of harmonic oscillator wave functions below some cutoff energy.

    **p (Parabola)**
    The parabola which will be used as the harmonic oscillator's potential surface.

    **cutoff_energy (float)**
    The energy at which to stop creating basis functions. That is, all basis functions created
    will have energy less than or equal to `cutoff_energy`


    """
    hos = (HarmonicOscillator.from_parabola(p, n=i) for i in count())
    return list(takewhile(lambda ho: ho.energy <= cutoff_energy, hos))


@attrs(frozen=True)
class EigenBasis:
    """A class for representing an eigenbasis for a hamiltonian"""

    # TODO: validate shapes, properties
    states: Tuple[Callable, ...] = attrib(converter=tuple, desc="A set of eigen states")
    energies: Tuple[float, ...] = attrib(converter=tuple, desc="The energies of the eigen states")
    ao_S: np.ndarray = attrib(desc="The overlap matrix in the original basis")
    eigvecs: np.ndarray = attrib(
        desc="The eigenvectors of the hamiltonian. Each column is a vector."
    )
    ao_basis: List[Callable] = attrib(desc="The original basis")

    @document_me
    def __len__(self):
        """The size of the eigenbasis"""
        return len(self.states)

    @staticmethod
    def from_basis(basis: List[Callable], H: np.ndarray, S: np.ndarray) -> EigenBasis:
        """Create an eigenbasis from another basis, given a hamiltonian and overlap matrix

        **H (np.ndarray)**
        The hamiltonian matrix in the basis

        **S (np.ndarray)**
        The overlap matrix in the basis

        """
        # Sort vals/vecs
        eigvals, eigvecs = eigh(H, S)
        idx = np.argsort(eigvals)
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]

        states = [LinearComb(c, basis) for c in eigvecs.T]
        return EigenBasis(
            states=tuple(states), energies=tuple(eigvals), eigvecs=eigvecs, ao_S=S, ao_basis=basis
        )

    def transformed(self, matrix: np.ndarray) -> np.ndarray:
        """Given a matrix in the original basis, return the matrix in the Eigen basis."""
        return inv(self.eigvecs) @ inv(self.ao_S) @ matrix @ self.eigvecs


@cache
def get_expansion_coeffs(state: Callable, basis: List[Callable]) -> List[float]:
    """Given a state and a basis, return the expansion coefficents for that state"""
    return [Overlap()(state, basis_func) for basis_func in basis]
