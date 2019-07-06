from dataclasses import dataclass

from typing import Tuple

import numba
import math
import numpy as np
from scipy import special, integrate
from scipy.misc import derivative

from . import utils

___all__ = [
    'HarmonicOscillatorWaveFunction',
    'overlap1d',
]


class Basis:

    def __init__(self, factory, x0, y0, z0, alpha):
        """Create a basis object from some parameters and a function-factory function.

        params
        ======
            factory: function g(x0, y0, z0, alpha) that returns a function f(x, y, z)
            x0, y0, z0: centers of basis function
            alpha: parameter for basis (decay coefficient or some other such nonsense).
        """
        self.f = factory(x0, y0, z0, alpha)
        self.x0 = x0
        self.y0 = y0
        self.z0 = z0
        self.alpha = alpha

    def __call__(self, x, y, z):
        return self.f(x, y, z)

    def overlap(self, other, distance_cutoff=5.0, interval=5.0, xmin=None, xmax=None, ymin=None, ymax=None,
                zmin=None, zmax=None):
        """Return the overlap of this with another

        This is to reduce the calculations performed by naively integrating over
        all space for all combinations of basis functions.
        Return 0.0 if the basis functions are separated by more than distance_cutoff.
        Integrate from interval away from the centers. If bounds are not given.

        If the midpoint of the centers of the functions is more than distance_cutoff from the boundaries,
        0.0 is returned.

        """
        assert isinstance(other, __class__)
        r = np.sqrt((self.x0 - other.x0)**2 + (self.y0 - other.y0)**2 + (self.z0 - other.z0)**2)
        if r > distance_cutoff:
            return 0.0
        else:
            if xmin is None:
                xmin = min(self.x0, other.x0) - interval
            if ymin is None:
                ymin = min(self.y0, other.y0) - interval
            if zmin is None:
                zmin = min(self.z0, other.z0) - interval
            if xmax is None:
                xmax = max(self.x0, other.x0) + interval
            if ymax is None:
                ymax = max(self.y0, other.y0) + interval
            if zmax is None:
                zmax = max(self.z0, other.z0) + interval

            #  Check to see if the midpoint of the overlapped functions falls outside integration range.
            #  The maximum value between the centers, and if it's sufficiently far from the region of integration
            #  we can return 0.0 and not lose anything.

            midpoint = (0.5 * (self.x0 + other.x0), 0.5 * (self.y0 + other.y0), 0.5 * (self.z0 + other.z0))
            if midpoint[0] > (xmax + interval) or midpoint[0] < (xmin - interval):
                return 0.0
            if midpoint[1] > (ymax + interval) or midpoint[1] < (ymin - interval):
                return 0.0
            if midpoint[2] > (zmax + interval) or midpoint[2] < (zmin - interval):
                return 0.0

            # Finally actually do the overlap if none of the zero-conditions apply
            ovlp = overlap_factory(self.f, other.f)
            integrator = utils.integrate_region_gen_3D([[xmin, xmax], [ymin, ymax], [zmin, zmax]])
            return integrator(ovlp)


def overlap_factory(f1, f2):
    """Generate a jitted function of f1(x, y, z) * f2(x, y, z). Only works on f1, f2 that are numba jitted."""
    @numba.jit('float64(float64, float64, float64)', nopython=True)
    def overlap(x, y, z):
        return f1(x, y, z) * f2(x, y, z)
    return overlap


def pz_gaussian_orbital_factory(x0, y0, z0, alpha):
    """Return a function for a gaussian pz orbital function centered at (x0, y0, z0) and decay coeff alpha.

    f = N (z - z0) * exp(-(ar)^2)

    The pz orbital is compiled at runtime using numba's jit LLVM compiler. As such, it is fast!"""

    #TODO: FIX NORMALIZATION CONSTANT.
    N = alpha ** 2.5

    @numba.jit('float64(float64, float64, float64)', nopython=True)
    def gaussian_orbital(x, y, z):

        #TODO: MAKE THIS FUNCTION WHATEVER IS DESIRED.
        r = math.sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)
        return N * (z - z0) * math.exp(-(alpha * r)**2)
    return gaussian_orbital


def pz_orbital_factory(x0, y0, z0, alpha):
    """Return a function for pz orbital function centered at (x0, y0, z0) and decay coeff alpha.

    The pz orbital is compiled at runtime using numba's jit LLVM compiler. As such, it is fast!"""
    N = alpha ** 2.5 / math.sqrt(math.pi)

    @numba.jit('float64(float64, float64, float64)', nopython=True)
    def pz_orbital(x, y, z):
        r = math.sqrt((x - x0)**2 + (y - y0)**2 + (z - z0)**2)
        return N * (z - z0) * math.exp(-alpha * r)
    return pz_orbital


def overlap1d(
        first: callable,
        second: callable,
        *args,
        lower_limit: float=-np.inf,
        upper_limit: float=np.inf
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


def kinetic_integral(
        first: callable,
        second: callable,
        *args,
        lower_limit: float=-np.inf,
        upper_limit: float=np.inf
) -> float:
    """Return kinetic energy integral of two functions

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
    float
        The value of the kinetic energy integral.

    """
    def kinetic_func(x, *args_):
        return -0.5 * derivative(second, x, dx=1e-8, n=2, args=args_)
    return overlap1d(first, kinetic_func, *args, lower_limit=lower_limit, upper_limit=upper_limit)


def make_potential_integral(potential: callable):
    """

    Parameters
    -----------

    potential
        A function V(x) that takes x and returns the potential energy.
:
    """

    def potential_integral(
            first: callable,
            second: callable,
            *args,
            lower_limit: float=-np.inf,
            upper_limit: float=np.inf
    ) -> float:
        """Return kinetic energy integral of two functions

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
        float
            The value of the kinetic energy integral.

        """
        def potential_func(x, *args_):
            return potential(x) * second(x, *args_)
        return overlap1d(first, potential_func, *args, lower_limit=lower_limit, upper_limit=upper_limit)

    return potential_integral


@dataclass
class HarmonicOscillator:
    """A 1D quantum harmonic oscillator wave function.

    Parameters
    -----------
    n
        The quantum number
    center
        The center of the potential well
    omega
        The angular frequency of the oscillator
    mass
        The mass of the particle

    Examples
    ---------

    >>> ho = HarmonicOscillator(n=2, center=0.5)
    >>> ho
    HarmonicOscillator(n=2, center=0.5, omega=1, mass=1.0)
    >>> ho(0.5)
    -0.5311259660135984
    >>> ho(1000)
    0.0
    >>> ho(-1000)
    0.0

    """

    n: int
    center: float
    mass: int = 1
    omega: float = 1.0

    @staticmethod
    def from_potential_points(
            point1: Tuple[float, float],
            point2: Tuple[float, float],
            point3: Tuple[float, float],
            n: int,
            mass: float=1.0
    ):
        """Create a harmonic oscillator wave function from 3 samples of the potential.

        The three points are fit to a parabola (harmonic potential), then the parameters
        for the harmonic oscillator are determined, and the corresponding wave function
        generated and returned.

        Parameters
        -----------
        point1
            A sample point of the potential
        point2
            A sample point of the potential
        point3
            A sample point of the potential
        n
            The quantum number of the resulting wave function
        mass
            The mass of the particle

        Returns
        --------
        HarmonicOscillator
            The harmonic oscillator wave function resulting from the potential at the 3
            points.

        Examples
        ---------

        >>> ho = HarmonicOscillator.from_potential_points(
        ...     point1=(0.5, 1),
        ...     point2=(2.0, 0.5),
        ...     point3=(3.0, 1.5),
        ...     n=0
        ... )
        >>> ho
        HarmonicOscillator(n=0, center=1.5624999999999998, omega=1.0327955589886444, mass=1.0)


        """
        a, b, c = utils.parabola_from_points(point1, point2, point3)
        center = -b / (2*a)

        # a = m/2 * w**2
        # 2a / m = w**2
        # sqrt(2a/m) = w
        w = np.sqrt(2*a / mass)
        return HarmonicOscillator(n=n, center=center, omega=w, mass=mass)



    @property
    def N(self):
        return (1.0 / math.sqrt(2 ** n * math.factorial(n)) * ((mass * omega) / math.pi) ** 0.25)

    @property
    def _hermite(self):
        return special.hermite(self.n)

    def __call__(self, x):
        y = (np.sqrt(self.mass * self.omega)) * (x - self.center)
        return self.N * np.exp(-y**2 / 2.0) * self._hermite(y)
