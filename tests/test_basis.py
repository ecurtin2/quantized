from math import isclose

from transit_chem.basis import HarmonicOscillator
from transit_chem.config import HARMONIC_OSCILLATOR_MAX_N, SMALL_NUMBER, LARGE_NUMBER, FLOAT_TOL

from hypothesis import given
from hypothesis.strategies import floats, integers
from pytest import raises
import pytest


@given(
    integers(min_value=0, max_value=HARMONIC_OSCILLATOR_MAX_N),
    floats(allow_nan=False, allow_infinity=False),
    floats(min_value=SMALL_NUMBER, max_value=LARGE_NUMBER, allow_nan=False, allow_infinity=False),
    floats(min_value=SMALL_NUMBER, max_value=LARGE_NUMBER, allow_nan=False, allow_infinity=False),
)
def test_harmonic_oscillator_properties(
    n: int, center: float, mass: float, omega: float
):
    ho = HarmonicOscillator(n=n, center=center, mass=mass, omega=omega)

    if n % 2 == 0:
        # Even n are symmetric about center
        assert isclose(ho(center + 3.1), ho(center - 3.1), abs_tol=FLOAT_TOL)
    else:
        # While odd n are antisymmetric
        assert isclose(ho(center + 4.7), -ho(center - 4.7), abs_tol=FLOAT_TOL)


@pytest.mark.parametrize(
    "n,center,mass,omega", [(0, 1.3, -1.0, 1.0), (-1, 0.0, 1.0, 0.0)]
)
def test_harmonic_oscillator_raises_on_invalid(n, center, mass, omega):
    with raises(ValueError):
        HarmonicOscillator(n=n, center=center, mass=mass, omega=omega)
