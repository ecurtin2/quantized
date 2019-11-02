from math import isclose

import pytest
from hypothesis import given
from hypothesis.strategies import floats, integers
from pytest import raises

from transit_chem.basis import (
    HarmonicOscillator,
    overlap1d,
    kinetic_integral,
    make_potential_integral,
)
from transit_chem.config import (
    FLOAT_TOL,
    HARMONIC_OSCILLATOR_MAX_N,
    LARGE_NUMBER,
    SMALL_NUMBER,
)


@given(
    integers(min_value=0, max_value=HARMONIC_OSCILLATOR_MAX_N),
    floats(
        min_value=-LARGE_NUMBER,
        max_value=LARGE_NUMBER,
        allow_nan=False,
        allow_infinity=False,
    ),
    floats(
        min_value=SMALL_NUMBER,
        max_value=LARGE_NUMBER,
        allow_nan=False,
        allow_infinity=False,
    ),
    floats(
        min_value=SMALL_NUMBER,
        max_value=LARGE_NUMBER,
        allow_nan=False,
        allow_infinity=False,
    ),
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


def test_overlap1d_harmonic_oscillators_is_orthonormal():
    ho0 = HarmonicOscillator(n=0, center=0.0)
    ho1 = HarmonicOscillator(n=1, center=0.0)
    assert isclose(overlap1d(ho0, ho0), 1.0, abs_tol=SMALL_NUMBER)
    assert isclose(overlap1d(ho0, ho1), 0.0, abs_tol=SMALL_NUMBER)


def test_kinetic_integral_of_harmonic_oscillator():
    ho0 = HarmonicOscillator(n=0, center=0.0)

    ke = kinetic_integral(ho0, ho0)
    pe_integral = make_potential_integral(ho0.potential)
    pe = pe_integral(ho0, ho0)

    assert isclose(pe + ke, ho0.energy, abs_tol=SMALL_NUMBER)
