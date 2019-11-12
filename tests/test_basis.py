from math import isclose

import pytest
from hypothesis import given
from hypothesis.strategies import integers
from pytest import raises

from transit_chem.basis import HarmonicOscillator
from transit_chem.operators import overlap
from transit_chem.config import (
    FLOAT_TOL,
    HARMONIC_OSCILLATOR_MAX_N,
    SMALL_NUMBER,
)


from utils import reasonable_floats, reasonable_pos_floats


@given(
    integers(min_value=0, max_value=HARMONIC_OSCILLATOR_MAX_N),
    reasonable_floats,
    reasonable_pos_floats,
    reasonable_pos_floats,
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
    assert isclose(overlap(ho0, ho0), 1.0, abs_tol=SMALL_NUMBER)
    assert isclose(overlap(ho0, ho1), 0.0, abs_tol=SMALL_NUMBER)


def test_total_energy_integral_of_harmonic_oscillator():

    for n in range(10):
        ho0 = HarmonicOscillator(n=n, center=0.0)

        ke = overlap(ho0, ho0.__kinetic__())

        def v(x):
            return ho0.potential(x) * ho0(x)

        pe = overlap(ho0, v)

        assert isclose(pe + ke, ho0.energy, abs_tol=SMALL_NUMBER)
