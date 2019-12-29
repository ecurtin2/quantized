from math import isclose

import numpy as np
import pytest
from hypothesis import given
from hypothesis.strategies import integers
from pytest import raises

from transit_chem import Harmonic
from transit_chem.basis import EigenBasis, HarmonicOscillator, TimeEvolvingState
from transit_chem.config import FLOAT_TOL, HARMONIC_OSCILLATOR_MAX_N, SMALL_NUMBER
from transit_chem.operators import Hamiltonian, overlap
from transit_chem.utils import pairwise_array_from_func
from utils import is_diagonal, reasonable_floats, reasonable_pos_floats


@given(
    integers(min_value=0, max_value=HARMONIC_OSCILLATOR_MAX_N),
    reasonable_floats,
    reasonable_pos_floats,
    reasonable_pos_floats,
)
def test_harmonic_oscillator_properties(n: int, center: float, mass: float, omega: float):
    ho = HarmonicOscillator(n=n, center=center, mass=mass, omega=omega)

    if n % 2 == 0:
        # Even n are symmetric about center
        assert isclose(ho(center + 3.1), ho(center - 3.1), abs_tol=FLOAT_TOL)
    else:
        # While odd n are antisymmetric
        assert isclose(ho(center + 4.7), -ho(center - 4.7), abs_tol=FLOAT_TOL)


@pytest.mark.parametrize("n,center,mass,omega", [(0, 1.3, -1.0, 1.0), (-1, 0.0, 1.0, 0.0)])
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


def test_eigen_basis_on_diagonals():
    H = np.array([[1.0, 0.0], [0.0, 2.0]])
    S = np.eye(2)

    basis = [lambda x: 2 for _ in range(2)]
    eigb = EigenBasis.from_basis(basis, H, S)

    expected_vecs = np.array([[1.0, 0.0], [0.0, 1.0]])
    expected_energies = [1.0, 2.0]

    assert np.allclose(expected_vecs, eigb.eigvecs, atol=SMALL_NUMBER)
    assert all(
        isclose(e1, e2, abs_tol=SMALL_NUMBER) for e1, e2 in zip(expected_energies, eigb.energies)
    )


def test_eigen_basis_with_offdiagonal():
    e = 1.0
    d = 0.1

    H = np.array([[e, d], [d, e]])
    S = np.eye(2)
    basis = [lambda _: 2 for _ in range(2)]
    eigb = EigenBasis.from_basis(basis, H, S)

    # The eigenvectors for the symmetric problem
    # are:
    #    (-1, 1) with eigenvalue  e - d
    #    (1, 1) with eigenvalue e + d
    # normalized, it's 1/sqrt(2) instead of 1
    x = 1 / np.sqrt(2)
    expected_vecs = np.array([[-x, x], [x, x]])
    expected_energies = [e - d, e + d]

    assert np.allclose(expected_vecs, eigb.eigvecs, atol=SMALL_NUMBER)
    assert all(
        isclose(e1, e2, abs_tol=SMALL_NUMBER) for e1, e2 in zip(expected_energies, eigb.energies)
    )


def test_eigen_basis_non_orthogonal():
    """
    This is the best thing I could think of to test this.
    Given a set of harmonic oscillators that aren't centered,
    see if we can reproduce the solution to a QHO
    """

    basis = [HarmonicOscillator(n=i, center=0.25) for i in range(2)] + [
        HarmonicOscillator(n=i, center=-0.25) for i in range(2)
    ]

    S = pairwise_array_from_func(basis, overlap)
    H = pairwise_array_from_func(basis, Hamiltonian(Harmonic(center=0.0)))
    eigb = EigenBasis.from_basis(basis, H, S)

    # check the first 3 energy levels, we won't have converged
    # the higher ones wrt basis set size
    expected_energies = [(n + 0.5) for n in range(3)]

    diffs = [e1 - e2 for e1, e2 in zip(sorted(eigb.energies), expected_energies)]

    # a little lenient due to convergence of basis to keep test fast
    assert all(isclose(d, 0.0, abs_tol=1e-3) for d in diffs)


def test_time_evolution_of_stationary_state():
    # HO is the eigenbasis
    basis = [HarmonicOscillator(n=1, center=0.0), HarmonicOscillator(n=0, center=0.0)]
    eigbasis = EigenBasis(
        states=basis,
        energies=[b.energy for b in basis],
        ao_S=np.eye(2),
        eigvecs=np.eye(2),
        ao_basis=basis,
    )

    initial_state = basis[0]

    time_evolving = TimeEvolvingState(initial_state, eigbasis)

    x = np.linspace(-4, 4, 500)

    initial_density = np.conj(initial_state(x)) * initial_state(x)
    d1 = np.conj(time_evolving(x, t=0)) * time_evolving(x, t=0)
    d2 = np.conj(time_evolving(x, t=0.3)) * time_evolving(x, t=0.3)

    # They have different complex components
    assert not np.allclose(time_evolving(x, t=0.0), time_evolving(x, t=0.3))

    # But the density is constant over time
    assert np.allclose(initial_density, d1, atol=SMALL_NUMBER)
    assert np.allclose(initial_density, d2, atol=SMALL_NUMBER)


def test_eigenbasis_transformation_orthogonal():
    basis = [lambda _: 2 for _ in range(2)]

    d = 0.1
    e = 1.0

    H = np.array([[e, d], [d, e]])
    S = np.eye(2)

    eigb = EigenBasis.from_basis(basis, H, S)

    xform_H = eigb.transformed(H)

    assert is_diagonal(xform_H)
    assert np.allclose(np.diag(xform_H), eigb.energies, atol=SMALL_NUMBER)


def test_eigenbasis_transformation_non_orthogonal():

    basis = [HarmonicOscillator(n=i, center=0.25) for i in range(2)] + [
        HarmonicOscillator(n=i, center=-0.25) for i in range(2)
    ]

    S = pairwise_array_from_func(basis, overlap)
    H = pairwise_array_from_func(basis, Hamiltonian(Harmonic(center=0.0)))
    eigb = EigenBasis.from_basis(basis, H, S)
    xform_H = eigb.transformed(H)

    assert is_diagonal(xform_H)
    assert np.allclose(np.diag(xform_H), eigb.energies, atol=SMALL_NUMBER)
