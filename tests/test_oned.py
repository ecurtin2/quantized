from math import isclose

import numpy as np
import pytest

from transit_chem import operators as op
from transit_chem.basis import HarmonicOscillator
from transit_chem.config import conf
from transit_chem.potentials import TripleWell
from utils import is_diagonal, is_hermitian, is_identity


@pytest.fixture
def ho_eigen_basis():
    return [
        HarmonicOscillator(n=0, center=0),
        HarmonicOscillator(n=1, center=0),
        HarmonicOscillator(n=2, center=0),
        HarmonicOscillator(n=3, center=0),
    ]


def test_ho_eigen_basis_overlap_is_diagonal(ho_eigen_basis):
    S = op.Overlap().matrix(ho_eigen_basis)
    assert is_diagonal(S)
    assert is_identity(S)


def test_ho_eigen_basis_kinetic(ho_eigen_basis):
    K = op.Kinetic().matrix(ho_eigen_basis)
    energies = np.asarray([b.energy for b in ho_eigen_basis])
    expected = energies / 2.0
    assert np.allclose(np.diag(K), expected)
    assert is_hermitian(K)


def test_ho_eigen_basis_potential(ho_eigen_basis):
    V = op.Potential(ho_eigen_basis[0].potential).matrix(ho_eigen_basis)

    energies = np.asarray([b.energy for b in ho_eigen_basis])
    expected = energies / 2.0
    assert np.allclose(np.diag(V), expected)
    assert is_hermitian(V)


def test_ho_eigen_basis_hamiltonian(ho_eigen_basis):
    H = op.Hamiltonian(ho_eigen_basis[0].potential).matrix(ho_eigen_basis)
    energies = np.asarray([b.energy for b in ho_eigen_basis])

    assert np.allclose(np.diag(H), energies)
    assert is_diagonal(H)


def test_triple_well():
    w1d = 1
    w1h = 2
    bl = 5
    bd = 1
    w3w = 1.5
    w3d = 0.5
    v = TripleWell.from_params(
        well1_depth=w1d,
        well1_halfwidth=w1h,
        bridge_length=bl,
        bridge_depth=bd,
        well3_halfwidth=w3w,
        well3_depth=w3d,
    )

    assert isclose(v(0), 0.0, abs_tol=conf.small_number)
    assert isclose(v(w1h), w1d, abs_tol=conf.small_number)
    assert isclose(v(w1h + bl / 2.0), w1d - bd, abs_tol=conf.small_number)
    assert isclose(v(w1h + bl), w1d, abs_tol=conf.small_number)
    assert isclose(v(w1h + bl + w3w), w1d - w3d, abs_tol=conf.small_number)

    # Potential should go as x^2, so bigger than x at large x
    assert v(conf.large_number) > conf.large_number
    assert v(-conf.large_number) > conf.large_number


def test_triple_well_runs_on_numpy_array():
    w1d = 1
    w1h = 2
    bl = 5
    bd = 1
    w3w = 1.5
    w3d = 0.5
    x = np.linspace(-10, 10, 100)
    v = TripleWell.from_params(
        well1_depth=w1d,
        well1_halfwidth=w1h,
        bridge_length=bl,
        bridge_depth=bd,
        well3_halfwidth=w3w,
        well3_depth=w3d,
    )
    v(x)
