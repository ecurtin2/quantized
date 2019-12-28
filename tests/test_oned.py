from math import isclose

import numpy as np
import pytest

from transit_chem import operators as op
from transit_chem.basis import HarmonicOscillator
from transit_chem.config import LARGE_NUMBER, SMALL_NUMBER
from transit_chem.oneD import TripleWellPotential
from transit_chem.utils import pairwise_array_from_func


@pytest.fixture
def ho_eigen_basis():
    return [
        HarmonicOscillator(n=0, center=0),
        HarmonicOscillator(n=1, center=0),
        HarmonicOscillator(n=2, center=0),
        HarmonicOscillator(n=3, center=0),
    ]


def is_diagonal(x: np.ndarray) -> bool:
    n, m = x.shape
    off_diags = x[~np.eye(n, m, dtype=bool)]
    return np.allclose(off_diags, 0)


def is_identity(x: np.ndarray) -> bool:
    n, m = x.shape
    return np.allclose(np.eye(n, m), x)


def is_hermitian(x: np.ndarray) -> bool:
    return np.allclose(x, np.conj(x).T)


def test_ho_eigen_basis_overlap_is_diagonal(ho_eigen_basis):
    S = pairwise_array_from_func(ho_eigen_basis, op.overlap)
    assert is_diagonal(S)
    assert is_identity(S)


def test_ho_eigen_basis_kinetic(ho_eigen_basis):
    K = pairwise_array_from_func(ho_eigen_basis, op.kinetic, symmetric=True)
    energies = np.asarray([b.energy for b in ho_eigen_basis])
    expected = energies / 2.0
    assert np.allclose(np.diag(K), expected)
    assert is_hermitian(K)


def test_ho_eigen_basis_potential(ho_eigen_basis):
    V = pairwise_array_from_func(
        ho_eigen_basis, op.Potential(ho_eigen_basis[0].potential), symmetric=True
    )
    energies = np.asarray([b.energy for b in ho_eigen_basis])
    expected = energies / 2.0
    assert np.allclose(np.diag(V), expected)
    assert is_hermitian(V)


def test_ho_eigen_basis_hamiltonian(ho_eigen_basis):
    H = pairwise_array_from_func(
        ho_eigen_basis, op.Hamiltonian(ho_eigen_basis[0].potential), symmetric=True
    )
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
    v = TripleWellPotential.from_params(
        well1_depth=w1d,
        well1_halfwidth=w1h,
        bridge_length=bl,
        bridge_depth=bd,
        well3_halfwidth=w3w,
        well3_depth=w3d,
    )

    assert isclose(v(0), 0.0, abs_tol=SMALL_NUMBER)
    assert isclose(v(w1h), w1d, abs_tol=SMALL_NUMBER)
    assert isclose(v(w1h + bl / 2.0), w1d - bd, abs_tol=SMALL_NUMBER)
    assert isclose(v(w1h + bl), w1d, abs_tol=SMALL_NUMBER)
    assert isclose(v(w1h + bl + w3w), w1d - w3d, abs_tol=SMALL_NUMBER)

    # Potential should go as x^2, so bigger than x at large x
    assert v(LARGE_NUMBER) > LARGE_NUMBER
    assert v(-LARGE_NUMBER) > LARGE_NUMBER


def test_triple_well_runs_on_numpy_array():
    x = np.linspace(-10, 10, 100)
    triple_well = TripleWellPotential(
        center1=(0, 0),
        center2=(1, 0.2),
        center3=(2, 0.1),
        barrier12=(0.5, 0.4),
        barrier23=(1.4, 0.3),
    )
    triple_well(x)
