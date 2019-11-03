import numpy as np
import pytest

from transit_chem.basis import HarmonicOscillator
from transit_chem.utils import pairwise_array_from_func
from transit_chem.oneD import TripleWellPotential
from transit_chem import operators as op


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


def test_ho_eigen_basis_overlap_is_diagonal(ho_eigen_basis):
    S = pairwise_array_from_func(ho_eigen_basis, op.overlap)
    assert is_diagonal(S)
    assert is_identity(S)


def test_ho_eigen_basis_kinetic(ho_eigen_basis):
    K = pairwise_array_from_func(
        ho_eigen_basis, op.kinetic, symmetric=True
    )
    energies = np.asarray([b.energy for b in ho_eigen_basis])
    expected = energies / 2.0
    assert np.allclose(np.diag(K), expected)


def test_ho_eigen_basis_potential(ho_eigen_basis):
    V = pairwise_array_from_func(
        ho_eigen_basis, op.Potential(ho_eigen_basis[0].potential), symmetric=True
    )
    energies = np.asarray([b.energy for b in ho_eigen_basis])
    expected = energies / 2.0
    assert np.allclose(np.diag(V), expected)


def test_ho_eigen_basis_hamiltonian(ho_eigen_basis):
    H = pairwise_array_from_func(
        ho_eigen_basis, op.Hamiltonian(ho_eigen_basis[0].potential), symmetric=True
    )
    energies = np.asarray([b.energy for b in ho_eigen_basis])

    assert np.allclose(np.diag(H), energies)
    assert is_diagonal(H)

