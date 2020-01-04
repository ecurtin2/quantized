import numpy as np
from pytest import raises

from transit_chem.operators import ExtendedHuckelHamiltonian
from transit_chem.molecule import Molecule
from transit_chem.atom import Atom
from transit_chem.elements import C, N, Uuo
from transit_chem.utils import isclose

from utils import is_diagonal


def test_extended_huckel_hamiltonian_if_s_orthogonal():
    S = np.eye(2)
    m = Molecule([Atom(C, 0, 0, 0), Atom(N, 0, 0, 0)])

    h = ExtendedHuckelHamiltonian(S=S, molecule=m).matrix()

    assert is_diagonal(h)

    assert h[0, 0] == -C.voie
    assert h[1, 1] == -N.voie


def test_extended_huckel_hamiltonian_s_ones():
    S = np.ones((2, 2))

    m = Molecule([Atom(C, 0, 0, 0), Atom(N, 0, 0, 0)])

    h = ExtendedHuckelHamiltonian(S=S, molecule=m).matrix()

    assert h[0, 0] == -C.voie
    assert h[1, 1] == -N.voie
    expected = 1.75 * (-C.voie + -N.voie) / 2.0
    assert isclose(h[1, 0], expected)
    assert isclose(h[0, 1], expected)


def test_extended_huckel_hamiltonian_raises_if_missing_voie():
    m = Molecule([Atom(Uuo, 0, 0, 0)])
    S = np.eye(1)
    with raises(ValueError):
        ExtendedHuckelHamiltonian(S, m)
