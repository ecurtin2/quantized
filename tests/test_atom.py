from transit_chem.elements import H
from transit_chem.atom import Atom


def test_atom__init__():
    h1 = Atom("H", 0, 0, 0)
    h2 = Atom(H, 0, 0, 0)
    assert h1 == h2


def test_atom_mass():
    assert Atom(H, 0, 0, 0).mass == 1
    for i in range(1, 100):
        assert Atom(i, 0, 0, 0).mass == i


def test_atom_translated():
    a = Atom(H, 1, 2, 3)
    assert a.translated(x=2.0, y=3.0, z=4.0) == Atom(H, 3.0, 5.0, 7.0)


def test_atom_scaled():
    a = Atom(H, 1, 2, 4)
    assert a.scaled(2.0) == Atom(H, 2.0, 4.0, 8.0)


def test_atom_rotated():
    a = Atom(H, 0.0, 2.0, 0.0)
    r = a.rotation_matrix_to(1, 0, 0)
    assert a.rotated(r) == Atom(H, 2.0, 0.0, 0.0)
