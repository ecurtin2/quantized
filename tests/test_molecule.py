import numpy as np

from quantized.elements import H, Li
from quantized.atom import Atom
from quantized.molecule import Molecule

from utils import allclose

from pytest import raises


def test_rotate():
    atom = Atom(H, 0, 2, 0)
    r = atom.rotation_matrix_to(x=4, y=0, z=0)
    rotated = atom.rotated(r)

    expected = Atom(H, 2, 0, 0)
    assert rotated == expected

    m = Molecule([atom])
    assert m.with_atom_aligned_to(atom, x=1.0, y=0.0, z=0.0) == Molecule([expected])
    assert m.rotated_about_z(np.pi / 2) == Molecule(atoms=[Atom(H, -2.0, 0.0, 0.0)])


def test_center_of_mass():
    a1 = Atom(H, -1, 0, 0)
    a2 = Atom(H, 1, 0, 0)
    m = Molecule([a1, a2])
    assert allclose(m.center_of_mass, (0, 0, 0))


def test_molecule___len__():
    a1 = Atom(H, -1, 0, 0)
    a2 = Atom(H, 1, 0, 0)
    m = Molecule([a1, a2])
    assert len(m) == 2


def test_molecule_coords():
    a1 = Atom(H, -1, 0, 0)
    a2 = Atom(H, 1, 0, 0)
    m = Molecule([a1, a2])

    expected = [[-1, 0, 0], [1, 0, 0]]

    assert np.allclose(m.coords, expected, atol=1e-10)


def test_molecule_R():
    a1 = Atom(H, 0, 0, 0)
    a2 = Atom(H, 1, 0, 0)
    m = Molecule([a1, a2])

    expected = [[0, 1], [1, 0]]
    assert np.allclose(m.R, expected, atol=1e-10)


def test_molecule_mass():
    a1 = Atom(3, 0, 0, 0)
    a2 = Atom(25, 1, 0, 0)
    m = Molecule([a1, a2])
    assert m.mass == 28


def test_molecule_translated():
    a1 = Atom(H, 20, 0, 0)
    a2 = Atom(H, 10, 0, 0)
    m = Molecule([a1, a2])

    e1 = Atom(H, 20, 10, -3)
    e2 = Atom(H, 10, 10, -3)
    expected = Molecule([e1, e2])

    assert m.translated(y=10, z=-3) == expected


def test_molecule_from_xyz():
    expected = Molecule([Atom(Li, 3.0, 4.0, 5.0), Atom(H, 2.3, 1.3, -1.0)])
    xyz = """first line not needed
    second is comment and ignored
    3 3.0 4.0 5.0
    1 2.3 1.3 -1.0
    """
    assert Molecule.from_xyz(xyz) == expected


def test_molecule_com_as_origin():
    m = Molecule([Atom(Li, 23, 1, 2), Atom(H, 24, 8, -23)])
    translated = m.com_as_origin()
    assert allclose(translated.center_of_mass, (0, 0, 0))


def test_molecule_rotated_raises_if_wrong_shape():
    with raises(ValueError):
        Molecule([Atom(H, 1, 1, 1)]).rotated(np.array([1]))


def test_molecule_scaled():
    m = Molecule([Atom(Li, 1, 1, 1), Atom(H, -1, 10, -30)])
    expected = Molecule([Atom(Li, 10, 10, 10), Atom(H, -10, 100, -300)])
    assert m.scaled(10) == expected


def test_molecule_flipped_x():
    m = Molecule([Atom(Li, 1, 1, 1), Atom(H, -1, 10, -30)])
    expected = Molecule([Atom(Li, -1, 1, 1), Atom(H, 1, 10, -30)])
    assert m.flipped_x() == expected


def test_molecule_sorted():
    m = Molecule([Atom(Li, 1, 1, 1), Atom(H, -1, 10, -30)])
    expected = Molecule([Atom(H, -1, 10, -30), Atom(Li, 1, 1, 1)])
    assert m.sorted(lambda a: a.z) == expected
    assert m.sorted(lambda a: a.mass) == expected


def test_molecule_rotated_about_x():
    a1 = Atom(H, 5, 1, 0)
    a2 = Atom(H, 2, 0, 1)
    m = Molecule([a1, a2])
    expected = Molecule([Atom(H, 5, 0, -1), Atom(H, 2, 1, 0)])
    assert m.rotated_about_x(-np.pi / 2) == expected


def test_molecule_rotated_about_y():
    a1 = Atom(H, 1, 2, 0)
    a2 = Atom(H, 0, 6, -1)
    m = Molecule([a1, a2])
    expected = Molecule([Atom(H, 0, 2, 1), Atom(H, 1, 6, 0)])
    assert m.rotated_about_y(np.pi / 2) == expected


def test_molecule_rotated_about_z():
    a1 = Atom(H, -1, 0, 100)
    a2 = Atom(H, 0, -1, 200)
    m = Molecule([a1, a2])
    expected = Molecule([Atom(H, 0, -1, 100), Atom(H, 1, 0, 200)])
    assert m.rotated_about_z(np.pi / 2) == expected
