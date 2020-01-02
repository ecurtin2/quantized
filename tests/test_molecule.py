import numpy as np

from transit_chem.atom import Atom
from transit_chem.molecule import Molecule


def test_rotate():
    atom = Atom(1, 0, 2, 0)
    r = atom.rotation_matrix_to(x=4, y=0, z=0)
    rotated = atom.rotated(r)

    expected = Atom(1, 2, 0, 0)
    assert rotated == expected

    m = Molecule([atom])
    assert m.with_atom_aligned_to(atom, x=1.0, y=0.0, z=0.0) == Molecule([expected])
    assert m.rotated_about_z(np.pi / 2) == Molecule(atoms=[Atom(1, -2.0, 0.0, 0.0)])
