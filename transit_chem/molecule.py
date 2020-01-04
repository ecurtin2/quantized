from functools import partial
from typing import Callable, List

import attr
import numpy as np

from transit_chem.atom import Atom
from transit_chem.utils import pairwise_array_from_func


@attr.s(frozen=True, cmp=False)
class Molecule:
    atoms: List[Atom] = attr.ib()

    def __len__(self):
        return len(self.atoms)

    @property
    def coords(self):
        return np.asarray([a.coords for a in self.atoms])

    @property
    def R(self):
        """Calculate the distances for each atom-atom pair."""
        return pairwise_array_from_func(self.atoms, Atom.distance)

    @property
    def mass(self):
        return sum(a.mass for a in self.atoms)

    @property
    def center_of_mass(self):
        """Determine the center of mass of the molecule."""
        mx = sum(a.x * a.mass for a in self.atoms)
        my = sum(a.y * a.mass for a in self.atoms)
        mz = sum(a.z * a.mass for a in self.atoms)
        return mx / self.mass, my / self.mass, mz / self.mass

    def translated(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> "Molecule":
        f = partial(Atom.translated, x=x, y=y, z=z)
        return self.map(f)

    @staticmethod
    def from_xyz(xyz: str) -> "Molecule":
        """Convert from xyzfile and basistype into list of atoms.

        :param xyz_file: Name of xyzfile containing atoms and coordinates
        :type xyz_file: str
        :param basis_type: Descriptor for basis function. Currently supports
        only one type across all atom types.
        :type basis_type: str
        """
        text = (line.split() for line in xyz.splitlines()[2:] if line.strip())

        atoms = [Atom(element=e, x=x, y=y, z=z) for e, x, y, z in text]
        return Molecule(atoms)

    def com_as_origin(self) -> "Molecule":
        x, y, z = self.center_of_mass
        return self.translated(x=-x, y=-y, z=-z)

    def with_atom_aligned_to(self, atom: Atom, x: float, y: float, z: float) -> "Molecule":
        r = atom.rotation_matrix_to(x=x, y=y, z=z)
        return self.rotated(r)

    def rotated(self, r: np.array) -> "Molecule":
        if not r.shape == (3, 3):
            raise ValueError(f"Rotation matrix has invalid shape: {r.shape}")
        rotate = partial(Atom.rotated, r=r)
        return self.map(rotate)

    def map(self, f: Callable) -> "Molecule":
        return attr.evolve(self, atoms=[f(a) for a in self.atoms])

    def scaled(self, factor: float) -> "Molecule":
        return self.map(partial(Atom.scaled, factor=factor))

    def flipped_x(self):
        return self.map(Atom.flipped_x)

    def sorted(self, atomic_key: Callable) -> "Molecule":
        return attr.evolve(self, atoms=sorted(self.atoms, key=atomic_key))

    def rotated_about_x(self, angle: float) -> "Molecule":
        f = partial(Atom.rotated_about_x, angle=angle)
        return self.map(f)

    def rotated_about_y(self, angle: float) -> "Molecule":
        f = partial(Atom.rotated_about_y, angle=angle)
        return self.map(f)

    def rotated_about_z(self, angle: float) -> "Molecule":
        f = partial(Atom.rotated_about_z, angle=angle)
        return self.map(f)

    def __iter__(self):
        return iter(self.atoms)

    def __eq__(self, other):
        return all(a1 == a2 for a1, a2 in zip(self, other))
