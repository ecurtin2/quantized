import copy
import glob
import sys

import numpy as np

from . import atom


class Molecule(object):
    """Molecule class. Contains list of atoms as well as collective properties.

    Contains functions for rotating, aligning and scaling. Iterating over the molecule gives each atom.
    """

    def __init__(self, xyz_file, basis_type='minimalpz'):
        """Constructor for Molecule from an xyzfile.

        The xyzfile can be either of these forms:

        C 1.0  1.0  1.0
        H 2.0 2.0 2.0
        ...

        or

        6 1.0 1.0 1.0
        1 1.0 1.0 1.0
        ...

        :param xyz_file: Name of xyzfile containing atoms and coordinates
        :type xyz_file: str
        :param basis_type: Descriptor for basis function. Currently supports
        only one type across all atom types.
        :type basis_type: str
        """
        self.xyzfile = xyz_file
        self.comment = ''
        self._atoms_saved, self.comment = self.xyz_to_atoms(self.xyzfile, basis_type)
        self.atoms = copy.deepcopy(self._atoms_saved)

    @property
    def Natoms(self):
        return len(self.atoms)

    @property
    def coords(self):
        return np.asarray([atom.coords for atom in self.atoms])

    def trim_hydrogens(self):
        self.atoms = [atom for atom in self._atoms_saved if not atom.symbol == 'H']

    def show_hydrogens(self):
        self.atoms = self._atoms_saved

    def __iter__(self):
        self._gen = iter(self.atoms)
        return self

    def __next__(self):
        return next(self._gen)

    @property
    def R(self):
        """Calculate the distances for each atom-atom pair."""
        r = np.zeros((self.Natoms, self.Natoms))
        for i in range(self.Natoms):
            for j in range(self.Natoms):
                diffs = np.asarray(self.atoms[i].coords) - np.asarray(self.atoms[j].coords)
                r[i, j] = np.sqrt(np.sum(diffs**2))
        return r

    @property
    def center_of_mass(self):
        """Determine the center of mass of the molecule."""
        total_mass = 0.0
        mx = 0.0
        my = 0.0
        mz = 0.0
        for atom in self.atoms:
            total_mass += atom.Z
            mx += atom.Z * atom.coords[0]
            my += atom.Z * atom.coords[1]
            mz += atom.Z * atom.coords[2]
        return mx / total_mass, my / total_mass, mz / total_mass

    def translate(self, x, y, z):
        for atom in self.atoms:
            atom.translate(x, y, z)

    @staticmethod
    def xyz_to_atoms(xyz_file, basis_type):
        """Convert from xyzfile and basistype into list of atoms.

        :param xyz_file: Name of xyzfile containing atoms and coordinates
        :type xyz_file: str
        :param basis_type: Descriptor for basis function. Currently supports
        only one type across all atom types.
        :type basis_type: str
        """
        with open(xyz_file) as f:
            text = f.readlines()
            comment = [line.rstrip() for line in text[:2]]
            text = text[2:] # First two are number and comment lines
            atomlist = [line.rstrip().split() for line in text]
            atomtypes, coords = zip(*[(line[0], line[1:]) for line in atomlist])
            coords = [[float(c) for c in line] for line in coords]

        atoms = []
        for i, atomtype in enumerate(atomtypes):
            atoms.append(atom.Atom(atomtype, coords[i], basis_type))
        return atoms, comment

    def align_to_x(self, atom1, atom2, xy_anchor=None):
        """Align the line connecting atom1 and atom2 to the x axis and set the midpoint to the origin.

        after aligning atom1 and atom2 to x axis, rotate the molecule about the x axis such that
        the xy_anchor atom lies in the xy plane.
        """
        midpoint = [0.5 * (coord2 + coord1) for coord1, coord2 in zip(atom1.coords, atom2.coords)]
        self.translate(-midpoint[0], -midpoint[1], -midpoint[2])

        #  Rotating atom1 to be in the xz plane
        x, y, z = atom1.coords
        theta = - np.arctan(y / x)
        Rz = np.asarray(
            [
                [np.cos(theta), -np.sin(theta), 0],
                [np.sin(theta), np.cos(theta), 0],
                [     0       ,       0      , 1]
             ]
        )
        for atom in self.atoms:
            atom.rotate(Rz)

        # Rotate to put atom1 in on the +x axis from within xz plane
        x, y, z = atom1.coords
        theta = - np.arctan(z / x)
        Ry = np.asarray(
            [
                [np.cos(theta), 0, -np.sin(theta)],
                [     0       , 1,        0      ],
                [np.sin(theta), 0,  np.cos(theta)]
             ]
        )
        for atom in self.atoms:
            atom.rotate(Ry)


        # Rotate to put xy anchor on the xy plane
        if xy_anchor is not None:
            x, y, z = xy_anchor.coords
            theta = - np.arctan(z / y)
            Rx = np.asarray(
                [
                    [1,      0       ,       0       ],
                    [0, np.cos(theta), -np.sin(theta)],
                    [0, np.sin(theta),  np.cos(theta)]
                 ]
            )
            for atom in self.atoms:
                atom.rotate(Rx)



    def scale(self, factor):
        for atom in self.atoms:
            atom.scale(factor)

    def flip_x(self):
        for atom in self.atoms:
            atom.flip_x()

    def sort(self, atomic_key):
        self.atoms = sorted(self.atoms, key=atomic_key)

    def __str__(self):
        my_list = self.comment + [str(atom) for atom in self.atoms]
        s = '\n'.join(my_list)
        return s


def main(argv):
    for f in glob.glob("../analysis/xyz_files/*.xyz"):
        M = Molecule(f)
        M.trim_hydrogens()
        nitrogens = [atom for atom in M.atoms if atom.Z == 7]
        if nitrogens:
            M.align_to_x(*nitrogens)

        #  Align the nitro (acceptor) group to be +x
        oxygens = [atom for atom in M.atoms if atom.Z == 8]
        if oxygens:
            if all(oxygen.coords[0] < 0 for oxygen in oxygens):
                M.flip_x()
        M.sort(atomic_key=lambda atom: atom.coords[0])
        print(M)
        print(M.coords)


if __name__ == '__main__':
    main(sys.argv)
