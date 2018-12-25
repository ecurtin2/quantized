import os
import numpy as np

from . import units


class Atom(object):
    """Atom class containing coordinates, basis and mass."""
    dir = os.path.dirname(__file__)
    fname = os.path.join(dir, "elements.dat")
    with open(fname) as f:
        splitlines = [line.split(' ') for line in f.readlines()]
        Z, symbols = zip(*splitlines)
        Z = [int(i) for i in Z]
        symbols = [s.strip() for s in symbols]

    ZtoSymbol = {z: symbol for z, symbol in zip(Z, symbols)}
    SymboltoZ = {symbol: z for z, symbol in zip(Z, symbols)}
    minimal_pz_orbital_coeff_dict = {
                          1: None,
                          6: 1.5679,
                          7: 1.9170,
                          8: 2.2266,
                         }

    # Atomic number to Valence Orbital Ionization Energy in eV
    VOIEdict_eV = {
                1: None,
                6: 10.77,
                7: 13.19,
                8: 15.80,
               }

    def __init__(self, atom_type, coords, basistype='minimalpz'):
        """Constructor for Atom class. Currently supports only minimal pz basis

        :param atom_type: Atom symbol or atomic number.
        :type atom_type: str or int
        :param coords: (x, y, z) of atom.
        :type coords: iterable
        :param basistype: String describing the basis. Only 'minimalpz' implemented
        so far.
        :type basistype: string
        """

        self.Z = self.interprete_atom_type(atom_type)
        self.coords = coords
        self.alpha = self.minimal_pz_orbital_coeff_dict[self.Z]
        self.basis_type = basistype
        self.VOIE = self.get_voie()

    @property
    def symbol(self):
        return __class__.ZtoSymbol[self.Z]

    @property
    def x(self):
        return self.coords[0]

    @property
    def y(self):
        return self.coords[1]

    @property
    def z(self):
        return self.coords[2]

    def translate(self, x, y, z):
        self.coords = (self.coords[0] + x, self.coords[1] + y, self.coords[2] + z)

    def scale(self, factor):
        self.coords = tuple(factor * x for x in self.coords)

    def rotate(self, R):
        coords = np.asarray(self.coords)
        coords = R @ coords
        self.coords = tuple(coords)

    def flip_x(self):
        self.coords = (-self.coords[0], self.coords[1], self.coords[2])

    def interprete_atom_type(self, z):
        """Convert atom_type to atomic number. Accepts symbol or number.

        :param z: Atom symbol or atomic number.
        :type z: str or int

        """
        try:
            z = int(z)
        except ValueError:
            try:
                z = self.SymboltoZ[z]
            except:
                raise ValueError('Atom type not understood.')
        return z

    def get_voie(self):
        """Return the Valence Orbital Ionization Energy in Hartrees."""
        voie = self.VOIEdict_eV[self.Z]
        if voie is not None:
            return units.conversion_factors['EV_TO_HARTREE'] * voie
        else:
            return None

    def __str__(self):
        s = "{} {:10.5f} {:10.5f} {:10.5f}".format(self.symbol, *self.coords)
        return s


def main():
    a = Atom('H', coords=(4, 5, 6))
    print(a)


if __name__ == "__main__":
    main()
