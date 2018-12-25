import numpy as np

from . import basis
from . import transittimeanalyzer
from . import utils


class ThreeD(transittimeanalyzer.TransitTime):

    def __init__(self, nregions, boundaries, molecule, basis, id='unlabeled'):
        """
        :param boundaries: A list of the x coordinates of planes dividing the regions.
        :type boundaries: iterable

        """
        super(ThreeD, self).__init__(nregions)  # Call base class __init__
        self.Ndim = 3
        self.molecule = molecule

        #  These are the limits of integration in x for 'all space'. The basis class takes care of the other
        #  dimensions, but has no easy way to do it in x due to the integrals over the various regions.
        #
        self.neg_infty = -3 + min(atom.x for atom in self.molecule)
        self.pos_infty =  3 + max(atom.x for atom in self.molecule)

        bounds = sorted([self.neg_infty, self.pos_infty] + list(boundaries))
        self.boundaries = list(utils.pairwise(bounds))

        self.basis = basis
        self.Nbasis = len(self.basis)
        self.basis_funcs = self.basis
        self.initial_state = None
        self.id = id

    def get_H(self):
        """Create the Hamiltonian under the Extended Hueckel Approximation."""
        self.get_S()
        self.H = np.zeros((self.Nbasis, self.Nbasis))

        for i in range(self.Nbasis):
            self.H[i, i] = - self.molecule.atoms[i].VOIE

        for i in range(self.Nbasis):
            for j in range(self.Nbasis):
                if i != j:
                    self.H[i, j] = (1.75 * (self.H[i, i] + self.H[j, j]) * self.S[i, j] / 2.0)

    def get_S(self):
        """Create the overlap matrix"""
        self.S = np.zeros((self.Nbasis, self.Nbasis))
        for i in range(self.Nbasis):
            for j in range(i + 1):
                self.S[i, j] = self.S[j, i] = self.basis[i].overlap(self.basis[j])

        if not np.all(np.isclose(np.diag(self.S), 1.0, atol=1e-2)):
            print(self.S)
            raise ValueError('Basis is not normalized!')

    def create_mo_func(self, i):
        """Create the i'th MO function, f(x, y, z)."""
        def mo_func(x, y, z):
            val = 0.0
            for j in range(self.Nbasis):
                val += self.basis_funcs[j](x, y, z) * self.eigvecs[j, i]
            return val
        return mo_func

    def get_eigen_basis(self):
        """Create the eigenbasis as a list of functions."""
        self.eigenbasis = []
        for i in range(self.Nbasis):
            mo_func = self.create_mo_func(i)
            self.eigenbasis.append(mo_func)

    def calc_partial_overlap(self):
        """Calculate the partial overlap integrals in the MO basis"""
        self.S_partial = np.zeros((self.Nregions, self.Nbasis, self.Nbasis))
        S_partial_ao = np.zeros((self.Nregions, self.Nbasis, self.Nbasis))

        # Integrate the overlap in the AO basis over each region.
        # S is symmetric so only calculate the upper triangular form and copy.
        for n_region in range(self.Nregions):
            for i in range(self.Nbasis):
                for j in range(i + 1):
                    xmin, xmax = self.boundaries[n_region]
                    Sij = self.basis[i].overlap(self.basis[j], xmin=xmin, xmax=xmax)
                    S_partial_ao[n_region, i, j] = S_partial_ao[n_region, j, i] = Sij

        assert np.allclose(np.sum(S_partial_ao, axis=0), self.S, atol=1e-2)

        # Transform from ao -> mo basis
        for n_region in range(self.Nregions):
            self.S_partial[n_region] = self.matrix_ao_to_mo(S_partial_ao[n_region])

        #  Adding up partial overlaps in MO basis should give the identity matrix.
        #  Since eigenfunctions of hermitian operators form an orthogonal set.
        S = np.sum(self.S_partial, axis=0)
        assert np.allclose(S, np.eye(*S.shape), atol=1e-2)

    def func_to_ao(self, func):
        """Given a vector of values at each grid point, return coefficients in AO basis.

            Currently only works for PzBasis classes.
        """
        ao_expansion_coeffs = np.zeros(self.Nbasis)

        for i, basis in enumerate(self.basis):
            ao_expansion_coeffs[i] = basis.overlap(func)
        ao_expansion_coeffs /= np.linalg.norm(ao_expansion_coeffs)
        return ao_expansion_coeffs

    def func_to_mo(self, func):
        """Given a 3D function f(x, y, z), return its expansion coefficients in MO basis"""
        ao_expansion_coeffs = self.func_to_ao(func)
        mo_expansion_coeffs = self.ao_to_mo(ao_expansion_coeffs)
        return mo_expansion_coeffs

    def set_initial_state(self, initial_state):
        """Set the initial state of the time-dependent calculation.

        :param initial_state: 3D function of cartesian coordinates f(x, y, z)
        :type initial_state: function
        """
        self.initial_state = initial_state
        self.state_mo_coefficients = self.func_to_mo(self.initial_state)

