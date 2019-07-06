from abc import abstractmethod, abstractproperty, ABCMeta
import itertools
import json
import logging

import numpy as np
from scipy import linalg

from transittime import utils
from transittime.log_helper import log_block, TimeLevel


mod_log = logging.getLogger(__name__)


class TransitTime(metaclass=ABCMeta):
    """The transit time Class. This is the base class from which OneD and
    ThreeD inherit. The time-dependence is all in this class. The details
    of spatially dependent quantities, like potentials and integrators are
    defined in the child classes.

    Metaclass=ABCMeta defines this class as an Abstract base class. This means it can't be directly
    instantiated and only one of the subclasses can be instantiated. This is to protect you from
    yourself.

    """
    @log_block('Transittime Init', mod_log, logging.DEBUG)
    def __init__(self, nregions):
        """Initializer for Transit time object.

        :param nregions: Number of spatially separated regions to consider.
        :type nregions: int
        """
        self.Nregions = nregions
        self.times = np.zeros(0)
        self.OccProb = np.zeros((0, 0))
        self.ChangeMat = np.zeros((0, 0, 0))
        self.S_partial = np.zeros((0, 0, 0))
        self.S = np.zeros((0, 0))
        self.T = np.zeros((0, 0))
        self.V = np.zeros((0, 0))
        self.H = np.zeros((0, 0))
        self.eigvals = np.zeros((0, 0))
        self.eigvecs = np.zeros((0, 0))
        self.eigenbasis = np.zeros((0, 0))
        self.P_not = np.zeros(0)
        self.Po = np.zeros(0)
        self.basis_vec = []
        self.Nbasis = 0
        self.state_mo_coefficients = None

    @log_block("Calculating both time independent and time-dependent quantities.",
               mod_log, logging.INFO, TimeLevel.ALWAYS)
    def calc_all(self, initialstate):
        """Calculate both time independent and time-dependent quantities.

        :param initialstate: The initial state of the system.
        :type initialstate: (3D) function of cartesian coordinates f(x, y, z).
                            (1D) Array containing the values of initial state on a grid.
        :param times: NumPy array of time steps to be used.
        :type times: ndarray
        """
        self.calc_time_independent()
        self.set_initial_state(initialstate)
        self.calc_time_dependent()

    @log_block("Calculate all the time independent quantities", mod_log, logging.INFO, TimeLevel.OFTEN)
    def calc_time_independent(self):
        """Calculate all the time independent quantities.

        Calculates the eigenbasis as well as partial overlap matrix.
        """
        self.diagonalize_H()
        self.get_eigen_basis()
        self.calc_partial_overlap()

    @log_block("Calculate all the time dependent quantities", mod_log, logging.INFO, TimeLevel.OFTEN)
    def calc_time_dependent(self, guess_t_max=500, p_not_threshold=0.05, stepsize=1.0, maxsteps=10**4):
        """Calculate the time dependent quantities. Must be after time-independent is done.

        Calculates occupancy probabilities, the change matrix and pnot.

        """

        # Recalculate until P_not < threshold or  n_steps = t_max // stepsize >
        p_not = 1.0
        t_max = guess_t_max
        n_steps = 1
        while (p_not >= p_not_threshold) and (n_steps < maxsteps):
            n_steps = round(t_max / stepsize)
            times = np.linspace(0, t_max, n_steps)
            self.times = times
            self.calc_occupancy_probability()
            self.get_change_mat()
            self.get_p_not(3)
            p_not = self.P_not[-2]
            t_max *= 2

    def set_initial_state(self, initialstate): raise NotImplementedError

    def calc_partial_overlap(self): raise NotImplementedError

    def normalize(self, vec): raise NotImplementedError

    def get_S(self): raise NotImplementedError

    def get_H(self): raise NotImplementedError

    def vec_to_ao(self, vec): raise NotImplementedError

    def vec_to_mo(self, vec): raise NotImplementedError

    def propagate_state(self, state): raise NotImplementedError

    def get_eigen_basis(self): raise NotImplementedError

    @log_block("Diagonalize the Hamiltonian", mod_log, level=logging.INFO, time_level=TimeLevel.RARELY)
    def diagonalize_H(self):
        """Calculate H, S, V, T and diagonalize. Form eigenbasis.

        """
        # Get matrices
        self.get_S()
        self.get_H()

        # Diagonalize and sort
        self.eigvals, self.eigvecs = linalg.eigh(self.H, self.S)
        idx = np.argsort(self.eigvals)
        self.eigvals = self.eigvals[idx]
        self.eigvecs = self.eigvecs[:, idx]

        #  Normalize eigenvectors
        norm = np.linalg.norm(self.eigvecs, axis=0)
        self.eigvecs /= norm[None, :]

    @log_block('ao to mo transformation', mod_log, level=logging.DEBUG)
    def ao_to_mo(self, ao_vec):
        """Given a vector of coefficients in AO basis, return coefficients in MO basis."""
        mo = np.dot(self.eigvecs.T, ao_vec)
        mo /= np.linalg.norm(mo)
        return mo

    @log_block('Mo to ao transformation', mod_log, level=logging.DEBUG)
    def mo_to_ao(self, mo_vec):
        """Given a vector of coefficients in MO basis, return coefficients in AO basis."""
        ao = np.dot(self.eigvecs, mo_vec)
        ao /= np.linalg.norm(ao)
        return ao

    @log_block('Matrix ao to mo transformation', mod_log, level=logging.DEBUG)
    def matrix_ao_to_mo(self, matrix):
        inv = np.linalg.inv
        return inv(self.eigvecs) @ inv(self.S) @ matrix @ self.eigvecs

    @log_block('Determine the occupancy probability', mod_log, level=logging.INFO)
    def calc_occupancy_probability(self):
        """Calculate the time-dependent occupancy probabilities"""
        c_conjugate = []
        c = []
        for i in range(self.Nbasis):
            c.append(self.state_mo_coefficients[i] * np.exp(-1j * self.eigvals[i] * self.times))
            c_conjugate.append(self.state_mo_coefficients[i] * np.exp(1j * self.eigvals[i] * self.times))

        occ_prob = np.zeros((self.Nregions, len(self.times)), dtype=np.complex128)

        for i in range(self.Nregions):
            for n in range(self.Nbasis):
                for m in range(self.Nbasis):
                    occ_prob[i] += c_conjugate[n] * c[m] * self.S_partial[i, n, m]
        self.OccProb = np.abs(np.real(occ_prob))
        self.OccProb /= np.sum(self.OccProb, axis=0)  # Normalize
        return self.OccProb

    @log_block("Determine the change/transition matrix", mod_log, level=logging.INFO)
    def get_change_mat(self):
        """calculate the change matrix over time, currently supports 3 regions"""
        A = np.zeros((len(self.times) - 1, self.Nregions, self.Nregions))  # change matrix
        dOccProb = [np.diff(i) for i in self.OccProb]  # change for each region
        sign_of_change = [np.sign(i) for i in dOccProb]
        two_states_increasing = (sign_of_change[0]
                                 * sign_of_change[1]
                                 * sign_of_change[2] < 0)

        for t in range(len(self.times) - 1):
            if two_states_increasing[t]:
                # If the state is decreasing in probability, the diagonal is
                # determined
                for j in range(3):
                    if sign_of_change[j][t] < 0:
                        A[t, j, j] = self.OccProb[j, t + 1] / self.OccProb[j, t]
                        for k in range(3):
                            if k != j:
                                A[t, k, j] = dOccProb[k][t] / self.OccProb[j, t]

                    # If the state is increasing, diagonal is 1. Off diagonal
                    # columns = 0
                    elif sign_of_change[j][t] > 0:
                        A[t, j, j] = 1
                        for k in range(3):
                            if k != j:
                                A[t, k, j] = 0

            elif not two_states_increasing[t]:
                for j in range(3):
                    # If the state is decreasing in probability, the diagonal
                    # is determined
                    if sign_of_change[j][t] < 0:
                        A[t, j, j] = self.OccProb[j, t + 1] / self.OccProb[j, t]
                        for k in range(3):
                            if k != j:
                                A[t, j, k] = 0

                # If the state is increasing, diagonal is 1. Off diagonal
                # column element is 1 - diagonal
                for j in range(3):
                    if sign_of_change[j][t] > 0:
                        A[t, j, j] = 1
                        for k in range(3):
                            if k != j:
                                A[t, j, k] = 1 - A[t, k, k]

            else:  # just in case weird equivalence happens, repeat last time step
                A[t, :, :] = A[t - 1, :, :]
        self.ChangeMat = A
        return A

    @log_block('Determine p_not - this gets called multiple times until it finds desired threshold',
               mod_log, level=logging.INFO, time_level=TimeLevel.OFTEN)
    def get_p_not(self, acceptor):
        """Determine P_not, the probability than acceptor state has never been occupied.

        :param acceptor: index of the state to consider.
        :type acceptor: int
        TODO: GET ACCEPTOR/DONOR
        """
        self.Po = self.OccProb[:, 0]

        #  Make mask to remove the acceptor rows and columns
        mask = np.full_like(self.ChangeMat, True, dtype=bool)
        mask[:, acceptor - 1, :] = False
        mask[:, :, acceptor - 1] = False
        newshape = self.ChangeMat.shape
        newshape = newshape[0], newshape[1] - 1, newshape[2] - 1
        A_tilde = self.ChangeMat[mask].reshape(newshape)

        P_o_tilde = self.Po[:2]
        #  Taking [a, b, c, d, ...] -> [(b @ a), (c @ b @ a), (d @ c @ b @ a), ...] and making into an array.
        a_prod = np.asarray(list(itertools.accumulate(A_tilde, lambda a, b: b @ a)))
        #  Broadcasting rules means that each 2x2 submatrix of a_prod matrix multiplies P_o_tilde
        self.P_not = np.sum(a_prod @ P_o_tilde, axis=1)
