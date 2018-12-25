import numpy as np

from transit_chem.utils import parabola_from_points


def triple_well_potential(center1, barrier12, center2, barrier23, center3):
    """Construct a well potential from the points.

    x                                                         x
    x                                                         x
    x                                                         x
    x                                                        x
    x                                                        x
     x          barrier12                                   x
     x         x         x            barrier23            x
      x       x          x           x         x          x
       center1            x         x           x       x
                           x       x             center3
                            center2

    Returns
    --------
    callable
        A function f(x) that returns the y value of the potential at coordinate x.

    Examples
    ---------

    >>> f = triple_well_potential(
    ...     center1=(0, 0),
    ...     barrier12=(0.5, 1),
    ...     center2=(2, 0.5),
    ...     barrier23=(3, 1.5),
    ...     center3=(4, 1.0),
    ... )
    >>> f(0)
    array(0.)
    >>> f(3)
    array(1.5)
    >>> f([1, 2, 3])
    array([0.56666667, 0.5       , 1.5       ])

    """

    # Wells are in expected left to right ordering.
    if not center1[0] < barrier12[0] < center2[0] < barrier23[0] < center3[0]:
        raise ValueError('Points are not in ascending x-value order.')

    # Wells are below the barriers
    assert center1[1] < barrier12[1]
    assert center2[1] < barrier12[1]
    assert center2[1] < barrier23[1]
    assert center3[1] < barrier23[1]

    def fit_well(center, barrier):
        center_x, center_y = center
        barrier_x, barrier_y = barrier
        # Third point is reflecting barrier about center
        x = - barrier_x + 2*center_x
        y = barrier_y
        return np.poly1d(parabola_from_points(center, barrier, (x, y)))

    well1 = fit_well(center1, barrier12)
    well2 = np.poly1d(parabola_from_points(barrier12, center2, barrier23))
    well3 = fit_well(center3, barrier23)

    # Pre fetch to avoid repeated lookup
    barrier12_x = barrier12[0]
    barrier23_x = barrier23[0]

    def potential(x):
        if x < barrier12_x:
            return well1(x)
        elif barrier12_x <= x <= barrier23_x:
            return well2(x)
        elif x > barrier23_x:
            return well3(x)

    # vectorize so it can work on numpy arrays
    return np.vectorize(potential)


def centers_barriers_from_params(well1_depth, well1_halfwidth,
                                 bridge_length, bridge_depth,
                                 well3_halfwidth, well3_depth):
    center1 = [0, 0]
    barrier1 = [center1[0] + well1_halfwidth, center1[1] + well1_depth]
    center2 = [barrier1[0] + 0.5 * bridge_length, barrier1[1] - bridge_depth]
    barrier2 = [barrier1[0] + bridge_length, barrier1[1]]
    center3 = [barrier2[0] + well3_halfwidth, barrier2[1] - well3_depth]

    centers = np.asarray([center1, center2, center3])
    barriers = np.asarray([barrier1, barrier2])
    return centers, barriers

#
#     def get_eigen_basis(self):
#         """Create the normalized eigen basis from the AO basis and the eigenvectors."""
#         self.eigenbasis = [np.zeros_like(self.basis_vec[0]) for _ in range(self.Nbasis)]
#         for j in range(self.Nbasis):
#             for i in range(self.Nbasis):
#                 self.eigenbasis[j] += self.eigvecs[i, j] * self.basis_vec[i]
#         self.eigenbasis = [i * self.normalize(i) for i in self.eigenbasis]
#
#     def get_H(self):
#         """Calculate the hamiltonian"""
#         self.V = np.zeros((self.Nbasis, self.Nbasis))  # Potential
#         for i in range(self.Nbasis):
#             for j in range(i + 1):
#                 Vij = self.integrate(self.basis_vec[j] * self.potential_array
#                                     *self.basis_vec[i])
#                 self.V[i, j] = self.V[j, i] = Vij
#
#         self.T = np.zeros((self.Nbasis, self.Nbasis))
#         for i in range(self.Nbasis):
#             for j in range(self.Nbasis):
#                 # < i | T | j > = - 0.5 p*p, apply p in each direction since it is Hermitian
#                 # the finite difference behaves better this way
#                 di = utils.np_finite_diff(self.basis_vec[i], self.coords, 1)
#                 dj = utils.np_finite_diff(self.basis_vec[j], self.coords, 1)
#                 # Duplicate the first value so the length is the same.
#                 # This is good approximation with smooth function and
#                 # dense enough grid.
#                 di = np.insert(di, 0, di[0])
#                 dj = np.insert(dj, 0, dj[0])
#                 # negative sign not here since - 0.5 * p * p = -0.5 * i(dxi)*i(dxj) = 0.5(dxi)(dxj)
#                 self.T[i, j] = (0.5 / self.mass) * self.integrate(di * dj)
#
#         self.H = self.V + self.T
#
#     def get_S(self):
#         """Calculate the overlap matrix."""
#         self.S = np.zeros((self.Nbasis, self.Nbasis))
#         for i in range(self.Nbasis):
#             for j in range(i + 1):
#                 Sij = self.overlap(self.basis_vec[j], self.basis_vec[i])
#                 self.S[i, j] = self.S[j, i] = Sij
#
#     def calc_partial_overlap(self):
#         """Calculate the partial overlap integrals in the MO basis."""
#         self.S_partial = np.zeros((self.Nregions, self.Nbasis, self.Nbasis))
#         for n_region in range(self.Nregions):
#             for i in range(self.Nbasis):
#                 for j in range(i + 1):
#                     Sij = self.region_integrators[n_region](
#                         self.eigenbasis[i] * self.eigenbasis[j]
#                     )
#                     self.S_partial[n_region, i, j] = self.S_partial[n_region, j, i] = Sij
#         return self.S_partial
#
#     def vec_to_ao(self, vec):
#         """Given a vector of values at each grid point, return coefficients in AO basis."""
#         vec *= self.normalize(vec)
#         c = np.zeros(self.Nbasis)
#         for i in range(self.Nbasis):
#             c[i] = np.dot(vec, self.basis_vec[i])
#         c /= np.linalg.norm(c)
#         vec_in_ao_basis = np.zeros(len(self.coords))
#         for i in range(self.Nbasis):
#             vec_in_ao_basis += c[i] * self.basis_vec[i]
#         return c, vec_in_ao_basis
#
#     def vec_to_mo(self, vec):
#         """Given a vector of values at each grid point, return coefficients in MO basis.
#
#         :param vec: 1D Vector of values along the self.coords.
#         :type vec: np.ndarray
#         """
#         vec *= self.normalize(vec)
#         c = np.zeros(self.Nbasis)
#         for i in range(self.Nbasis):
#             c[i] = np.dot(vec, self.eigenbasis[i])
#         c /= np.linalg.norm(c)
#         vec_in_mo_basis = np.zeros(len(self.coords))
#         for i in range(self.Nbasis):
#             vec_in_mo_basis += c[i] * self.eigenbasis[i]
#         return c, vec_in_mo_basis * self.normalize(vec_in_mo_basis)

#
#     def normalized(self, vec):
#         """Return a normalized version of given vector.
#
#         Vec is not modified in place.
#
#         :param vec: 1D Vector of values along the self.coords.
#         :type vec: np.ndarray
#         """
#         return 1.0 / math.sqrt(self.overlap(np.conjugate(vec), vec))
#
#     def time_evolve_vec(self, vec):
#         """Given a vector return the time-dependence of the vector as a complex matrix. Time is an input vector
#         of times. Each row of the matrix is the wave function's position vector at the corresponding
#         point in time in the t vector given."""
#         vec *= self.normalize(vec)
#         c, _ = self.vec_to_mo(vec)
#         td = np.zeros((len(self.times), len(self.coords)), dtype=np.complex128)
#         for i in range(self.Nbasis):
#             time_part = np.exp(-1j * self.eigvals[i] * self.times)
#             td += c[i] * np.outer(time_part, self.eigenbasis[i])
#         return td
#
#     def set_initial_state(self, initialstate):
#         """Set the initial state for the time dependent solution.
#
#         :param initialstate: 1D Vector of values along the self.coords.
#         :type initialstate: np.ndarray
#         """
#         initialstate *= self.normalize(initialstate)
#         self.initial_state = initialstate
#         self.state_mo_coefficients = self.vec_to_mo(initialstate)[0]
#         self.initial_state_energy = np.dot(self.state_mo_coefficients**2, self.eigvals)
#
#     def time_evolve_initialstate(self, state):
#         self.time_evolved_state = self.time_evolve_vec(state)
#         self.time_evolved_prob_density = self.time_evolved_state.conj() * self.time_evolved_state
