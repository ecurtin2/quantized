from typing import Tuple

import attr

from transit_chem.utils import Parabola


@attr.s()
class TripleWellPotential:
    center1: Tuple[float, float] = attr.ib()
    barrier12: Tuple[float, float] = attr.ib()
    center2: Tuple[float, float] = attr.ib()
    barrier23: Tuple[float, float] = attr.ib()
    center3: Tuple[float, float] = attr.ib()
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
    """

    def __attrs_post_init__(self):
        # Wells are in expected left to right ordering.
        if (
            not self.center1[0]
            < self.barrier12[0]
            < self.center2[0]
            < self.barrier23[0]
            < self.center3[0]
        ):
            raise ValueError("Points are not in ascending x-value order.")

        # Wells are below the barriers
        assert self.center1[1] < self.barrier12[1]
        assert self.center2[1] < self.barrier12[1]
        assert self.center2[1] < self.barrier23[1]
        assert self.center3[1] < self.barrier23[1]

        def fit_well(center, barrier):
            center_x, center_y = center
            barrier_x, barrier_y = barrier
            # Third point is reflecting barrier about center
            x = -barrier_x + 2 * center_x
            y = barrier_y
            return Parabola.from_points(center, barrier, (x, y))

        self.well1 = fit_well(self.center1, self.barrier12)
        self.well2 = Parabola.from_points(self.barrier12, self.center2, self.barrier23)
        self.well3 = fit_well(self.center3, self.barrier23)

        # Pre fetch to avoid repeated lookup
        self.barrier12_x = self.barrier12[0]
        self.barrier23_x = self.barrier23[0]

    def __call__(self, x):
        if x < self.barrier12_x:
            return self.well1(x)
        elif self.barrier12_x <= x <= self.barrier23_x:
            return self.well2(x)
        elif x > self.barrier23_x:
            return self.well3(x)

    @staticmethod
    def from_params(
        well1_depth: float,
        well1_halfwidth: float,
        bridge_length: float,
        bridge_depth: float,
        well3_halfwidth: float,
        well3_depth: float,
    ):
        center1 = (0, 0)
        barrier1 = (center1[0] + well1_halfwidth, center1[1] + well1_depth)
        center2 = (barrier1[0] + 0.5 * bridge_length, barrier1[1] - bridge_depth)
        barrier2 = (barrier1[0] + bridge_length, barrier1[1])
        center3 = (barrier2[0] + well3_halfwidth, barrier2[1] - well3_depth)
        return TripleWellPotential(center1, barrier1, center2, barrier2, center3)


#
#     def get_eigen_basis(self):
#         """Create the normalized eigen basis from the AO basis and the eigenvectors."""
#         self.eigenbasis = [np.zeros_like(self.basis_vec[0]) for _ in range(self.Nbasis)]
#         for j in range(self.Nbasis):
#             for i in range(self.Nbasis):
#                 self.eigenbasis[j] += self.eigvecs[i, j] * self.basis_vec[i]
#         self.eigenbasis = [i * self.normalize(i) for i in self.eigenbasis]
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
