from functools import partial

import numpy as np
from loguru import logger
from matplotlib import pyplot as plt

from transit_chem.basis import (
    EigenBasis,
    HarmonicOscillator,
    TimeEvolvingState,
    harmonic_basis_from_parabola,
)
from transit_chem.oneD import TripleWellPotential
from transit_chem.operators import Hamiltonian, overlap
from transit_chem.utils import pairwise_array_from_func

logger.info("Running transit time calculation")

V = TripleWellPotential.from_params(
    well1_depth=0.7,
    well1_halfwidth=1.0,
    bridge_length=8.0,
    bridge_depth=0.5,
    well3_halfwidth=1.0,
    well3_depth=0.5,
)

cutoff_energy = 4
basis = (
    harmonic_basis_from_parabola(V.well1, cutoff_energy=cutoff_energy)
    + harmonic_basis_from_parabola(V.well2, cutoff_energy=cutoff_energy)[:3]
    + harmonic_basis_from_parabola(V.well3, cutoff_energy=cutoff_energy)
)


logger.info("Calculating the Hamiltonian...")
H = pairwise_array_from_func(basis, Hamiltonian(V))

logger.info("Calculating the overlap matrix...")
S = pairwise_array_from_func(basis, overlap)

logger.info("Calculating the eigen basis...")
eig_basis = EigenBasis.from_basis(basis, H, S)

initial_state = HarmonicOscillator(n=1, center=V.center1[0])

logger.info("Time evolving the initial state")
time_evolving = TimeEvolvingState(eigen_basis=eig_basis, initial_state=initial_state)

overlap_well1 = partial(overlap, lower_limit=-np.inf, upper_limit=V.barrier12[0])
overlap_well2 = partial(overlap, lower_limit=V.barrier12[0], upper_limit=V.barrier23[0])
overlap_well3 = partial(overlap, lower_limit=V.barrier23[0], upper_limit=np.inf)

s1t = time_evolving.observable(overlap_well1)
s2t = time_evolving.observable(overlap_well2)
s3t = time_evolving.observable(overlap_well3)


t = np.linspace(0, 1000, 10000)

fig = plt.figure()
plt.plot(t, s1t(t), label="Well 1")
plt.plot(t, s2t(t), label="Well 2")
plt.plot(t, s3t(t), label="Well 3")
plt.plot(t, s1t(t) + s2t(t) + s3t(t), "k--", label="Sum")
plt.legend()
plt.savefig("occupancies_over_time.png")


def hopping_matrix(s1, s2, s3, delta_t: float):
    occ_probs = [s1, s2, s3]

    def f(t: float) -> np.array:
        d1 = s1(t) - s1(t - delta_t)
        d2 = s2(t) - s2(t - delta_t)
        d3 = s3(t) - s3(t - delta_t)
        increasing = [d > 0 for d in (d1, d2, d3)]

        d_occ_probs = [d1, d2, d3]

        two_states_increasing = d1 * d2 * d3 < 0

        A = np.zeros((3, 3))
        if two_states_increasing:
            # If the state is decreasing in probability, the diagonal is
            # determined
            for j in range(3):
                if not increasing[j]:
                    A[j, j] = occ_probs[j](t) / occ_probs[j](t - delta_t)
                    for k in range(3):
                        if k != j:
                            A[k, j] = d_occ_probs[k] / occ_probs[j](t)

                # If the state is increasing, diagonal is 1. Off diagonal
                # columns = 0
                elif increasing[j] > 0:
                    A[j, j] = 1
                    for k in range(3):
                        if k != j:
                            A[k, j] = 0

        else:
            for j in range(3):
                # If the state is decreasing in probability, the diagonal
                # is determined
                if not increasing[j]:
                    A[j, j] = occ_probs[j](t) / occ_probs[j](t - delta_t)
                    for k in range(3):
                        if k != j:
                            A[j, k] = 0

            # If the state is increasing, diagonal is 1. Off diagonal
            # column element is 1 - diagonal
            for j in range(3):
                if increasing[j]:
                    A[j, j] = 1
                    for k in range(3):
                        if k != j:
                            A[j, k] = 1 - A[k, k]

        return A

    return f


At = hopping_matrix(s1t, s2t, s3t, delta_t=1e-2)
A = np.stack([At(t) for t in range(300)])
print(A)
