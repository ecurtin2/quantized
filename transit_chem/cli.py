"""Console script for transit_chem."""
from functools import partial
from itertools import takewhile
from pathlib import Path

import click
import numpy as np
from loguru import logger

from transit_chem.basis import (
    EigenBasis,
    HarmonicOscillator,
    TimeEvolvingState,
    harmonic_basis_from_parabola,
)
from transit_chem.hopping_matrix import hopping_matrix, p_not_generator
from transit_chem.operators import Hamiltonian, overlap
from transit_chem.plotting import plot_occupancies_over_time
from transit_chem.potentials import TripleWell
from transit_chem.utils import pairwise_array_from_func


@click.command()
def main():
    logger.info("Running transit time calculation")

    V = TripleWell.from_params(
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

    logger.info("Calculating the Hamiltonian")
    H = pairwise_array_from_func(basis, Hamiltonian(V))

    logger.info("Calculating the overlap matrix")
    S = pairwise_array_from_func(basis, overlap)

    logger.info("Calculating the eigen basis")
    eig_basis = EigenBasis.from_basis(basis, H, S)

    logger.info("Time evolving the initial state")
    time_evolving = TimeEvolvingState(
        eigen_basis=eig_basis, initial_state=HarmonicOscillator(n=0, center=V.center1[0])
    )

    overlap_well1 = partial(overlap, lower_limit=-np.inf, upper_limit=V.barrier12[0])
    overlap_well2 = partial(overlap, lower_limit=V.barrier12[0], upper_limit=V.barrier23[0])
    overlap_well3 = partial(overlap, lower_limit=V.barrier23[0], upper_limit=np.inf)

    s1t = time_evolving.observable(overlap_well1)
    s2t = time_evolving.observable(overlap_well2)
    s3t = time_evolving.observable(overlap_well3)

    At = hopping_matrix(s1t, s2t, s3t)

    delta_t = 1e-2
    p0 = np.array([s1t(0), s2t(0), s3t(0)])
    p_not = p_not_generator(acceptor=2, hopping_matrix=At, p0=p0, delta_t=delta_t)

    from itertools import count
    from tqdm import tqdm

    logger.info("Starting p_not generator computation")
    val = None
    n_iters = None
    tau = 0.9
    for val, n_iters in tqdm(takewhile(lambda x: x[0] >= tau, zip(p_not, count()))):
        pass
    val = next(p_not)
    n_iters += 1
    logger.info(
        f"P_not ({val}) reached tau ({tau}) after {n_iters} iterations. t = {n_iters * delta_t}"
    )

    times = np.linspace(0, 10, 10)
    plot_occupancies_over_time(times, s1t, s2t, s3t, save_to=Path() / "occupancies_over_time.png")

    return 0


if __name__ == "__main__":
    main()
