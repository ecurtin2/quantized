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
from transit_chem.plotting import plot_occupancies_over_time, plot_p_not
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
        + harmonic_basis_from_parabola(V.well2, cutoff_energy=cutoff_energy)
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
    tau = 0.1

    time_gen = (delta_t * i for i in count())
    result = tqdm(takewhile(lambda x: x[0] >= tau, zip(p_not, time_gen)))
    p_not_list, times = zip(*result)

    logger.info(
        f"P_not ({p_not_list[-1]}) reached tau ({tau}) after {len(p_not_list)} iterations. t = {times[-1]}"
    )

    plot_occupancies_over_time(
        np.array(times), s1t, s2t, s3t, save_to=Path() / "occupancies_over_time.png"
    )

    plot_p_not(p_not_list, s3t, times, save_to=Path() / "p_not.png")

    return 0


if __name__ == "__main__":
    main()
