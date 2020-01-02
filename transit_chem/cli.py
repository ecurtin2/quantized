"""Console script for transit_chem."""
from pathlib import Path

import click
import numpy as np
from loguru import logger

from transit_chem.basis import (
    EigenBasis,
    HarmonicOscillator,
    harmonic_basis_from_parabola,
)
from transit_chem.time_evolution import TimeEvolvingState
from transit_chem.hopping_matrix import HoppingMatrix, Pnot, OccupancyProbabilites
from transit_chem.operators import Hamiltonian, Overlap
from transit_chem.plotting import plot_occupancies_over_time, plot_p_not, plot_hopping_matrix
from transit_chem.potentials import TripleWell


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
    well1_basis = harmonic_basis_from_parabola(V.well1, cutoff_energy=cutoff_energy)
    well2_basis = harmonic_basis_from_parabola(V.well2, cutoff_energy=cutoff_energy)
    well3_basis = harmonic_basis_from_parabola(V.well3, cutoff_energy=cutoff_energy)
    basis = well1_basis + well2_basis + well3_basis
    logger.info(
        f"Basis set size is well1: {len(well1_basis)} well2: {len(well2_basis)}"
        f" well3: {len(well3_basis)} a total size of {len(basis)}"
    )

    logger.info("Calculating the Hamiltonian")
    H = Hamiltonian(V).matrix(basis)

    logger.info("Calculating the overlap matrix")
    S = Overlap().matrix(basis)

    logger.info("Calculating the eigen basis")
    eig_basis = EigenBasis.from_basis(basis, H, S)

    logger.info("Time evolving the initial state")
    time_evolving = TimeEvolvingState(
        eigen_basis=eig_basis, initial_state=HarmonicOscillator(n=0, center=V.center1[0])
    )
    occ_probs = OccupancyProbabilites.from_1d_state(
        time_evolving, borders=(-np.inf, V.barrier12[0], V.barrier23[0], np.inf)
    )
    At = HoppingMatrix(occ_probs)

    logger.info("Starting p_not generator computation")
    tau = 0.1
    times, p_not_list, p_not = Pnot.converged_with_timestep(At, tau=tau, tolerance=2.0)
    logger.info(
        f"P_not ({p_not_list[-1]}) reached tau ({tau}) after {len(p_not_list)} iterations. t = {times[-1]}"
    )

    plot_occupancies_over_time(
        np.array(times), occ_probs, save_to=Path() / "occupancies_over_time.png"
    )
    plot_p_not(p_not_list, times, p_not.acceptor_occ_prob, save_to=Path() / "p_not.png")
    plot_hopping_matrix(
        At, times=times, delta_t=p_not.delta_t, save_to=Path() / "hopping_matrix.png"
    )


if __name__ == "__main__":
    main()
