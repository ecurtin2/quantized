"""Console script for transit_chem."""
import sys
from functools import partial
from pathlib import Path

import click
import numpy as np
from loguru import logger

from transit_chem import config
from transit_chem.basis import EigenBasis, HarmonicOscillator, harmonic_basis_from_parabola
from transit_chem.hopping_matrix import HoppingMatrix, OccupancyProbabilites, Pnot
from transit_chem.operators import Hamiltonian, Overlap
from transit_chem.plotting import plot_matrix, plot_occupancies_over_time, plot_p_not
from transit_chem.potentials import TripleWell
from transit_chem.time_evolution import TimeEvolvingState


@click.group()
@click.option(
    "-l",
    "--log-level",
    type=click.Choice(["TRACE", "DEBUG", "INFO", "SUCCESS", "WARNING", "ERROR", "CRITICAL"]),
    default="INFO",
)
def cli(log_level: str):
    logger.remove()
    logger.add(sys.stdout, level=log_level)


@cli.command("run-oned")
@click.option("-o", "--output-path", type=str, default=".", show_default=True)
@click.option(
    "--cutoff-energy",
    type=float,
    help="The cutoff energy for the basis set. Basis functions will all be below this energy",
    default=4.0,
    show_default=True,
)
@click.option(
    "--until-pnot",
    type=float,
    help="The cutoff value of Pnot to calculate to. For tau90, use --until-pnot=0.1",
    default=0.9,
    show_default=True,
)
@click.option(
    "--pnot_conv_tol",
    type=float,
    default=1.0,
    help="Tolerance for converging until-pnot value with respect to time step.",
    show_default=True,
)
@click.option(
    "--acceptor",
    type=click.Choice(["1", "2", "3"]),
    default="3",
    help="Which well is the acceptor. ",
    show_default=True,
)
@click.option("--well1_depth", type=float, default=0.7, show_default=True)
@click.option("--well1_halfwidth", type=float, default=1.0, show_default=True)
@click.option("--bridge_length", type=float, default=8.0, show_default=True)
@click.option("--bridge_depth", type=float, default=0.5, show_default=True)
@click.option("--well3_halfwidth", type=float, default=1.0, show_default=True)
@click.option("--well3_depth", type=float, default=0.5, show_default=True)
def oned(
    output_path: str = ".",
    log_level: str = "INFO",
    cutoff_energy: float = 4.0,
    until_pnot: float = 0.9,
    pnot_conv_tol: float = 1.0,
    acceptor: str = 3,
    well1_depth=0.7,
    well1_halfwidth=1.0,
    bridge_length=8.0,
    bridge_depth=0.5,
    well3_halfwidth=1.0,
    well3_depth=0.5,
):
    logger.info("Running transit time calculation")

    triple_well = TripleWell.from_params(
        well1_depth=well1_depth,
        well1_halfwidth=well1_halfwidth,
        bridge_length=bridge_length,
        bridge_depth=bridge_depth,
        well3_halfwidth=well3_halfwidth,
        well3_depth=well3_depth,
    )

    well1_basis = harmonic_basis_from_parabola(triple_well.well1, cutoff_energy=cutoff_energy)
    well2_basis = harmonic_basis_from_parabola(triple_well.well2, cutoff_energy=cutoff_energy)
    well3_basis = harmonic_basis_from_parabola(triple_well.well3, cutoff_energy=cutoff_energy)
    basis = well1_basis + well2_basis + well3_basis
    logger.info(
        f"Basis set size is well1: {len(well1_basis)} well2: {len(well2_basis)}"
        f" well3: {len(well3_basis)} a total size of {len(basis)}"
    )

    logger.info("Calculating the Hamiltonian")
    H = Hamiltonian(triple_well).matrix(basis)

    logger.info("Calculating the overlap matrix")
    S = Overlap().matrix(basis)

    logger.info("Calculating the eigen basis")
    eig_basis = EigenBasis.from_basis(basis, H, S)

    logger.info("Time evolving the initial state")
    time_evolving = TimeEvolvingState(
        eigen_basis=eig_basis, initial_state=HarmonicOscillator(n=0, center=triple_well.center1[0])
    )
    occ_probs = OccupancyProbabilites.from_1d_state(
        time_evolving, borders=(-np.inf, triple_well.barrier12[0], triple_well.barrier23[0], np.inf)
    )

    hopping_matrix = HoppingMatrix(occ_probs)

    logger.info("Starting p_not generator computation")
    p_not = Pnot.converged_with_timestep(
        hopping_matrix, until_equals=until_pnot, tolerance=pnot_conv_tol, acceptor=int(acceptor) - 1
    )

    output_path = Path(output_path)
    plot_occupancies_over_time(
        np.array(p_not.times), occ_probs, save_to=output_path / "occupancies_over_time.png"
    )
    plot_p_not(p_not, p_not.times, p_not.acceptor_occ_prob, save_to=output_path / "p_not.png")

    plot_matrix(
        partial(hopping_matrix, delta_t=p_not.delta_t),
        p_not.times,
        save_to=output_path / "hopping_matrix.png",
    )


@cli.group("config")
def config_():
    pass


config_.command("show")(config.show)
config_.command("restore-defaults")(config.restore_defaults)


@config_.command("set")
@click.argument("args", nargs=-1)
def set_items(args):
    kwargs = dict(s.split("=") for s in args)
    config.set_items(**kwargs)


@cli.command("clear-cache")
def clear_cache():
    click.confirm(
        f"This will remove the contents of {config.conf.cache_dir}, recursively. Continue?",
        abort=True,
    )
    config.clear_cache()


if __name__ == "__main__":
    cli()
