from pathlib import Path
from typing import List, Iterable, Union, Callable

import numpy as np
from matplotlib import pyplot as plt

from transit_chem.hopping_matrix import OccupancyProbabilites, HoppingMatrix


def plot_occupancies_over_time(
    times: Iterable[float],
    occ_probs: OccupancyProbabilites,
    save_to: Union[str, Path, None] = None,
):
    fig = plt.figure()
    t = np.array(times)

    s_list = [np.array([s(t) for t in times]) for s in occ_probs]

    for i, s in enumerate(s_list):
        plt.plot(t, s, label=f"Well {i+1}")
    plt.plot(t, sum(s_list), "k--", label="Sum")
    plt.legend()
    if save_to is not None:
        plt.savefig(save_to)
    return fig


def plot_p_not(
    p_not: List[float],
    times: Iterable[float],
    acceptor_occ_prob: Callable,
    save_to: Union[str, Path, None] = None,
):
    fig = plt.figure()

    occprob = np.array([1.0 - acceptor_occ_prob(t) for t in times])

    plt.plot(times, p_not, label="pnot")
    plt.plot(times, occprob, label="occprob")
    plt.legend()
    if save_to is not None:
        plt.savefig(save_to)
    return fig


def plot_hopping_matrix(
    m: HoppingMatrix, times: List[float], delta_t: float, save_to: Union[str, Path, None] = None,
):
    fig, axes = plt.subplots(m.N, m.N)

    a_mat = m(np.array(times), delta_t=delta_t)

    for i in range(m.N):
        for j in range(m.N):
            axes[i, j].plot(times, a_mat[:, i, j], label=f"{i}{j}")

    plt.legend()
    if save_to is not None:
        plt.savefig("hopping_matrix.png")
    return fig, axes
