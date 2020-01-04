from pathlib import Path
from typing import Callable, Iterable, List, Union

import numpy as np
from matplotlib import pyplot as plt


def plot_occupancies_over_time(
    times: Iterable[float], occ_probs: Iterable[Callable], save_to: Union[str, Path, None] = None,
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
    p_not: Callable,
    times: List[float],
    acceptor_occ_prob: Callable,
    save_to: Union[str, Path, None] = None,
):
    fig = plt.figure()

    occprob = np.array([1.0 - acceptor_occ_prob(t) for t in times])
    pnot_over_time = [p_not(t) for t in times]

    plt.plot(times, pnot_over_time, label="pnot")
    plt.plot(times, occprob, label="occprob")
    plt.legend()
    if save_to is not None:
        plt.savefig(save_to)
    return fig


def plot_matrix(
    f: Callable[[np.array], np.array],
    x: List[float],
    save_to: Union[str, Path, None] = None,
    axis: int = 0,
):
    allowed = (0, 1, 2)
    if axis not in allowed:
        raise ValueError(f"Axis must be in {allowed}, got {axis}")

    mat = f(np.array(x))
    if not len(mat.shape) == 3:
        raise ValueError(f"Array must be 3d")

    M, N = [x for i, x in enumerate(mat.shape) if i != axis]
    fig, axes = plt.subplots(M, N)

    for i in range(M):
        for j in range(N):
            if axis == 0:
                axes[i, j].plot(x, mat[:, i, j], label=f"{i}{j}")
            elif axis == 1:
                axes[i, j].plot(x, mat[i, :, j], label=f"{i}{j}")
            elif axis == 2:
                axes[i, j].plot(x, mat[i, j, :], label=f"{i}{j}")

    plt.legend()
    if save_to is not None:
        plt.savefig(save_to)
    return fig, axes
