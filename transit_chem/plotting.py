from pathlib import Path
from typing import Callable, Iterable, Union

import numpy as np
from matplotlib import pyplot as plt


def plot_occupancies_over_time(
    times: Iterable[float],
    s1t: Callable,
    s2t: Callable,
    s3t: Callable,
    save_to: Union[str, Path, None] = None,
):
    fig = plt.figure()
    t = np.array(times)

    s1 = np.array([s1t(t) for t in times])
    s2 = np.array([s2t(t) for t in times])
    s3 = np.array([s3t(t) for t in times])

    plt.plot(t, s1, label="Well 1")
    plt.plot(t, s2, label="Well 2")
    plt.plot(t, s3, label="Well 3")
    plt.plot(t, s1 + s2 + s3, "k--", label="Sum")
    plt.legend()
    if save_to is not None:
        plt.savefig(save_to)
    return fig


def plot_p_not(
    p_not: Iterable[float],
    acceptor_occprob: callable,
    times: Iterable[float],
    save_to: Union[str, Path, None] = None,
):
    fig = plt.figure()
    t = np.array(list(times))
    p = np.array(list(p_not))

    occprob = np.array([1.0 - acceptor_occprob(t) for t in times])

    plt.plot(t, p)
    plt.plot(t, occprob)
    if save_to is not None:
        plt.savefig(save_to)
    return fig
