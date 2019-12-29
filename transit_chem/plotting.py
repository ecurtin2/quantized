from pathlib import Path
from typing import Callable, Iterable, Union

import numpy as np
from matplotlib import pyplot as plt


def plot_occupancies_over_time(
    times: Iterable[float],
    s1t: Callable,
    s2t: Callable,
    s3t: Callable,
    save_to: Union[str, Path, None],
):
    fig = plt.figure()
    t = np.array(times)
    plt.plot(t, s1t(t), label="Well 1")
    plt.plot(t, s2t(t), label="Well 2")
    plt.plot(t, s3t(t), label="Well 3")
    plt.plot(t, s1t(t) + s2t(t) + s3t(t), "k--", label="Sum")
    plt.legend()
    if save_to is not None:
        plt.savefig(save_to)
    return fig
