from pathlib import Path
from typing import Callable, Iterable, List, Union

from matplotlib.lines import Line2D
import matplotlib.animation as animation

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation

from quantized.time_evolution import TimeEvolvingState
from quantized.hopping_matrix import Pnot, OccupancyProbabilites


def plot_occupancies_over_time(
    times: Iterable[float], occ_probs: Iterable[Callable], save_to: Union[str, Path, None] = None
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
    name: str = "M",
    **kwargs,
):
    allowed = (0, 1, 2)
    if axis not in allowed:
        raise ValueError(f"Axis must be in {allowed}, got {axis}")

    mat = f(np.array(x))
    if not len(mat.shape) == 3:
        raise ValueError(f"Array must be 3d")

    M, N = [x for i, x in enumerate(mat.shape) if i != axis]
    fig, axes = plt.subplots(M, N, **kwargs)

    for i in range(M):
        for j in range(N):
            label = f"{name}[{i},{j}]"
            if axis == 0:
                axes[i, j].plot(x, mat[:, i, j])
                axes[i, j].set_title(label)
            elif axis == 1:
                axes[i, j].plot(x, mat[i, :, j])
                axes[i, j].set_title(label)
            elif axis == 2:
                axes[i, j].plot(x, mat[i, j, :])
                axes[i, j].set_title(label)

    if save_to is not None:
        plt.savefig(save_to)
    return fig, axes


def animate_state(
    fig,
    ax,
    initial_state: Callable,
    time_dependent_state: TimeEvolvingState,
    potential: Callable,
    nframes: int,
    interval: int,
) -> FuncAnimation:

    x = np.linspace(*ax.get_xlim(), 1000)
    ax.plot(x, potential(x), "k--")
    ax.plot(x, initial_state(x) * np.conj(initial_state(x)), color="C1", label="$\phi(0)$")

    line, = ax.plot([], [], color="C0", label="$\phi(t)$")
    ax.legend(loc="upper center")

    def density(x, t):
        return np.real(time_dependent_state(x, t) * np.conj(time_dependent_state(x, t)))

    def init():
        line.set_data([], [])
        return (line,)

    def animate(t):
        y = density(x, t / 10)
        line.set_data(x, y)
        return (line,)

    anim = FuncAnimation(fig, animate, init_func=init, frames=nframes, interval=interval, blit=True)

    return anim


def animate_pnot_and_matrix(
    potential: Callable,
    initial: Callable,
    initial_energy,
    time_state: TimeEvolvingState,
    pnot: Pnot,
    occ_probs: OccupancyProbabilites,
):
    class SubplotAnimation(animation.TimedAnimation):
        def __init__(self):
            fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(10, 4))
            fig.tight_layout()
            self.nframes = 200
            self.x = np.linspace(-4, 4, 1000)

            ax1.set_xlabel("x")
            ax1.set_ylabel("$\phi(t)$")
            self.line1 = Line2D([], [], color="C0", label="$\phi(t)$")
            ax1.add_line(self.line1)
            ax1.set_xlim(-4, 4)
            ax1.set_ylim(0, 4)
            ax1.plot(self.x, potential(self.x), "k--")
            ax1.plot(
                self.x,
                initial(self.x) * np.conj(initial(self.x)) + initial_energy,
                color="C1",
                label="$\phi(0)$",
            )
            ax1.legend()

            ax2.set_xlabel("time")
            ax2.set_ylabel("Occupancy Probability")
            self.donor_occprob = Line2D([], [], color="C0", label="Donor Density")
            self.acceptor_occprob = Line2D([], [], color="C1", label="Acceptor Density")

            ax2.add_line(self.donor_occprob)
            ax2.add_line(self.acceptor_occprob)
            ax2.legend()
            ax2.set_xlim(0, pnot.tau90)
            ax2.set_ylim(0, 1)

            ax3.set_xlabel("time")
            ax3.set_ylabel("Pnot(t)")
            self.pnot = Line2D([], [], color="C0", label="Pnot")
            self.pn = Line2D([], [], color="C1", label="Non-Acceptor Density")
            ax3.add_line(self.pn)
            ax3.add_line(self.pnot)
            ax3.legend()
            ax3.set_xlim(0, pnot.tau90)
            ax3.set_ylim(0, 1)

            animation.TimedAnimation.__init__(self, fig, interval=50, blit=True)

            self.t = []

        def _draw_frame(self, i):
            t = np.linspace(0, i, 1000)

            y = time_state(self.x, i)
            self.line1.set_data(self.x, np.abs(np.conj(y) * y) + initial_energy)
            self.donor_occprob.set_data(t, occ_probs(np.array(t))[0])
            self.acceptor_occprob.set_data(t, occ_probs(np.array(t))[2])
            self.pnot.set_data(t, pnot(t))
            self.pn.set_data(t, 1 - occ_probs(np.array(t))[2])
            self._drawn_artists = [
                self.line1,
                self.donor_occprob,
                self.acceptor_occprob,
                self.pnot,
                self.pn,
            ]

        def new_frame_seq(self):
            return iter(np.linspace(0, pnot.tau90, self.nframes))

        def _init_draw(self):
            lines = [self.line1, self.donor_occprob, self.acceptor_occprob, self.pnot, self.pn]
            for l in lines:
                l.set_data([], [])

    return SubplotAnimation()
