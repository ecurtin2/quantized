import numpy as np

from transit_chem import plotting as p
from pytest import raises


def test_plot_occupancies(tmp_path):
    times = list(range(5))
    fake_occ_probs = tuple([lambda x: x, lambda y: 2 * y])

    p.plot_occupancies_over_time(times, fake_occ_probs, save_to=tmp_path)


def test_plot_p_not(tmp_path):
    times = list(range(5))
    p.plot_p_not(lambda x: x, times, lambda x: x, save_to=tmp_path)


def test_plot_matrix(tmp_path):

    n = 3
    nrows = ncols = 3
    times = list(range(n))

    def matrix(x: np.array):
        return np.arange(n * nrows * ncols).reshape((n, nrows, ncols))

    for i in (0, 1, 2):
        p.plot_matrix(matrix, times, save_to=tmp_path, axis=i)

    def bad_f(x: np.array):
        return np.array([1, 2, 3])

    with raises(ValueError):
        p.plot_matrix(bad_f, times, save_to=tmp_path)

    with raises(ValueError):
        p.plot_matrix(matrix, times, save_to=tmp_path, axis=3)
