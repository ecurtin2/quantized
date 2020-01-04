from math import isclose

import numpy as np
from hypothesis import given

from transit_chem.config import conf
from transit_chem.utils import Parabola, pairwise_array_from_func
from utils import reasonable_floats


@given(reasonable_floats, reasonable_floats, reasonable_floats)
def test_parabola_properties(a: float, b: float, c: float):

    p = Parabola(a, b, c)
    assert isclose(p(0), p.c, abs_tol=conf.small_number)

    if p.has_vertex:
        x1 = p.vertex + 3.4
        x2 = p.vertex - 3.4
        assert isclose(p(x1), p(x2), abs_tol=conf.small_number)


def f(x, y):
    return x + y


def test_pairwise_array_from_func():

    mylist = [1, 2, 3]

    expected = np.array([[2, 3, 4], [3, 4, 5], [4, 5, 6]])

    result = pairwise_array_from_func(mylist, f)

    assert np.allclose(result, expected)
