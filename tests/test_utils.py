from math import isclose

from hypothesis import given

from transit_chem.config import SMALL_NUMBER
from transit_chem.utils import Parabola
from utils import reasonable_floats


@given(reasonable_floats, reasonable_floats, reasonable_floats)
def test_parabola_properties(a: float, b: float, c: float):

    p = Parabola(a, b, c)
    assert isclose(p(0), p.c, abs_tol=SMALL_NUMBER)

    if p.has_vertex:
        x1 = p.vertex + 3.4
        x2 = p.vertex - 3.4
        assert isclose(p(x1), p(x2), abs_tol=SMALL_NUMBER)
