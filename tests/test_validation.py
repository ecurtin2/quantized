import attr
from pytest import raises

from transit_chem.validation import Range


class DummyAttrib:
    @property
    def name(self):
        return "DummyAttrib"


def test_range_throws_if_min_gt_max():

    with raises(ValueError):
        Range(10, 1)


def test_range_on_examples():
    inst = None
    attrib = DummyAttrib()

    Range(1, 3)(inst, attrib, 2)
    Range("a", "z")(inst, attrib, "m")
    Range(-10.0, 100.0)(inst, attrib, 75.3)
    Range(1, 1)(inst, attrib, 1)
    Range(3.4, 11.3)(inst, attrib, 11.3)
    Range("b", "x")(inst, attrib, "b")

    with raises(ValueError):
        Range(1, 3)(inst, attrib, 5)
        Range("l", "m")(inst, attrib, "n")
        Range(-10.0, -3.0)(inst, attrib, 0.0)


@attr.s
class A:
    x: int = attr.ib(validator=[Range(0, 10)])
    y: float = attr.ib(validator=[Range(-3.0, 3.0)])


def test_range_attrs_instance():
    A(3, 0.0)
    A(0, 1.5)
    A(5, -3.0)

    with raises(ValueError):
        A(-1, 0.0)
        A(5, -10.0)
        A(2.3, 0.0)
        A("A", 2.3)
