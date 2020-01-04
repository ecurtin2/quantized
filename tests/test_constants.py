from math import isclose

from transit_chem import constants


def test_constants():
    assert isclose(constants.ANGSTROM_TO_BOHR * constants.BOHR_TO_ANGSTROM, 1.0, abs_tol=1e-5)
    _ = constants.ATU_TO_FEMTOSECOND
    assert isclose(constants.HARTREE_TO_EV * constants.EV_TO_HARTREE, 1.0, abs_tol=1e-5)
