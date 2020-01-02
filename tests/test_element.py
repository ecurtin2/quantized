from transit_chem.elements import element_from_string, H, Uup


def test_element_from_string():

    assert element_from_string("Uup") is Uup
    assert element_from_string(1) is H
    assert element_from_string("1") is H
    assert element_from_string(H) is H
