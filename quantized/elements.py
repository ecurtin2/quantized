"""This is a module docstring

It spans multiple lines
"""


from attr.validators import optional
from typing import Optional, Union

from quantized.attr_wrapped import attrs, attrib
from quantized.validation import Range, positive


@attrs(frozen=True)
class Element:
    """Class to hold element information"""

    z: int = attrib(repr=False, validator=[Range(1, 118)], desc="The atomic number of the element")
    name: str = attrib(desc="The symbol of the element, e.g. H, He, O")
    voie: Optional[float] = attrib(
        default=None,
        repr=False,
        desc="The valence orbital ionization energy in eV",
        validator=optional(positive),
    )
    alpha: Optional[float] = attrib(
        default=None, repr=False, desc="The orbital decay coefficient", validator=optional(positive)
    )


H = Element(z=1, name="H")
He = Element(z=2, name="He")
Li = Element(z=3, name="Li")
Be = Element(z=4, name="Be")
B = Element(z=5, name="B")
C = Element(z=6, name="C", voie=10.77, alpha=1.5679)
N = Element(z=7, name="N", voie=13.9, alpha=1.9170)
O = Element(z=8, name="O", voie=15.8, alpha=2.2266)  # noqa: E741
F = Element(z=9, name="F")
Ne = Element(z=10, name="Ne")
Na = Element(z=11, name="Na")
Mg = Element(z=12, name="Mg")
Al = Element(z=13, name="Al")
Si = Element(z=14, name="Si")
P = Element(z=15, name="P")
S = Element(z=16, name="S")
Cl = Element(z=17, name="Cl")
Ar = Element(z=18, name="Ar")
K = Element(z=19, name="K")
Ca = Element(z=20, name="Ca")
Sc = Element(z=21, name="Sc")
Ti = Element(z=22, name="Ti")
V = Element(z=23, name="V")
Cr = Element(z=24, name="Cr")
Mn = Element(z=25, name="Mn")
Fe = Element(z=26, name="Fe")
Co = Element(z=27, name="Co")
Ni = Element(z=28, name="Ni")
Cu = Element(z=29, name="Cu")
Zn = Element(z=30, name="Zn")
Ga = Element(z=31, name="Ga")
Ge = Element(z=32, name="Ge")
As = Element(z=33, name="As")
Se = Element(z=34, name="Se")
Br = Element(z=35, name="Br")
Kr = Element(z=36, name="Kr")
Rb = Element(z=37, name="Rb")
Sr = Element(z=38, name="Sr")
Y = Element(z=39, name="Y")
Zr = Element(z=40, name="Zr")
Nb = Element(z=41, name="Nb")
Mo = Element(z=42, name="Mo")
Tc = Element(z=43, name="Tc")
Ru = Element(z=44, name="Ru")
Rh = Element(z=45, name="Rh")
Pd = Element(z=46, name="Pd")
Ag = Element(z=47, name="Ag")
Cd = Element(z=48, name="Cd")
In = Element(z=49, name="In")
Sn = Element(z=50, name="Sn")
Sb = Element(z=51, name="Sb")
Te = Element(z=52, name="Te")
I = Element(z=53, name="I")  # noqa: E741
Xe = Element(z=54, name="Xe")
Cs = Element(z=55, name="Cs")
Ba = Element(z=56, name="Ba")
La = Element(z=57, name="La")
Ce = Element(z=58, name="Ce")
Pr = Element(z=59, name="Pr")
Nd = Element(z=60, name="Nd")
Pm = Element(z=61, name="Pm")
Sm = Element(z=62, name="Sm")
Eu = Element(z=63, name="Eu")
Gd = Element(z=64, name="Gd")
Tb = Element(z=65, name="Tb")
Dy = Element(z=66, name="Dy")
Ho = Element(z=67, name="Ho")
Er = Element(z=68, name="Er")
Tm = Element(z=69, name="Tm")
Yb = Element(z=70, name="Yb")
Lu = Element(z=71, name="Lu")
Hf = Element(z=72, name="Hf")
Ta = Element(z=73, name="Ta")
W = Element(z=74, name="W")
Re = Element(z=75, name="Re")
Os = Element(z=76, name="Os")
Ir = Element(z=77, name="Ir")
Pt = Element(z=78, name="Pt")
Au = Element(z=79, name="Au")
Hg = Element(z=80, name="Hg")
Tl = Element(z=81, name="Tl")
Pb = Element(z=82, name="Pb")
Bi = Element(z=83, name="Bi")
Po = Element(z=84, name="Po")
At = Element(z=85, name="At")
Rn = Element(z=86, name="Rn")
Fr = Element(z=87, name="Fr")
Ra = Element(z=88, name="Ra")
Ac = Element(z=89, name="Ac")
Th = Element(z=90, name="Th")
Pa = Element(z=91, name="Pa")
U = Element(z=92, name="U")
Np = Element(z=93, name="Np")
Pu = Element(z=94, name="Pu")
Am = Element(z=95, name="Am")
Cm = Element(z=96, name="Cm")
Bk = Element(z=97, name="Bk")
Cf = Element(z=98, name="Cf")
Es = Element(z=99, name="Es")
Fm = Element(z=100, name="Fm")
Md = Element(z=101, name="Md")
No = Element(z=102, name="No")
Lr = Element(z=103, name="Lr")
Rf = Element(z=104, name="Rf")
Db = Element(z=105, name="Db")
Sg = Element(z=106, name="Sg")
Bh = Element(z=107, name="Bh")
Hs = Element(z=108, name="Hs")
Mt = Element(z=109, name="Mt")
Ds = Element(z=110, name="Ds")
Rg = Element(z=111, name="Rg")
Cn = Element(z=112, name="Cn")
Uut = Element(z=113, name="Uut")
Fl = Element(z=114, name="Fl")
Uup = Element(z=115, name="Uup")
Lv = Element(z=116, name="Lv")
Uus = Element(z=117, name="Uus")
Uuo = Element(z=118, name="Uuo")


all_elements = [e for e in globals().values() if isinstance(e, Element)]


def element_from_string(s: Union[str, int, Element]) -> Element:
    if isinstance(s, Element):
        return s

    try:
        z = int(s)
        return next(e for e in all_elements if e.z == z)
    except ValueError:
        return next(e for e in all_elements if e.name.lower() == s.lower())


elements_dict = {a: [getattr(e, a) for e in all_elements] for a in ["name", "z", "voie", "alpha"]}