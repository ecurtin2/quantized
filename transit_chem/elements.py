from typing import Optional

import attr


@attr.s
class Element:
    z: int = attr.ib()
    voie: Optional[float] = attr.ib(default=None)
    alpha: Optional[float] = attr.ib(default=None)


H = Element(z=1)
He = Element(z=2)
Li = Element(z=3)
Be = Element(z=4)
B = Element(z=5)
C = Element(z=6, voie=10.77, alpha=1.5679)
N = Element(z=7, voie=13.9, alpha=1.9170)
O = Element(z=8, voie=15.8, alpha=2.2266)  # noqa: E741
F = Element(z=9)
Ne = Element(z=10)
Na = Element(z=11)
Mg = Element(z=12)
Al = Element(z=13)
Si = Element(z=14)
P = Element(z=15)
S = Element(z=16)
Cl = Element(z=17)
Ar = Element(z=18)
K = Element(z=19)
Ca = Element(z=20)
Sc = Element(z=21)
Ti = Element(z=22)
V = Element(z=23)
Cr = Element(z=24)
Mn = Element(z=25)
Fe = Element(z=26)
Co = Element(z=27)
Ni = Element(z=28)
Cu = Element(z=29)
Zn = Element(z=30)
Ga = Element(z=31)
Ge = Element(z=32)
As = Element(z=33)
Se = Element(z=34)
Br = Element(z=35)
Kr = Element(z=36)
Rb = Element(z=37)
Sr = Element(z=38)
Y = Element(z=39)
Zr = Element(z=40)
Nb = Element(z=41)
Mo = Element(z=42)
Tc = Element(z=43)
Ru = Element(z=44)
Rh = Element(z=45)
Pd = Element(z=46)
Ag = Element(z=47)
Cd = Element(z=48)
In = Element(z=49)
Sn = Element(z=50)
Sb = Element(z=51)
Te = Element(z=52)
I = Element(z=53)  # noqa: E741
Xe = Element(z=54)
Cs = Element(z=55)
Ba = Element(z=56)
La = Element(z=57)
Ce = Element(z=58)
Pr = Element(z=59)
Nd = Element(z=60)
Pm = Element(z=61)
Sm = Element(z=62)
Eu = Element(z=63)
Gd = Element(z=64)
Tb = Element(z=65)
Dy = Element(z=66)
Ho = Element(z=67)
Er = Element(z=68)
Tm = Element(z=69)
Yb = Element(z=70)
Lu = Element(z=71)
Hf = Element(z=72)
Ta = Element(z=73)
W = Element(z=74)
Re = Element(z=75)
Os = Element(z=76)
Ir = Element(z=77)
Pt = Element(z=78)
Au = Element(z=79)
Hg = Element(z=80)
Tl = Element(z=81)
Pb = Element(z=82)
Bi = Element(z=83)
Po = Element(z=84)
At = Element(z=85)
Rn = Element(z=86)
Fr = Element(z=87)
Ra = Element(z=88)
Ac = Element(z=89)
Th = Element(z=90)
Pa = Element(z=91)
U = Element(z=92)
Np = Element(z=93)
Pu = Element(z=94)
Am = Element(z=95)
Cm = Element(z=96)
Bk = Element(z=97)
Cf = Element(z=98)
Es = Element(z=99)
Fm = Element(z=100)
Md = Element(z=101)
No = Element(z=102)
Lr = Element(z=103)
Rf = Element(z=104)
Db = Element(z=105)
Sg = Element(z=106)
Bh = Element(z=107)
Hs = Element(z=108)
Mt = Element(z=109)
Ds = Element(z=110)
Rg = Element(z=111)
Cn = Element(z=112)
Uut = Element(z=113)
Fl = Element(z=114)
Uup = Element(z=115)
Lv = Element(z=116)
Uus = Element(z=117)
Uuo = Element(z=118)
