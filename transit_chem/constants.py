from scipy.constants import physical_constants as c

HARTREE_TO_EV = c["Hartree energy in eV"][0]
EV_TO_HARTREE = 1.0 / c["Hartree energy in eV"][0]
BOHR_TO_METER = c["Bohr radius"][0]
BOHR_TO_ANGSTROM = 0.529177
ANGSTROM_TO_BOHR = 1.88973
ATU_TO_FEMTOSECOND = c["atomic unit of time"][0] / 1e-15
