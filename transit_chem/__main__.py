import numpy as np

import transittime


def main(argv=None):

    x = np.linspace(-10, 10, 500)
    one = transittime.OneD(coordinates=x, nregions=3)
    t = np.linspace(0, 1500, 1500)
    initial_state = one.basis_funcs[0](one.coords, 0, 0.0)
    one.calc_all(initialstate=initial_state)

    #fname = '../analysis/xyz_files/2C.xyz'
    #mol = transittime.Molecule(fname)
    #mol.trim_hydrogens()
    #mol.scale(transittime.units.conversion_factors['ANGSTROM_TO_BOHR'])

    ##  Aligns the line joining the two nitrogens with the x axis.
    #nitrogens = [atom for atom in mol if atom.symbol == 'N']
    #if nitrogens:
    #    mol.align_to_x(*nitrogens)
    #else:
    #    raise AttributeError('Did not find nitrogens!')

    ##  Align the nitro (acceptor) group to be +x
    #oxygens = [atom for atom in mol if atom.symbol == 'O']
    #if oxygens:
    #    if all(oxygen.coords[0] < 0 for oxygen in oxygens):
    #        mol.flip_x()

    ##  Sort into ascending order of x coordinates for convenience
    #mol.sort(atomic_key=lambda atom: atom.x)

    ##  This works since the molecule is aligned/shifted as expected.
    #donor_n = min(nitrogens, key=lambda atom: atom.x)
    #acceptor_n = max(nitrogens, key=lambda atom: atom.x)
    #donor_index = mol.atoms.index(donor_n)
    #acceptor_index = mol.atoms.index(acceptor_n)

    ##  Since they are aligned on the x axis already and atoms are sorted by x value, the adjacent atoms in
    ##  the list will be the relevant carbons for the boundary

    #donor_c = mol.atoms[donor_index + 1]
    #acceptor_c = mol.atoms[acceptor_index - 1]

    #donor_bridge_bound = 0.5 * (donor_n.x + donor_c.x)
    #acceptor_bridge_bound = 0.5 * (acceptor_n.x + acceptor_c.x)
    #bounds = (donor_bridge_bound, acceptor_bridge_bound)

    #ThreeDCalc = transittime.threeD.ThreeD(nregions=3, boundaries=bounds, molecule=mol, id='testing')
    #ThreeDCalc.calc_time_independent()

    ##  Initial state is Pz orbital on donor nitrogen. This was most convenient way to do it.
    #initialstate = transittime.basis.PzBasis(*donor_n.coords, donor_n.alpha)
    #ThreeDCalc.set_initial_state(initialstate)  # does overlap integrals, can be slow
    #ThreeDCalc.calc_time_dependent(maxsteps=10**6)


if __name__ == '__main__':
    import sys
    main(sys.argv)
