
from ase.optimize import BFGS, FIRE
from gpaw import GPAW, PW
import ase.io, sys

curr_structure = ase.io.read('/home/alg/installs/pyrovskite/tests/data/testvcopt_fcopt.cif')

calc = GPAW(xc='PBE',
            mode=PW(450, dedecut='estimate'),
            ("kpts={'density': %d, 'gamma': True},", 6.0)
            parallel={'sl_auto': True},
            txt= '/home/alg/installs/pyrovskite/tests/data/testvcopt' + '_fixopt.txt')

curr_structure.calc = calc
fix_relax = FIRE(curr_structure)
fix_relax.run(fmax=0.02)

curr_structure.write('/home/alg/installs/pyrovskite/tests/data/testvcopt' + '_OPT.cif')
                      