#!/usr/bin/env python
#
# Author: Oliver J. Backhouse <olbackhouse@gmail.com>
#

'''
An example of ways to restart an AGF2 calculation.
'''

import numpy
from pyscf import gto, scf, agf2, lib

mol = gto.M(atom='O 0 0 0; H 0 0 1; H 0 1 0', basis='cc-pvdz', verbose=5)

mf = scf.RHF(mol)
mf.conv_tol = 1e-12
mf.chkfile = 'agf2.chk'
mf.run()

# Run an AGF2 calculation:
gf2 = agf2.AGF2(mf)
gf2.conv_tol = 1e-7
gf2.max_cycle = 3
gf2.run()

# Restore the Mole and SCF first
mol = lib.chkfile.load_mol('agf2.chk')
mf = scf.RHF(mol)
mf.__dict__.update(lib.chkfile.load('agf2.chk', 'scf'))

# Restore the AGF2 calculation
dic = agf2.chkfile.load_agf2('agf2.chk')
gf2a = agf2.AGF2(mf)
gf2a.__dict__.update(agf2.chkfile.load_agf2('agf2.chk')[1])
gf2a.max_cycle = 50
gf2a.run()
