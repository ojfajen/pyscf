#!/usr/bin/env python

import numpy as np
from pyscf import lib, gto, scf, dft
import sys
from keywords import keywords
from thc import thc


def run_thc():
    keys = keywords(sys.argv[1])
    THC = thc(keys)
    print("Fitting Z from AO integrals")
    THC.do_thc()
    print("Reconstructing g from Z")
    g = THC.get_g_reconstruct(THC.X,THC.Z)
    print("Max reconstructed integral: ",np.max(g))
    print("Max original integral: ",np.max(THC.ao_eri))
    print("Max integral difference: ",np.max(THC.ao_eri-g))

run_thc()