#!/usr/bin/env python

import numpy as np
from pyscf import lib, gto, scf, dft
import sys
from keywords import keywords, parse_geom, parse_input
from thc import thc


def run_thc():
    keys = keywords(sys.argv[1])
    THC = thc(keys)
    THC.do_thc()
    g = THC.get_g_reconstruct()
    print(np.max(g))
    print(np.max(THC.ao_eri))
    print(np.max(THC.ao_eri-g))

run_thc()