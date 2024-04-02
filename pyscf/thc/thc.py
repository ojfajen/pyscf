#!/usr/bin/env python

import numpy as np
from pyscf import lib, gto, scf, dft
import sys
from keywords import keywords

class thc:
    def __init__(self,keywords):
        self.keywords = keywords
        self.setup()
        print("Starting RHF calculation")
        self.mf = self.mol.RHF().run()
        self.c = self.mf.mo_coeff
        print("Fetching MO integrals")
        self.mo_eri = self.mol.ao2mo(self.c)


    def setup(self):
        self.mol = gto.M(atom=self.keywords.coords,basis=self.keywords.basis,unit=self.keywords.units)
        self.grid = self.get_grid()
        self.ao_eri = self.get_ao_eri()
        return

    def get_grid(self):
        grids = dft.gen_grid.Grids(self.mol)
        grids.level = self.keywords.gridlevel
        grids.build()
        return grids

    def get_ao_eri(self):
        print("Fetching AO integrals")
        if self.mol.cart:
            intor = 'int2e_cart'
        else:
            intor = 'int2e_sph'    
        return self.mol.intor(intor,aosym='s1')

    def do_thc(self):
        self.get_X()
        self.get_Z(self.X,self.ao_eri)
        return

    def get_X(self):
        self.X = self.mol.eval_gto('GTOval_sph', self.grid.coords)
        self.X_mo = self.X.dot(self.mf.mo_coeff)
        # need to deal with weights here
        return 

    def get_Z(self,X,eri):
        Y = np.linalg.inv(X.T @ X)
        Zhalf = Y @ eri
        self.Z = Zhalf @ Y.T
        return 

    def get_g_reconstruct(self,X,Z):
        return X.T @ X @ Z @ X.T @ X

