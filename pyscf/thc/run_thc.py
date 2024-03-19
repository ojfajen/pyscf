#!/usr/bin/env python

import numpy as np
from pyscf import lib, gto, scf, dft
import sys

def do_thc():
    geom, basis, units = init(sys.argv[1])
    mol, grid, ao_eri = setup(geom,basis,units)
    mf = scf.RHF(mol).run()
    coords, weights = grid.coords, grid.weights
    X = get_X(coords,weights,mol)
    Z = get_Z(X,ao_eri)
    g = get_g_reconstruct(X,Z)
    print(np.max(g))
    print(np.max(ao_eri))
    print(np.max(ao_eri-g))

def init(infile):
    geomfile, basis, units = parse_input(infile)
    geom = parse_geom(geomfile)
    return geom, basis, units

def parse_input(infile):
    geomfile = ''
    basis = ''
    units = ''
    with open(infile,'r') as file:
        for line in file:
            if 'coords' in line:
                geomfile = line.split()[1]
            elif 'basis' in line:
                basis = line.split()[1]
            elif 'units' in line:
                units = line.split()[1]
    return geomfile, basis, units

def parse_geom(geomfile):
    geom = ''
    f = open(geomfile,'r')
    lines = f.readlines()[2:]
    f.close()
    for line in lines:
        geom += (line[:-2])
        geom += ('; ')
    return geom
            

def setup(coords,basis,units='Angstrom'):
    mol = gto.M(atom=coords,basis=basis,unit=units)
    grid = get_grid(mol)
    ao_eri = get_ao_eri(mol)
    return mol, grid, ao_eri

def get_grid(mol):
    grids = dft.gen_grid.Grids(mol)
    grids.level = 4
    grids.build()
    return grids

def get_ao_eri(mol):
    if mol.cart:
        intor = 'int2e_cart'
    else:
        intor = 'int2e_sph'    
    return mol.intor(intor,aosym='s1')

def get_X(coords,weights,mol):
    ao = mol.eval_gto('GTOval_sph', coords)
    #mo = ao.dot(mf.mo_coeff)
    # need to deal with weights here
    return ao

def get_Z(X,ao_eri):
    Y = np.linalg.inv(X.T @ X)
    Zhalf = Y @ ao_eri
    Z = Zhalf @ Y.T
    return Z

def get_g_reconstruct(X,Z):
    return X.T @ X @ Z @ X.T @ X

do_thc()