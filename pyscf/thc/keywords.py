import numpy as np
import sys

class keywords:
    def __init__(self,infile):
        self.geomfile, self.basis, self.units, self.gridlevel = parse_input(infile)
        self.coords = parse_geom(self.geomfile)

def parse_input(infile):
    geomfile = ''
    basis = ''
    units = ''
    gridlevel = 4
    with open(infile,'r') as file:
        for line in file:
            if 'coordinates' in line:
                geomfile = line.split()[1]
            elif 'basis' in line:
                basis = line.split()[1]
            elif 'units' in line:
                units = line.split()[1]
            elif 'thcgridlevel' in line:
                gridlevel = int(line.split()[1])
    return geomfile, basis, units, gridlevel

def parse_geom(geomfile):
    geom = ''
    f = open(geomfile,'r')
    lines = f.readlines()[2:]
    f.close()
    for line in lines:
        geom += (line[:-2])
        geom += ('; ')
    return geom


