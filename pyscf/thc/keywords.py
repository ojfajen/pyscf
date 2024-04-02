import numpy as np
import sys

class keywords:
    def __init__(self,infile):
        self.parse_input(infile)
        self.parse_geom(self.geomfile)
        print("Running PySCF/THC for input file: ",infile)
        print(f"We have {self.natom} atoms")
        print(f"Using basis: {self.basis}")
        print(f"And grid level: {self.gridlevel}")

    def parse_input(self,infile):
        self.geomfile = ''
        self.basis = ''
        self.units = ''
        self.gridlevel = 4
        with open(infile,'r') as file:
            for line in file:
                if 'coordinates' in line:
                    self.geomfile = line.split()[1]
                elif 'basis' in line:
                    self.basis = line.split()[1]
                elif 'units' in line:
                    self.units = line.split()[1]
                elif 'thcgridlevel' in line:
                    self.gridlevel = int(line.split()[1])
        return 

    def parse_geom(self,geomfile):
        self.coords = ''
        f = open(geomfile,'r')
        lines = f.readlines() 
        self.natom = int(lines[:1][0])
        f.close()
        for line in lines[2:]:
            self.coords += (line[:-2])
            self.coords += ('; ')
        return 


