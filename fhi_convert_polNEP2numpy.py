#!/usr/bin/env python

import numpy as np
import argparse

"""
Reads a MD file (will try to infer the format) and writes an .extxyz file. Different outputfiles can be given if supported by ASE
"""
from ase.io import read, write

parser = argparse.ArgumentParser(
    prog='fhitk_convert_polarizabilities2numpy.py',
    description='Converts polarizabilites from NEP (GPUMD)  to numpy array',
)

# inputfile = fhi_aims_outputfile
# outputfile = extxyz_outputfile

parser.add_argument('-i', '--inputfile', default='polarizability_train.out',
                    help='input file, an AIMS ouptut file to be parsed')
parser.add_argument('-o', '--outputfile', default='polarizabilities',
                    help='outputfile file where the geoemtries are writen in an .extxyz file')

args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile

pol = np.loadtxt(infile)

# For loop method:
# pol_33 = []
#
# for i in range(len(pol)):
#     pol_33.append([pol[i][0], pol[i][3], pol[i][4], pol[i][3], pol[i][1], pol[i][5], pol[i][4], pol[i][5], pol[i][2]])
#
# dt = np.dtype([
#     ('polarizability', (float, (3, 3))),
# ])
# data_array = np.array(pol_33, dtype=dt)

# List comprehension method
data_array = np.array([(p[0], p[3], p[4], p[3], p[1], p[5], p[4], p[5], p[2])
                       for p in pol]).astype(np.dtype([
                           ('polarizability', (float, (3, 3))),
                       ]))

np.save(outfile, data_array)
