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

# --------------------------------------------
# Step 1: Load the data
# pol = np.loadtxt(infile)  # shape (N, 6); each row: [xx, yy, zz, xy, xz, yz]
#
# # Step 2: Build the 3x3 tensor for each row
# N = pol.shape[0]
# pol_tensor = np.empty((N, 3, 3))
#
# # Fill in the symmetric tensor:
# pol_tensor[:, 0, 0] = pol[:, 0]  # xx
# pol_tensor[:, 0, 1] = pol[:, 3]  # xy
# pol_tensor[:, 0, 2] = pol[:, 4]  # xz
#
# pol_tensor[:, 1, 0] = pol[:, 3]  # yx = xy
# pol_tensor[:, 1, 1] = pol[:, 1]  # yy
# pol_tensor[:, 1, 2] = pol[:, 5]  # yz
#
# pol_tensor[:, 2, 0] = pol[:, 4]  # zx = xz
# pol_tensor[:, 2, 1] = pol[:, 5]  # zy = yz
# pol_tensor[:, 2, 2] = pol[:, 2]  # zz
#
# # Step 3: Create a structured array with a field for the tensor
# dt = np.dtype([('polarizability', float, (3, 3))])
# data_array = np.empty(N, dtype=dt)
# data_array['polarizability'] = pol_tensor
#
# # Save the structured array to a file
# np.save(outfile, data_array)
# --------------------------------------------


# --------------------------------------------
# Alternatively, making it more compact
pol = np.loadtxt(infile)  # Load data

# Reshape data into (N, 3, 3) symmetric tensors
pol_33 = np.stack([pol[:, [0, 3, 4]], pol[:, [3, 1, 5]], pol[:, [4, 5, 2]]], axis=1)

# Create structured array and save
np.save(outfile, np.array(list(zip(pol_33)), dtype=[('polarizability', float, (3, 3))]))
# --------------------------------------------
