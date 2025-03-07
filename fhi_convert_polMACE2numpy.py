#!/usr/bin/env python

import numpy as np
import argparse

"""
Reads a MD file (will try to infer the format) and writes an .extxyz file. Different outputfiles can be given if supported by ASE
"""
from ase.io import read, write

parser = argparse.ArgumentParser(
    prog="fhitk_convert_polarizabilities2numpy.py",
    description="Converts polarizabilites from NEP (GPUMD)  to numpy array",
)

# inputfile = fhi_aims_outputfile
# outputfile = extxyz_outputfile

parser.add_argument(
    "-i",
    "--inputfile",
    default="polarizability_train.out",
    help="input file, an AIMS ouptut file to be parsed",
)
parser.add_argument(
    "-o",
    "--outputfile",
    default="polarizabilities",
    help="outputfile file where the geoemtries are writen in an .extxyz file",
)

args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile

pol = np.loadtxt(infile)  # Load data

# Reshape data into (N, 3, 3) symmetric tensors
pol_33 = np.stack([pol[:, [0, 1, 2]], pol[:, [3, 4, 5]], pol[:, [6, 7, 8]]], axis=1)

# Create structured array and save
np.save(outfile, np.array(list(zip(pol_33)), dtype=[("polarizability", float, (3, 3))]))
# --------------------------------------------
