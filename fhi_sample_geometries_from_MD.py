#!/usr/bin/env python3

import os

import argparse
from ase.io import read
from ase.io.aims import write_aims
import numpy as np

'''
sample num_samples unique geometries from a MD, and creates the corresponding geometry.in file of each geometry in a separate folder
'''

# Read the arguments
parser = argparse.ArgumentParser(
    prog='sample_geometries.py',
    description='sample unique geometries from a MD, and creates the corresponding geometry.in file of each geometry in a separate folder: prefix + number',
)

inputfile = 'md.extxyz'
outputfile = 'geometry.in'

parser.add_argument('-i', '--inputfile', default=inputfile,
                    help='input file: MD file to take the samples from [md.extxyz]')
parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='output file: geometry.in file in AIMS format There is no reason to change the default value [geometry.in]')
parser.add_argument('-p', '--prefix', default='md_sample_',
                    help='Prefix for the created folders. [md_sample_]')
parser.add_argument('-n', '--samples', default=1,
                    help='Total number of geometries to sample. [1]')
parser.add_argument('-c', '--write_control', action='store_true',
                    help='Write template control file in the current folder')

args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile
num_samples = int(args.samples)
prefix_folder = args.prefix
write_control = args.write_control

# Read the MD
print('Reading Molecular dynamics file...', end='')
md = read(infile, index=':')
print('Done!')

# Sampling from the Md
print('Selecting random geometries...', end='')
random_indices = np.random.choice(len(md), num_samples, replace=False)
print('Done!')

# Wrtite the aims geometry.in files
print('Writting selected geometries...', end='')
for i in random_indices:
    # write(outputfile, md[i], append=True)
    ifolder = prefix_folder+str(i)
    os.makedirs(ifolder, exist_ok=True)
    write_aims(ifolder + '/geometry.in', md[i])
print('Done!')

if write_control:
    # TODO: write the control file
    print('Control file writting not yet implemented. sory!')
    print('Bring your own control file')
