#!/usr/bin/env python3

import os

import argparse
from ase.io import read, write
from ase.io.aims import write_aims
import numpy as np

'''
Reads the geometry and forces from an aims.out files in the prefix folders, or folders in the current dir, and generates a train.xyz file that contains the polarizabilities
'''

# Read the arguments
parser = argparse.ArgumentParser(
    prog='fhi_convert_polarizabilities_aims2train.py',
    description='Reads polarizabilities from aims.out files in the prefix folders and generates a train.xyz file',
)

inputfile = 'aims.out'
outputfile = 'train.xyz'
# outputfile_extxyz = 'sampled.extxyz'

parser.add_argument('-i', '--inputfile', default=inputfile,
                    help='input file: aims.out . Read in the aims outputfile to extract the forces (and maybe polarizabilities)')
parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='output file: train.xyz file.')
parser.add_argument('-p', '--prefix', default='',
                    help='Prefix of the folders. []')
parser.add_argument('-n', '--samples', default=0,
                    help='Total number of geometries to sample. If == 0, all will be included [0]')

args = parser.parse_args()
inputfile = args.inputfile
outfile = args.outputfile
num_samples = int(args.samples)
prefix_folder = args.prefix

print(f'Reading {inputfile} files in folders with prefix {prefix_folder} ...', end='', flush=True)
# Get a list of subdirectories (folders) within the current directory
subdirs = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith(prefix_folder)]
print('Done!')

indices = range(len(subdirs))
if num_samples > 0:
    indices = np.random.choice(len(subdirs), num_samples, replace=False)

final_samples = len(indices)
print(f'Writing {outputfile} with a total of {final_samples} samples ...', end='', flush=True)
# Loop over each subdir and read a files inside it
for i in indices:
    path = os.path.join(subdirs[i], inputfile)
    if os.path.exists(path):
        # This is done because ASE does not read the polarizabilities. But for forces is not needed
        # maybe read dielectric_tensor?
        # polarizability_tensor = []

        # with open(path) as f:
        #     lines = f.readlines()
        #     for i, line in enumerate(lines):
        #         if 'Polarizability (Bohr^3) :' in line:    # Line when using DFPT_centralised (new module)
        #             for n in range(3):
        #                 polarizability_tensor.append(lines[i+1+n].split()[0:3])
        #             polarizability_tensor = np.array(polarizability_tensor, dtype=float)
        #         # Full Polarizability tensor for crystals
        #         if 'DFPT for dielectric_constant:' in line:
        #             for n in range(3):
        #                 polarizability_tensor.append(lines[i+1+n].split()[0:3])
        #             polarizability_tensor = np.array(polarizability_tensor, dtype=float)
        #     if len(polarizability_tensor) == 0:
        #         print('File ', f.name, ' contains no polarizability')
        #         continue
        iaims = read(path)
        # iaims.info['REF_polarizability'] = polarizability_tensor.flatten()
        write(outputfile, iaims, append=True)
print('Done!')
