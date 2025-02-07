#!/usr/bin/env python3

import numpy as np
import argparse

from config import dataset_all_outfile
from utils import merge_datasets, write_extxyz_file

"""
Combine dipoles and polarizabilities numpy arrays into a merged .extxyz file
It compares coordinates before combining data
Output file is ready for training MACE ML potentials
"""


parser = argparse.ArgumentParser(
    prog='fhi_merge_all_data.py',
    description='Merge dipoles and polarizabilities numpy arrays into an .extxyz file',
)

outputfile = dataset_all_outfile

parser.add_argument('-d', '--dipoles', default=[], nargs='+', required=True,
                    help='dipole file or files, as numpy array containing keys: step, coordinates dipoles and forces')

parser.add_argument('-p', '--polarizabilities', default=[], nargs='+', required=True,
                    help='polarizability file or files as numpy array containing step, coordinates, polarizabilities and energies')

parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='outputfile file where the geoemtries are writen in an .extxyz file')

parser.add_argument('--pbc', default=False, action='store_true',
                    help='Add this flag if the system is periodic')

parser.add_argument('--sort', default=False, action='store_true',
                    help='sort .extxyz by time step')

parser.add_argument('--convert_dipoles', default=False, action='store_true',
                    help='convert dipoles from eAng to Debye')

parser.add_argument('--convert_alphas', default=False, action='store_true',
                    help='convert polarizabilities from Bohr^3 to me A^2/V ')

# parser.add_argument('-a', '--all', action='store_true', default=True,
#                     help='all steps will be read/writen')

args = parser.parse_args()
dipole_files = args.dipoles
polarizability_files = args.polarizabilities
is_periodic = args.pbc
do_sort = args.sort
convert_dipoles = args.convert_dipoles
convert_alphas = args.convert_alphas
outfile = args.outputfile

# Handle errors
if len(dipole_files) == 0:
    print('ERROR: Give me my dipoles!')
if len(polarizability_files) == 0:
    print('ERROR: Give me my polarizabilities!')

# Reading dipoles
print('Joining dipole files:', dipole_files)
dipoles = np.load(dipole_files[0])
if len(dipole_files) > 1:
    for i in range(1, len(dipole_files)):
        dipoles = np.concatenate((dipoles, np.load(dipole_files[i])), axis=0)
print('Dipole shape:', dipoles.shape)

# Reading polarizabilities
print('Joining polarizabilities files:', polarizability_files)
polarizabilities = np.load(polarizability_files[0])
if len(polarizability_files) > 1:
    for i in range(1, len(polarizability_files)):
        polarizabilities = np.concatenate((polarizabilities, np.load(polarizability_files[i])), axis=0)
print('polarizabilities shape:', polarizabilities.shape)

# Merging arrays
print('Merging dipoles and polarizabilities...', end='', flush=True)
merged_dataset = merge_datasets(dipoles, polarizabilities,
        search_range=2, periodic=is_periodic,
        convert_dipoles=convert_dipoles, convert_alphas=convert_alphas)
print('OK')
print('Merged data shape:', merged_dataset.shape)

if do_sort:
    merged_dataset = np.sort(merged_dataset, order='step')

# Write .extxyz file
print('Writting extxyz file...', end='', flush=True)
write_extxyz_file(merged_dataset)
print('OK')
