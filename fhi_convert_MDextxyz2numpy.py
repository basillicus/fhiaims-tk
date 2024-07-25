#!/usr/bin/env python3

import os
import sys
import numpy as np
from tqdm import tqdm

"""
Parse  extxyz outfile  and saves information to a numpy dataset. Options not implemented
Generates a numpy binary file with the following structure:
[steps (int,1), species(str, n_atoms), coordinates (float, (n_atoms x 3), total_energy(float, 1), forces(float, (n_atoms x 3), velocities (float, (ntaoms x 3)]

can be loaded with np.load('md.npy', allow_pickle=True)
subarrays can be created like:
    (d['step'], d['polarizability'])
"""

import argparse
from config import fhi_aims_outputfile, information_outputfile

parser = argparse.ArgumentParser(
    prog='fhi_convert_MDextxyz2nump.py',
    description='Converts a MD in extxyz format into a numpy array with a dataset strcuture',
)

inputfile = fhi_aims_outputfile
outputfile = information_outputfile

parser.add_argument('-i', '--inputfile', default=inputfile, help='This is the input file for the script, an AIMS ouptut file to be parsed')
parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='This is the outputfile file for the script, where the geoemtries are writen in and .extxyz file')
parser.add_argument('--pbc', default=False, action='store_true',
                    help='Set the system has Periodic Boundary Condition')

args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile
pbc = args.pbc

if pbc:
    search_pattern = 'Latti'
else:
    search_pattern = 'Prop'

# Find file to be parsed
def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result


parsing_file = infile

for_the_array = []
lattice_vector = []
atoms = []
energy = []
polarizability_tensor = []
polarizability_elements = []
n_atoms = None
step = 0

# Parse the output file
with open(parsing_file) as f:
    n_atoms = int(f.readline())
    lines = f.readlines()
    completion = len(lines)
    print('reading file...')
    for i, line in tqdm(enumerate(lines), total=completion, colour='blue'):
        # sys.stdout.write('\r')
        # sys.stdout.write("[%-20s] %d%%" % ('='*int(i/completion*20), i/completion*100+1))
        # sys.stdout.flush()
        if line.startswith(search_pattern):
            atomic_coordinates = []
            species = []
            forces = []
            velocities = []
            step += 1
            for n in range(n_atoms):
                split_line = lines[i+1+n].split()
                species.append(split_line[0])
                atomic_coordinates.append(split_line[1:4])
                velocities.append(split_line[4:7])
                forces.append(split_line[7:10])
            species = np.array(species, dtype=str)
            atomic_coordinates = np.array(atomic_coordinates, dtype=float)
            velocities = np.array(velocities, dtype=float)
            forces = np.array(forces, dtype=float)
            for_the_array.append((step, species, atomic_coordinates, velocities, forces))

print()
print('Parsing file completed')
print('Writing numpy dataset')
data_type = np.dtype([
    ('step', int),
    ('species', 'S2', (n_atoms,)),
    ('coordinates', (float, (n_atoms, 3))),  # n_atoms is known at this point
    ('velocities', (float, (n_atoms, 3))),
    ('forces', (float, (n_atoms, 3))),  # n_atoms is known at this point
    # ('energy', float),
    # ('lattice_vector', (float, (3, 3))),  # n_atoms is known at this point
])

data_array = np.array(for_the_array, dtype=data_type)
#
# data_array.sort()
# np.save('md_numpy', data_array, allow_pickle=True)
# np.save(outfile, data_array, allow_pickle=True)
np.save(outfile, data_array, allow_pickle=False)
print('Done!')

# STEP BY STEP
# data_array = np.empty(len(for_the_array), dtype=data_type)
#
# # Populate the structured array with the collected data
# for i, item in enumerate(for_the_array):
#     step, species, atomic_coordinates, total_energy, forces, polarizability = item
#     data_array[i] = (step, species, atomic_coordinates, total_energy, forces, polarizability)
