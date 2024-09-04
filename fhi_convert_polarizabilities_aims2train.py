#!/usr/bin/env python3

import os

import argparse
from ase.io import read, write
from ase.io.aims import write_aims
import numpy as np

'''
Reads the geometry and polarizability from an aims.out files in the prefix folders and generates a train.xyz file that contains the polarizabilities
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
parser.add_argument('-p', '--prefix', default='md_sample_',
                    help='Prefix of the folders. [md_sample_]')
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
# Loop over each subdir and read a file named "example.txt" inside it
for i in indices:
    path = os.path.join(subdirs[i], inputfile)
    if os.path.exists(path):
        # lattice_vector = []
        # atoms = []
        # atomic_coordinates = []
        # species = []
        # forces = []
        polarizability_tensor = []
        # polarizability_elements = []
        # polarizability = None
        # n_atoms = None

        with open(path) as f:
            lines = f.readlines()
            for i, line in enumerate(lines):
                # if 'Number of atoms' in line:
                #     n_atoms = int(line.split()[5])
                # if 'Number of lattice vectors' in line:
                #     n_lattice_vectors = int(line.split()[6])
                # if 'lattice_vector' in line:
                #     lattice_vector.append(line.split()[1:4])
                # if '  Atom  ' in line:
                #     atomic_coordinates = []
                #     species = []
                #     for n in range(n_atoms):
                #         species.append(lines[i+1+n].split()[3])
                #         atomic_coordinates.append(lines[i+1+n].split()[4:7])
                #     atomic_coordinates = np.array(atomic_coordinates, dtype=float)
                # if 'Total atomic forces' in line:
                #     for n in range(n_atoms):
                #         forces.append(lines[i+1+n].split()[2:5])
                #     forces = np.array(forces, dtype=float)
                # # Full Polarizability tensor for molecules
                # if 'DFPT for polarizability:' in line:    # Line when not using DFPT_centralised (old module)
                #     for n in range(3):
                #         polarizability_tensor.append(lines[i+1+n].split()[0:3])
                #     polarizability_tensor = np.array(polarizability_tensor, dtype=float)
                if 'Polarizability (Bohr^3) :' in line:    # Line when using DFPT_centralised (new module)
                    for n in range(3):
                        polarizability_tensor.append(lines[i+1+n].split()[0:3])
                    polarizability_tensor = np.array(polarizability_tensor, dtype=float)
                # Full Polarizability tensor for crystals
                if 'DFPT for dielectric_constant:' in line:
                    for n in range(3):
                        polarizability_tensor.append(lines[i+1+n].split()[0:3])
                    polarizability_tensor = np.array(polarizability_tensor, dtype=float)
                # if '| Total energy of the ' in line:
                #     total_energy = float(line.split()[11])
            if len(polarizability_tensor) == 0:
                print('File ', f.name, ' contains no polarizability')
                continue
            iaims = read(path)
            iaims.info['polarizability'] = polarizability_tensor.flatten()
            write(outputfile, iaims, append=True)
print('Done!')
