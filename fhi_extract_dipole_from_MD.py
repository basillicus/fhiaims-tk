#!/usr/bin/env python3

import os

import sys
import argparse
from ase.io import read, write
from ase.io.aims import write_aims
import numpy as np

'''
Reads the geometry and dipole from an MD aims.out file  and generates a numpy array
'''

# Read the arguments
parser = argparse.ArgumentParser(
    prog='fhi_extract_dipole_from_MD.py',
    description='Reads dipoles from MD aims.out files to generate dipoles.npy file',
)

inputfile = 'aims.out'
outputfile = 'dipoles'
# outputfile_extxyz = 'sampled.extxyz'

parser.add_argument('-i', '--inputfile', default=inputfile,
                    help='input file: aims.out . Read in the aims outputfile to extract the forces (and maybe polarizabilities)')
parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='output file: train.xyz file.')
# parser.add_argument('-p', '--prefix', default='',
#                     help='Prefix of the folders. [md_sample_]')
# parser.add_argument('-n', '--samples', default=0,
#                     help='Total number of geometries to sample. If == 0, all will be included [0]')

args = parser.parse_args()
inputfile = args.inputfile
outfile = args.outputfile
# prefix_folder = args.prefix
# num_samples = int(args.samples)

# Parse the output file
def parse_MD(infile):
    md_time_step = None
    md_step = []
    lattice_vector = []
    n_atoms = None
    atoms = []
    species = []
    velocities = []
    forces = []
    temperature = []
    dipole = []

    array_atoms = []
    array_species = []
    array_lattice_vector = []
    array_velocities = []
    array_forces = []
    for_the_array = []

    print('Reading:', infile)
    with open(infile) as f:
        lines = f.readlines()
        completion = len(lines)
        for i, line in enumerate(lines):
            #             sys.stdout.write('\r')
            # sys.stdout.write("[%-20s] %d%%" % ('='*int(i/completion*100), i/completion*100+1))
            # sys.stdout.flush()
            if 'Molecular dynamics time step' in line:
                md_time_step = float(line.split()[5])
            if 'Number of atoms' in line:
                n_atoms = int(line.split()[5])
            #  if 'Number of lattice vectors' in line:
            #      n_lattice_vectors = int(line.split()[6])
            #  if 'lattice_vector' in line:
            #      lattice_vector.append(line)
            if line.startswith('  | Time step number'):
                md_step.append(int(line.split()[5]))
            if line.startswith('  | Total dipole moment'):
                # dipole.append(line.split()[6:9])
                dipole.append(np.array(np.fromstring(' '.join(line.split()[6:9]), dtype=float, sep=' ')))
            if 'tomic structure (and velocities) ' in line:
                atomic_coordinates = []
                species = []
                velocities = []
                for n in range(n_atoms):
                    species.append(lines[i+2+2*n].split()[4])
                    atomic_coordinates.append(lines[i+2+2*n].split()[1:4])
                    velocities.append(lines[i+3+2*n].split()[1:4])
                atomic_coordinates = np.array(atomic_coordinates, dtype=float)
                velocities = np.array(velocities, dtype=float)
                array_species.append(species)
                array_atoms.append(atomic_coordinates)
                array_velocities.append(velocities)
            if 'Total atomic forces' in line:
                forces = []
                for n in range(n_atoms):
                    forces.append(lines[i+1+n].split()[2:5])
                forces = np.array(forces, dtype=float)
                array_forces.append(forces)
            if '| Total time     ' in line:
                total_time = line.split()[-2:-1]
            if 'Temperature (nuclei) ' in line:
                temperature.append(float(line.split()[4]))

    for i in range(len(md_step)):
        for_the_array.append((md_step[i], array_species[i], array_atoms[i], array_velocities[i], dipole[i], temperature[i], array_forces[i]))

# Define the data type for the structured array
    data_type = np.dtype([
        ('step', int),
        ('species', 'S2', (n_atoms,)),
        ('coordinates', (float, (n_atoms, 3))),  # n_atoms is known at this point
        ('velocities', (float, (n_atoms, 3))),
        ('dipole', (float, (3,))),  # n_atoms is known at this point
        ('temperature', float),
        ('forces', (float, (n_atoms, 3))),
        # ('energy', float),
        # ('polarizability', (float, (3, 3))),
    ])

    data_array = np.array(for_the_array, dtype=data_type)
    data_array.sort()
    np.save(outputfile, data_array, allow_pickle=True)
    return data_array


print(f'Parsing {inputfile} file ...', end='', flush=True)
data_array = parse_MD(inputfile)


# Get a list of subdirectories (folders) within the current directory
# subdirs = [d for d in os.listdir('.') if os.path.isdir(d) and d.startswith(prefix_folder)]

# for i in indices:
#     path = os.path.join(subdirs[i], inputfile)
#     if os.path.exists(path):
#         polarizability_tensor = []
#
#         with open(path) as f:
#             lines = f.readlines()
#             for i, line in enumerate(lines):
#                 if 'Polarizability (Bohr^3) :' in line:    # Line when using DFPT_centralised (new module)
#                     for n in range(3):
#                         polarizability_tensor.append(lines[i+1+n].split()[0:3])
#                     polarizability_tensor = np.array(polarizability_tensor, dtype=float)
#                 # Full Polarizability tensor for crystals
#                 if 'DFPT for dielectric_constant:' in line:
#                     for n in range(3):
#                         polarizability_tensor.append(lines[i+1+n].split()[0:3])
#                     polarizability_tensor = np.array(polarizability_tensor, dtype=float)
#             if len(polarizability_tensor) == 0:
#                 print('File ', f.name, ' contains no polarizability')
#                 continue
#             iaims = read(path)
#             iaims.info['REF_polarizability'] = polarizability_tensor.flatten()
#             write(outputfile, iaims, append=True)
print('Done!')
