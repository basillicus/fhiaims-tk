#!/usr/bin/env python3

import os
import sys
import numpy as np

"""
Parse  aims outfiles recursively and extract information. Options not implemented
Generates a numpy binary file with the following structure:
[steps (int,1), species(str, n_atoms), coordinates (float, (n_atoms x 3), total_energy(float, 1), forces(float, (n_atoms x 3), polarizability (float, 6)]

can be loaded with np.load('polarizabilities.npy', allow_pickle=True)
subarrays can be created like:
    (d['step'], d['polarizability'])
"""

import argparse
from config import fhi_aims_outputfile, information_outputfile

parser = argparse.ArgumentParser(
    prog='fhi_get_output_info.py',
    description='Extracts information from an aims ouptput (aims.out) file',
)

inputfile = fhi_aims_outputfile
outputfile = information_outputfile

parser.add_argument('-i', '--inputfile', default=inputfile, help='This is the input file for the script, an AIMS ouptut file to be parsed')
parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='This is the outputfile file for the script, where the geoemtries are writen in and .extxyz file')

# TODO: Implement this options?
#
parser.add_argument('-k', '--kpoints', action='store_true', help='get k-points grid')
parser.add_argument('-T', '--totaltime', action='store_true', help='get total time (Wall time)')
parser.add_argument('-n', '--iterations', action='store_true', help='get total iterations')

parser.add_argument('-d', '--dielectric', action='store_true', help='get polarizability/dielectric')

parser.add_argument('-t', '--steptime', action='store_true', help='get  a single step time')
parser.add_argument('-v', '--volume', action='store_true', help='get cell volume (not yet implemeted)')
parser.add_argument('-p', '--pressure', action='store_true', help='get pressure (not yet implemeted)')

args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile
get_kpoints = args.kpoints
get_total_time = args.totaltime
get_iterations = args.iterations
get_polarizabilty = args.dielectric

# Find file to be parsed
def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result


files = find_all(infile, './')

for_the_array = []
completion = len(files)
for i, parsing_file in enumerate(files):
    sys.stdout.write('\r')
    sys.stdout.write("[%-20s] %d%%" % ('='*int(i/completion*20), i/completion*100+1))
    sys.stdout.flush()
    lattice_vector = []
    atoms = []
    atomic_coordinates = []
    species = []
    forces = []
    polarizability_tensor = []
    polarizability_elements = []
    kpoints = None
    total_time = None
    total_scf_iterations = None
    polarizability = None
    n_atoms = None
    # Parse the output file
    with open(parsing_file) as f:
        lines = f.readlines()
        # Get the step number from folder's name
        # step = parsing_file.split('_').split('/')[0]
        step = parsing_file.split('_')[-1].split('/')[0]
        for i, line in enumerate(lines):
            if 'Number of atoms' in line:
                n_atoms = int(line.split()[5])
            if 'Number of lattice vectors' in line:
                n_lattice_vectors = int(line.split()[6])
            if 'lattice_vector' in line:
                lattice_vector.append(line.split()[1:4])
            if '  Atom  ' in line:
                atomic_coordinates = []
                species = []
                for n in range(n_atoms):
                    species.append(lines[i+1+n].split()[3])
                    atomic_coordinates.append(lines[i+1+n].split()[4:7])
                atomic_coordinates = np.array(atomic_coordinates, dtype=float)
            if 'Total atomic forces' in line:
                for n in range(n_atoms):
                    forces.append(lines[i+1+n].split()[2:5])
                forces = np.array(forces, dtype=float)
            if 'k_grid      ' in line:
                kpoints = line.split()[-3:]
            if '| Total time     ' in line:
                total_time = line.split()[-2:-1]
            if 'Self-consistency cycle converged' in line:
                total_scf_iterations = int(lines[i+4].split()[4])
            # Full Polarizability tensor for molecules
            if 'DFPT for polarizability:' in line:    # Line when not using DFPT_centralised (old module)
                for n in range(3):
                    polarizability_tensor.append(lines[i+1+n].split()[0:3])
                polarizability_tensor = np.array(polarizability_tensor, dtype=float)
            if 'Polarizability (Bohr^3) :' in line:    # Line when using DFPT_centralised (new module)
                for n in range(3):
                    polarizability_tensor.append(lines[i+1+n].split()[0:3])
                polarizability_tensor = np.array(polarizability_tensor, dtype=float)
            # Full Polarizability tensor for crystals
            if 'DFPT for dielectric_constant:' in line:
                for n in range(3):
                    polarizability_tensor.append(lines[i+1+n].split()[0:3])
                polarizability_tensor = np.array(polarizability_tensor, dtype=float)
            # Polarizability components
            # NOTE: not used, and we avoid the formating error '*****'
                # if '| Polarizability' in line:
                #     polarizability = line.split()[2:]
                #     polarizability = np.array(polarizability, dtype=float)
                # if 'Polarizability (Bohr)' in line:
                #     for n in range(3):
                #         polarizability_elements.append(lines[i+1+n].split()[0:3])
                #         pe = polarizability_elements
                #     polarizability = np.array([pe[0], pe[4], pe[8], pe[1], pe[2], pe[5]], dtype=float)
            if '| Total energy of the ' in line:
                total_energy = float(line.split()[11])

        # Check the calculation is finished properly:
        if len(atomic_coordinates) == 0:
            print('File ', f.name, ' contains no coordinates')
            continue
        if len(polarizability_tensor) == 0:
            print('File ', f.name, ' contains no polarizability')
            continue

        if len(lattice_vector) == 0:     # is a molecule
            for_the_array.append((step, species, atomic_coordinates, total_energy, polarizability_tensor))
        else:   # is a  a crystal
            lattice_vector = np.array(lattice_vector)
            for_the_array.append((step, species, lattice_vector, atomic_coordinates, total_energy, polarizability_tensor))

print()
# Define the data type for the structured array
if len(lattice_vector) == 0:     # is a molecule
    data_type = np.dtype([
        ('step', int),
        ('species', 'S2', (n_atoms,)),
        # ('lattice_vector', (float, (3, 3))),  # n_atoms is known at this point
        ('coordinates', (float, (n_atoms, 3))),  # n_atoms is known at this point
        ('energy', float),
        # # ('forces', (float, (n_atoms, 3))),  # n_atoms is known at this point
        ('polarizability', (float, (3, 3))),
    ])
else:
    data_type = np.dtype([
        ('step', int),
        ('species', 'S2', (n_atoms,)),
        ('lattice_vector', (float, (3, 3))),  # n_atoms is known at this point
        ('coordinates', (float, (n_atoms, 3))),  # n_atoms is known at this point
        ('energy', float),
        # # ('forces', (float, (n_atoms, 3))),  # n_atoms is known at this point
        ('polarizability', (float, (3, 3))),
    ])

data_array = np.array(for_the_array, dtype=data_type)
data_array.sort()
# np.save('polarizabilities', data_array, allow_pickle=True)
np.save('polarizabilities', data_array)

# STEP BY STEP
# data_array = np.empty(len(for_the_array), dtype=data_type)
#
# # Populate the structured array with the collected data
# for i, item in enumerate(for_the_array):
#     step, species, atomic_coordinates, total_energy, forces, polarizability = item
#     data_array[i] = (step, species, atomic_coordinates, total_energy, forces, polarizability)
