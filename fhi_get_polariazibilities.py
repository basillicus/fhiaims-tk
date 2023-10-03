#!python

import os
import numpy as np

"""
Parse  aims outfilesi recursively and extract information. Options not implemented
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
for parsing_file in files:
    lattice_vector = []
    atoms = []
    species = []
    forces = []
    polarizability_tensor = []
    kpoints = None
    total_time = None
    total_scf_iterations = None
    polarizability = None
    n_atoms = None
    # Parse the output file
    with open(parsing_file) as f:
        lines = f.readlines()
        # Get the step number from folder's name
        step = parsing_file.split('_')[2]
        for i, line in enumerate(lines):
            # IMPROVEMENT Check the n_atoms does not change
            if 'Number of atoms' in line:
                n_atoms = int(line.split()[5])
            if 'Number of lattice vectors' in line:
                n_lattice_vectors = int(line.split()[6])
            if 'lattice_vector' in line:
                lattice_vector.append(line)
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
            # Full Polarizability tensor
            if 'DFPT for polarizability:' in line:
                for n in range(3):
                    polarizability_tensor.append(lines[i+1+n].split()[0:3])
                polarizability_tensor = np.array(polarizability_tensor, dtype=float)
            # Polarizability components
            if '| Polarizability' in line:
                polarizability = line.split()[2:]
                polarizability = np.array(polarizability, dtype=float)
            if '| Total energy of the ' in line:
                total_energy = float(line.split()[11])

        for_the_array.append((step, species, atomic_coordinates, total_energy, forces, polarizability_tensor))

# Define the data type for the structured array
data_type = np.dtype([
    ('step', int),
    ('species', 'S2', (n_atoms,)),
    ('coordinates', (float, (n_atoms, 3))),  # n_atoms is known at this point
    ('energy', float),
    ('forces', (float, (n_atoms, 3))),  # n_atoms is known at this point
    ('polarizability', (float, (3, 3))),
])

data_array = np.array(for_the_array, dtype=data_type)
np.save('polarizabilities', data_array, allow_pickle=True)

# STEP BY STEP
# data_array = np.empty(len(for_the_array), dtype=data_type)
#
# # Populate the structured array with the collected data
# for i, item in enumerate(for_the_array):
#     step, species, atomic_coordinates, total_energy, forces, polarizability = item
#     data_array[i] = (step, species, atomic_coordinates, total_energy, forces, polarizability)
