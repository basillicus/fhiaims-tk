#!/usr/bin/env python3

import os
import sys
import numpy as np

"""
Parse  aims recursively and extract geometry and velocity from the requested step. Options not yet implemented
Internally creates a numpy array with the following structure:
[steps (int,1), species(str, n_atoms), coordinates (float, (n_atoms x 3), velocities (n_atoms x 3), forces(float, (n_atoms x 3)]

and saves the binary file that can be loaded with np.load('MD.npy')
subarrays can be created like:
    (d['step'], d['polarizability'])

Generates a geoemtry.in file including the velocities for each requested step
"""
import argparse
from config import fhi_aims_outputfile, numpy_outputfile

parser = argparse.ArgumentParser(
    prog='fhi_get_output_info.py',
    description='Extracts information from an aims ouptput (aims.out) file',
)

inputfile = fhi_aims_outputfile
outputfile = numpy_outputfile

parser.add_argument('-i', '--inputfile', default=inputfile, help='Input file for the script: an AIMS MD ouptut file to be parsed')
parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='outputfile file  where the information will be writen as an Numpy binary array')
parser.add_argument('-l', '--loadfile', default=None,
                    help='Load a previously processed loadfile .npy file')

parser.add_argument('-n', '--steps', default=1, nargs='*',
                    help='Space separated integer values. Extract information from the requested steps. Positive values get the requested time step, negative values get the -step from the last (-1 implies the last step)(NEGATIVE VALUES NOT YET IMPLEMENTED)')


# TODO: Implement this options:
# parser.add_argument('-n', '--steps', default=1, nargs='*',
#                     help='Space separated integer values. Extract information from the requested steps. Positive values get the requested time step, negative values get the -step from the last (-1 implies the last step)(NEGATIVE VALUES NOT YET IMPLEMENTED)')
# 
# parser.add_argument('-r', '--random', default=0, help='Extracts as many as random geometries as requested')
# parser.add_argument('-s', '--startingstep', default=0, help='Starting step from which start the random sampling')
# parser.add_argument('-f', '--finalstep', default=-1, help='Last step from which take the random sampling')


args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile
loadfile = args.loadfile
steps = args.steps
# random = args.random

# Parse the output file
def parse_MD():
    md_time_step = None
    md_step = []
    lattice_vector = []
    n_atoms = None
    atoms = []
    species = []
    velocities = []
    forces = []
    temperature = []

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
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%%" % ('='*int(i/completion*100), i/completion*100+1))
            sys.stdout.flush()
            if 'Molecular dynamics time step' in line:
                md_time_step = float(line.split()[5])
            if 'Number of atoms' in line:
                n_atoms = int(line.split()[5])
            if 'Number of lattice vectors' in line:
                n_lattice_vectors = int(line.split()[6])
            if 'lattice_vector' in line:
                lattice_vector.append(line)
            if '| Time step number' in line:
                md_step.append(int(line.split()[5]))
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
        for_the_array.append((md_step[i], array_species[i], array_atoms[i], array_velocities[i], temperature[i], array_forces[i]))

# Define the data type for the structured array
    data_type = np.dtype([
        ('step', int),
        ('species', 'S2', (n_atoms,)),
        ('coordinates', (float, (n_atoms, 3))),  # n_atoms is known at this point
        ('velocities', (float, (n_atoms, 3))),
        ('temperature', float),
        ('forces', (float, (n_atoms, 3))),
        # ('energy', float),
        # ('polarizability', (float, (3, 3))),
    ])

    data_array = np.array(for_the_array, dtype=data_type)
    data_array.sort()
    np.save(outputfile, data_array, allow_pickle=True)
    return data_array

def write_step(dd):
    geom_file_name = 'geometry.in.'

    # for s in steps:
    output_file = geom_file_name + str(dd[0][0])

    # Open the output file for writing
    with open(output_file, 'w') as f:
        for row in dd:
            species = row[1]
            coordinates = row[2]
            velocities = row[3]
            for i in range(len(species)):
                # Write atom coordinates
                f.write(f'atom  {coordinates[i][0]:.6f} {coordinates[i][1]:.6f} {coordinates[i][2]:.6f} {species[i].decode()} \n')

                # Write velocities
                f.write(f'velocity  {velocities[i][0]:.6f} {velocities[i][1]:.6f} {velocities[i][2]:.6f}\n')

    # Print a message indicating the file has been written
    print(f'Data has been written to {output_file}')


if not loadfile:
    data_array = parse_MD()
else:
    try:
        print('Loading file:', loadfile)
        data_array = np.load(loadfile)
    except:
        print('File not found:', loadfile)
        print('Parsing file ', infile)
        data_array = parse_MD()

print(steps)
for s in steps:
    write_step(data_array[data_array['step'] == int(s)])
