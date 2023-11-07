#!/usr/bin/env python3

import os
import sys
import numpy as np

"""
Parse  aims MD outfile and plots T or pseudo-Hamiltoninan (conserverd quantity). Options not implemented
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
parser.add_argument('-l', '--loadfile', default=None,
                    help='Load a previously processed loadfile .npy file')

# TODO: Implement this options?
parser.add_argument('-n', '--steps', default=1, nargs='*',
                    help='Space separated integer values. Extract information from the requested steps. Positive values get the requested time step, negative values get the -step from the last (-1 implies the last step)(NEGATIVE VALUES NOT YET IMPLEMENTED)')

parser.add_argument('-r', '--random', default=0, help='Extracts as many as random geometries as requested')
parser.add_argument('-s', '--startingstep', default=0, help='Starting step from which start the random sampling')
parser.add_argument('-f', '--finalstep', default=0, help='Steps to remove from the end')
parser.add_argument('-w', '--window', default=None, help='window of moving average')
parser.add_argument('-a', '--averaged', action="store_true", help='Plot cumulative average')

args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile
loadfile = args.loadfile
starting_step = int(args.startingstep)
final_step = int(args.finalstep)
steps = int(args.steps)
random = args.random
window = args.window
averaged = args.averaged


# Parse the output file
def parse_MD():
    md_time_step = None
    md_step = []
    temperature = []

    for_the_array = []

    print('Reading:', infile)
    with open(infile) as f:
        lines = f.readlines()
        completion = len(lines)
        for i, line in enumerate(lines):
            sys.stdout.write('\r')
            sys.stdout.write("[%-20s] %d%%" % ('='*int(i/completion*20), i/completion*100+1))
            sys.stdout.flush()
#             if 'Molecular dynamics time step' in line:
#                 md_time_step = float(line.split()[5])
#             if 'Number of atoms' in line:
#                 n_atoms = int(line.split()[5])
#             if 'Number of lattice vectors' in line:
#                 n_lattice_vectors = int(line.split()[6])
#             if 'lattice_vector' in line:
#                 lattice_vector.append(line)
            if '| Time step number' in line:
                md_step.append(int(line.split()[5]))
#             if 'tomic structure (and velocities) ' in line:
#                 atomic_coordinates = []
#                 species = []
#                 velocities = []
#                 for n in range(n_atoms):
#                     species.append(lines[i+2+2*n].split()[4])
#                     atomic_coordinates.append(lines[i+2+2*n].split()[1:4])
#                     velocities.append(lines[i+3+2*n].split()[1:4])
#                 atomic_coordinates = np.array(atomic_coordinates, dtype=float)
#                 velocities = np.array(velocities, dtype=float)
#                 array_species.append(species)
#                 array_atoms.append(atomic_coordinates)
#                 array_velocities.append(velocities)
#             if 'Total atomic forces' in line:
#                 forces = []
#                 for n in range(n_atoms):
#                     forces.append(lines[i+1+n].split()[2:5])
#                 forces = np.array(forces, dtype=float)
#                 array_forces.append(forces)
#             if '| Total time     ' in line:
#                 total_time = line.split()[-2:-1]
            if 'Temperature (nuclei) ' in line:
                temperature.append(float(line.split()[4]))

    for i in range(len(md_step)):
        # for_the_array.append((md_step[i], array_species[i], array_atoms[i], array_velocities[i], temperature[i], array_forces[i]))
        for_the_array.append((md_step[i], temperature[i]))

# Define the data type for the structured array
    data_type = np.dtype([
        ('step', int),
        # ('species', 'S2', (n_atoms,)),
        # ('coordinates', (float, (n_atoms, 3))),  # n_atoms is known at this point
        # ('velocities', (float, (n_atoms, 3))),
        ('temperature', float),
        # ('forces', (float, (n_atoms, 3))),
        # ('energy', float),
        # ('polarizability', (float, (3, 3))),
    ])

    data_array = np.array(for_the_array, dtype=data_type)
    data_array.sort()
    np.save('step_vs_T', data_array, allow_pickle=True)
    return data_array

def write_step(dd):
    geom_file_name = 'geometry.in.'

    # for s in steps:
    output_file = geom_file_name + str(steps)

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

def plot_temperature(d):

    import matplotlib.pyplot as plt

    def moving_average(x, w):
        return np.convolve(x, np.ones(w), 'valid') / w

    ts = d['step']
    T = d['temperature']
    timestep = ts[starting_step:len(ts) - final_step]
    temperature = T[starting_step:len(T) - final_step]
    print(final_step)
    print(len(timestep))
    print(len(temperature))

    plt.plot(timestep, temperature, '--', lw=0.3, label='Instant T', color='grey')
    if averaged:
        avg = np.cumsum(temperature)/np.arange(1, len(temperature)+1)
        print(len(avg))
        plt.plot(timestep, avg, label='cumulative avg')
    if window:
        ravg = moving_average(temperature[::-1], int(window))
        x = np.array(range(len(ravg)))
        x += int(window) + starting_step
        plt.plot(x, ravg[::-1], label=f'rolling avg {window}')
    plt.title('FHI-AIMS MD')
    plt.legend()
    plt.show()


if not loadfile:
    data_array = parse_MD()
else:
    try:
        print('Loading file:', loadfile)
        data_array = np.load(loadfile)
    except FileNotFoundError as e:
        print('File not found:', loadfile)
        print('Parsing file ', infile)
        data_array = parse_MD()

# write_step(data_array[data_array['step'] == steps])
plot_temperature(data_array)
