#!/usr/bin/env python3

""" Generates individual geometries.in files, one from each MD step, from an CP2K file containing a MD
Returns one folder per geometry where the polarazibility will be calculated
"""

import os
import sys
import shutil

import argparse
from config import fhi_aims_outputfile, visualization_outputfile

parser = argparse.ArgumentParser(
    prog='fhi_prepare_polarizabilities_from_MD.py',
    description='Prepares polarizability calculations from a MD',
)

inputfile = cp2k_outputfile
outputfile = 'polar'

parser.add_argument('-i', '--inputfile', default=None, required=True, help='This is the input file for the script, an AIMS ouptut file to be parsed')
parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='This is the prefix for the outputfile folder, where the geoemtries will be writen in aims format')
parser.add_argument('-c', '--copycontrol', action='store_true', default=False,
                    help='If control.in file is in the folder, it will be copied to the folders')
parser.add_argument('-s', '--steps', default=1, type=int,
                    help='Every how many steps the polarizability will be calculated')
parser.add_argument('-n', '--initialstep', default=0, type=int,
                    help='start form this step')
parser.add_argument('-f', '--finalstep', default=None, type=int,
                    help='end at this step')
parser.add_argument('-d', '--idstart', default=0, type=int,
                    help='starting number for calculation ID')

# Global variables initialization
args = parser.parse_args()
init_geom = args.inputfile
outfile = args.outputfile
step_interval = args.steps
do_cp_control = args.copycontrol
initial_step = args.initialstep
final_step = args.finalstep
calc_ID = args.idstart

if not init_geom:
    print('Plase, provide the cp2k (.xyz) file with the MD geometries. Use -i inputfile')
    sys.exit(0)

if do_cp_control:
    if not os.path.exists('control.in'):
        print('WARNING: control.in not found. Will not be copied to folders!')
        do_cp_control = False

# Per step variables initialization
lattice_vector = []
try:
    with open('cell.inc') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('A'):
                lattice_vector.append(line.split()[1:4])
            if line.startswith('B'):
                lattice_vector.append(line.split()[1:4])
            if line.startswith('C'):
                lattice_vector.append(line.split()[1:4])
except FileNotFoundError as e:
    print('WARNING: cell.inc file not found! If your system is periodic, lattice parameters are missing. Please bring here the file cell.inc')
n_lattice_vectors = len(lattice_vector)

with open(init_geom) as f:
    # The first line of the xyz file is the number of atoms
    n_atoms = int(f.readline())

# forces = []

print(f'Parsing {init_geom} ...', end="", flush=True)

# Parse the output file
with open(init_geom) as f:
    lines = f.readlines()
    nlines = len(lines)
    for i, line in enumerate(lines):
        if 'Number of atoms' in line:
            n_atoms = int(line.split()[5])
        if 'Number of lattice vectors' in line:
            n_lattice_vectors = int(line.split()[6])
        if 'lattice_vector' in line:
            lattice_vector.append(line)
        if '  atom  ' in line:
            atoms.append(line)
        if 'Total atomic forces' in line:
            forces.append(lines[i+1:i+1+n_atoms])

print("OK!")

# Work out the interval
steps = len(forces)
if not final_step:
    final_step = steps
elif final_step > steps:
    final_steps = steps
step = initial_step


print(f'Writing {outputfile} folders...', end=' ', flush=True)
while step < final_step:
    # dirname = outfile + '_step_' + str(step) + '_ID_' + str(calc_ID)
    dirname = outfile + '_step_' + str(step)
    os.makedirs(dirname, exist_ok=True)
    if do_cp_control:
        shutil.copy('control.in', dirname)
    os.chdir(dirname)
    fout = open('geometry.in', 'w')
    if n_lattice_vectors > 0:
        # Work out the lattice parameters
        idx_lattice_vector = step * n_lattice_vectors
        for vec in range(n_lattice_vectors):
            fout.write(lattice_vector[idx_lattice_vector + vec])

    # Write atoms with their forces
    idx_atoms = step * n_atoms
    for atom in range(n_atoms):
        iatom = atoms[idx_atoms + atom]
        fout.write(iatom)
    step += step_interval
    calc_ID += 1
    fout.close()
    os.chdir('../')
print("OK!")
print("Done!")
