#!/usr/bin/env python3

""" Generates individual geometries.in files, one from each MD step, from an aims.out file containing a MD 
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

inputfile = fhi_aims_outputfile
outputfile = 'polar'

parser.add_argument('-i', '--inputfile', default=inputfile, help='This is the input file for the script, an AIMS ouptut file to be parsed')
parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='This is the prefix for the outputfile folder, where the geometries will be writen in aims format')
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

if do_cp_control:
    if not os.path.exists('control.in'):
        print('WARNING: control.in not found. Will not be copied to folders!')
        do_cp_control = False

# Per step variables initialization
lattice_vector = []
n_lattice_vectors = 0
atoms = []
forces = []
isteps = []

print(f'Parsing {init_geom} ...', end="", flush=True)

# Parse the output file
with open(init_geom) as f:
    lines = f.readlines()
    for i, line in enumerate(lines):
        if 'Time step number' in line:
            isteps.append(int(line.split()[5]))
        if 'Number of atoms' in line:
            n_atoms = int(line.split()[5])
        if 'Number of lattice vectors' in line:
            n_lattice_vectors = int(line.split()[6])
        if 'lattice_vector' in line:
            lattice_vector.append(line)
        if '   atom   ' in line:
            atoms.append(line)
        if 'Total atomic forces' in line:
            forces.append(lines[i+1:i+1+n_atoms])

print("OK!")

# Work out the interval
nsteps = len(isteps)
first_step = isteps[0]
last_step = isteps[-1]
if not initial_step:
    initial_step = first_step
elif initial_step < first_step:
    print('WARNING. Initial step is smaller than available first step')
    print('setting initial step to ', first_step)
    initial_step = first_step
elif initial_step > last_step:
    print('ERROR. Initial step is larger than available last step')
    print('Reduce inital step or extend your MD')
    sys.exit()
elif initial_step > last_step:
    print('ERROR. inital step is larger than available last step')
    print('reduce inital step ')
    sys.exit()

if not final_step:
    final_step = last_step
elif final_step > last_step:
    print('WARNING. Final step is larger than available last step')
    print('setting final step to ', last_step)
    final_step = last_step
elif final_step < first_step:
    print('ERROR. final step is smaller than available first step')
    print('Increase final step ')
    sys.exit()
if final_step < initial_step:
    print('ERROR. final step is smaller than initial step. Invert values?')
    sys.exit()


print(f'Writing {outputfile} folders, from {initial_step} to {final_step}...', end=' ', flush=True)
n = 0
while n < nsteps:
    # dirname = outfile + '_step_' + str(step) + '_ID_' + str(calc_ID)
    step = isteps[n]
    update_step = 1
    if step >= initial_step and step <= final_step:
        update_step = step_interval
        dirname = outfile + '_step_' + str(step)
        os.makedirs(dirname, exist_ok=True)
        if do_cp_control:
            shutil.copy('control.in', dirname)
        os.chdir(dirname)
        fout = open('geometry.in', 'w')
        if n_lattice_vectors > 0:
            # Work out the lattice parameters
            # NOTE: I leave this here just in case in the future the lattice vectors
            # change with the dynamics and need to be updated at each step
            # idx_lattice_vector = n * n_lattice_vectors
            idx_lattice_vector = 0 * n_lattice_vectors
            for vec in range(n_lattice_vectors):
                fout.write(lattice_vector[idx_lattice_vector + vec])

        # Write atoms with their forces
        idx_atoms = n * n_atoms
        for atom in range(n_atoms):
            iatom = atoms[idx_atoms + atom]
            fout.write(iatom)
        calc_ID += 1
        fout.close()
        os.chdir('../')
    n += update_step
print("OK!")
print("Done!")
