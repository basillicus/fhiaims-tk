#!python

import argparse
import os
import sys

from config import fhi_basis_set_path

# TODO: Work out the tier

parser = argparse.ArgumentParser(
    prog='fhi_add_basisSet',
    description='Add basis sets to the control.in file',
    # epilog="""The edited file will be printed out. If you want it to overwrite your current control.in execute the script as:\n

    #     fhi_add_basisSet -s X Y Z > control.in
    # """
)

personal_path = fhi_basis_set_path

parser.add_argument('-p', '--path', default=personal_path,
                    help='set basis set path. Change the personal_path variable in the script so you do not have to set it everytime')
parser.add_argument('-l', '--level', default='light', choices=['light', 'intermediate', 'tight', 'very-tight'], help='set basis set level')
parser.add_argument('-t', '--tier', default=1, choices=[1, 2, 3, 4], help='set basis set tier level (Not yet implemented. Do it manually)')
parser.add_argument('-s', '--species', default=None, nargs='*',
                    help='set species to be added. If none given, will guess it from the geometry.in file')
parser.add_argument('-f', '--force-overwrite', default=True, help='if True, will overwrite the current control.in (Not yet implemented)')

args = parser.parse_args()

basis_set_path = args.path
basis_set_level = args.level
basis_set_tier = args.tier
species = args.species

# Work out required species if not given
if not species:
    species = []
    try:
        with open('geometry.in', 'r') as f:
            geoemtry = f.readlines()
    except FileNotFoundError:
        print("""!ERROR: Trying to guess the species you need, but the file geometry.in is not here.
        Try to call the script wiht -s and give the species you need""")
        sys.exit()
    for line in geoemtry:
        if 'atom' in line:
            add_atom = line.split()[-1]
            if not add_atom in species:
                species.append(add_atom)

# Load control.in
control = open('control.in', 'r')

# First Remove Existing Basis sets
lines = control.readlines()

for i, line in enumerate(lines):
    if '##########' in line:
        remove_from = i
        break
lines = lines[0:i]

if not basis_set_path.endswith('/'):
    basis_set_path += '/'
work_path = os.path.join(basis_set_path, basis_set_level)


warnings = []
# Get the file for each species
basis_set_files = os.listdir(work_path)
for s in species:
    pattern = f'_{s}_'
    found = False

    for bsf in basis_set_files:
        if pattern in bsf:
            found = True
            add_file = os.path.join(work_path, bsf)
            with open(add_file, 'r') as f:
                lines.extend(f.readlines())
            break

    # If you made it all the way down here is because the basis set is not where it is supoused to be
    if not found:
        warnings.append(f'## WARNING! Basis set for {s} not found for level {basis_set_level}!')

for line in lines:
    print(line[:-1])

for w in warnings:
    print(w)
