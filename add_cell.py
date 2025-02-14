#!/usr/bin/env python3
import numpy as np

from ase.io import read
from ase.io.extxyz import write_extxyz

import argparse

from config import extxyz_outputfile

parser = argparse.ArgumentParser(
    prog='add_cell.py',
    description='Add  a new cell to an .extxyz file which is no periodic. Useful for NEP potentials, that require a lattice')

inputfile = 'train.xyz'
outputfile = extxyz_outputfile
energy_key = 'REF_energy'
force_key = 'REF_forces'
virial_key = 'REF_virial'
dipole_key = 'REF_mu'
polarizability_key = 'REF_alpha'

parser.add_argument('-i', '--inputfile', default=inputfile,
                    help='input file, an extxyz file to be parsed')

parser.add_argument('-c', '--add_cell', default=[], nargs='+', type=float, required=True,
                    help='Add a lattice cell when the system is not periodic. 9 componets: [x1 x2 .... z2 z3]')

parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='outputfile file where the geoemtries are writen in an .extxyz file')

parser.add_argument('--pbc', default=False, action='store_true',
                    help='Add this flag if the system is periodic')

parser.add_argument('--energy_key', default=energy_key, help='Energy keyword for the extyz file')
parser.add_argument('--forces_key', default=force_key, help='Force keyword for the extyz file')
parser.add_argument('--virial_key', default=virial_key, help='Virial keyword for the extyz file')
parser.add_argument('--dipole_key', default=dipole_key, help='Dipole keyword for the extyz file')
parser.add_argument('--polarizability_key', default=polarizability_key, help='Polarizability keyword for the extyz file')
#
# parser.add_argument('-a', '--all', action='store_true', default=True,
#                     help='all steps will be read/writen')


args = parser.parse_args()

inputfile = args.inputfile
outputfile = args.outputfile

is_periodic = args.pbc
add_cell = args.add_cell
energy_key = args.energy_key
forces_key = args.forces_key
virial_key = args.virial_key
dipole_key = args.dipole_key
polarizability_key = args.polarizability_key

print(f'Reading {inputfile} file...', flush=True, end=' ')
geometries = read(inputfile, index=':')
print('Done!')

print(f'Adding cell: {add_cell} to geometries...', flush=True, end= ' ')
list_geoms = []
for geom in geometries:
    geom.set_pbc(True)
    geom.set_cell(np.array(add_cell).reshape(3,3))
    list_geoms.append(geom)
print('Done!')

print(f'Writing {outputfile} file...', flush=True, end=' ')
# write(list_geoms, outputfile)
with open(outputfile, "w") as f:
    write_extxyz(f, list_geoms)
print('Done!')
