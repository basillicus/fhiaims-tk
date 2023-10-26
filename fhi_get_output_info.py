#!/usr/bin/env python3

import os

"""
Parse an aims outfile and extract requested information.
Prints a table with the requested information. By default gets the total energy.
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

output = ''
# Find file to be parsed
def find_all(name, path):
    result = []
    for root, dirs, files in os.walk(path):
        if name in files:
            result.append(os.path.join(root, name))
    return result


files = find_all(infile, './')

of = open(outfile, 'a')

of.write('# file  \t\tTotal E \t n_atoms \t K Points \t  SCF iter \t Total Time \t Polarizability\n')
for parsing_file in files:
    lattice_vector = []
    atoms = []
    forces = []
    kpoints = None
    total_time = None
    total_scf_iterations = None
    polarizability = None
    n_atoms = None
    # Parse the output file
    with open(parsing_file) as f:
        lines = f.readlines()
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
            if 'k_grid      ' in line:
                kpoints = line.split()[-3:]
            if '| Total time     ' in line:
                total_time = line.split()[-2:-1]
            if 'Self-consistency cycle converged' in line:
                total_scf_iterations = lines[i+4].split()[4]
            if '| Polarizability' in line:
                polarizability = line.split()[2:]
            if '| Total energy of the ' in line:
                total_energy = line.split()[11]

        line = [str(parsing_file), '\t', str(total_energy), '\t', str(n_atoms), '\t', str(kpoints), '\t',
                str(total_scf_iterations), '\t', str(total_time), '\t', str(polarizability), '\n']

# TODO: Print it according to the calling options
        of.write(' '.join(line))

# for line in output:
#     print(output)
#
# # Create the extxyz outputfile
# fout = open(outfile, 'w')
# steps = len(forces)
# step = 0
#
# #               step      atom     force components
# # print(forces[steps-1][n_atoms-1].split()[2:5])
#
# while step < steps:
#     fout.write(f'{n_atoms}\n')
#
#     # Work out the comment line of the extxyz file
#     idx_lattice_vector = step * n_lattice_vectors
#     tmp = []
#     for vec in range(n_lattice_vectors):
#         tmp.append(lattice_vector[idx_lattice_vector + vec].split()[1:4])
#
#     ilattice = ''
#     for lat in tmp:
#         ilattice += ' '.join(lat)
#         ilattice += ' '
#
#     # Write atoms with their forces
#     idx_atoms = step * n_atoms
#     for atom in range(n_atoms):
#         # print(comment.join(lattice_vector[idx_lattice_vector + vec].split()[1:]))
#         iatom = atoms[idx_atoms + atom].split()[1:4]
#         ispecies = str(atoms[idx_atoms + atom].split()[4])
#         iforces = forces[step][atom].split()[2:5]
#         fout.write(ispecies + '\t' + '\t'.join(iatom) + '\t' + '\t'.join(iforces) + '\n')
#     step += 1
#
# print('Geometry optimization steps: ', steps)
