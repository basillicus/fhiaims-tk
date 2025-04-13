#!/usr/bin/env python3

import os

"""
Convert an MD simulation output file (AIMS style) to an extxyz file.

The script reads an AIMS output file containing MD data,
and writes an extended XYZ file for visualization or further processing.
It optionally centers the geometries if the corresponding flag is set.

Command-line arguments:
    - -i, --inputfile: Input file name (default from config, e.g. fhi_aims_outputfile).
    - -o, --outputfile: Output filename for the extxyz file (default from config). Outpts any file format  suported by ASE
    - -a, --all: Flag to read and write all steps (default True).
    - -c, --center: Flag to disable centering of atoms (default False).
"""

import argparse
from config import fhi_aims_outputfile, extxyz_outputfile

from ase.io import read, write
# Note on argparse: The '-c/--center' flag is a bit counter-intuitive: 
# passing -c disables centering. Consider renaming it (perhaps --no-center) for clarity.
parser = argparse.ArgumentParser(
    prog='fhi_convert_MD2extxyz.py',
    description='Converts MD files to extxyz (or others if supported by ASE)',
)

inputfile = fhi_aims_outputfile
outputfile = extxyz_outputfile

parser.add_argument('-i', '--inputfile', default=inputfile,
                    help='input file, an AIMS ouptut file to be parsed')
parser.add_argument('-o', '--outputfile', default=outputfile,
                    help='outputfile file where the geoemtries are writen in an .extxyz file')
parser.add_argument('-a', '--all', action='store_true', default=True,
                    help='all steps will be read/writen')
parser.add_argument('-c', '--center', action='store_false', default=True,
                    help='Controls centering atoms within the unit cell. If -c is given, atoms will NOT be centered')

args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile
do_center = args.center
all = args.all

print('Input file: ', infile)
print('Output file: ', outfile)

if all:
    print('Reading input file. This may take a while...')
    md = read(infile, index=':')
    if do_center:
        for i in range(len(md)):
            md[i].center()
    write(outfile, md)
