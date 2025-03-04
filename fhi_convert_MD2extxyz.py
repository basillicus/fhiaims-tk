#!/usr/bin/env python3

import os

"""
Reads a MD file (will try to infer the format) and writes an .extxyz file. Different outputfiles can be given if supported by ASE
"""

import argparse
from config import fhi_aims_outputfile, extxyz_outputfile

from ase.io import read, write

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
parser.add_argument('-c', '--center', action='store_false', default=False,
                    help='center atoms within the unit cell')

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
