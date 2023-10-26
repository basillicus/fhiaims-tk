#!/usr/bin/env python3

"""
Convert geometry files from and to fhi-aims format
"""
from ase.io import read, write
from ase.io.aims import read_aims, write_aims

import argparse
from config import visualization_outputfile

parser = argparse.ArgumentParser(
    prog='fhi_convert_xyz2fhi.py',
    description='Interconvert geometry files',
)

inputfile = 'geometry.in'
outputfile = visualization_outputfile

parser.add_argument('-i', '--inputfile', default=inputfile, help='Input geometry file')
parser.add_argument('-o', '--outputfile', default=None, help='output geometry file. Extension sensitive')
parser.add_argument('-f', '--format', default=None, help='output geometry file format. Outpuf file  extension will be ignored')


args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile
format = args.format

if infile == "geometry.in":
    geom = read_aims(infile)
else:
    geom = read(infile)

if outfile != 'geometry.in':
    if outfile == None:
        outfile = outputfile
    write(outfile, geom)
elif outfile == 'geometry.in' or format == 'aims':
    if not outfile:
        outfile = 'geometry.in'
    write_aims(outfile, geom)

# if not format:
