#!python

"""
Convert geometry files from and to fhi-aims format
"""
from ase.io import read, write
from ase.io.aims import read_aims, write_aims

import argparse
from config import visualization_outputfile



def read_cp2k(infile, cellfile):
    cell = []
    atoms = []
    with open(cellfile) as f:
        c = f.readlines()
    with open(infile) as f:
         a = f.readlines()
    for i in c:
        cell.append(i[1:4])
    for i in a:
        atoms.append(i[1:4])

parser = argparse.ArgumentParser(
    prog='fhi_convert_cp2k2fhi.py',
    description='Interconvert geometry files',
)


inputfile = 'coordinates.inc'
cellfile = 'cell.inc'
outputfile = visualization_outputfile

parser.add_argument('-i', '--inputfile', default=inputfile, help='Input geometry file')
parser.add_argument('-c', '--inputfile', default=cellfile, help='Input cell file')
parser.add_argument('-o', '--outputfile', default=None, help='output geometry file. Extension sensitive')
parser.add_argument('-f', '--format', default=None, help='output geometry file format. Outpuf file  extension will be ignored')


args = parser.parse_args()
infile = args.inputfile
outfile = args.outputfile
format = args.format

if infile == "coordinates.inc":
    geom = read_cp2k(infile, cellfile)
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
