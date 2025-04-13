#!/usr/bin/env python3
"""
Convert geometry files between CP2K and FHI-AIMS formats.

This script reads geometry and cell information from CP2K-style input files (coordinates.inc and cell.inc)
and writes out a geometry file in the FHI-AIMS format. The output file format may be inferred from its
extension, supporting both standard and AIMS-specific formats.

Command-line arguments:
    - -i, --inputfile: Input geometry file (default "coordinates.inc").
    - -c, --inputfile: Input cell file (default "cell.inc").  <-- Note: there’s a naming conflict here!
    - -o, --outputfile: Output filename; if not provided, a default from config is used.
    - -f, --format: Format specifier for the output file.
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

# Note on argparse: The cell file argument re-uses "--inputfile" with -c. This is a problem;
# each argument should have a unique destination. Consider renaming one of them.
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
