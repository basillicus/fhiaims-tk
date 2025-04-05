#!/usr/bin/env python3

import os
import sys

import argparse
from ase.io import read, write
from ase.io.aims import write_aims
import numpy as np

"""
Reads the aims.out files in the prefix folders and generates a train.xyz file
"""

# Read the arguments
parser = argparse.ArgumentParser(
    prog="fhi_convert_aims2train.py",
    description="Reads the aims.out files in the prefix folders and generates a train.xyz file",
)

inputfile = "aims.out"
outputfile = "train.xyz"
# outputfile_extxyz = 'sampled.extxyz'

parser.add_argument(
    "-i",
    "--inputfile",
    default=inputfile,
    help="input file: aims.out . Read in the aims outputfile to extract the forces (and maybe polarizabilities)",
)
parser.add_argument(
    "-o", "--outputfile", default=outputfile, help="output file: train.xyz file."
)
parser.add_argument(
    "-p", "--prefix", default="md_sample_", help="Prefix of the folders. [md_sample_]"
)
parser.add_argument(
    "-l", "--list", default=None, help="File with  a list of files to be read]"
)
parser.add_argument(
    "-n",
    "--samples",
    default=0,
    help="Total number of geometries to sample. If == 0, all will be included [0]",
)

args = parser.parse_args()
inputfile = args.inputfile
outfile = args.outputfile
num_samples = int(args.samples)
prefix_folder = args.prefix
list_of_files = args.list

if list_of_files:
    with open(list_of_files, "r") as f:
        files = f.readlines()
        print(
            f"Reading {inputfile} files in path given in {list_of_files} ...",
            end="",
            flush=True,
        )
    subdirs = [os.path.dirname(file) for file in files]
else:
    print(
        f"Reading {inputfile} files in folders in  current folder with prefix {prefix_folder} ...",
        end="",
        flush=True,
    )
    # Get a list of subdirectories (folders) within the current directory
    subdirs = [
        d for d in os.listdir(".") if os.path.isdir(d) and d.startswith(prefix_folder)
    ]

    print("Done!")


indices = range(len(subdirs))
if num_samples > 0:
    indices = np.random.choice(len(subdirs), num_samples, replace=False)

final_samples = len(indices)
print(
    f"Writing {outputfile} with a total of {final_samples} samples ...",
    end="",
    flush=True,
)
# Loop over each subdir and read a file named "example.txt" inside it
for i in indices:
    path = os.path.join(subdirs[i], inputfile)
    if os.path.exists(path):
        iaims = read(path)
        write(outputfile, iaims, append=True)
print("Done!")
