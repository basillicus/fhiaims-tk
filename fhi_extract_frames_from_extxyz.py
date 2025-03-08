#!/usr/bin/env python3

"""Generates individual geometries.in files, one from each frame in the input file.
Input file could be any file format supported by ASE
Returns one folder per geometry where the polarazibility will be calculated
"""

import argparse
import os
import shutil
import sys

from ase.io import read, write

from config import aims_species_directory, fhi_aims_outputfile, visualization_outputfile

parser = argparse.ArgumentParser(
    prog="fhi_prepare_polarizabilities_from_MD.py",
    description="Prepares polarizability calculations from a MD",
)

inputfile = fhi_aims_outputfile
outputfile = "geometry.in"
prefix = "frame"

parser.add_argument(
    "-i",
    "--inputfile",
    default=inputfile,
    help="This is the input file for the script, an AIMS ouptut file to be parsed",
)
parser.add_argument(
    "-o",
    "--outputfile",
    default=outputfile,
    help="This is the the outputfile to write in each folder, where the geometries will be writen in the requested format",
)
parser.add_argument(
    "-p",
    "--prefix",
    default=prefix,
    help="This is the prefix for the outputfile folder, where the geometries will be writen in aims format",
)
parser.add_argument(
    "-g",
    "--gencontrol",
    action="store_true",
    default=False,
    help="It will generate a control.in file for AIMS calculation",
)
parser.add_argument(
    "-c",
    "--copycontrol",
    action="store_true",
    default=False,
    help="If control.in file is in the folder, it will be copied to the folders",
)
parser.add_argument(
    "-s",
    "--steps",
    default=1,
    type=int,
    help="Every how many steps the polarizability will be calculated",
)
parser.add_argument(
    "-n", "--initialstep", default=0, type=int, help="start form this step"
)
parser.add_argument(
    "-f", "--finalstep", default=None, type=int, help="end at this step"
)
parser.add_argument(
    "-d", "--idstart", default=0, type=int, help="starting number for calculation ID"
)

# Global variables initialization
args = parser.parse_args()
init_geom = args.inputfile
outfile = args.outputfile
step_interval = args.steps
do_cp_control = args.copycontrol
initial_step = args.initialstep
final_step = args.finalstep
calc_ID = args.idstart
gencontrol = args.gencontrol

# Per step variables initialization
lattice_vector = []
n_lattice_vectors = 0
atoms = []
forces = []
isteps = []

print(f"Parsing {init_geom} ...", end="", flush=True)

# Parse the output file
geom = read(init_geom, index=":")
print("OK!")


def gen_control_file(geom):
    # from ase.io import write
    from ase.calculators.aims import Aims

    print("species_dir", aims_species_directory)
    input_parameters = {
        "species_dir": aims_species_directory,
        "xc": "pw-lda",
        "occupation_type": "gaussian 0.01",
        "sc_accuracy_etot": 1e-4,
        "sc_accuracy_eev": 1e-3,
        "compute_forces": True,
    }

    if geom.pbc.any():
        input_parameters.update({
            "k_grid": [1, 1, 8],
            "DFPT": "dielectric",
        })
    else:
        input_parameters.update({
            "DFPT": "polarizability",
            "output": ["dipole"],
        })

    calculator = Aims(**input_parameters)
    print(calculator.parameters)
    geom.calc = calculator
    geom.calc.write_input(geom)


if gencontrol:
    gen_control_file(geom[0])

if do_cp_control:
    if not os.path.exists("control.in"):
        print("WARNING: control.in not found. Will not be copied to folders!")
        print(
            "You can use -g --gencontrol option to generate a reference control.in file"
        )
        do_cp_control = False

# Work out the interval
nsteps = len(geom)
first_step = 0
last_step = len(geom)
if not initial_step:
    initial_step = first_step
elif initial_step < first_step:
    print("WARNING. Initial step is smaller than available first step")
    print("setting initial step to ", first_step)
    initial_step = first_step
elif initial_step > last_step:
    print("ERROR. Initial step is larger than available last step")
    print("Reduce inital step or extend your MD")
    sys.exit()
elif initial_step > last_step:
    print("ERROR. inital step is larger than available last step")
    print("reduce inital step ")
    sys.exit()

if not final_step:
    final_step = last_step
elif final_step > last_step:
    print("WARNING. Final step is larger than available last step")
    print("setting final step to ", last_step)
    final_step = last_step
elif final_step < first_step:
    print("ERROR. final step is smaller than available first step")
    print("Increase final step ")
    sys.exit()
if final_step < initial_step:
    print("ERROR. final step is smaller than initial step. Invert values?")
    sys.exit()


print(
    f"Writing {outputfile} folders, from {initial_step} to {final_step}...",
    end=" ",
    flush=True,
)
n = 0
while n < nsteps:
    # dirname = outfile + '_step_' + str(step) + '_ID_' + str(calc_ID)
    step = n
    update_step = 1
    if step >= initial_step and step <= final_step:
        update_step = step_interval
        dirname = prefix + "_step_" + str(step)
        os.makedirs(dirname, exist_ok=True)
        if do_cp_control:
            shutil.copy("control.in", dirname)
        os.chdir(dirname)
        fout = open("geometry.in", "w")
        write(outputfile, geom[step])
        os.chdir("../")
    n += update_step
print("OK!")
print("Done!")
