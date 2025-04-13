#!/usr/bin/env python

import sys
import numpy as np
from ase import Atoms
from ase.io import write

# Load your structured NumPy array.
# This assumes you saved the array with np.save. Adjust if using a different format.
data = np.load(sys.argv[1])

atoms_list = []

# Loop over each MD frame in the structured array.
for frame in data:
    # The 'species' field is stored as fixed-length bytes (e.g. S2). Decode each entry.
    # If your MD always has 32 atoms, this produces a list of 32 element symbols.
    species = [s.decode('utf-8').strip() for s in frame['species']]

    # Extract the lattice vectors (cell) and atomic coordinates.
    cell = frame['lattice_vector']
    positions = frame['coordinates']

    # Optionally, store extra properties in the info dictionary.
    # Here we add the energy and polarizability for this frame.
    info = {
        'energy': frame['energy'],
        'polarizability': frame['polarizability']
    }

    # Create the ASE Atoms object.
    # pbc=True enables periodic boundary conditions (adjust if needed).
    atoms = Atoms(symbols=species, positions=positions, cell=cell, pbc=True, info=info)
    atoms_list.append(atoms)

# Write all frames to an extended XYZ file.
# The extxyz format can store extra information from the info dictionary.
write("trajectory.extxyz", atoms_list, format="extxyz")
print("Trajectory successfully written to trajectory.extxyz")
