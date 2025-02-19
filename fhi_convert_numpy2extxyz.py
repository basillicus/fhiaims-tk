#!/usr/bin/env python

import sys
import numpy as np
from ase import Atoms
from ase.io import write

# Load the structured NumPy array
data = np.load(sys.argv[1])

atoms_list = []

# Loop over each MD frame
for frame in data:
    # Extract atomic species (decode if present)
    species = (
        [s.decode("utf-8").strip() for s in frame["species"]]
        if "species" in frame.dtype.names
        else None
    )

    # Extract atomic positions
    positions = frame["coordinates"] if "coordinates" in frame.dtype.names else None

    # Extract lattice vectors
    cell = frame["lattice_vector"] if "lattice_vector" in frame.dtype.names else None

    # Build the ASE Atoms object
    atoms = Atoms(symbols=species, positions=positions, cell=cell, pbc=True)

    # Collect extra per-frame properties (e.g., energy, dipole, etc.)
    info = {}
    for field in frame.dtype.names:
        if field not in {"species", "coordinates", "lattice_vector"}:  # Skip core fields
            info[field] = frame[field]

    # Attach additional properties
    atoms.info.update(info)
    atoms_list.append(atoms)

# Write all frames to an extended XYZ file
write("trajectory.extxyz", atoms_list, format="extxyz")
print("Trajectory successfully written to trajectory.extxyz")

