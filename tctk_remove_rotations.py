#!/usr/bin/env python3
"""
Remove rotational motion from an MD trajectory while also handling forces.

This script reads a trajectory file (supported by ASE), computes the optimal
rotation for each frame using quaternion minimization (from a common utility),
and rotates forces (if available) along with positions, momenta, etc.
It also offers a utility function to write out a custom XYZ file that includes forces.

There is no explicit command-line interface here apart from providing input/output filenames inside the code.
"""


import numpy as np
from ase import io
from scipy.spatial.transform import Rotation
from ase.atoms import Atoms
import os


def get_optimal_rotation(P, Q):
    """
    Calculate the optimal rotation to align P to Q using quaternions.

    Parameters:
    P (numpy.ndarray): Coordinates to rotate (n_atoms, 3)
    Q (numpy.ndarray): Target coordinates (n_atoms, 3)

    Returns:
    numpy.ndarray: Rotation matrix (3, 3)
    """
    # Initialize correlation matrix
    H = np.zeros((3, 3))

    # Calculate correlation matrix
    for p, q in zip(P, Q):
        H += np.outer(p, q)

    # Create the key matrix
    K = np.array(
        [
            [
                H[0, 0] + H[1, 1] + H[2, 2],
                H[1, 2] - H[2, 1],
                H[2, 0] - H[0, 2],
                H[0, 1] - H[1, 0],
            ],
            [
                H[1, 2] - H[2, 1],
                H[0, 0] - H[1, 1] - H[2, 2],
                H[0, 1] + H[1, 0],
                H[2, 0] + H[0, 2],
            ],
            [
                H[2, 0] - H[0, 2],
                H[0, 1] + H[1, 0],
                H[1, 1] - H[0, 0] - H[2, 2],
                H[1, 2] + H[2, 1],
            ],
            [
                H[0, 1] - H[1, 0],
                H[2, 0] + H[0, 2],
                H[1, 2] + H[2, 1],
                H[2, 2] - H[0, 0] - H[1, 1],
            ],
        ]
    )

    # Get eigenvector corresponding to maximum eigenvalue
    eigenvalues, eigenvectors = np.linalg.eigh(K)
    quaternion = eigenvectors[:, np.argmax(eigenvalues)]

    # Normalize quaternion and create rotation object
    try:
        quat = np.array([quaternion[0], quaternion[1], quaternion[2], quaternion[3]])
        quat_norm = np.linalg.norm(quat)

        if quat_norm < 1e-10:  # Check for zero norm
            return np.eye(3)  # Return identity matrix if no rotation needed

        quat_normalized = quat / quat_norm
        rot = Rotation.from_quat(
            [
                quat_normalized[1],
                quat_normalized[2],
                quat_normalized[3],
                quat_normalized[0],
            ]
        )
        return rot.as_matrix()
    except:
        print("Warning: Failed to compute rotation matrix, returning identity")
        return np.eye(3)


def check_and_extract_forces(atoms_list):
    """
    Check if forces are available in the trajectory.
    ASE doesn't store forces directly, so we need to extract them.

    Parameters:
    atoms_list (list): List of ASE Atoms objects

    Returns:
    list or None: List of force arrays if available, None otherwise
    """
    # First, check if forces are available in any special arrays
    for atoms in atoms_list:
        if hasattr(atoms, "arrays") and "forces" in atoms.arrays:
            return [atoms.arrays["forces"] for atoms in atoms_list]

    # Check if we can extract forces from the calculator
    try:
        forces = []
        for atoms in atoms_list:
            if atoms.calc is not None:
                forces.append(atoms.get_forces())
            else:
                # If no calculator, try to get forces from info dict
                if "forces" in atoms.info:
                    forces.append(atoms.info["forces"])
                else:
                    return None  # Forces not available
        return forces
    except:
        return None


def remove_rotations(atoms_list):
    """
    Remove rotational motion by aligning all frames to the first frame
    using quaternion-based optimal rotation.

    Parameters:
    atoms_list (list): List of ASE Atoms objects

    Returns:
    list: Aligned Atoms objects
    """
    if not atoms_list:
        return []

    # Get reference structure (first frame)
    ref_atoms = atoms_list[0]
    ref_com = ref_atoms.get_center_of_mass()
    ref_pos = ref_atoms.get_positions() - ref_com

    # Extract forces if available
    forces_list = check_and_extract_forces(atoms_list)
    has_forces = forces_list is not None

    if has_forces:
        print(f"Forces detected in trajectory. Will rotate them accordingly.")
    else:
        print(f"No forces found in trajectory.")

    aligned_atoms_list = []

    for i, atoms in enumerate(atoms_list):
        new_atoms = atoms.copy()

        # Center at COM
        com = new_atoms.get_center_of_mass()
        pos = new_atoms.get_positions() - com

        # Get optimal rotation matrix
        R = get_optimal_rotation(pos, ref_pos)

        # Apply rotation to positions
        rotated_pos = (pos @ R.T) + com
        new_atoms.set_positions(rotated_pos)

        # If momenta exist, rotate them too
        if new_atoms.has("momenta"):
            momenta = new_atoms.get_momenta()
            rotated_momenta = momenta @ R.T
            new_atoms.set_momenta(rotated_momenta)

        # Handle forces if they exist
        if has_forces:
            # Rotate forces
            rotated_forces = forces_list[i] @ R.T

            # # Store rotated forces in the atoms object's info dictionary
            # if "info" not in new_atoms.__dict__:
            #     new_atoms.info = {}
            # new_atoms.info["forces"] = rotated_forces

            # If arrays exists, also store there for compatibility
            if hasattr(new_atoms, "arrays"):
                new_atoms.arrays["forces"] = rotated_forces

        if i % 100 == 0:  # Progress indicator
            print(f"Processed frame {i}")

        aligned_atoms_list.append(new_atoms)

    return aligned_atoms_list


def write_xyz_with_forces(filename, atoms_list):
    """
    Write a custom XYZ file that includes forces.

    Parameters:
    filename (str): Output filename
    atoms_list (list): List of ASE Atoms objects with forces stored in info
    """
    with open(filename, "w") as f:
        for atoms in atoms_list:
            n_atoms = len(atoms)
            f.write(f"{n_atoms}\n")
            f.write(
                'Lattice="0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0" Properties=species:S:1:pos:R:3:forces:R:3\n'
            )

            positions = atoms.get_positions()
            symbols = atoms.get_chemical_symbols()

            # Get forces from info dict
            if "forces" in atoms.info:
                forces = atoms.info["forces"]
            elif hasattr(atoms, "arrays") and "forces" in atoms.arrays:
                forces = atoms.arrays["forces"]
            else:
                # If no forces, use zeros
                forces = np.zeros_like(positions)

            for i in range(n_atoms):
                f.write(
                    f"{symbols[i]} {positions[i, 0]:.8f} {positions[i, 1]:.8f} {positions[i, 2]:.8f} "
                    f"{forces[i, 0]:.8f} {forces[i, 1]:.8f} {forces[i, 2]:.8f}\n"
                )


def main():
    try:
        input_file = "trajectory.xyz"
        output_file = "trajectory_aligned.xyz"

        # Check if command line arguments were provided
        import sys

        if len(sys.argv) > 1:
            input_file = sys.argv[1]
        if len(sys.argv) > 2:
            output_file = sys.argv[2]

        # Read trajectory
        print(f"Reading trajectory from '{input_file}'...")
        traj = io.read(input_file, index=":")

        if not traj:
            print("Error: No frames read from trajectory")
            return

        print(f"Successfully read {len(traj)} frames")

        # Remove rotations
        print("Removing rotations...")
        aligned_traj = remove_rotations(traj)

        # Write aligned trajectory
        print(f"Writing aligned trajectory to '{output_file}'...")

        # # Check if we need to preserve forces
        # if (hasattr(aligned_traj[0], "info") and "forces" in aligned_traj[0].info) or (
        #     hasattr(aligned_traj[0], "arrays") and "forces" in aligned_traj[0].arrays
        # ):
        #     write_xyz_with_forces(output_file, aligned_traj)
        # else:
        #     io.write(output_file, aligned_traj)

        io.write(output_file, aligned_traj)

        print("Done!")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        import traceback


# Note on argparse: This script does not use argparse for CLI input – consider adding one
        traceback.print_exc()


if __name__ == "__main__":
    main()
