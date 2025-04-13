#!/usr/bin/env python
"""
Align an MD trajectory by removing rotations via a quaternion-based optimal rotation.

This script reads a trajectory file (supported by ASE),
selects a reference frame (default first frame), and aligns all frames
by applying an optimal rotation (computed with quaternions) so that the rotational
motion is removed. Optionally, momenta and forces are rotated along.

Command-line arguments:
    - -i, --input: Input trajectory file (default "trajectory.xyz").
    - -o, --output: Output aligned trajectory file (default "traj_aligned.xyz").
    - --ref: Index of the reference frame to align to (default 0).

"""

import numpy as np
import argparse
from ase import io
from scipy.spatial.transform import Rotation
from ase.atoms import Atoms


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


def remove_rotations(atoms_list, ref=0):
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
    ref_atoms = atoms_list[ref]
    ref_com = ref_atoms.get_center_of_mass()
    ref_pos = ref_atoms.get_positions() - ref_com

    # ref_forces = ref_atoms.get_forces()
    # print("Forces:", ref_forces)

    aligned_atoms_list = []

    for i, atoms in enumerate(atoms_list):
        new_atoms = atoms.copy()

        # print("ATOMS:", atoms)
        # print("ATOMS LIST [i]:", atoms_list[i])
        #
        # print("atoms_list[i]", atoms_list[i].get_forces())
        # print("new atoms", new_atoms.get_forces())
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
            momenta = atoms_list[i].get_momenta()
            rotated_momenta = momenta @ R.T
            new_atoms.set_momenta(rotated_momenta)

        # If forces exist, rotate them too
        if new_atoms.has("forces"):
            forces = atoms_list[i].get_forces()
            rotated_forces = forces @ R.T
            new_atoms.set_forces(rotated_forces)

        if i % 100 == 0:  # Progress indicator
            print(f"Processed frame {i}")

        aligned_atoms_list.append(new_atoms)

    return aligned_atoms_list


def main():

# Note on argparse: Could add more specifics about accepted file formats.
    parser = argparse.ArgumentParser(
        description="Align an MD trajectory to remove rotations using quaternion-based alignment."
    )
    parser.add_argument(
        "-i",
        "--input",
        default="trajectory.xyz",
        help="Input trajectory file (format supported by ASE)",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="traj_aligned.xyz",
        help="Output trajectory file (format supported by ASE)",
    )
    parser.add_argument(
        "--ref",
        type=int,
        default=0,
        help="Index of reference frame to align to (default: 0)",
    )
    args = parser.parse_args()
    ref = args.ref
    try:
        # Read trajectory
        print("Reading trajectory...")
        traj = io.read(args.input, index=":")

        if not traj:
            print("Error: No frames read from trajectory")
            return

        print(f"Successfully read {len(traj)} frames")

        # Remove rotations
        print("Removing rotations...")
        aligned_traj = remove_rotations(traj, ref)

        # Write aligned trajectory
        print("Writing aligned trajectory...")
        io.write(args.output, aligned_traj)
        print("Done!")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        raise


if __name__ == "__main__":
    main()
