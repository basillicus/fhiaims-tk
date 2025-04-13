#!/usr/bin/env python3
"""
Remove the angular momentum (rotational motion) from an MD trajectory.

This script ensures that each frame’s rotational (angular) contribution is removed
by computing the inertia tensor, total angular momentum, and the corresponding
angular velocity. It then subtracts the rotational velocity from the atomic velocities.

Command-line interface:
    - Reads a trajectory from "trajectory.xyz" and writes the corrected trajectory to "trajectory_norot.xyz".

"""

import numpy as np
from ase import io
from ase.atoms import Atoms

def calculate_inertia_tensor(positions, masses):
    """
    Calculate the inertia tensor for a set of atoms.

    Parameters:
    positions (numpy.ndarray): Atomic positions (n_atoms, 3)
    masses (numpy.ndarray): Atomic masses (n_atoms,)

    Returns:
    numpy.ndarray: Inertia tensor (3, 3)
    """
    I = np.zeros((3, 3))

    for pos, mass in zip(positions, masses):
        r2 = np.sum(pos**2)
        I += mass * (r2 * np.eye(3) - np.outer(pos, pos))

    return I

def calculate_angular_momentum(positions, velocities, masses):
    """
    Calculate total angular momentum.

    Parameters:
    positions (numpy.ndarray): Atomic positions (n_atoms, 3)
    velocities (numpy.ndarray): Atomic velocities (n_atoms, 3)
    masses (numpy.ndarray): Atomic masses (n_atoms,)

    Returns:
    numpy.ndarray: Angular momentum vector (3,)
    """
    L = np.zeros(3)
    for pos, vel, mass in zip(positions, velocities, masses):
        L += mass * np.cross(pos, vel)
    return L

def remove_angular_momentum(atoms_list, dt=1.0):
    """
    Remove rotational motion from a trajectory by calculating and subtracting
    angular velocity components.

    Parameters:
    atoms_list (list): List of ASE Atoms objects representing the trajectory
    dt (float): Time step between frames

    Returns:
    list: New list of Atoms objects with rotations removed
    """
    # Make sure we have velocities
    if not atoms_list[0].has('momenta'):
        # Calculate velocities from positions if not present
        for i in range(len(atoms_list)-1):
            vel = (atoms_list[i+1].positions - atoms_list[i].positions) / dt
            atoms_list[i].set_velocities(vel)
        # Set last frame velocities equal to second-to-last frame
        atoms_list[-1].set_velocities(atoms_list[-2].get_velocities())

    new_atoms_list = []

    # Get masses
    masses = atoms_list[0].get_masses()

    for atoms in atoms_list:
        # Center system at center of mass
        com = atoms.get_center_of_mass()
        positions = atoms.positions - com
        velocities = atoms.get_velocities()

        # Calculate inertia tensor and angular momentum
        I = calculate_inertia_tensor(positions, masses)
        L = calculate_angular_momentum(positions, velocities, masses)

        # Calculate angular velocity
        omega = np.linalg.solve(I, L)

        # Subtract rotational velocity component
        for i, (pos, mass) in enumerate(zip(positions, masses)):
            v_rot = np.cross(omega, pos)
            velocities[i] -= v_rot

        # Create new Atoms object with updated velocities
        new_atoms = atoms.copy()
        new_atoms.set_velocities(velocities)
        new_atoms_list.append(new_atoms)

    return new_atoms_list

def main():
# Note on argparse: This script does not use argparse for CLI input – consider adding one 
    # Read trajectory
    traj = io.read('trajectory.xyz', index=':')

    # Remove rotations
    aligned_traj = remove_angular_momentum(traj)

    # Write new trajectory
    io.write('trajectory_norot.xyz', aligned_traj)


if __name__ == '__main__':
    main()
