import numpy as np
from ase import Atoms
from ase.io.extxyz import write_extxyz, read_extxyz


def merge_datasets(dipoles, polarizabilities, search_range=5, periodic=False, convert_dipoles=False, convert_alphas=False):
    """
    Merges dipole and polarizability datasets based on matching coordinates.
    Uses step values from dipoles as the reference.

    Parameters:
    dipoles (numpy.ndarray): Structured array containing dipole-related data.
    polarizabilities (numpy.ndarray): Structured array containing polarizability-related data.
    search_range (int): Number of steps before/after to search for coordinate matches.
    periodic (bool): Whether the system is periodic, including lattice and stress if True.
    convert_dipoles (bool): Convert dipoles from eAng to Debye (as per MACE output dipole)
    convert_alphas (bool): Convert alphas from Bohr**3 to meA^2/V (as per MACE RMSE polarizability output)


    Returns:
    numpy.ndarray: Merged dataset with matched dipole and polarizability information.
    """
    dipole_dict = {d['step']: d for d in dipoles}
    merged_data = []

    num_atoms = dipoles[0]['coordinates'].shape[0]  # Infer number of atoms

    for polar_entry in polarizabilities:
        step = polar_entry['step']
        matching_dipole = dipole_dict.get(step)

        if matching_dipole is not None and np.allclose(polar_entry['coordinates'], matching_dipole['coordinates']):
            correct_step = matching_dipole['step']
        else:
            correct_step = None
            for offset in range(1, search_range + 1):
                for test_step in (step - offset, step + offset):
                    if test_step in dipole_dict and np.allclose(polar_entry['coordinates'], dipole_dict[test_step]['coordinates']):
                        correct_step = test_step
                        break
                if correct_step is not None:
                    break

        if correct_step is not None:
            matching_dipole = dipole_dict[correct_step]
            merged_entry = (
                matching_dipole['step'],
                matching_dipole['species'],
                matching_dipole['coordinates'],
                matching_dipole['velocities'],
                matching_dipole['dipole'],
                matching_dipole['forces'],
                polar_entry['energy'],
                polar_entry['polarizability']
            )
            if periodic:
                merged_entry += (matching_dipole['lattice'], matching_dipole['stress'])
            merged_data.append(merged_entry)

    merged_dtype = [
        ('step', '<i8'),
        ('species', 'S2', (num_atoms,)),
        ('coordinates', '<f8', (num_atoms, 3)),
        ('velocities', '<f8', (num_atoms, 3)),
        ('dipole', '<f8', (3,)),
        ('forces', '<f8', (num_atoms, 3)),
        ('energy', '<f8'),
        ('polarizability', '<f8', (3, 3))
    ]
    if periodic:
        merged_dtype.extend([('lattice', '<f8', (3, 3)), ('stress', '<f8', (3, 3))])


    merged_array = np.array(merged_data, dtype=np.dtype(merged_dtype))
    if convert_dipoles:
        merged_array['dipole'] = merged_array['dipole'] / 0.20819434   # eAng to Debye

    if convert_dipoles:
        # Bohr^3 to Angstrom^3
        merged_array['polarizability'] = merged_array['polarizability'] * 0.14818471

    return merged_array


def write_extxyz_file(merged_dataset, filename="merged_data.extxyz", periodic=False, lattice=None,
                      force_key="REF_forces", energy_key="REF_energy", dipole_key="REF_dipole",
                      polarizability_key="REF_polarizability", stress_key="REF_stress"):
    """
    Writes the merged dataset to an .extxyz file in ASE format.

    Parameters:
    merged_dataset (numpy.ndarray): The structured array containing merged data.
    filename (str): Name of the output .extxyz file.
    periodic (bool): If True, includes lattice and periodic boundary conditions.
    lattice (numpy.ndarray): 3x3 lattice vectors for periodic systems, default is None.
    """
    atoms_list = []

    for entry in merged_dataset:
        atoms = Atoms(symbols=[s.decode() for s in entry['species']], positions=entry['coordinates'])
        atoms.info[energy_key] = entry['energy']
        atoms.info[dipole_key] = " ".join(map(str, entry['dipole']))
        atoms.info[polarizability_key] = " ".join(map(str, entry['polarizability'].flatten()))

        if periodic and lattice is not None:
            atoms.set_pbc(True)
            atoms.set_cell(entry['lattice'])
            atoms.info["Lattice"] = " ".join(map(str, entry['lattice'].flatten()))
            atoms.info[stress_key] = " ".join(map(str, entry['stress'].flatten()))
        else:
            atoms.info["pbc"] = "F F F"

        atoms.set_momenta(entry['velocities'])
        atoms.set_array(force_key, entry['forces'])
        atoms_list.append(atoms)

    with open(filename, "w") as f:
        write_extxyz(f, atoms_list)


def read_extxyz_file(filename="merged_data.extxyz", force_key="REF_forces", energy_key="REF_energy",
                     dipole_key="REF_dipole", polarizability_key="REF_polarizability", stress_key="REF_stress"):
    """
    Reads an .extxyz file and returns a structured numpy array.

    Parameters:
    filename (str): Name of the input .extxyz file.

    Returns:
    numpy.ndarray: Structured array containing the extracted data.
    """
    atoms_list = list(read_extxyz(filename, index=slice(None)))

    num_atoms = len(atoms_list[0].get_chemical_symbols())
    periodic = "Lattice" in atoms_list[0].info

    merged_dtype = [
        ('step', '<i8'),
        ('species', 'S2', (num_atoms,)),
        ('coordinates', '<f8', (num_atoms, 3)),
        ('velocities', '<f8', (num_atoms, 3)),
        ('dipole', '<f8', (3,)),
        ('forces', '<f8', (num_atoms, 3)),
        ('energy', '<f8'),
        ('polarizability', '<f8', (3, 3))
    ]
    if periodic:
        merged_dtype.extend([('lattice', '<f8', (3, 3)), ('stress', '<f8', (3, 3))])

    merged_data = []
    for atoms in atoms_list:
        species = np.array(atoms.get_chemical_symbols(), dtype='S2')
        coordinates = atoms.get_positions()
        velocities = atoms.get_momenta() if "momenta" in atoms.arrays else np.zeros((num_atoms, 3))
        forces = atoms.get_array(force_key) if force_key in atoms.arrays else np.zeros((num_atoms, 3))
        energy = atoms.info.get(energy_key, 0.0)
        dipole = np.array([float(x) for x in atoms.info.get(dipole_key, "0 0 0").split()])
        polarizability = np.array([float(x) for x in atoms.info.get(polarizability_key, "0 0 0 0 0 0 0 0 0").split()]).reshape(3, 3)

        lattice = np.array(atoms.get_cell()) if periodic else np.zeros((3, 3))
        stress = np.array([float(x) for x in atoms.info.get(stress_key, "0 0 0 0 0 0 0 0 0").split()]).reshape(3, 3) if periodic else np.zeros((3, 3))

        merged_entry = (0, species, coordinates, velocities, dipole, forces, energy, polarizability, lattice,
                        stress) if periodic else (0, species, coordinates, velocities, dipole, forces, energy, polarizability)
        merged_data.append(merged_entry)

    return np.array(merged_data, dtype=np.dtype(merged_dtype))
