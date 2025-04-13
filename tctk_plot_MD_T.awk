#!/usr/bin/awk -f

BEGIN {
    md_time_step = 0
    n_atoms = 0
    n_lattice_vectors = 0
    md_step_index = 0
    in_atomic_structure = 0
    in_forces = 0
    atomic_coordinates = ""
    species = ""
    velocities = ""
    total_time = 0
}

{
    if ($0 ~ /Molecular dynamics time step/) {
        md_time_step = $6
    }
    if ($0 ~ /Number of atoms/) {
        n_atoms = $6
    }
    if ($0 ~ /Number of lattice vectors/) {
        n_lattice_vectors = $6
    }
    if ($0 ~ /lattice_vector/) {
        lattice_vector[n_lattice_vectors] = $0
        n_lattice_vectors++
    }
    if ($0 ~ /\| Time step number/) {
        md_step[md_step_index] = $6
        md_step_index++
    }
    if ($0 ~ /atomic structure \(and velocities\)/) {
        in_atomic_structure = 1
        atomic_coordinates = ""
        species = ""
        velocities = ""
    }
    if (in_atomic_structure && $0 !~ /atomic structure \(and velocities\)/ && $0 !~ /\|/) {
        species = species $5 " "
        atomic_coordinates = atomic_coordinates $2 " " $3 " " $4 "\n"
        velocities = velocities $2 " " $3 " " $4 "\n"
    }
    if (in_atomic_structure && $0 ~ /\| Total time/) {
        in_atomic_structure = 0
        atomic_coordinates = atomic_coordinates substr(atomic_coordinates, 1, length(atomic_coordinates) - 1)
        species = species substr(species, 1, length(species) - 1)
        velocities = velocities substr(velocities, 1, length(velocities) - 1)
        atomic_structure[md_step[md_step_index - 1]] = atomic_coordinates
        species_array[md_step[md_step_index - 1]] = species
        velocities_array[md_step[md_step_index - 1]] = velocities
    }
    if ($0 ~ /Total atomic forces/) {
        in_forces = 1
        forces = ""
    }
    if (in_forces && $0 !~ /Total atomic forces/) {
        forces = forces $3 " " $4 " " $5 "\n"
    }
    if (in_forces && $0 ~ /\| Total time/) {
        in_forces = 0
        forces = forces substr(forces, 1, length(forces) - 1)
        forces_array[md_step[md_step_index - 1]] = forces
    }
    if ($0 ~ /Temperature \(nuclei\)/) {
        temperature[md_step[md_step_index - 1]] = $5
    }
}

END {
    for (i = 0; i < md_step_index; i++) {
        print "MD Step:", md_step[i]
        print "Species:", species_array[md_step[i]]
        print "Atomic Coordinates:", atomic_structure[md_step[i]]
        print "Velocities:", velocities_array[md_step[i]]
        print "Forces:", forces_array[md_step[i]]
        print "Temperature:", temperature[md_step[i]]
        print ""
    }
}

