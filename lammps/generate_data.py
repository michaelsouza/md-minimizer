"""
generate_data.py

This script generates a LAMMPS data file for a 9x10 spring network,
ensuring the simulation box is large enough for a tensile test.
It now assigns unique atom types to the top and bottom bars for robust grouping.
It also randomly assigns one of two types to each spring (bond).
"""

import numpy as np
import random


def generate_lammps_data():
    """Generates a LAMMPS data file for a 9x10 spring network."""

    # --- System Parameters ---
    Nx = 9
    Ny = 10
    spacing = 1.0
    buffer = 0.1  # A small buffer added to the box boundaries

    # --- Simulation Parameters (should match in.springs_loop) ---
    pull_steps = 10
    displacement_per_step_factor = 0.1  # This is 0.1 * D

    # --- Calculations ---
    L = (Nx - 1) * spacing
    D = (Ny - 1) * spacing
    N_particles = Nx * Ny
    total_displacement = pull_steps * (displacement_per_step_factor * D)
    final_y_height = D + total_displacement

    x = np.linspace(0, L, Nx)
    y = np.linspace(0, D, Ny)
    xv, yv = np.meshgrid(x, y)
    positions = np.stack([xv.flatten(), yv.flatten(), np.zeros(N_particles)]).T

    # Generate bonds
    bonds = []
    # Horizontal bonds
    for i in range(Ny):
        for j in range(Nx - 1):
            p1 = i * Nx + j
            p2 = i * Nx + j + 1
            bonds.append((p1 + 1, p2 + 1))
    # Vertical bonds
    for i in range(Ny - 1):
        for j in range(Nx):
            p1 = i * Nx + j
            p2 = (i + 1) * Nx + j
            bonds.append((p1 + 1, p2 + 1))

    # --- Write all data to the file ---
    with open("data.springs", "w", encoding="utf-8") as f:
        f.write("# LAMMPS data file for a 9x10 spring network\n")
        f.write("# Atom types: 1=mobile, 2=bottom_bar, 3=top_bar\n")
        f.write("# Bond types: 1=stiff, 2=soft (randomly assigned)\n\n")
        f.write(f"{N_particles} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write("3 atom types\n")
        f.write("2 bond types\n\n")

        f.write(f"{-buffer} {L + buffer} xlo xhi\n")
        f.write(f"{-buffer} {final_y_height + buffer} ylo yhi\n")
        f.write("-0.6 0.6 zlo zhi\n\n")

        # Masses for each atom type
        f.write("Masses\n\n")
        f.write("1 1.0\n")
        f.write("2 1.0\n")
        f.write("3 1.0\n\n")

        # Atoms section
        f.write("Atoms\n\n")
        for i in range(N_particles):
            atom_id = i + 1
            molecule_id = 1
            if i < Nx:
                atom_type = 2
            elif i >= N_particles - Nx:
                atom_type = 3
            else:
                atom_type = 1
            px, py, pz = positions[i]
            f.write(f"{atom_id} {molecule_id} {atom_type} {px:.2f} {py:.2f} {pz:.2f}\n")

        # Bonds section
        f.write("\nBonds\n\n")
        for i, bond in enumerate(bonds):
            bond_id = i + 1
            # --- MODIFIED: Randomly assign bond type 1 or 2 ---
            bond_type = random.choice([1, 2])
            p1, p2 = bond
            f.write(f"{bond_id} {bond_type} {p1} {p2}\n")

    print(
        "'data.springs' file with typed atoms and bonds has been successfully generated."
    )


if __name__ == "__main__":
    generate_lammps_data()
