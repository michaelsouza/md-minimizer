"""
generate_data.py

This script generates a LAMMPS data file for a 9x10 spring network,
ensuring standard whitespace characters are used for compatibility.
"""

import numpy as np


def generate_lammps_data():
    """Generates a LAMMPS data file for a 9x10 spring network."""

    # Parameters
    Nx = 9
    Ny = 10
    spacing = 1.0
    L = (Nx - 1) * spacing
    D = (Ny - 1) * spacing
    N_particles = Nx * Ny

    # Generate positions on a grid
    x = np.linspace(0, L, Nx)
    y = np.linspace(0, D, Ny)
    xv, yv = np.meshgrid(x, y)
    positions = np.stack([xv.flatten(), yv.flatten(), np.zeros(N_particles)]).T

    # Generate bonds between nearest neighbors
    bonds = []
    # Horizontal bonds
    for i in range(Ny):
        for j in range(Nx - 1):
            p1 = i * Nx + j
            p2 = i * Nx + j + 1
            bonds.append((p1+1, p2+1))
    # Vertical bonds
    for i in range(Ny - 1):
        for j in range(Nx):
            p1 = i * Nx + j
            p2 = (i + 1) * Nx + j
            bonds.append((p1+1, p2+1))

    # Write all data to the file
    with open("data.springs", "w", encoding="utf-8") as f:
        f.write("# LAMMPS data file for a 9x10 spring network\n\n")
        f.write(f"{N_particles} atoms\n")
        f.write(f"{len(bonds)} bonds\n")
        f.write("1 atom types\n")
        f.write("1 bond types\n\n")

        # Simulation box dimensions
        f.write(f"-0.1 {L+0.1} xlo xhi\n")
        f.write(f"-0.1 {D+0.1} ylo yhi\n")
        f.write("-0.6 0.6 zlo zhi\n\n")

        # Masses
        f.write("Masses\n\n")
        f.write("1 1.0\n\n")

        # Atoms section with corrected format
        f.write("Atoms\n\n")
        for i in range(N_particles):
            atom_id = i + 1
            molecule_id = 1
            atom_type = 1
            px, py, pz = positions[i]
            # Format: atom-ID molecule-ID atom-type x y z
            f.write(f"{atom_id} {molecule_id} {atom_type} {px:.2f} {py:.2f} {pz:.2f}\n")

        # Bonds section
        f.write("\nBonds\n\n")
        for i, bond in enumerate(bonds):
            bond_id = i + 1
            bond_type = 1
            p1, p2 = bond
            f.write(f"{bond_id} {bond_type} {p1} {p2}\n")

    print("'data.springs' file has been successfully generated.")


if __name__ == "__main__":
    generate_lammps_data()
