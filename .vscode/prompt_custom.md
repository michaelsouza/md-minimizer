I've used the generate_data.py file to create the data.springs file.

===
# generate_data.py
import numpy as np

# Parameters
Nx = 10
Ny = 10
spacing = 1.0
L = (Nx - 1) * spacing
D = (Ny - 1) * spacing
N_particles = Nx * Ny

# Generate positions
x = np.linspace(0, L, Nx)
y = np.linspace(0, D, Ny)
xv, yv = np.meshgrid(x, y)
positions = np.stack([xv.flatten(), yv.flatten(), np.zeros(N_particles)]).T

# Generate bonds
bonds = []
# Horizontal
for i in range(Ny):
    for j in range(Nx - 1):
        p1 = i * Nx + j + 1
        p2 = i * Nx + j + 2
        bonds.append((p1, p2))
# Vertical
for i in range(Ny - 1):
    for j in range(Nx):
        p1 = i * Nx + j + 1
        p2 = (i + 1) * Nx + j + 1
        bonds.append((p1, p2))

# Write to file
with open("data.springs", "w") as f:
    f.write("# LAMMPS data file for a 10x10 spring network\n\n")
    f.write(f"{N_particles}  atoms\n")
    f.write(f"{len(bonds)}  bonds\n")
    f.write("1  atom types\n")
    f.write("1  bond types\n\n")
    f.write(f"0.0 {L} xlo xhi\n")
    f.write(f"0.0 {D} ylo yhi\n")
    f.write("-0.5 0.5 zlo zhi\n\n")
    f.write("Masses\n\n")
    f.write("1 1.0\n\n")

    f.write("Atoms\n\n")
    for i in range(N_particles):
        atom_id = i + 1
        atom_type = 1
        px, py, pz = positions[i]
        f.write(f"{atom_id} {atom_type} {px:.2f} {py:.2f} {pz:.2f}\n")

    f.write("\nBonds\n\n")
    for i, bond in enumerate(bonds):
        bond_id = i + 1
        bond_type = 1
        p1, p2 = bond
        f.write(f"{bond_id} {bond_type} {p1} {p2}\n")

print("data.springs file has been generated.")

===

lmp -in in.springs
LAMMPS (29 Aug 2024 - Update 1)
OMP_NUM_THREADS environment is not set. Defaulting to 1 thread. (src/src/comm.cpp:98)
  using 1 OpenMP thread(s) per MPI task
Reading data file ...
  orthogonal box = (0 0 -0.5) to (9 9 0.5)
  1 by 1 by 1 MPI processor grid
  reading atoms ...
ERROR: Incorrect format in Atoms section of data file: 1 1 0.00 0.00 0.00
For more information see https://docs.lammps.org/err0002 (src/src/atom.cpp:1089)
Last command: read_data       data.springs # Read particle positions and bonds from data file
