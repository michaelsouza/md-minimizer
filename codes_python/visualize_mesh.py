"""
visualize_mesh.py

This script reads a LAMMPS data file ('data.springs') and visualizes the
initial 2D mesh of atoms and bonds using matplotlib.
It now color-codes bonds based on their type.

This is useful for verifying the output of 'generate_data.py'.
"""

import sys
from typing import Dict, List, Tuple
import matplotlib.pyplot as plt


def parse_lammps_data(
    filename: str,
) -> Tuple[Dict[int, Tuple[float, float]], List[Tuple[int, int, int]]]:
    """
    Parses a LAMMPS data file to extract atom positions and typed bonds.

    Args:
        filename: The path to the LAMMPS data file.

    Returns:
        A tuple containing:
        - A dictionary mapping atom IDs to their (x, y) coordinates.
        - A list of tuples, where each tuple represents a bond
          (bond_type, atom1_id, atom2_id).
    """
    atoms = {}
    bonds = []

    # Flags to indicate which section of the file we are in
    reading_atoms = False
    reading_bonds = False

    try:
        with open(filename, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()

                # Skip empty lines or comments
                if not line or line.startswith("#"):
                    continue

                if line.startswith("Atoms"):
                    reading_atoms = True
                    reading_bonds = False
                    continue
                elif line.startswith("Bonds"):
                    reading_atoms = False
                    reading_bonds = True
                    continue

                if reading_atoms:
                    parts = line.split()
                    if len(parts) >= 5:
                        atom_id = int(parts[0])
                        x, y = float(parts[3]), float(parts[4])
                        atoms[atom_id] = (x, y)

                elif reading_bonds:
                    parts = line.split()
                    if len(parts) >= 4:
                        # --- MODIFIED: Read bond type (column 2) ---
                        bond_type = int(parts[1])
                        p1, p2 = int(parts[2]), int(parts[3])
                        bonds.append((bond_type, p1, p2))

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        print("Please run 'generate_data.py' first to create it.")
        return {}, []

    return atoms, bonds


def visualize_mesh(
    atoms: Dict[int, Tuple[float, float]], bonds: List[Tuple[int, int, int]]
):
    """
    Generates and displays a plot of the mesh from atom and typed bond data.

    Args:
        atoms: A dictionary of atom positions (ID -> (x, y)).
        bonds: A list of bonds (bond_type, p1, p2).
    """
    if not atoms:
        print("No atom data to visualize.")
        return

    plt.style.use("seaborn-v0_8-whitegrid")
    _, ax = plt.subplots(figsize=(10, 10))

    # --- MODIFIED: Plot Bonds with different colors ---
    bond_colors = {1: "darkblue", 2: "firebrick"}
    bond_labels = {1: "Stiff Spring (Type 1)", 2: "Soft Spring (Type 2)"}
    plotted_labels = set()  # To avoid duplicate legend entries

    for bond_type, p1_id, p2_id in bonds:
        if p1_id in atoms and p2_id in atoms:
            x_coords = [atoms[p1_id][0], atoms[p2_id][0]]
            y_coords = [atoms[p1_id][1], atoms[p2_id][1]]
            color = bond_colors.get(bond_type, "grey")
            label = bond_labels.get(bond_type)

            # Plot with a label only once to keep the legend clean
            if label and label not in plotted_labels:
                ax.plot(
                    x_coords,
                    y_coords,
                    "-",
                    linewidth=1.5,
                    color=color,
                    zorder=1,
                    label=label,
                )
                plotted_labels.add(label)
            else:
                ax.plot(x_coords, y_coords, "-", linewidth=1.5, color=color, zorder=1)

    # --- Plot Atoms ---
    x_vals = [pos[0] for pos in atoms.values()]
    y_vals = [pos[1] for pos in atoms.values()]
    ax.scatter(x_vals, y_vals, color="orangered", s=50, zorder=2, ec="black")

    # --- Formatting ---
    ax.set_xlabel("X-coordinate", fontsize=12)
    ax.set_ylabel("Y-coordinate", fontsize=12)
    ax.set_title("Initial Spring Network Mesh", fontsize=16, fontweight="bold")
    ax.legend(title="Bond Types")
    ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()
    plt.savefig("output/initial_mesh.png")

    print("Mesh visualization saved to 'output/initial_mesh.png'")


if __name__ == "__main__":
    DATA_FILE = sys.argv[1]
    atom_positions, bond_connections = parse_lammps_data(DATA_FILE)
    visualize_mesh(atom_positions, bond_connections)
