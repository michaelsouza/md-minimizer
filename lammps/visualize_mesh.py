"""
visualize_mesh.py

This script reads a LAMMPS data file ('data.springs') and visualizes the 
initial 2D mesh of atoms and bonds using matplotlib.

This is useful for verifying the output of 'generate_data.py'.
"""

from typing import Dict, List, Tuple
import matplotlib.pyplot as plt


def parse_lammps_data(
    filename: str,
) -> Tuple[Dict[int, Tuple[float, float]], List[Tuple[int, int]]]:
    """
    Parses a LAMMPS data file to extract atom positions and bonds.

    Args:
        filename: The path to the LAMMPS data file.

    Returns:
        A tuple containing:
        - A dictionary mapping atom IDs to their (x, y) coordinates.
        - A list of tuples, where each tuple represents a bond between two atom IDs.
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
                        p1, p2 = int(parts[2]), int(parts[3])
                        bonds.append((p1, p2))

    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        print("Please run 'generate_data.py' first to create it.")
        return {}, []

    return atoms, bonds


def visualize_mesh(atoms: Dict[int, Tuple[float, float]], bonds: List[Tuple[int, int]]):
    """
    Generates and displays a plot of the mesh from atom and bond data.

    Args:
        atoms: A dictionary of atom positions (ID -> (x, y)).
        bonds: A list of bonds connecting atom IDs.
    """
    if not atoms:
        print("No atom data to visualize.")
        return

    plt.style.use("seaborn-v0_8-whitegrid")
    _, ax = plt.subplots(figsize=(10, 10))

    # --- Plot Bonds ---
    # Draw the lines first so they appear behind the atom markers.
    for p1_id, p2_id in bonds:
        # Check if both atoms in the bond exist in the atoms dictionary
        if p1_id in atoms and p2_id in atoms:
            x_coords = [atoms[p1_id][0], atoms[p2_id][0]]
            y_coords = [atoms[p1_id][1], atoms[p2_id][1]]
            ax.plot(
                x_coords, y_coords, "k-", linewidth=1.5, color="royalblue", zorder=1
            )

    # --- Plot Atoms ---
    x_vals = [pos[0] for pos in atoms.values()]
    y_vals = [pos[1] for pos in atoms.values()]
    ax.scatter(x_vals, y_vals, color="orangered", s=50, zorder=2, ec="black")

    # --- Formatting ---
    ax.set_xlabel("X-coordinate", fontsize=12)
    ax.set_ylabel("Y-coordinate", fontsize=12)
    ax.set_title("Initial Spring Network Mesh", fontsize=16, fontweight="bold")

    # Ensure the scaling is equal in x and y directions to see the true geometry
    ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()
    plt.savefig("output/initial_mesh.png")

    print("Mesh visualization saved to 'output/initial_mesh.png'")


if __name__ == "__main__":
    DATA_FILE = "data.springs"
    atom_positions, bond_connections = parse_lammps_data(DATA_FILE)
    visualize_mesh(atom_positions, bond_connections)
