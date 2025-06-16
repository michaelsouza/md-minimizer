"""
visualize_simulation.py

This script reads a series of LAMMPS dump files from a looped simulation,
plots the state of the spring network for each step, and compiles the
resulting images into a movie. It now color-codes bonds by type.

The visualization style is matched to 'visualize_mesh.py' for consistency.
"""

import glob
from typing import Dict, List, Tuple

import imageio.v2 as imageio
import matplotlib.pyplot as plt


def parse_data_for_viz(
    data_filename: str,
) -> Tuple[
    Dict[int, Tuple[float, float]], List[Tuple[int, int, int]], Dict[str, float]
]:
    """
    Parses a LAMMPS data file to extract bonds (with types) and box dimensions.

    Args:
        data_filename: The path to the LAMMPS data file ('data.springs').

    Returns:
        A tuple containing:
        - A dictionary mapping atom IDs to their initial (x, y) coordinates.
        - A list of tuples, representing bonds (bond_type, atom1_id, atom2_id).
        - A dictionary with the simulation box boundaries ('xlo', 'xhi', 'ylo', 'yhi').
    """
    atoms = {}
    bonds = []
    box = {}

    reading_atoms = False
    reading_bonds = False

    with open(data_filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            if "xlo xhi" in line:
                parts = line.split()
                box["xlo"], box["xhi"] = float(parts[0]), float(parts[1])
            elif "ylo yhi" in line:
                parts = line.split()
                box["ylo"], box["yhi"] = float(parts[0]), float(parts[1])
            elif line.startswith("Atoms"):
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
    return atoms, bonds, box


def read_lammps_dump(dump_filename: str) -> Dict[int, Tuple[float, float]]:
    """
    Reads a single LAMMPS dump file and returns atom positions.

    Args:
        dump_filename: Path to the dump file for a single frame.

    Returns:
        A dictionary mapping atom IDs to their (x, y) coordinates for that frame.
    """
    positions = {}
    with open(dump_filename, "r", encoding="utf-8") as f:
        # Skip header lines
        for _ in range(9):
            next(f)
        # Read atom data
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 4:
                atom_id = int(parts[0])
                x, y = float(parts[2]), float(parts[3])
                positions[atom_id] = (x, y)
    return positions


def visualize_frame(
    step_num: int,
    positions: Dict[int, Tuple[float, float]],
    bonds: List[Tuple[int, int, int]],
    box: Dict[str, float],
) -> str:
    """
    Generates and saves a single plot of the mesh for a given simulation step.

    Args:
        step_num: The current simulation step number (for the title).
        positions: A dictionary of atom positions for the current step.
        bonds: A list of bonds connecting atom IDs (bond_type, p1, p2).
        box: A dictionary with the simulation box boundaries for setting axis limits.

    Returns:
        The filename of the saved image.
    """
    plt.style.use("seaborn-v0_8-whitegrid")
    _, ax = plt.subplots(figsize=(8, 10))

    # --- MODIFIED: Plot bonds with different colors ---
    bond_colors = {1: "darkblue", 2: "firebrick"}

    for bond_type, p1_id, p2_id in bonds:
        if p1_id in positions and p2_id in positions:
            x_coords = [positions[p1_id][0], positions[p2_id][0]]
            y_coords = [positions[p1_id][1], positions[p2_id][1]]
            color = bond_colors.get(bond_type, "grey")
            ax.plot(x_coords, y_coords, "-", linewidth=1.5, color=color, zorder=1)

    # Plot atoms
    x_vals = [pos[0] for pos in positions.values()]
    y_vals = [pos[1] for pos in positions.values()]
    ax.scatter(x_vals, y_vals, color="orangered", s=50, zorder=2, ec="black")

    # Formatting
    ax.set_xlabel("X-coordinate", fontsize=12)
    ax.set_ylabel("Y-coordinate", fontsize=12)
    ax.set_title(f"Deformation at Step {step_num+1}", fontsize=16, fontweight="bold")

    # Set fixed axis limits based on the full box size to prevent jitter
    ax.set_xlim(box["xlo"], box["xhi"])
    ax.set_ylim(box["ylo"], box["yhi"])
    ax.set_aspect("equal", adjustable="box")

    plt.tight_layout()
    output_filename = f"output/frame_{step_num:03d}.png"
    plt.savefig(output_filename)
    plt.close()
    return output_filename


def create_movie(image_files: List[str], output_filename: str):
    """
    Creates an animated GIF from a list of image files.

    Args:
        image_files: A sorted list of filenames for the animation frames.
        output_filename: The name of the output movie file (e.g., 'movie.gif').
    """
    with imageio.get_writer(
        f"output/{output_filename}", mode="I", duration=0.2, loop=0
    ) as writer:
        for filename in image_files:
            image = imageio.imread(filename)
            writer.append_data(image)
    print(f"\nMovie saved to 'output/{output_filename}'")


if __name__ == "__main__":
    DATA_FILE = "data.springs"
    DUMP_PATTERN = "output/dump.springs_loop.*.lammpstrj"

    # 1. Get bond connectivity and full box size from the data file
    try:
        _, bond_list, box_dims = parse_data_for_viz(DATA_FILE)
        if not bond_list or not box_dims:
            raise FileNotFoundError  # Trigger error if data is incomplete
    except FileNotFoundError:
        print(f"Error: Could not properly read '{DATA_FILE}'.")
        print("Please ensure it exists and is correctly formatted.")
        exit()

    # 2. Find all dump files from the simulation
    dump_files = sorted(glob.glob(DUMP_PATTERN))
    if not dump_files:
        print(f"Error: No dump files found matching the pattern '{DUMP_PATTERN}'.")
        print("Please run the LAMMPS simulation with the loop script first.")
        exit()

    # 3. Generate a plotted frame for each dump file
    frame_filenames = []
    print("Generating frames for movie...")
    for i, dump_file in enumerate(dump_files):
        atom_positions = read_lammps_dump(dump_file)
        filename = visualize_frame(i, atom_positions, bond_list, box_dims)
        frame_filenames.append(filename)
        print(f"  ... saved {filename}")

    # 4. Create the movie from the generated frames
    create_movie(frame_filenames, "simulation_movie.gif")
