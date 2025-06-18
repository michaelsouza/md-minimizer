"""
create_network.py

This single script creates a spring network model and directly generates a
LAMMPS data file ready for simulations with avalanche dynamics.

It merges the functionalities of:
1. create_network.py: Generates the network topology.
2. json2lammps.py: Converts network data to LAMMPS format.
3. add_breaking_thresholds.py: Assigns unique, random breaking
   thresholds to each breakable spring.

The output is a single LAMMPS data file that includes:
- Atom positions and types (mobile, fixed, pull).
- Bond topology.
- A unique bond type for each breakable spring.
- A 'Bond Coeffs' section defining a random breaking length for each of
  these unique bond types, enabling the simulation of fracture avalanches
  as described by Noguchi et al. (2024).

Usage:
    python lammps/create_network.py [OPTIONS]
"""

import os
import math
import random
import argparse
import networkx as nx
import matplotlib.pyplot as plt


# --- Core Network Generation Functions ---


def get_node_indices(node_id, N):
    """
    Convert a 1D node ID to 2D indices (j, i).
    """
    if not (0 <= node_id < N * N):
        raise ValueError(f"node_id {node_id} out of range for N={N}")
    j = node_id // N
    i = node_id % N
    return j, i


def get_node_id(j, i, N):
    """
    Convert 2D indices (j, i) to a 1D node ID.
    """
    if i == -1:
        i = N - 1
    if i == N:
        i = 0
    return j * N + i


def get_neighbors(node_id, N):
    """
    Get the IDs of all neighbors of a given node.
    """
    j, i = get_node_indices(node_id, N)
    neighbor_ids = []
    if j < (N - 1):
        neighbor_ids.append(get_node_id(j + 1, i, N))
        if j % 2 == 0:
            neighbor_ids.append(get_node_id(j + 1, i - 1, N))
        else:
            neighbor_ids.append(get_node_id(j + 1, i + 1, N))
    neighbor_ids.append(get_node_id(j, i + 1, N))
    return list(set(neighbor_ids))


def is_unbreakable(node1_id, node2_id, N, L_matrix):
    """
    Determines if the edge between two nodes is unbreakable based on the
    matrix structure defined by L_matrix.

    Args:
        node1_id: ID of the first node.
        node2_id: ID of the second node.
        N: Grid size.
        L_matrix: Matrix size parameter.

    Returns:
        True if the edge should be unbreakable, False otherwise.
    """
    if L_matrix <= 0:  # No matrix structure
        return False

    j1, i1 = get_node_indices(node1_id, N)
    j2, i2 = get_node_indices(node2_id, N)

    # Ensure j1 <= j2 for easier checking
    if j1 > j2:
        j1, j2 = j2, j1
        i1, i2 = i2, i1
        node1_id, node2_id = node2_id, node1_id  # Keep IDs consistent if needed

    # --- Check for Horizontal Unbreakable Spring ---
    # Springs are horizontal if they connect nodes in the same row
    # These exist *on* rows that are multiples of L_matrix
    if j1 == j2 and j1 % L_matrix == 0:
        # Check if they are horizontal neighbors (i differs by 1, wrapping around)
        if abs(i1 - i2) == 1 or abs(i1 - i2) == N - 1:
            return True

    # --- Check for Zigzag Unbreakable Spring ---
    # Springs connect row j to row j+1, where j is a multiple of L_matrix
    if (
        i1 % L_matrix == 0 and i2 == i1
    ):  # Must have same column index for this type of connection
        return True

    # If none of the above conditions are met, it's breakable
    return False


def create_spring_network(N, L_matrix, l0=1.0):
    """
    Create a spring network with the specified parameters.
    """
    print(f"Generating {N}x{N} spring network with L_matrix={L_matrix}...")
    G = nx.Graph()
    for j in range(N):
        for i in range(N):
            node_id = j * N + i
            y_coord = j * l0 * math.sqrt(3) / 2
            x_coord = (i + 0.5 * (j % 2)) * l0
            G.add_node(node_id, pos=(x_coord, y_coord), row=j, col=i)

    for node_id in range(N * N):
        neighbors = get_neighbors(node_id, N)
        for neighbor_id in neighbors:
            if not G.has_edge(node_id, neighbor_id):
                is_unbreak = is_unbreakable(node_id, neighbor_id, N, L_matrix)
                G.add_edge(node_id, neighbor_id, is_unbreak=is_unbreak)

    print(f"Generated {G.number_of_nodes()} nodes and {G.number_of_edges()} edges.")
    return G


# --- LAMMPS File Generation Function ---


def write_lammps_data_file(
    G, data_filename, thresholds_filename, l0=1.0, strain_mean=0.1, strain_std=0.02
):
    """
    Write the spring network to a LAMMPS data file and a separate file
    for breaking thresholds.
    """
    print(f"Generating LAMMPS data file: {data_filename}")
    nodes = list(G.nodes(data=True))
    edges = list(G.edges(data=True))
    num_atoms = len(nodes)
    num_bonds = len(edges)
    max_row = max(data["row"] for _, data in nodes) if nodes else 0

    atom_lines, x_coords, y_coords = [], [], []
    for node_id, data in nodes:
        atom_type = 1  # Default mobile atom
        if data["row"] == 0:
            atom_type = 2  # Fixed bottom row
        elif data["row"] == max_row:
            atom_type = 3  # Pull top row
        x, y = data["pos"]
        # For atom_style bond, we need: atom-ID molecule-ID atom-type x y z
        atom_lines.append(f"{node_id+1} 1 {atom_type} {x:.6f} {y:.6f} 0.0")
        x_coords.append(x)
        y_coords.append(y)

    bond_lines = []
    bond_coeff_lines = []
    breaking_threshold_lines = []

    # Bond type 1 is for unbreakable bonds, types 2+ are for breakable bonds.
    bond_coeff_lines.append("1 1.0 1.0")
    bond_id, breakable_type = 1, 2

    for u, v, data in edges:
        if data.get("is_unbreak"):
            bond_lines.append(f"{bond_id} 1 {u+1} {v+1}")
        else:
            # Each breakable bond gets a unique type to allow individual breaking
            bond_coeff_lines.append(f"{breakable_type} 1.0 1.0")

            # Determine breaking length and write to separate file
            breaking_strain = random.uniform(0.0, 1.0)
            breaking_length = l0 * (1 + breaking_strain)
            breaking_threshold_lines.append(f"{breakable_type} {breaking_length:.6f}")

            bond_lines.append(f"{bond_id} {breakable_type} {u+1} {v+1}")
            breakable_type += 1
        bond_id += 1

    num_atom_types, num_bond_types = 3, breakable_type - 1

    # Write the main data file
    with open(data_filename, "w", encoding="utf-8", newline="\n") as f:
        f.write(
            f"LAMMPS data file for spring network (N={int(math.sqrt(num_atoms))})\n\n"
        )
        f.write(f"{num_atoms} atoms\n{num_bonds} bonds\n\n")
        f.write(f"{num_atom_types} atom types\n{num_bond_types} bond types\n\n")
        buffer = 1.0
        f.write(f"{min(x_coords)-buffer:.6f} {max(x_coords)+buffer:.6f} xlo xhi\n")
        f.write(f"{min(y_coords)-buffer:.6f} {max(y_coords)+buffer:.6f} ylo yhi\n")
        f.write(f"{-1.0:.6f} {1.0:.6f} zlo zhi\n\n")
        f.write("Masses\n\n1 1.0\n2 1.0\n3 1.0\n\n")
        f.write("Bond Coeffs\n\n")
        f.write("\n".join(bond_coeff_lines))
        f.write("\n\n")
        f.write("Atoms # id molecule-id type x y z\n\n")
        f.write("\n".join(atom_lines) + "\n\n")
        f.write("Bonds # id type p1 p2\n\n")
        f.write("\n".join(bond_lines) + "\n")
    print(f"Successfully wrote LAMMPS data file with {num_bond_types} bond types.")

    # Write the breaking thresholds file
    print(f"Generating breaking thresholds file: {thresholds_filename}")
    with open(thresholds_filename, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Bond Type, Breaking Length\n")
        f.write("\n".join(breaking_threshold_lines) + "\n")
    print("Successfully wrote breaking thresholds file.")


def display_network(G, save_path):
    """
    Display the spring network as a PNG image, EXCLUDING periodic boundary
    connections for better visualization.
    """
    print(f"Saving network visualization to {save_path}...")
    plt.figure(figsize=(12, 12))
    pos = nx.get_node_attributes(G, "pos")

    # Determinar o tamanho da rede N para identificar as bordas
    # O número de nós é N*N, então N é a raiz quadrada do número de nós.
    N = int(math.sqrt(G.number_of_nodes()))

    # --- Filtro para remover as ligações de contorno periódico ---
    # Uma ligação (u, v) é periódica se conecta a coluna 0 à coluna N-1.
    # Isso acontece quando a diferença absoluta dos índices de coluna é N-1.

    # Lista de arestas quebráveis, excluindo as periódicas
    breakable = [
        (u, v)
        for u, v, d in G.edges(data=True)
        if not d.get("is_unbreak")
        and abs(G.nodes[u]["col"] - G.nodes[v]["col"]) != N - 1
    ]

    # Lista de arestas inquebráveis, excluindo as periódicas
    unbreakable = [
        (u, v)
        for u, v, d in G.edges(data=True)
        if d.get("is_unbreak") and abs(G.nodes[u]["col"] - G.nodes[v]["col"]) != N - 1
    ]

    # --- Desenha os componentes da rede com as listas filtradas ---
    nx.draw_networkx_nodes(G, pos, node_size=15, node_color="black")
    nx.draw_networkx_edges(
        G, pos, edgelist=breakable, width=0.3, style="solid", edge_color="black"
    )
    nx.draw_networkx_edges(
        G, pos, edgelist=unbreakable, width=1.5, style="solid", edge_color="red"
    )

    # O bloco para desenhar os labels dos nós já deve estar comentado ou removido.

    plt.axis("equal")
    plt.margins(0.1)
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate a LAMMPS data file for a spring network simulation.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--N", type=int, default=12, help="Grid size (NxN nodes).")
    parser.add_argument(
        "--L_matrix",
        type=int,
        default=4,
        help="Size of the unbreakable matrix structure.",
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        default=".",
        help="Directory to save the output files.",
    )
    parser.add_argument(
        "--skip_png",
        action="store_true",
        help="Do not generate the network visualization PNG file.",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    network_graph = create_spring_network(args.N, args.L_matrix)

    base_filename = f"N{args.N}_Lmat{args.L_matrix}"
    lammps_filepath = os.path.join(args.output_dir, f"{base_filename}.data")
    thresholds_filepath = os.path.join(
        args.output_dir, f"{base_filename}_breaking_thresholds.dat"
    )

    write_lammps_data_file(network_graph, lammps_filepath, thresholds_filepath)

    if not args.skip_png:
        png_filepath = os.path.join(args.output_dir, f"{base_filename}.png")
        display_network(network_graph, png_filepath)

    print("\nScript finished successfully.")
    print("To run the simulation, use the generated files with LAMMPS:")
    print(
        f"python spring_network.py "
        f"-v data_file {lammps_filepath} "
        f"-v thresholds_filename {thresholds_filepath}"
    )
