"""
visualize_simulation.py

Este script lê uma série de arquivos de dump do LAMMPS, plota o estado da
rede para cada passo e compila as imagens em um filme. As ligações
inquebráveis (tipo 1) são coloridas de vermelho.
"""

import glob
import os
from typing import Dict, List, Tuple, Set

import imageio.v2 as imageio
import matplotlib.pyplot as plt


def parse_initial_unbreakable_bonds(data_filename: str) -> Set[Tuple[int, int]]:
    """
    Lê o arquivo de dados inicial do LAMMPS para encontrar todas as
    ligações inquebráveis (tipo 1).

    Args:
        data_filename: Caminho para o arquivo .data inicial.

    Returns:
        Um conjunto de tuplas, cada uma representando uma ligação inquebrável com
        os IDs dos átomos ordenados (menor, maior).
    """
    unbreakable_bonds = set()
    reading_bonds = False
    try:
        with open(data_filename, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                if line.startswith("Bonds"):
                    reading_bonds = True
                    continue
                if line.startswith("Atoms"):  # Seção de átomos vem antes
                    reading_bonds = False

                if reading_bonds:
                    parts = line.split()
                    if len(parts) >= 4:
                        bond_type = int(parts[1])
                        if bond_type == 1:  # Tipo 1 é inquebrável
                            p1, p2 = int(parts[2]), int(parts[3])
                            # Ordena os IDs para consistência na busca
                            unbreakable_bonds.add(tuple(sorted((p1, p2))))
    except FileNotFoundError:
        print(f"Aviso: Arquivo de dados inicial '{data_filename}' não encontrado.")
        print("Não será possível colorir as ligações inquebráveis.")
    return unbreakable_bonds


def read_atom_dump(dump_filename: str) -> Dict[int, Tuple[float, float]]:
    """Lê um único arquivo de dump de átomos do LAMMPS e retorna as posições."""
    positions = {}
    with open(dump_filename, "r", encoding="utf-8") as f:
        lines = f.readlines()
        item_line_index = -1
        for i, line in enumerate(lines):
            if line.strip().startswith("ITEM: ATOMS"):
                item_line_index = i
                break
        if item_line_index == -1:
            return positions
        for line in lines[item_line_index + 1 :]:
            parts = line.strip().split()
            if len(parts) >= 4:
                atom_id = int(parts[0])
                x, y = float(parts[2]), float(parts[3])
                positions[atom_id] = (x, y)
    return positions


def read_bond_dump(bond_dump_filename: str) -> List[Tuple[int, int]]:
    """Lê um único arquivo de dump de ligações do LAMMPS."""
    bonds = []
    if not os.path.exists(bond_dump_filename):
        return bonds
    with open(bond_dump_filename, "r", encoding="utf-8") as f:
        lines = f.readlines()
        item_line_index = -1
        for i, line in enumerate(lines):
            if line.strip().startswith("ITEM: ENTRIES"):
                item_line_index = i
                break
        if item_line_index == -1:
            return bonds
        for line in lines[item_line_index + 1 :]:
            parts = line.strip().split()
            if len(parts) >= 2:
                p1, p2 = int(parts[0]), int(parts[1])
                bonds.append((p1, p2))
    return bonds


def visualize_frame(
    step_num: int,
    positions: Dict[int, Tuple[float, float]],
    active_bonds: List[Tuple[int, int]],
    unbreakable_bonds: Set[Tuple[int, int]],  # <-- NOVO ARGUMENTO
    box: Dict[str, float],
) -> str:
    """Gera e salva um único quadro da simulação."""
    plt.style.use("seaborn-v0_8-whitegrid")
    _, ax = plt.subplots(figsize=(8, 10))

    for p1_id, p2_id in active_bonds:
        if p1_id in positions and p2_id in positions:
            pos1 = positions[p1_id]
            pos2 = positions[p2_id]

            # --- LÓGICA DE COR MODIFICADA ---
            # Verifica se a ligação pertence ao conjunto de inquebráveis
            bond_pair = tuple(sorted((p1_id, p2_id)))
            color = "red" if bond_pair in unbreakable_bonds else "darkblue"

            ax.plot(
                [pos1[0], pos2[0]],
                [pos1[1], pos2[1]],
                "-",
                linewidth=1.5,
                color=color,
                zorder=1,
            )

    if positions:
        x_vals = [pos[0] for pos in positions.values()]
        y_vals = [pos[1] for pos in positions.values()]
        ax.scatter(x_vals, y_vals, color="orangered", s=50, zorder=2, ec="black")

    ax.set_xlabel("Coordenada X", fontsize=12)
    ax.set_ylabel("Coordenada Y", fontsize=12)
    ax.set_title(f"Deformação no Passo {step_num}", fontsize=16, fontweight="bold")
    ax.set_xlim(box["xlo"], box["xhi"])
    ax.set_ylim(box["ylo"], box["yhi"])
    ax.set_aspect("equal", adjustable="box")
    plt.tight_layout()
    output_filename = f"frame_{step_num:04d}.png"
    plt.savefig(output_filename)
    plt.close()
    return output_filename


def create_movie(image_files: List[str], output_filename: str):
    """Cria um GIF animado a partir de uma lista de arquivos de imagem."""
    with imageio.get_writer(output_filename, mode="I", duration=0.1, loop=0) as writer:
        for filename in image_files:
            image = imageio.imread(filename)
            writer.append_data(image)
    print(f"\nFilme salvo como '{output_filename}'")


if __name__ == "__main__":
    # --- ARQUIVOS DE ENTRADA ---
    # Altere DATA_FILE se o nome do seu arquivo for diferente
    DATA_FILE = "N12_Lmat4.data"
    ATOM_DUMP_PATTERN = "dump.atoms.step.*.lammpstrj"

    # 1. Carrega o conjunto inicial de ligações inquebráveis (tipo 1)
    print(f"Lendo ligações inquebráveis de '{DATA_FILE}'...")
    unbreakable_bonds_set = parse_initial_unbreakable_bonds(DATA_FILE)
    print(f"Encontradas {len(unbreakable_bonds_set)} ligações inquebráveis.")

    # 2. Ordena numericamente os arquivos de dump de átomos
    atom_dump_files = sorted(
        glob.glob(ATOM_DUMP_PATTERN),
        key=lambda f: int(os.path.basename(f).split(".")[3]),
    )
    if not atom_dump_files:
        print(f"Erro: Nenhum arquivo encontrado para o padrão '{ATOM_DUMP_PATTERN}'.")
        exit()

    # 3. Pré-carrega posições e calcula a caixa de plotagem global
    print("Lendo todos os passos para determinar a caixa de plotagem...")
    all_steps_positions = [read_atom_dump(df) for df in atom_dump_files]
    all_x = [p[0] for s in all_steps_positions for p in s.values()]
    all_y = [p[1] for s in all_steps_positions for p in s.values()]
    if not all_x or not all_y:
        print("Erro: Nenhuma posição de átomo encontrada.")
        exit()
    min_x, max_x = min(all_x), max(all_x)
    min_y, max_y = min(all_y), max(all_y)
    x_pad = (max_x - min_x) * 0.05 or 1
    y_pad = (max_y - min_y) * 0.05 or 1
    global_box_dims = {
        "xlo": min_x - x_pad,
        "xhi": max_x + x_pad,
        "ylo": min_y - y_pad,
        "yhi": max_y + y_pad,
    }
    print("Caixa de plotagem global calculada.")

    # 4. Gera um quadro para cada passo
    frame_filenames = []
    print("Gerando quadros para o filme...")
    for i, atom_positions in enumerate(all_steps_positions):
        current_atom_file = atom_dump_files[i]
        step_num = int(os.path.basename(current_atom_file).split(".")[3])
        bond_file = f"dump.bonds.step.{step_num}.txt"

        active_bonds_for_frame = read_bond_dump(bond_file)

        # Passa o conjunto de ligações inquebráveis para a função de plotagem
        filename = visualize_frame(
            step_num,
            atom_positions,
            active_bonds_for_frame,
            unbreakable_bonds_set,  # <-- PASSA A INFORMAÇÃO AQUI
            global_box_dims,
        )
        frame_filenames.append(filename)
        print(f"   ... salvo {filename}")

    # 5. Cria o filme
    create_movie(frame_filenames, "simulation_movie.gif")
