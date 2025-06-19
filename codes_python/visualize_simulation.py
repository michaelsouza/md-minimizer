"""
visualize_simulation_parallel.py

Este script lê uma série de arquivos de dump do LAMMPS, plota o estado da
rede para cada passo e compila as imagens em um filme. As ligações
inquebráveis (tipo 1) são coloridas de vermelho.
Utiliza multiprocessing para acelerar a leitura de arquivos e geração de frames.
"""

import glob
import os
from typing import Dict, List, Tuple, Set
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
import multiprocessing as mp

import imageio.v2 as imageio
import matplotlib.pyplot as plt
import matplotlib

# Use non-interactive backend for parallel processing
matplotlib.use("Agg")


def parse_initial_unbreakable_bonds(data_filename: str) -> Set[Tuple[int, int]]:
    """
    Lê o arquivo de dados inicial do LAMMPS para encontrar todas as
    ligações inquebráveis (tipo 1).
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
                if line.startswith("Atoms"):
                    reading_bonds = False

                if reading_bonds:
                    parts = line.split()
                    if len(parts) >= 4:
                        bond_type = int(parts[1])
                        if bond_type == 1:
                            p1, p2 = int(parts[2]), int(parts[3])
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


def process_single_step(
    args: Tuple,
) -> Tuple[int, Dict[int, Tuple[float, float]], List[Tuple[int, int]]]:
    """
    Processa um único passo: lê tanto átomos quanto ligações.
    Esta função é usada para paralelização.
    """
    atom_file, step_num = args

    # Lê posições dos átomos
    atom_positions = read_atom_dump(atom_file)

    # Lê ligações correspondentes
    bond_file = f"dump.bonds.step.{step_num}.txt"
    active_bonds = read_bond_dump(bond_file)

    return step_num, atom_positions, active_bonds


def visualize_frame_parallel(args: Tuple) -> str:
    """
    Versão da função visualize_frame adaptada para paralelização.
    Recebe todos os argumentos em uma tupla.
    """
    step_num, positions, active_bonds, unbreakable_bonds, box = args

    plt.style.use("seaborn-v0_8-whitegrid")
    _, ax = plt.subplots(figsize=(8, 10))

    for p1_id, p2_id in active_bonds:
        if p1_id in positions and p2_id in positions:
            pos1 = positions[p1_id]
            pos2 = positions[p2_id]

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
    plt.savefig(output_filename, dpi=100, bbox_inches="tight")
    plt.close()
    return output_filename


def create_movie(image_files: List[str], output_filename: str):
    """Cria um GIF animado a partir de uma lista de arquivos de imagem."""
    with imageio.get_writer(output_filename, mode="I", duration=0.1, loop=0) as writer:
        for filename in image_files:
            image = imageio.imread(filename)
            writer.append_data(image)
    print(f"\nFilme salvo como '{output_filename}'")


def calculate_global_box(
    all_positions: List[Dict[int, Tuple[float, float]]]
) -> Dict[str, float]:
    """Calcula as dimensões globais da caixa de plotagem."""
    all_x = [p[0] for positions in all_positions for p in positions.values()]
    all_y = [p[1] for positions in all_positions for p in positions.values()]

    if not all_x or not all_y:
        raise ValueError("Nenhuma posição de átomo encontrada.")

    min_x, max_x = min(all_x), max(all_x)
    min_y, max_y = min(all_y), max(all_y)
    x_pad = (max_x - min_x) * 0.05 or 1
    y_pad = (max_y - min_y) * 0.05 or 1

    return {
        "xlo": min_x - x_pad,
        "xhi": max_x + x_pad,
        "ylo": min_y - y_pad,
        "yhi": max_y + y_pad,
    }


if __name__ == "__main__":
    # Configuração do número de workers (ajuste conforme seu hardware)
    n_workers = min(mp.cpu_count(), 24)  # Limita a 8 para evitar sobrecarga

    print(f"Usando {n_workers} workers para paralelização")

    # --- ARQUIVOS DE ENTRADA ---
    DATA_FILE = "N12_Lmat4.data"
    ATOM_DUMP_PATTERN = "dump.atoms.step.*.lammpstrj"

    # 1. Carrega o conjunto inicial de ligações inquebráveis
    print(f"Lendo ligações inquebráveis de '{DATA_FILE}'...")
    unbreakable_bonds_set = parse_initial_unbreakable_bonds(DATA_FILE)
    print(f"Encontradas {len(unbreakable_bonds_set)} ligações inquebráveis.")

    # 2. Ordena arquivos de dump
    atom_dump_files = sorted(
        glob.glob(ATOM_DUMP_PATTERN),
        key=lambda f: int(os.path.basename(f).split(".")[3]),
    )
    if not atom_dump_files:
        print(f"Erro: Nenhum arquivo encontrado para o padrão '{ATOM_DUMP_PATTERN}'.")
        exit()

    print(f"Encontrados {len(atom_dump_files)} arquivos de dump.")

    # 3. PARALELIZAÇÃO: Lê todos os passos simultaneamente
    print("Lendo todos os passos em paralelo...")
    step_args = []
    for atom_file in atom_dump_files:
        step_num = int(os.path.basename(atom_file).split(".")[3])
        step_args.append((atom_file, step_num))

    # Usa ThreadPoolExecutor para I/O (leitura de arquivos)
    step_data = {}
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        future_to_step = {
            executor.submit(process_single_step, args): args[1] for args in step_args
        }

        for future in as_completed(future_to_step):
            step_num = future_to_step[future]
            try:
                step_num, positions, bonds = future.result()
                step_data[step_num] = (positions, bonds)
                print(f"   ... passo {step_num} carregado")
            except Exception as exc:
                print(f"Erro ao processar passo {step_num}: {exc}")

    # 4. Calcula caixa global
    print("Calculando caixa de plotagem global...")
    all_positions = [data[0] for data in step_data.values()]
    global_box_dims = calculate_global_box(all_positions)
    print("Caixa de plotagem global calculada.")

    # 5. PARALELIZAÇÃO: Gera frames simultaneamente
    print("Gerando quadros em paralelo...")
    frame_args = []
    for step_num in sorted(step_data.keys()):
        positions, bonds = step_data[step_num]
        frame_args.append(
            (step_num, positions, bonds, unbreakable_bonds_set, global_box_dims)
        )

    # Usa ProcessPoolExecutor para CPU-intensive tasks (geração de plots)
    frame_results = {}
    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        future_to_step = {
            executor.submit(visualize_frame_parallel, args): args[0]
            for args in frame_args
        }

        for future in as_completed(future_to_step):
            step_num = future_to_step[future]
            try:
                filename = future.result()
                frame_results[step_num] = filename
                print(f"   ... frame {step_num} gerado: {filename}")
            except Exception as exc:
                print(f"Erro ao gerar frame {step_num}: {exc}")

    # 6. Ordena os arquivos de frame e cria o filme
    frame_filenames = [frame_results[step] for step in sorted(frame_results.keys())]
    create_movie(frame_filenames, "simulation_movie.gif")

    print(f"\nProcessamento concluído! {len(frame_filenames)} frames gerados.")
    print("Filme salvo como 'simulation_movie.gif'")
