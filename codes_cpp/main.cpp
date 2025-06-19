// main.cpp
//
// Versão em C++ do script spring_network.py para simulação de avalanches
// em redes de molas.
//
// Otimizações:
// - O loop de avalanche acessa diretamente as estruturas de dados internas do LAMMPS
//   para verificar o comprimento das ligações e quebrá-las.
// - As ligações quebradas têm seu tipo alterado para 0, removendo-as
//   efetivamente da minimização de energia sem o custo de comandos de script.
//
// Compilação (usando CMake):
// mkdir build && cd build
// cmake ..
// make

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <map>
#include <cmath>
#include <chrono>
#include <mpi.h>
#include "lammps.h"     // Official LAMMPS header
#include "library.h"    // Provides library function prototypes
#include "lmptype.h"
#include "atom.h"

// The manual extern "C" block is removed.
// All function prototypes are now correctly sourced from library.h

using LAMMPS_NS::tagint;

// Função para parsear o arquivo de limiares de quebra
std::map<int, double> parse_thresholds(const std::string& filename) {
    std::map<int, double> thresholds;
    std::ifstream file(filename);

    if (!file.is_open()) {
        throw std::runtime_error("Erro: Não foi possível abrir o arquivo de limiares: " + filename);
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        int bond_type;
        double break_len;
        if (sscanf(line.c_str(), "%d %lf", &bond_type, &break_len) == 2) {
            thresholds[bond_type] = break_len;
        }
    }

    if (thresholds.empty()) {
        throw std::runtime_error("Nenhum limiar válido encontrado em " + filename);
    }
    
    std::cout << "Info: Limiares de quebra lidos para " << thresholds.size() << " tipos de ligação." << std::endl;
    return thresholds;
}


int main(int argc, char* argv[]) {
    // --- Initialize MPI ---
    MPI_Init(&argc, &argv);
    
    // --- Version Message ---
    std::cout << "md-minimizer C++ v1.8" << std::endl;

    // --- Argumentos da Linha de Comando ---
    if (argc < 4) {
        std::cerr << "Uso: " << argv[0] << " <config_file> <data_file> <thresholds_file> [total_steps] [strain_inc]" << std::endl;
        MPI_Finalize();
        return 1;
    }
    std::string config_file = argv[1];
    std::string data_file = argv[2];
    std::string thresholds_file = argv[3];
    int total_steps = (argc > 4) ? std::stoi(argv[4]) : 10;
    double strain_inc = (argc > 5) ? std::stod(argv[5]) : 0.1;
    
    // --- Carregar Limiares de Quebra ---
    auto thresholds = parse_thresholds(thresholds_file);

    // --- Inicialização do LAMMPS ---
    void *lammps = lammps_open(0, NULL, MPI_COMM_WORLD, NULL);
    if (!lammps) {
        MPI_Finalize();
        throw std::runtime_error("Failed to initialize LAMMPS");
    }

    // --- Carregar Configuração Estática do Arquivo ---
    std::string setup_cmds = "variable data_file string " + data_file + "\n" +
                             "include " + config_file;
    lammps_commands_string(lammps, setup_cmds.c_str());

    // --- Loop Principal de Deformação (Lógica Dinâmica) ---
    long long num_broken_total = 0;
    for (int step_id = 0; step_id < total_steps; ++step_id) {
        auto step_start_time = std::chrono::high_resolution_clock::now();
        std::cout << "--- Strain Step " << step_id + 1 << "/" << total_steps << " ---" << std::endl;

        // Aplica o deslocamento
        std::string displace_cmd = "displace_atoms top_atoms move 0 " + std::to_string(strain_inc) + " 0";
        lammps_command(lammps, displace_cmd.c_str());

        // Fixa os átomos do topo durante a relaxação
        lammps_command(lammps, "fix 2 top_atoms setforce 0.0 0.0 0.0");

        // --- Loop da Avalanche ---
        while (true) {
            auto minimize_start_time = std::chrono::high_resolution_clock::now();
            lammps_command(lammps, "min_style cg");
            lammps_command(lammps, "minimize 1.0e-5 1.0e-7 1000 10000");
            auto minimize_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> minimize_duration = minimize_end_time - minimize_start_time;
            std::cout << "   time (minimize): " << minimize_duration.count() << " s" << std::endl;

            int broken_this_iter = 0;

            // --- Acesso Direto aos Dados do LAMMPS ---
            auto access_start_time = std::chrono::high_resolution_clock::now();
            
            // Per-atom data
            double **x = (double **)lammps_extract_atom(lammps, "x");
            tagint *tag = (tagint *)lammps_extract_atom(lammps, "tag");
            
            // Global data
            int nbonds = *(int *)lammps_extract_global(lammps, "nbonds");
            int *bond_type = (int *)lammps_extract_global(lammps, "bond_type");
            tagint **bond_atom = (tagint **)lammps_extract_global(lammps, "bond_atom");

            if (!x || !tag || !bond_type || !bond_atom || !(&nbonds)) {
                std::cerr << "Error: Failed to extract required data pointers from LAMMPS." << std::endl;
                break;
            }
            
            // lammps_get_natoms() returns double, so we cast it.
            int nlocal = static_cast<int>(lammps_get_natoms(lammps));

            std::map<tagint, int> tag_to_local_idx;
            for(int i = 0; i < nlocal; ++i) {
                tag_to_local_idx[tag[i]] = i;
            }

            double boxlo[3], boxhi[3];
            lammps_extract_box(lammps, boxlo, boxhi, NULL, NULL, NULL, NULL, NULL);
            double x_period = boxhi[0] - boxlo[0];

            for (int i = 0; i < nbonds; ++i) {
                int current_type = bond_type[i];
                if (current_type <= 1) continue;
                if (thresholds.find(current_type) == thresholds.end()) continue;

                double break_len = thresholds.at(current_type);
                tagint atom1_tag = bond_atom[i][0];
                tagint atom2_tag = bond_atom[i][1];
                
                if (tag_to_local_idx.find(atom1_tag) == tag_to_local_idx.end() ||
                    tag_to_local_idx.find(atom2_tag) == tag_to_local_idx.end()) continue;

                int idx1 = tag_to_local_idx[atom1_tag];
                int idx2 = tag_to_local_idx[atom2_tag];
                
                double dx = x[idx1][0] - x[idx2][0];
                double dy = x[idx1][1] - x[idx2][1];
                dx = dx - x_period * round(dx / x_period);
                double dist_sq = dx * dx + dy * dy;

                if (sqrt(dist_sq) > break_len) {
                    bond_type[i] = 0; // Set bond type to 0 to "break" it
                    broken_this_iter++;
                }
            }
            auto access_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> access_duration = access_end_time - access_start_time;
            std::cout << "   time (breakage): " << access_duration.count() << " s" << std::endl;

            num_broken_total += broken_this_iter;
            std::cout << "   Avalanche iteration broke " << broken_this_iter << " bonds." << std::endl;

            if (broken_this_iter == 0) {
                break;
            }
        }

        lammps_command(lammps, "unfix 2");
        
        auto step_end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> step_duration = step_end_time - step_start_time;

        std::cout << "Finished strain step " << step_id + 1 << "; cumulative broken = " << num_broken_total << std::endl;
        std::cout << "Total time for step: " << step_duration.count() << " s\n" << std::endl;
    }

    // --- Finalização ---
    lammps_close(lammps);
    MPI_Finalize();
    std::cout << "Simulação finalizada." << std::endl;

    return 0;
}