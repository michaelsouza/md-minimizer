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
// cmake -DLAMMPS_LIB_PATH=/path/to/your/liblammps.so -DLAMMPS_INC_PATH=/path/to/your/lammps/src ..
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
#include "lammps.h"
#include "library.h"
#include "lmptype.h"
#include "atom.h"

extern "C" {
void *lammps_open(int, char**, MPI_Comm, void**);
void lammps_close(void*);
char *lammps_command(void*, const char*);
void lammps_commands_list(void*, int, const char**);
void *lammps_open_no_mpi(int, char**, void**);
}

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
    // --- Argumentos da Linha de Comando (versão simplificada) ---
    if (argc < 3) {
        std::cerr << "Uso: " << argv[0] << " <data_file> <thresholds_file> [total_steps] [strain_inc]" << std::endl;
        return 1;
    }
    std::string data_file = argv[1];
    std::string thresholds_file = argv[2];
    int total_steps = (argc > 3) ? std::stoi(argv[3]) : 10;
    double strain_inc = (argc > 4) ? std::stod(argv[4]) : 0.1;
    
    // --- Carregar Limiares de Quebra ---
    auto thresholds = parse_thresholds(thresholds_file);

    // --- Inicialização do LAMMPS ---
    // Initialize LAMMPS
    void *lammps_ptr = nullptr;
    void *lammps = lammps_open_no_mpi(0, nullptr, &lammps_ptr);
    if (!lammps) {
        throw std::runtime_error("Failed to initialize LAMMPS");
    }

    // Os argumentos da linha de comando para o LAMMPS podem ser passados aqui
    // Ex: "-log", "none", "-screen", "none"
    char *lmp_argv[] = {(char*)"-log", (char*)"none", (char*)"-screen", (char*)"none"};
    int lmp_argc = sizeof(lmp_argv)/sizeof(char*);

    // Individual LAMMPS commands with error checking
    char* error_msg;
    error_msg = lammps_command(lammps, "units lj");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "dimension 2");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "boundary p s p");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "atom_style bond");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "bond_style harmonic");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "pair_style none");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, ("read_data " + data_file).c_str());
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "group bottom_atoms type 2");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "group top_atoms type 3");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "group mobile_atoms type 1");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "fix 1 bottom_atoms setforce 0.0 0.0 0.0");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "thermo 1");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }
    error_msg = lammps_command(lammps, "thermo_style custom step pe press pyy");
    if (error_msg != nullptr) {
        std::cerr << "LAMMPS error: " << error_msg << std::endl;
        lammps_close(lammps);
        return 1;
    }

    // --- Loop Principal de Deformação ---
    long long num_broken_total = 0;
    for (int step_id = 0; step_id < total_steps; ++step_id) {
        auto step_start_time = std::chrono::high_resolution_clock::now();
        std::cout << "--- Strain Step " << step_id + 1 << "/" << total_steps << " ---" << std::endl;

        // Aplica o deslocamento (incremento de deformação)
        std::string displace_cmd = "displace_atoms top_atoms move 0 " + std::to_string(strain_inc) + " 0";
        lammps_command(lammps, displace_cmd.c_str());

        // Fixa os átomos do topo durante a relaxação/avalanche
        lammps_command(lammps, "fix 2 top_atoms setforce 0.0 0.0 0.0");

        // --- Loop da Avalanche (Otimizado) ---
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

            // Ponteiros para as estruturas de dados internas
            double **x = (double **)lammps_extract_atom(lammps, "x");
            int *bond_type = (int *)lammps_extract_atom(lammps, "bond_type");
            int **bond_atom = (int **)lammps_extract_atom(lammps, "bond_atom");
            tagint *tag = (tagint *)lammps_extract_atom(lammps, "tag");
            
            int nlocal = static_cast<int>(lammps_get_natoms(lammps));
            int nbonds = *(int *)lammps_extract_global(lammps, "nbond");

            // Mapeamento de tag global para índice local do átomo
            std::map<tagint, int> tag_to_local_idx;
            for(int i = 0; i < nlocal; ++i) {
                tag_to_local_idx[tag[i]] = i;
            }

            // Dimensões da caixa para tratar condições de contorno periódicas
            double boxlo[3], boxhi[3];
            lammps_extract_box(lammps, boxlo, boxhi, NULL, NULL, NULL, NULL, NULL);
            double x_period = boxhi[0] - boxlo[0];

            for (int i = 0; i < nbonds; ++i) {
                int current_type = bond_type[i];

                // Ignora ligações inquebráveis (tipo 1) ou já quebradas (tipo 0)
                if (current_type <= 1) continue;

                // Verifica se o tipo da ligação tem um limiar de quebra definido
                if (thresholds.find(current_type) == thresholds.end()) continue;

                double break_len = thresholds.at(current_type);

                tagint atom1_tag = bond_atom[i][0];
                tagint atom2_tag = bond_atom[i][1];

                int idx1 = tag_to_local_idx[atom1_tag];
                int idx2 = tag_to_local_idx[atom2_tag];
                
                double dx = x[idx1][0] - x[idx2][0];
                double dy = x[idx1][1] - x[idx2][1];
                
                // Aplica a condição de contorno periódica em x (Minimum Image Convention)
                dx = dx - x_period * round(dx / x_period);

                double dist_sq = dx * dx + dy * dy;
                
                if (sqrt(dist_sq) > break_len) {
                    bond_type[i] = 0;
                    broken_this_iter++;
                }
            }
            auto access_end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> access_duration = access_end_time - access_start_time;
            std::cout << "   time (breakage): " << access_duration.count() << " s" << std::endl;

            num_broken_total += broken_this_iter;
            std::cout << "   Avalanche iteration broke " << broken_this_iter << " bonds." << std::endl;

            if (broken_this_iter == 0) {
                break; // Fim da avalanche, sistema está estável
            }
        }

        // Libera os átomos do topo para o próximo passo de deslocamento
        lammps_command(lammps, "unfix 2");
        
        auto step_end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> step_duration = step_end_time - step_start_time;

        std::cout << "Finished strain step " << step_id + 1 << "; cumulative broken = " << num_broken_total << std::endl;
        std::cout << "Total time for step: " << step_duration.count() << " s\n" << std::endl;
    }

    // --- Finalização ---
    lammps_close(lammps);
    std::cout << "Simulação finalizada." << std::endl;

    return 0;
}