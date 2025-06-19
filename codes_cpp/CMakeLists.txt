cmake_minimum_required(VERSION 3.10)
project(spring_network_cpp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Explicit MPI detection
find_package(MPI REQUIRED)

add_executable(spring_network_cpp main.cpp)

# --- Define a macro para compilar com suporte a MPI do LAMMPS ---
# Isso garante que as funções corretas (como lammps_open) fiquem visíveis no library.h
target_compile_definitions(spring_network_cpp PRIVATE LAMMPS_LIB_MPI)

# --- Include directories ---
target_include_directories(spring_network_cpp PRIVATE 
    # Use o caminho absoluto para o diretório de includes do LAMMPS
    /home/michael/gitrepos/md-minimizer/codes_cpp/lammps-stable_29Aug2024_update3/build/includes/lammps
    ${MPI_CXX_INCLUDE_DIRS}  # Headers do MPI
)

# --- Link libraries ---
target_link_libraries(spring_network_cpp PRIVATE
    # Use o caminho absoluto para a biblioteca compartilhada do LAMMPS
    /home/michael/gitrepos/md-minimizer/codes_cpp/lammps-stable_29Aug2024_update3/build/liblammps.so
    ${MPI_CXX_LIBRARIES}  # Bibliotecas do MPI
)