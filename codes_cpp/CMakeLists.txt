cmake_minimum_required(VERSION 3.10)
project(spring_network_cpp CXX)  # Removed MPI from project()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Explicit MPI detection
find_package(MPI REQUIRED)

# Find pre-built LAMMPS library




add_executable(spring_network_cpp main.cpp)

# --- Include directories ---
target_include_directories(spring_network_cpp PRIVATE 
    /home/michael/gitrepos/md-minimizer/codes_cpp/lammps-stable_29Aug2024_update3/build/includes/lammps # Explicit path for lammps.h
    ${MPI_CXX_INCLUDE_DIRS}  # MPI headers
)

# --- Link libraries ---
target_link_libraries(spring_network_cpp PRIVATE
    /home/michael/gitrepos/md-minimizer/codes_cpp/lammps-stable_29Aug2024_update3/build/liblammps.so  # Explicitly link liblammps.so
    ${MPI_CXX_LIBRARIES}  # MPI libraries
)
