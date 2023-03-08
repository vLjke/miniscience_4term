cmake_minimum_required(VERSION 3.24)
project(gmsh_lab)

set(CMAKE_CXX_STANDARD 20)

add_executable(lab_gmsh lab_gmsh.cpp)

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
        IMPORTED_LOCATION "/home/vljke/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so"
        INTERFACE_INCLUDE_DIRECTORIES "/home/vljke/gmsh-4.11.1-Linux64-sdk/include"
        )

set(GMSH_LIBRARY gmsh)

target_link_libraries(lab_gmsh ${GMSH_LIBRARY})