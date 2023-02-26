cmake_minimum_required(VERSION 3.24)
project(science_programming_4term)

set(CMAKE_CXX_STANDARD 20)

add_executable(Lab1 Lab1.cpp)

add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
        IMPORTED_LOCATION "/home/vljke/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so"
        INTERFACE_INCLUDE_DIRECTORIES "/home/vljke/gmsh-4.11.1-Linux64-sdk/include"
        )

set(GMSH_LIBRARY gmsh)

target_link_libraries(Lab1 ${GMSH_LIBRARY})


