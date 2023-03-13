cmake_minimum_required(VERSION 3.24)
project(vtk_lab)

set(CMAKE_CXX_STANDARD 20)

add_executable(lab_vtk lab_vtk.cpp)
# add gmsh lib
add_library(gmsh SHARED IMPORTED)
set_target_properties(gmsh PROPERTIES
        IMPORTED_LOCATION "/home/vljke/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so"
        INTERFACE_INCLUDE_DIRECTORIES "/home/vljke/gmsh-4.11.1-Linux64-sdk/include"
        )

set(GMSH_LIBRARY gmsh)
# add vtk lib
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
# link vtk and gmsh
target_link_libraries(lab_vtk ${VTK_LIBRARIES})
target_link_libraries(lab_vtk ${GMSH_LIBRARY})