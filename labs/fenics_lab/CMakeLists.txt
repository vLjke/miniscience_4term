cmake_minimum_required(VERSION 3.24)
project(fenics_lab)

set(CMAKE_CXX_STANDARD 20)

# add dolphin lib
find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})
# link dolphin and mshr
add_executable(lab_fenics lab_fenics.cpp)
target_link_libraries(lab_fenics dolfin mshr)

add_executable(poisson examples/Poisson.cpp)
target_link_libraries(poisson dolfin)

add_executable(hyper_elasticity examples/HyperElasticity.cpp)
target_link_libraries(hyper_elasticity dolfin)