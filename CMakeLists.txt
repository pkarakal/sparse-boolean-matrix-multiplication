cmake_minimum_required(VERSION 3.16)
project(sparse-boolean-matrix-multiplication LANGUAGES CUDA CXX)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_STANDARD 17)

if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
    set(CMAKE_CUDA_ARCHITECTURES 52 61 75)
endif()

option(ENABLE_MMIO "Enable Matrix Market matrices support" ON)

if(ENABLE_MMIO)
    add_library(mmio OBJECT lib/mmio.h lib/mmio.cpp)
    link_libraries(mmio)
    set(USE_MMIO_MATRICES 1)
endif()

configure_file(lib/definitions.h.in lib/definitions.h)

add_library(coo2csc OBJECT lib/coo.cpp lib/coo.hpp)

add_executable(sbmm main.cu)
target_link_libraries(sbmm PUBLIC coo2csc)
set_target_properties(
        sbmm
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
