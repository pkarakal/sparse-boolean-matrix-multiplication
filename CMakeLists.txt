cmake_minimum_required(VERSION 3.16)
project(sparse-boolean-matrix-multiplication LANGUAGES CUDA CXX)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_STANDARD 17)

if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
    set(CMAKE_CUDA_ARCHITECTURES 52 61 75)
endif()

configure_file(lib/definitions.h.in lib/definitions.h)

add_library(bmm_lib OBJECT lib/coo.cpp lib/coo.hpp)

add_executable(sbmm main.cu)
target_link_libraries(sbmm bmm_lib)
set_target_properties(
        sbmm
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
