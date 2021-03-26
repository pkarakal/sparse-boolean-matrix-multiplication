cmake_minimum_required(VERSION 3.16)
project(sparse-boolean-matrix-multiplication LANGUAGES CUDA CXX)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_STANDARD 17)

if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
    set(CMAKE_CUDA_ARCHITECTURES 52 61 75)
endif()

add_executable(sbmm main.cu)
set_target_properties(
        sbmm
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
