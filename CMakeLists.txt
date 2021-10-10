cmake_minimum_required(VERSION 3.16)
project(sparse-boolean-matrix-multiplication LANGUAGES CUDA CXX)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CXX_STANDARD 17)

if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
    set(CMAKE_CUDA_ARCHITECTURES 52 61 75)
endif()

find_package(PkgConfig)

option(ENABLE_BOOST "Enable Boost support" ON)
option(ENABLE_MMIO "Enable Matrix Market matrices support" OFF)
option(ENABLE_RANDOM_MATRICES "Enable Random Matrices support" ON)

if(ENABLE_BOOST)
    include(cmake/FindBoost.cmake)
    if(BOOST_FOUND)
        link_libraries(${BOOST_LIBRARIES})
        set(USE_BOOST 1)
    endif()
endif()

if(ENABLE_MMIO)
    add_library(mmio OBJECT lib/mmio.h lib/mmio.cpp)
    link_libraries(mmio)
    set(USE_MMIO_MATRICES 1)
endif()

if(ENABLE_RANDOM_MATRICES)
    set(USE_RANDOM_MATRICES 1)
endif()

configure_file(lib/definitions.h.in lib/definitions.h)

add_library(coo2csc OBJECT lib/coo.cpp lib/coo.hpp)
add_library(bmm_lib OBJECT lib/bmm_lib.h lib/bmm_lib.cpp)

add_executable(sbmm main.cu)
add_executable(sbmm_seq src/sbmm_sequential.cpp)
target_link_libraries(sbmm PUBLIC bmm_lib coo2csc)
target_link_libraries(sbmm_seq PUBLIC bmm_lib coo2csc)
set_target_properties(
        sbmm
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
