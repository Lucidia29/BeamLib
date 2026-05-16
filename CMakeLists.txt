cmake_minimum_required(VERSION 3.16)
project(BeamLib LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find Eigen (header-only)
# Supports: installed Eigen3 (via find_package), or raw source tree (via EIGEN3_INCLUDE_DIR).
if(DEFINED EIGEN3_INCLUDE_DIR)
    if(NOT EXISTS "${EIGEN3_INCLUDE_DIR}/Eigen/Dense")
        message(FATAL_ERROR "EIGEN3_INCLUDE_DIR=${EIGEN3_INCLUDE_DIR} does not contain Eigen/Dense")
    endif()
    message(STATUS "Using Eigen from EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")
    add_library(Eigen3_target INTERFACE)
    target_include_directories(Eigen3_target INTERFACE "${EIGEN3_INCLUDE_DIR}")
    add_library(Eigen3::Eigen ALIAS Eigen3_target)
else()
    find_package(Eigen3 QUIET NO_MODULE)
    if(NOT Eigen3_FOUND)
        message(FATAL_ERROR
            "Eigen3 not found. Set -DEIGEN3_INCLUDE_DIR=<path-to-eigen-source> "
            "(directory containing Eigen/ subdirectory)")
    endif()
    message(STATUS "Using installed Eigen3: ${Eigen3_DIR}")
endif()

# BeamLib header-only INTERFACE library (becomes STATIC when .cpp files are added in Batch 6)
add_library(beamlib INTERFACE)
target_include_directories(beamlib INTERFACE "${CMAKE_SOURCE_DIR}/include")
target_link_libraries(beamlib INTERFACE Eigen3::Eigen)

# Tests
enable_testing()

add_executable(test_compile_check tests/test_compile_check.cpp)
target_link_libraries(test_compile_check beamlib)
add_test(NAME test_compile_check COMMAND test_compile_check)
