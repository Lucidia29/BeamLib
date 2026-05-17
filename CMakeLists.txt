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

# Make CTest portable across shells (Codex Round 4, issue 2):
# Windows + MinGW links against libgcc_s_seh-1.dll / libstdc++-6.dll. When
# CTest spawns the test executables from a shell that does not have the
# MinGW bin directory on PATH, Windows aborts with 0xc0000135 (missing
# runtime). Copying the needed DLLs next to the test binaries puts them
# in the same-directory step of the Windows DLL search order, so the
# tests run regardless of the caller's PATH.
if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    get_filename_component(BEAMLIB_MINGW_BIN "${CMAKE_CXX_COMPILER}" DIRECTORY)
    foreach(_dll libstdc++-6.dll libgcc_s_seh-1.dll libwinpthread-1.dll)
        if(EXISTS "${BEAMLIB_MINGW_BIN}/${_dll}")
            file(COPY "${BEAMLIB_MINGW_BIN}/${_dll}"
                 DESTINATION "${CMAKE_BINARY_DIR}")
        endif()
    endforeach()
endif()

# Tests
enable_testing()

add_executable(test_compile_check tests/test_compile_check.cpp)
target_link_libraries(test_compile_check beamlib)
add_test(NAME test_compile_check COMMAND test_compile_check)

add_executable(test_eb2d_cantilever tests/test_eb2d_cantilever.cpp)
target_link_libraries(test_eb2d_cantilever beamlib)
add_test(NAME test_eb2d_cantilever COMMAND test_eb2d_cantilever)

add_executable(test_eb2d_patch tests/test_eb2d_patch.cpp)
target_link_libraries(test_eb2d_patch beamlib)
add_test(NAME test_eb2d_patch COMMAND test_eb2d_patch)

add_executable(test_eb2d_moment tests/test_eb2d_moment.cpp)
target_link_libraries(test_eb2d_moment beamlib)
add_test(NAME test_eb2d_moment COMMAND test_eb2d_moment)

add_executable(test_eb2d_postprocess tests/test_eb2d_postprocess.cpp)
target_link_libraries(test_eb2d_postprocess beamlib)
add_test(NAME test_eb2d_postprocess COMMAND test_eb2d_postprocess)

add_executable(test_eb2d_simply_supported tests/test_eb2d_simply_supported.cpp)
target_link_libraries(test_eb2d_simply_supported beamlib)
add_test(NAME test_eb2d_simply_supported COMMAND test_eb2d_simply_supported)

add_executable(test_eb2d_rotated_beam tests/test_eb2d_rotated_beam.cpp)
target_link_libraries(test_eb2d_rotated_beam beamlib)
add_test(NAME test_eb2d_rotated_beam COMMAND test_eb2d_rotated_beam)

add_executable(test_eb2d_settlement tests/test_eb2d_settlement.cpp)
target_link_libraries(test_eb2d_settlement beamlib)
add_test(NAME test_eb2d_settlement COMMAND test_eb2d_settlement)

add_executable(test_eb2d_mass tests/test_eb2d_mass.cpp)
target_link_libraries(test_eb2d_mass beamlib)
add_test(NAME test_eb2d_mass COMMAND test_eb2d_mass)

# --- Batch 3: EB3D ---
add_executable(test_eb3d_cantilever tests/test_eb3d_cantilever.cpp)
target_link_libraries(test_eb3d_cantilever beamlib)
add_test(NAME test_eb3d_cantilever COMMAND test_eb3d_cantilever)

add_executable(test_eb3d_torsion tests/test_eb3d_torsion.cpp)
target_link_libraries(test_eb3d_torsion beamlib)
add_test(NAME test_eb3d_torsion COMMAND test_eb3d_torsion)

add_executable(test_eb3d_transform tests/test_eb3d_transform.cpp)
target_link_libraries(test_eb3d_transform beamlib)
add_test(NAME test_eb3d_transform COMMAND test_eb3d_transform)

add_executable(test_eb2d_eb3d_crosscheck tests/test_eb2d_eb3d_crosscheck.cpp)
target_link_libraries(test_eb2d_eb3d_crosscheck beamlib)
add_test(NAME test_eb2d_eb3d_crosscheck COMMAND test_eb2d_eb3d_crosscheck)

add_executable(test_eb3d_rigid_body tests/test_eb3d_rigid_body.cpp)
target_link_libraries(test_eb3d_rigid_body beamlib)
add_test(NAME test_eb3d_rigid_body COMMAND test_eb3d_rigid_body)

add_executable(test_eb3d_mass tests/test_eb3d_mass.cpp)
target_link_libraries(test_eb3d_mass beamlib)
add_test(NAME test_eb3d_mass COMMAND test_eb3d_mass)

# --- Batch 4: Timoshenko 2D ---
add_executable(test_timo2d_cantilever tests/test_timo2d_cantilever.cpp)
target_link_libraries(test_timo2d_cantilever beamlib)
add_test(NAME test_timo2d_cantilever COMMAND test_timo2d_cantilever)

add_executable(test_timo2d_deep_beam tests/test_timo2d_deep_beam.cpp)
target_link_libraries(test_timo2d_deep_beam beamlib)
add_test(NAME test_timo2d_deep_beam COMMAND test_timo2d_deep_beam)

add_executable(test_timo2d_slender_limit tests/test_timo2d_slender_limit.cpp)
target_link_libraries(test_timo2d_slender_limit beamlib)
add_test(NAME test_timo2d_slender_limit COMMAND test_timo2d_slender_limit)

add_executable(test_timo2d_shear_locking tests/test_timo2d_shear_locking.cpp)
target_link_libraries(test_timo2d_shear_locking beamlib)
add_test(NAME test_timo2d_shear_locking COMMAND test_timo2d_shear_locking)

add_executable(test_timo2d_phi_limit tests/test_timo2d_phi_limit.cpp)
target_link_libraries(test_timo2d_phi_limit beamlib)
add_test(NAME test_timo2d_phi_limit COMMAND test_timo2d_phi_limit)

add_executable(test_timo2d_mass tests/test_timo2d_mass.cpp)
target_link_libraries(test_timo2d_mass beamlib)
add_test(NAME test_timo2d_mass COMMAND test_timo2d_mass)
