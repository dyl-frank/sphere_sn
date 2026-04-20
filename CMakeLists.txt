# ===========================================================================
#  1D Spherical Sn Transport Code
#  CMake Build Configuration
# ===========================================================================
cmake_minimum_required(VERSION 3.14)

project(SphericalSnTransport
    VERSION   0.1.0
    LANGUAGES CXX
)

# ── C++ Standard ──────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Build type default ────────────────────────────────────────────────────
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# ── Compiler warnings ─────────────────────────────────────────────────────
add_compile_options(-Wall -Wextra -Wpedantic)

# ===========================================================================
#  Dependencies
# ===========================================================================

# ── Eigen3 ────────────────────────────────────────────────────────────────
# Install: spack install eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")

# ── Axom ──────────────────────────────────────────────────────────────────
# Install: spack install axom
find_package(axom REQUIRED)
message(STATUS "Found Axom: ${axom_DIR}")

# ── MFEM (Spack-compatible detection) ─────────────────────────────────────
if(NOT MFEM_DIR)
    message(FATAL_ERROR "MFEM_DIR must point to the MFEM install prefix (e.g. from spack location -i mfem)")
endif()

set(MFEM_INCLUDE_DIR "${MFEM_DIR}/include")
set(MFEM_LIBRARY_DIR "${MFEM_DIR}/lib")

find_library(MFEM_LIBRARY
    NAMES mfem
    PATHS ${MFEM_LIBRARY_DIR}
    REQUIRED
)
if(NOT MFEM_LIBRARY)
    message(FATAL_ERROR "Could not find libmfem in ${MFEM_LIBRARY_DIR}")
endif()

# ── METIS & zlib (dependencies of MFEM static library) ────────────────
# Look in Spack prefixes
execute_process(
    COMMAND spack location -i metis
    OUTPUT_VARIABLE METIS_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_library(METIS_LIBRARY
    NAMES metis
    PATHS "${METIS_PREFIX}/lib"
    REQUIRED
)

execute_process(
    COMMAND spack location -i zlib
    OUTPUT_VARIABLE ZLIB_PREFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_library(ZLIB_LIBRARY
    NAMES z
    PATHS "${ZLIB_PREFIX}/lib"
    REQUIRED
)

message(STATUS "MFEM include: ${MFEM_INCLUDE_DIR}")
message(STATUS "MFEM library: ${MFEM_LIBRARY}")
message(STATUS "METIS library: ${METIS_LIBRARY}")
message(STATUS "ZLIB library: ${ZLIB_LIBRARY}")

# ── Main executable ─────────────────────────────────────────────────
add_executable(sphere_sn
    src/main.cpp
    src/Functions.cpp
)

target_include_directories(sphere_sn PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${MFEM_INCLUDE_DIR}
)

target_link_libraries(sphere_sn PRIVATE
    Eigen3::Eigen
    axom::inlet
    ${MFEM_LIBRARY}  # MFEM must come before its dependencies
    ${METIS_LIBRARY}
    ${ZLIB_LIBRARY}
)
# ===========================================================================
#  Install (optional)
# ===========================================================================

install(TARGETS sphere_sn DESTINATION bin)