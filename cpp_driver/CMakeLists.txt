cmake_minimum_required(VERSION 3.15)
project(les_apriori LANGUAGES CXX)

# -----------------------------------------------------------------------
# C++17 standard — required for std::filesystem, structured bindings, etc.
# -----------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------
# Optimization flags
# -----------------------------------------------------------------------
add_compile_options(-O3 -march=native -ffast-math)

# -----------------------------------------------------------------------
# Find FFTW3 (double precision) and FFTW3f (single precision)
# This installation exposes variables, not imported targets.
# -----------------------------------------------------------------------
find_package(FFTW3  REQUIRED CONFIG
    PATHS /opt/homebrew/lib/cmake/fftw3
          /opt/homebrew/Cellar/fftw/3.3.10_3/lib/cmake/fftw3)

find_package(FFTW3f REQUIRED CONFIG
    PATHS /opt/homebrew/lib/cmake/fftw3
          /opt/homebrew/Cellar/fftw/3.3.10_3/lib/cmake/fftw3)

# Compile-time check: verify FFTW3 double-precision library and header are present.
# We use find_library + find_file rather than check_cxx_source_compiles to avoid
# linker path issues with the variable-based (non-imported-target) FFTW3 config.
find_library(FFTW3_DOUBLE_LIB
    NAMES fftw3
    HINTS ${FFTW3_LIBRARY_DIRS}
    REQUIRED)

find_file(FFTW3_HEADER
    NAMES fftw3.h
    HINTS ${FFTW3_INCLUDE_DIRS}
    REQUIRED)

if(NOT FFTW3_DOUBLE_LIB OR NOT FFTW3_HEADER)
    message(FATAL_ERROR
        "FFTW3 double-precision library (libfftw3) or header (fftw3.h) not found. "
        "On macOS: brew install fftw")
endif()

message(STATUS "FFTW3 double precision: OK")
message(STATUS "  library : ${FFTW3_DOUBLE_LIB}")
message(STATUS "  header  : ${FFTW3_HEADER}")
message(STATUS "FFTW3f single precision: OK (${FFTW3f_LIBRARIES})")

# -----------------------------------------------------------------------
# Find HDF5 (C interface only — project uses HDF5 C API)
# -----------------------------------------------------------------------
find_package(HDF5 REQUIRED COMPONENTS C)

message(STATUS "HDF5 version: ${HDF5_VERSION}")
message(STATUS "HDF5 includes: ${HDF5_INCLUDE_DIRS}")

# -----------------------------------------------------------------------
# Find Eigen3
# -----------------------------------------------------------------------
find_package(Eigen3 REQUIRED CONFIG
    PATHS /opt/homebrew/share/eigen3/cmake
          /opt/homebrew/Cellar/eigen/5.0.1/share/eigen3/cmake)

message(STATUS "Eigen3 version: ${Eigen3_VERSION}")

# -----------------------------------------------------------------------
# Collect all .cpp sources under src/cpp/
# NOTE: if you add a new .cpp file, re-run cmake to pick it up.
# -----------------------------------------------------------------------
file(GLOB_RECURSE CPP_SOURCES "src/cpp/*.cpp")

if(NOT CPP_SOURCES)
    # No .cpp files yet — create a placeholder so cmake doesn't error.
    message(STATUS "No .cpp sources found yet in src/cpp/ — add source files and re-run cmake.")
    # Provide a minimal stub so the target can be declared without sources
    set(CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/main.cpp")
endif()

add_executable(les_apriori ${CPP_SOURCES})

# -----------------------------------------------------------------------
# Include directories
# -----------------------------------------------------------------------
target_include_directories(les_apriori PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp   # for config.h and local headers
    ${FFTW3_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

# Also resolve full path to fftw3f library to avoid linker "library not found" errors.
# FFTW3's CMake config only sets the library *name* (e.g. "fftw3f"), not a full path.
# find_library resolves it to the actual .dylib/.a on this machine.
find_library(FFTW3F_LIB
    NAMES fftw3f
    HINTS ${FFTW3f_LIBRARY_DIRS}
    REQUIRED)

# Eigen3 exposes a proper imported target; FFTW3 does not — use resolved full paths.
target_link_libraries(les_apriori PRIVATE
    Eigen3::Eigen
    ${FFTW3_DOUBLE_LIB}
    ${FFTW3F_LIB}
    ${HDF5_C_LIBRARIES}
)

# -----------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------
message(STATUS "")
message(STATUS "=== les_apriori build configuration ===")
message(STATUS "  C++ standard : C++17")
message(STATUS "  Compiler     : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "  Flags        : -O3 -march=native -ffast-math")
message(STATUS "  FFTW3 (dbl)  : ${FFTW3_LIBRARIES}")
message(STATUS "  FFTW3 (flt)  : ${FFTW3f_LIBRARIES}")
message(STATUS "  HDF5         : ${HDF5_C_LIBRARIES}")
message(STATUS "  Eigen3       : ${Eigen3_VERSION}")
message(STATUS "  Sources      : ${CPP_SOURCES}")
message(STATUS "========================================")
message(STATUS "")
