# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(TIGHT_INCLUSION_TOPLEVEL_PROJECT OFF)
else()
    set(TIGHT_INCLUSION_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.14.0")
if(TIGHT_INCLUSION_TOPLEVEL_PROJECT)
    cmake_minimum_required(VERSION ${REQUIRED_CMAKE_VERSION})
else()
    # Don't use cmake_minimum_required here to avoid implicitly overriding parent policies
    if(${CMAKE_VERSION} VERSION_LESS ${REQUIRED_CMAKE_VERSION})
        message(FATAL_ERROR "CMake required version to build Tight Inclusion is ${REQUIRED_CMAKE_VERSION}")
    endif()
endif()

# Include user-provided default options if available. We do that before the main
# `project()` so that we can define the C/C++ compilers from the option file.
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/TightInclusionOptions.cmake)
    message(STATUS "Using local options file: ${CMAKE_CURRENT_SOURCE_DIR}/TightInclusionOptions.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/TightInclusionOptions.cmake)
endif()

################################################################################

project(TightInclusion
        DESCRIPTION "Tight Inclusion CCD"
        LANGUAGES CXX)

OPTION(TIGHT_INCLUSION_WITH_TESTS             "Enable test functions"                                   OFF)
OPTION(TIGHT_INCLUSION_WITH_GMP               "Enable rational based predicates, for debug"             OFF)
OPTION(TIGHT_INCLUSION_WITH_TIMER             "Enable profiling timers, for debug"                      OFF)
OPTION(TIGHT_INCLUSION_WITH_DOUBLE_PRECISION  "Enable double precision floating point numbers as input"  ON)
OPTION(TIGHT_INCLUSION_LIMIT_QUEUE_SIZE       "Enable limitation of maximal queue size"                 OFF)

# Option to supress progress output when on GH actions.
OPTION(TIGHT_INCLUSION_SUPPRESS_PROGRESS_OUTPUT "Enable limitation of maximal queue size"               OFF)
mark_as_advanced(TIGHT_INCLUSION_SUPPRESS_PROGRESS_OUTPUT)

include(CMakeDependentOption)
cmake_dependent_option(TIGHT_INCLUSION_FLOAT_WITH_DOUBLE_INPUT "Enable converting double queries to float" OFF "TIGHT_INCLUSION_WITH_DOUBLE_PRECISION" OFF)

# Set default minimum C++ standard
if(TIGHT_INCLUSION_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

### Configuration
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/tight_inclusion/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find/")

# Tight Inclusion utils
include(tight_inclusion_utils)

################################################################################
# Tight Inclusion Library
################################################################################

# Add an empty library and fill in the list of sources in `src/CMakeLists.txt`.
add_library(tight_inclusion)
add_library(tight_inclusion::tight_inclusion ALIAS tight_inclusion)

add_subdirectory(src/tight_inclusion)

# Public include directory for Tight Inclusion
target_include_directories(tight_inclusion PUBLIC src)

################################################################################
# Optional Definitions
################################################################################

# For MSVC, do not use the min and max macros.
target_compile_definitions(tight_inclusion PUBLIC NOMINMAX)

################################################################################
# Dependencies
################################################################################

# Extra warnings
include(tight_inclusion_warnings)
target_link_libraries(tight_inclusion PRIVATE tight_inclusion::warnings)

# libigl
include(eigen)
target_link_libraries(tight_inclusion PUBLIC Eigen3::Eigen)

# GMP (optional)
if(TIGHT_INCLUSION_WITH_GMP OR TIGHT_INCLUSION_WITH_TESTS)
    find_package(GMP REQUIRED)
    target_link_libraries(tight_inclusion PUBLIC gmp::gmp)
    target_compile_definitions(tight_inclusion PUBLIC TIGHT_INCLUSION_USE_GMP)
endif()

################################################################################
# Definitions
################################################################################

if(TIGHT_INCLUSION_WITH_TIMER)
    target_compile_definitions(tight_inclusion PRIVATE TIGHT_INCLUSION_USE_TIMER)
endif()

if (TIGHT_INCLUSION_WITH_DOUBLE_PRECISION)
    message(STATUS "Tight Inclusion: Using Double Precision Floating Points")
    target_compile_definitions(tight_inclusion PUBLIC TIGHT_INCLUSION_DOUBLE)
else()
    message(STATUS "Tight Inclusion: Using Single Precision Floating Points")
endif()

if(TIGHT_INCLUSION_LIMIT_QUEUE_SIZE)
    message(STATUS "Tight Inclusion: Limiting maximal queue size")
    target_compile_definitions(tight_inclusion PUBLIC TIGHT_INCLUSION_LIMIT_QUEUE_SIZE)
endif()

if(TIGHT_INCLUSION_FLOAT_WITH_DOUBLE_INPUT)
    message(STATUS "Tight Inclusion: Converting double inputs to float for tests")
    target_compile_definitions(tight_inclusion PUBLIC TIGHT_INCLUSION_FWDI)
endif()

################################################################################
# Compiler options
################################################################################

# Figure out AVX level support
message(STATUS "Searching for AVX...")
find_package(AVX)
string(REPLACE " " ";" SIMD_FLAGS "${AVX_FLAGS}")
target_compile_options(tight_inclusion PRIVATE ${SIMD_FLAGS})

# Use C++17
target_compile_features(tight_inclusion PUBLIC cxx_std_17)

################################################################################
# App
################################################################################

if(TIGHT_INCLUSION_TOPLEVEL_PROJECT)
    add_executable(Tight_Inclusion_bin "app/main.cpp")
    target_link_libraries(Tight_Inclusion_bin PUBLIC tight_inclusion)

    set(TIGHT_INCLUSION_SAMPLE_QUERIES_DIR "${PROJECT_SOURCE_DIR}/sample-queries")
    include(sample_queries)

    if(TIGHT_INCLUSION_SUPPRESS_PROGRESS_OUTPUT)
        target_compile_definitions(Tight_Inclusion_bin PUBLIC TIGHT_INCLUSION_SUPPRESS_PROGRESS_OUTPUT)
    endif()

    if (TIGHT_INCLUSION_WITH_TESTS)
        target_sources(Tight_Inclusion_bin PUBLIC "app/read_rational_csv.cpp")
        target_compile_definitions(Tight_Inclusion_bin PUBLIC TIGHT_INCLUSION_RUN_EXAMPLES)
        target_compile_definitions(Tight_Inclusion_bin PUBLIC
            TIGHT_INCLUSION_SAMPLE_QUERIES_DIR="${TIGHT_INCLUSION_SAMPLE_QUERIES_DIR}/")
    endif()
endif()
