# Detects whether this is a top-level project
get_directory_property(HAS_PARENT PARENT_DIRECTORY)
if(HAS_PARENT)
    set(TIGHT_INCLUSION_TOPLEVEL_PROJECT OFF)
else()
    set(TIGHT_INCLUSION_TOPLEVEL_PROJECT ON)
endif()

# Check required CMake version
set(REQUIRED_CMAKE_VERSION "3.18.0")
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

# Enable ccache if available
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    option(TIGHT_INCLUSION_WITH_CCACHE "Enable ccache when building Tight Inclusion" ${TIGHT_INCLUSION_TOPLEVEL_PROJECT})
else()
    option(TIGHT_INCLUSION_WITH_CCACHE "Enable ccache when building Tight Inclusion" OFF)
endif()
if(TIGHT_INCLUSION_WITH_CCACHE AND CCACHE_PROGRAM)
    message(STATUS "Enabling Ccache support")
    set(ccacheEnv
        CCACHE_BASEDIR=${CMAKE_BINARY_DIR}
        CCACHE_SLOPPINESS=clang_index_store,include_file_ctime,include_file_mtime,locale,pch_defines,time_macros
    )
    foreach(lang IN ITEMS C CXX)
        set(CMAKE_${lang}_COMPILER_LAUNCHER
            ${CMAKE_COMMAND} -E env ${ccacheEnv} ${CCACHE_PROGRAM}
        )
    endforeach()
endif()

################################################################################
# CMake Policies
################################################################################

cmake_policy(SET CMP0054 NEW) # Only interpret if() arguments as variables or keywords when unquoted.
cmake_policy(SET CMP0076 NEW) # target_sources() command converts relative paths to absolute.
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24")
    cmake_policy(SET CMP0135 NEW) # Set the timestamps of all extracted contents to the time of the extraction.
endif()

################################################################################

project(TightInclusion
        DESCRIPTION "Tight Inclusion CCD"
        LANGUAGES CXX
        VERSION "1.0.4")

option(TIGHT_INCLUSION_WITH_RATIONAL          "Enable rational based predicates (for debugging)"        OFF)
option(TIGHT_INCLUSION_WITH_TIMER             "Enable profiling timers (for debugging)"                 OFF)
option(TIGHT_INCLUSION_WITH_DOUBLE_PRECISION  "Enable double precision floating point numbers as input"  ON)
option(TIGHT_INCLUSION_LIMIT_QUEUE_SIZE       "Enable limitation of maximal queue size"                 OFF)

include(CMakeDependentOption)
cmake_dependent_option(TIGHT_INCLUSION_FLOAT_WITH_DOUBLE_INPUT "Enable converting double queries to float" OFF "TIGHT_INCLUSION_WITH_DOUBLE_PRECISION" ON)

# Set default minimum C++ standard
if(TIGHT_INCLUSION_TOPLEVEL_PROJECT)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)
endif()

### Configuration
set(TIGHT_INCLUSION_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src/tight_inclusion")
set(TIGHT_INCLUSION_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/src")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/tight_inclusion/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/recipes/")
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/find/")

# General CMake utils
include(tight_inclusion_cpm_cache)
include(tight_inclusion_use_colors)

# Generate position-independent code by default
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

################################################################################
# Tight Inclusion Library
################################################################################

# Add an empty library and fill in the list of sources in `src/CMakeLists.txt`.
add_library(tight_inclusion)
add_library(tight_inclusion::tight_inclusion ALIAS tight_inclusion)

# Fill in configuration options
configure_file(
  "${TIGHT_INCLUSION_SOURCE_DIR}/config.hpp.in"
  "${TIGHT_INCLUSION_SOURCE_DIR}/config.hpp")

# Add source and header files to tight_inclusion
add_subdirectory("${TIGHT_INCLUSION_SOURCE_DIR}")

# Public include directory for Tight Inclusion
target_include_directories(tight_inclusion PUBLIC "${TIGHT_INCLUSION_INCLUDE_DIR}")

################################################################################
# Optional Definitions
################################################################################

# For MSVC, do not use the min and max macros.
target_compile_definitions(tight_inclusion PUBLIC NOMINMAX)

################################################################################
# Dependencies
################################################################################

# Eigen
include(eigen)
target_link_libraries(tight_inclusion PUBLIC Eigen3::Eigen)

# Logger
include(spdlog)
target_link_libraries(tight_inclusion PUBLIC spdlog::spdlog)

# rational-cpp (optional)
if(TIGHT_INCLUSION_WITH_RATIONAL)
    include(rational_cpp)
    target_link_libraries(tight_inclusion PUBLIC rational::rational)
endif()

# Extra warnings (link last for highest priority)
include(tight_inclusion_warnings)
target_link_libraries(tight_inclusion PRIVATE tight_inclusion::warnings)

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
    add_subdirectory(app)
endif()
