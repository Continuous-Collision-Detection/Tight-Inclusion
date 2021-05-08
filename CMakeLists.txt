cmake_minimum_required(VERSION 3.8)
project(TI_CCD)

################################################################################
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(TIGHT_INCLUSION_TOPLEVEL_PROJECT ON)
else()
    set(TIGHT_INCLUSION_TOPLEVEL_PROJECT OFF)
endif()
### Configuration
set(TIGHT_INCLUSION_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)
include(TI_CCDUtils)


OPTION(TIGHT_INCLUSION_WITH_GMP   "Enable rational based predicates, for debug" OFF)
OPTION(TIGHT_INCLUSION_WITH_TESTS   "Enable test functions" OFF)
OPTION(TIGHT_INCLUSION_WITH_TIMER "Enable profiling timers, for debug"          OFF)
OPTION(TIGHT_INCLUSION_WITH_NO_ZERO_TOI 
"Continue to earch solution in higher resolution if toi is 0"         OFF)

include(${PROJECT_NAME}Dependencies)


# inclusion CCD
add_subdirectory(tight_inclusion)
prepend_current_path(INCLUSION_SOURCES)
inclusionCCD_copy_headers(${INCLUSION_SOURCES})
CCD_set_source_group(${INCLUSION_SOURCES})
add_library(tight_inclusion ${INCLUSION_SOURCES})
target_include_directories(tight_inclusion PUBLIC ${PROJECT_BINARY_DIR}/include)
target_link_libraries(tight_inclusion PUBLIC Eigen3::Eigen)
target_link_libraries(tight_inclusion PRIVATE warnings::all)

#Optional
#GMP
if(TIGHT_INCLUSION_WITH_GMP)
message(STATUS "TIGHT_INCLUSION_WITH_GMP is defined, now using rational root finder")
  target_include_directories(tight_inclusion PUBLIC ${GMP_INCLUDE_DIRS} src)
  target_link_libraries(tight_inclusion PUBLIC ${GMP_LIBRARIES})
  target_compile_definitions(tight_inclusion PUBLIC TIGHT_INCLUSION_USE_GMP)
endif()

if(TIGHT_INCLUSION_WITH_TIMER)
  target_compile_definitions(tight_inclusion PRIVATE TIGHT_INCLUSION_USE_TIMER)
endif()

# Figure out AVX level support
message(STATUS "Searching for AVX...")
find_package(AVX)
# Add SSE, AVX, and FMA flags to compiler flags
string(REPLACE " " ";" SIMD_FLAGS "${AVX_FLAGS}")
target_compile_options(tight_inclusion PRIVATE ${SIMD_FLAGS})
if (TIGHT_INCLUSION_WITH_NO_ZERO_TOI)
  target_compile_definitions(tight_inclusion PRIVATE TIGHT_INCLUSION_NO_ZERO_TOI)
endif()
target_compile_definitions(tight_inclusion PUBLIC NOMINMAX)

if(TIGHT_INCLUSION_TOPLEVEL_PROJECT)
if (TIGHT_INCLUSION_WITH_TESTS)
add_executable(Tight_Inclusion_bin
app/main.cpp
app/read_rational_csv.cpp
)
else()
add_executable(Tight_Inclusion_bin
app/main.cpp
)
endif()
target_link_libraries(Tight_Inclusion_bin PUBLIC tight_inclusion)

ticcd_download_sample_queries()
	if (NOT TIGHT_INCLUSION_WITH_GMP AND TIGHT_INCLUSION_WITH_TESTS)
		target_include_directories(tight_inclusion PUBLIC ${GMP_INCLUDE_DIRS} src)
		target_link_libraries(tight_inclusion PUBLIC ${GMP_LIBRARIES})
		target_compile_definitions(tight_inclusion PUBLIC TIGHT_INCLUSION_RUN_EXAMPLES)
    target_compile_definitions(Tight_Inclusion_bin PUBLIC
    TICCD_EXAMPLE_QUERIES_DIR="${CMAKE_CURRENT_SOURCE_DIR}/external/Sample-Queries/")
	endif()
endif()