cmake_minimum_required(VERSION 3.8)
project(TI_CCD)

################################################################################

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
OPTION(TIGHT_INCLUSION_WITH_TIMER "Enable profiling timers, for debug"          OFF)

include(${PROJECT_NAME}Dependencies)


# inclusion CCD
add_subdirectory(tight_inclusion)
prepend_current_path(INCLUSION_SOURCES)
inclusionCCD_copy_headers(${INCLUSION_SOURCES})
CCD_set_source_group(${INCLUSION_SOURCES})
add_library(tight_inclusion ${INCLUSION_SOURCES})
target_include_directories(tight_inclusion PUBLIC ${PROJECT_BINARY_DIR}/include)

target_link_libraries(tight_inclusion INTERFACE Eigen3::Eigen)
target_link_libraries(tight_inclusion PRIVATE warnings::all)

#Optional
#GMP
if(TIGHT_INCLUSION_WITH_GMP)
  target_include_directories(tight_inclusion PUBLIC ${GMP_INCLUDE_DIRS} src)
  target_link_libraries(tight_inclusion PUBLIC ${GMP_LIBRARIES})
endif()

if(TIGHT_INCLUSION_WITH_TIMER)
  target_compile_definitions(tight_inclusion PRIVATE TIGHT_INCLUSION_USE_TIMER)
endif()

# WARNING: This requires a check  as it will break compatability with non-AVX
#          processors.
target_compile_options(tight_inclusion PUBLIC "-mavx2")

################################################################################
# Executable
################################################################################
if(TIGHT_INCLUSION_TOPLEVEL_PROJECT)
  add_executable(${PROJECT_NAME}_bin src/main.cpp)
  target_link_libraries(${PROJECT_NAME}_bin PUBLIC tight_inclusion)

  set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/")
  target_compile_definitions(${PROJECT_NAME}_bin PUBLIC -DCCD_DATA_PATH=\"${DATA_PATH}\")
endif()
