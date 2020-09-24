cmake_minimum_required(VERSION 3.8)
project(TI_CCD)

################################################################################

# if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
#     set(CCD_TOPLEVEL_PROJECT ON)
# else()
#     set(CCD_TOPLEVEL_PROJECT OFF)
# endif()



# libigl
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" OFF)
option(LIBIGL_WITH_ANTTWEAKBAR       "Use AntTweakBar"              OFF)
option(LIBIGL_WITH_CGAL              "Use CGAL"                     OFF)
option(LIBIGL_WITH_COMISO            "Use CoMiso"                   OFF)
option(LIBIGL_WITH_CORK              "Use Cork"                     OFF)
option(LIBIGL_WITH_EMBREE            "Use Embree"                   OFF)
option(LIBIGL_WITH_LIM               "Use LIM"                      OFF)
option(LIBIGL_WITH_MATLAB            "Use Matlab"                   OFF)
option(LIBIGL_WITH_MOSEK             "Use MOSEK"                    OFF)
option(LIBIGL_WITH_OPENGL            "Use OpenGL"                   OFF)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"                     OFF)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"                    OFF)
option(LIBIGL_WITH_PNG               "Use PNG"                      OFF)
option(LIBIGL_WITH_PYTHON            "Use Python"                   OFF)
option(LIBIGL_WITH_TETGEN            "Use Tetgen"                   OFF)
option(LIBIGL_WITH_TRIANGLE          "Use Triangle"                 OFF)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"            OFF)
option(LIBIGL_WITH_XML               "Use XML"                      OFF)
option(LIBIGL_WITH_PREDICATES        "Use exact predicates"         ON )
################################################################################




### Configuration
set(CCD_EXTERNAL "${CMAKE_CURRENT_SOURCE_DIR}/external")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

include(Warnings)
include(UseColors)
include(${PROJECT_NAME}Dependencies)
include(TI_CCDUtils)


# inclusion CCD
add_subdirectory(tight_inclusion)
prepend_current_path(INCLUSION_SOURCES)
inclusionCCD_copy_headers(${INCLUSION_SOURCES})
CCD_set_source_group(${INCLUSION_SOURCES})
# message(STATUS "INCLUSION_SOURCE is ${INCLUSION_SOURCES}")
add_library(tight_inclusion ${INCLUSION_SOURCES})
target_include_directories(tight_inclusion PUBLIC ${PROJECT_BINARY_DIR}/include)
target_include_directories(tight_inclusion PUBLIC interval_ccd)

# boost
# find_package(Boost REQUIRED)
# target_link_libraries(tight_inclusion PUBLIC Boost::boost)

# libigl
target_link_libraries(tight_inclusion PUBLIC igl::core igl::predicates)

# Extra warnings
target_link_libraries(tight_inclusion PRIVATE warnings::all)
# #GMP
target_include_directories(tight_inclusion PUBLIC ${GMP_INCLUDE_DIRS} src)
target_link_libraries(tight_inclusion PUBLIC ${GMP_LIBRARIES})

################################################################################
# Executable
################################################################################
add_executable(${PROJECT_NAME}_bin
    # Add executable source here
    src/main.cpp
  )

  target_link_libraries(${PROJECT_NAME}_bin PUBLIC tight_inclusion)
  # WARNING: This requires a check  as it will break compatability with non-AVX
  #          processors.
  target_compile_options(tight_inclusion PUBLIC "-mavx2")
  set(DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data/")
  target_compile_definitions(${PROJECT_NAME}_bin PUBLIC -DCCD_DATA_PATH=\"${DATA_PATH}\")
