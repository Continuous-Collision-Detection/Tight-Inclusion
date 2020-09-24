# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(TI_CCDDownloadExternal)

################################################################################
# Required libraries
################################################################################

# libigl
if(NOT TARGET igl)
  ccd_download_libigl()
  add_subdirectory(${CCD_EXTERNAL}/libigl EXCLUDE_FROM_ALL)
  # Set Eigen directory to the one in libigl (needed for EVCTCD)
  set(ENV{EIGEN3_INCLUDE_DIR} "${CCD_EXTERNAL}/libigl/external/eigen/")
endif()

# # HDF5 Reader
# if(NOT TARGET HighFive::HighFive)
#   option(HIGHFIVE_UNIT_TESTS "Enable unit tests" OFF)
#   option(HIGHFIVE_EXAMPLES "Compile examples" OFF)
#   set(USE_EIGEN TRUE CACHE BOOL "Enable Eigen testing" FORCE)
#   ccd_download_high_five()
#   add_subdirectory(${CCD_EXTERNAL}/HighFive EXCLUDE_FROM_ALL)
#   add_library(HighFive::HighFive ALIAS HighFive)
# endif()


#GMP
find_package(GMPECCD)
IF(NOT ${GMP_FOUND})
        MESSAGE(FATAL_ERROR "Cannot find GMP")
ENDIF()