# Prepare dependencies
#
# For each third-party library, if the appropriate target doesn't exist yet,
# download it via external project, and add_subdirectory to build it alongside
# this project.


# Download and update 3rd_party libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
list(REMOVE_DUPLICATES CMAKE_MODULE_PATH)
include(TI_CCDDownloadExternal)


if(TIGHT_INCLUSION_WITH_GMP)
  #GMP
  find_package(GMPECCD)
  IF(NOT ${GMP_FOUND})
          MESSAGE(FATAL_ERROR "Cannot find GMP")
  ENDIF()
endif()


if(NOT TARGET Eigen3::Eigen)
  ccd_download_eigen()
  add_library(tccd_eigen INTERFACE)
  target_include_directories(tccd_eigen SYSTEM INTERFACE
    $<BUILD_INTERFACE:${TIGHT_INCLUSION_EXTERNAL}/eigen>
    $<INSTALL_INTERFACE:include>
  )
  set_property(TARGET tccd_eigen PROPERTY EXPORT_NAME Eigen3::Eigen)
  add_library(Eigen3::Eigen ALIAS tccd_eigen)
endif()