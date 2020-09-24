include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(CCD_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(CCD_EXTRA_OPTIONS "")
endif()

function(ccd_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${CCD_EXTERNAL}/${name}
        DOWNLOAD_DIR ${CCD_EXTERNAL}/.cache/${name}
        QUIET
        ${CCD_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

# libigl
function(ccd_download_libigl)
  ccd_download_project(libigl
    GIT_REPOSITORY https://github.com/libigl/libigl.git
    GIT_TAG        56f129e4403d3b7b04ddd786745fd3d574e95e04
  )
endfunction()

# Catch2 for testing
function(ccd_download_catch2)
    ccd_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.11.1
    )
endfunction()

# # HighFive - Header-only C++ HDF5 interface
# function(ccd_download_high_five)
#   ccd_download_project(HighFive
#     GIT_REPOSITORY https://github.com/BlueBrain/HighFive.git
#     GIT_TAG        67f5b3a67ec04e021f5f86b3993d2edf0d223fef
#   )
# endfunction()
