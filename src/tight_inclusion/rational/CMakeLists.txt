set(SOURCES
    ccd.cpp
    ccd.hpp
    interval_root_finder.cpp
    interval_root_finder.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(tight_inclusion PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################
