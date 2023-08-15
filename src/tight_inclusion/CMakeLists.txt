set(SOURCES
    avx.hpp
    avx.cpp
    ccd.cpp
    ccd.hpp
    interval_root_finder.cpp
    interval_root_finder.hpp
    interval.cpp
    interval.hpp
    logger.cpp
    logger.hpp
    timer.hpp
    types.hpp
)

source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" PREFIX "Source Files" FILES ${SOURCES})
target_sources(tight_inclusion PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

if(TIGHT_INCLUSION_WITH_GMP OR TIGHT_INCLUSION_WITH_TESTS)
    add_subdirectory(rational)
endif()
