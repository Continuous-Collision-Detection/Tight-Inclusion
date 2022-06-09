set(SOURCES
    avx.hpp
    avx.cpp
    ccd.cpp
    ccd.hpp
    interval_root_finder.cpp
    interval_root_finder.hpp
    interval.cpp
    interval.hpp
    timer.hpp
    types.hpp
)

tight_inclusion_prepend_current_path(SOURCES)
tight_inclusion_set_source_group(${SOURCES})
target_sources(tight_inclusion PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################

if(TIGHT_INCLUSION_WITH_GMP OR TIGHT_INCLUSION_WITH_TESTS)
    add_subdirectory(rational)
endif()
