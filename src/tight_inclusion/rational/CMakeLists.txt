set(SOURCES
    ccd.cpp
    ccd.hpp
    interval_root_finder.cpp
    interval_root_finder.hpp
    rational.hpp
)

tight_inclusion_prepend_current_path(SOURCES)
tight_inclusion_set_source_group(${SOURCES})
target_sources(tight_inclusion PRIVATE ${SOURCES})

################################################################################
# Subfolders
################################################################################
