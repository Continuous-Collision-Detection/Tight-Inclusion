# App options
option(TIGHT_INCLUSION_WITH_SAMPLE_QUERIES "Enable sample queries in the application" OFF)

# Add the executable
add_executable(Tight_Inclusion_bin "main.cpp")

# Link to the library
target_link_libraries(Tight_Inclusion_bin PUBLIC tight_inclusion)

# Fill in configuration options
configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in" "${CMAKE_CURRENT_SOURCE_DIR}/config.hpp")

# Link to the required libraries
if(TIGHT_INCLUSION_WITH_SAMPLE_QUERIES)
    include(ccd_query_io)
    target_link_libraries(Tight_Inclusion_bin PUBLIC ccd_io::ccd_io)
endif()

include(pbar)
target_link_libraries(Tight_Inclusion_bin PUBLIC pbar::pbar)

include(cli11)
target_link_libraries(Tight_Inclusion_bin PUBLIC CLI11::CLI11)