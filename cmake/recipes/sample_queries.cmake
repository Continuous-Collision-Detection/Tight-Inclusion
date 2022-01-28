message(STATUS "Downloading sample CCD queries")

include(FetchContent)
FetchContent_Declare(
    sample_queries
    GIT_REPOSITORY https://github.com/Continuous-Collision-Detection/Sample-Queries.git
    GIT_TAG 4d6cce33477d8d5c666c31c8ea23e1aea97be371
    GIT_SHALLOW FALSE
    SOURCE_DIR ${TIGHT_INCLUSION_SAMPLE_QUERIES_DIR}
)
FetchContent_GetProperties(sample_queries)
if(NOT sample_queries_POPULATED)
  FetchContent_Populate(sample_queries)
endif()
