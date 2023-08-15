if(TARGET pbar::pbar)
    return()
endif()

message(STATUS "Third-party: creating target 'pbar::pbar'")

include(CPM)
CPMAddPackage(
    NAME pbar
    GITHUB_REPOSITORY estshorter/pbar
    GIT_TAG 226b9ad291d72c6859456d5188f33c932e2e730b
    DOWNLOAD_ONLY ON
)

add_library(pbar INTERFACE)
add_library(pbar::pbar ALIAS pbar)

target_include_directories(pbar INTERFACE "${pbar_SOURCE_DIR}")