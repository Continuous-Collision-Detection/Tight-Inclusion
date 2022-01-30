# Copy header files into the build directory
function(tight_inclusion_copy_headers)
  foreach(filepath IN ITEMS ${ARGN})
    file(RELATIVE_PATH filename "${TIGHT_INCLUSION_SOURCE_DIR}" "${filepath}")
    if(${filename} MATCHES ".*\.(hpp|h|ipp|tpp)$")
      configure_file(${filepath} ${PROJECT_BINARY_DIR}/include/tight_inclusion/${filename})
    endif()
  endforeach()
endfunction()

function(tight_inclusion_prepend_current_path SOURCE_FILES)
  # Use recursive substitution to expand SOURCE_FILES
  unset(MODIFIED)
  foreach(SOURCE_FILE IN ITEMS ${${SOURCE_FILES}})
    list(APPEND MODIFIED "${CMAKE_CURRENT_SOURCE_DIR}/${SOURCE_FILE}")
  endforeach()
  set(${SOURCE_FILES} ${MODIFIED} PARENT_SCOPE)
endfunction()

# Set source group for IDE like Visual Studio or XCode
function(tight_inclusion_set_source_group)
  foreach(filepath IN ITEMS ${ARGN})
    get_filename_component(folderpath "${filepath}" DIRECTORY)
    get_filename_component(foldername "${folderpath}" NAME)
    source_group(foldername FILES "${filepath}")
  endforeach()
endfunction()
