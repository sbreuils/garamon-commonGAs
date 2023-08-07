find_package (Doxygen)

message("e4ga_SOURCE_DIR     " ${e4ga_SOURCE_DIR})
message("e4ga_BINARY_DIR     " ${e4ga_BINARY_DIR})
message("DOXYGEN_EXECUTABLE  " ${DOXYGEN_EXECUTABLE})

if(DOXYGEN_FOUND)
  configure_file (${e4ga_SOURCE_DIR}/doc/Doxyfile-html.cmake ${e4ga_BINARY_DIR}/doc/Doxyfile-html)
  add_custom_target (html ${DOXYGEN_EXECUTABLE} ${e4ga_BINARY_DIR}/doc/Doxyfile-html)
endif (DOXYGEN_FOUND)
