find_package (Doxygen)

message("e2ga_SOURCE_DIR     " ${e2ga_SOURCE_DIR})
message("e2ga_BINARY_DIR     " ${e2ga_BINARY_DIR})
message("DOXYGEN_EXECUTABLE  " ${DOXYGEN_EXECUTABLE})

if(DOXYGEN_FOUND)
  configure_file (${e2ga_SOURCE_DIR}/doc/Doxyfile-html.cmake ${e2ga_BINARY_DIR}/doc/Doxyfile-html)
  add_custom_target (html ${DOXYGEN_EXECUTABLE} ${e2ga_BINARY_DIR}/doc/Doxyfile-html)
endif (DOXYGEN_FOUND)
