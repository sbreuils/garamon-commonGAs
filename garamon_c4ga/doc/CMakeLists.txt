find_package (Doxygen)

message("c4ga_SOURCE_DIR     " ${c4ga_SOURCE_DIR})
message("c4ga_BINARY_DIR     " ${c4ga_BINARY_DIR})
message("DOXYGEN_EXECUTABLE  " ${DOXYGEN_EXECUTABLE})

if(DOXYGEN_FOUND)
  configure_file (${c4ga_SOURCE_DIR}/doc/Doxyfile-html.cmake ${c4ga_BINARY_DIR}/doc/Doxyfile-html)
  add_custom_target (html ${DOXYGEN_EXECUTABLE} ${c4ga_BINARY_DIR}/doc/Doxyfile-html)
endif (DOXYGEN_FOUND)
