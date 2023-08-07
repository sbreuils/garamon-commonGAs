find_package (Doxygen)

message("e3ga_SOURCE_DIR     " ${e3ga_SOURCE_DIR})
message("e3ga_BINARY_DIR     " ${e3ga_BINARY_DIR})
message("DOXYGEN_EXECUTABLE  " ${DOXYGEN_EXECUTABLE})

if(DOXYGEN_FOUND)
  configure_file (${e3ga_SOURCE_DIR}/doc/Doxyfile-html.cmake ${e3ga_BINARY_DIR}/doc/Doxyfile-html)
  add_custom_target (html ${DOXYGEN_EXECUTABLE} ${e3ga_BINARY_DIR}/doc/Doxyfile-html)
endif (DOXYGEN_FOUND)
