find_package(Doxygen REQUIRED)
find_package(Ghostscript REQUIRED)
find_package(PythonInterp REQUIRED)
find_package(WGET REQUIRED)

set(PACKAGE_NAME ${PROJECT_NAME})
set(PACKAGE_VERSION ${PROJECT_VERSION})
set(abs_top_srcdir ${PROJECT_SOURCE_DIR})
set(abs_top_builddir ${PROJECT_BINARY_DIR})

set(LEMON_DOC_SOURCE_BROWSER
    "NO"
    CACHE STRING "Include source into the doc (YES/NO).")
set(LEMON_DOC_USE_MATHJAX
    "NO"
    CACHE STRING "Use MathJax to display math formulae (YES/NO).")
set(LEMON_DOC_MATHJAX_RELPATH
    "http://www.mathjax.org/mathjax"
    CACHE STRING "MathJax library location.")

set(LEMON_DOC_LIBSTDC++_URL
    "http://gcc.gnu.org/onlinedocs/gcc-4.7.3/libstdc++/api"
    CACHE STRING "GCC libstdc++ doxygen doc url.")

configure_file(${PROJECT_SOURCE_DIR}/doc/Doxyfile.in
               ${PROJECT_BINARY_DIR}/doc/Doxyfile @ONLY)

configure_file(${PROJECT_SOURCE_DIR}/doc/mainpage.dox.in
               ${PROJECT_BINARY_DIR}/doc/mainpage.dox @ONLY)

# Copy doc from source (if exists)
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/html
   AND NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
  message(STATUS "Copy doc from source tree")
  execute_process(
    COMMAND cmake -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/html
            ${CMAKE_CURRENT_BINARY_DIR}/html)
endif()

if(DOXYGEN_EXECUTABLE
   AND PYTHONINTERP_FOUND
   AND GHOSTSCRIPT_EXECUTABLE)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/)
  set(GHOSTSCRIPT_OPTIONS
      -dNOPAUSE
      -dBATCH
      -q
      -dEPSCrop
      -dTextAlphaBits=4
      -dGraphicsAlphaBits=4
      -sDEVICE=pngalpha)
  add_custom_target(
    html
    COMMAND ${CMAKE_COMMAND} -E remove_directory gen-images
    COMMAND ${CMAKE_COMMAND} -E make_directory gen-images
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r20
      -sOutputFile=gen-images/grid_graph.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/grid_graph.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32
      -sOutputFile=gen-images/adaptors2.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/adaptors2.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32
      -sOutputFile=gen-images/connected_components.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/connected_components.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32
      -sOutputFile=gen-images/strongly_connected_components.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/strongly_connected_components.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32
      -sOutputFile=gen-images/node_biconnected_components.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/node_biconnected_components.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32
      -sOutputFile=gen-images/edge_biconnected_components.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/edge_biconnected_components.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r32
      -sOutputFile=gen-images/bipartite_partitions.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/bipartite_partitions.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r24
      -sOutputFile=gen-images/matching.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/matching.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r24
      -sOutputFile=gen-images/bipartite_matching.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/bipartite_matching.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r40
      -sOutputFile=gen-images/planar.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/planar.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r24
      -sOutputFile=gen-images/tsp.png ${CMAKE_CURRENT_SOURCE_DIR}/images/tsp.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8
      -sOutputFile=gen-images/nodeshape_0.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_0.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8
      -sOutputFile=gen-images/nodeshape_1.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_1.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8
      -sOutputFile=gen-images/nodeshape_2.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_2.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8
      -sOutputFile=gen-images/nodeshape_3.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_3.eps
    COMMAND
      ${GHOSTSCRIPT_EXECUTABLE} ${GHOSTSCRIPT_OPTIONS} -r8
      -sOutputFile=gen-images/nodeshape_4.png
      ${CMAKE_CURRENT_SOURCE_DIR}/images/nodeshape_4.eps
    COMMAND ${CMAKE_COMMAND} -E remove_directory html
    COMMAND ${DOXYGEN_EXECUTABLE} Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

  set_target_properties(html PROPERTIES PROJECT_LABEL BUILD_DOC)

  if(UNIX)
    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
      DESTINATION share/doc/lemon/html
      COMPONENT html_documentation)
  elseif(WIN32)
    install(
      DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/
      DESTINATION doc
      COMPONENT html_documentation)
  endif()

endif()

if(WGET_FOUND)
  add_custom_target(
    update-external-tags
    COMMAND ${WGET_EXECUTABLE} -N ${LEMON_DOC_LIBSTDC++_URL}/libstdc++.tag
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()
