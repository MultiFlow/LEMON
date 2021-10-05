find_path(SOPLEX_INCLUDE_DIR soplex.h HINTS ${SOPLEX_ROOT}/src)
find_library(SOPLEX_LIBRARY soplex HINTS ${SOPLEX_ROOT}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SOPLEX DEFAULT_MSG SOPLEX_LIBRARY
                                  SOPLEX_INCLUDE_DIR)

if(SOPLEX_FOUND)
  set(SOPLEX_INCLUDE_DIRS ${SOPLEX_INCLUDE_DIR})
  set(SOPLEX_LIBRARIES ${SOPLEX_LIBRARY})
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(SOPLEX_LIBRARIES "${SOPLEX_LIBRARIES};z")
  endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(SOPLEX_FOUND)

mark_as_advanced(SOPLEX_LIBRARY SOPLEX_INCLUDE_DIR)
