#[[
This file is a part of LEMON, a generic C++ optimization library.

Copyright (C) 2003-2021
Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
(Egervary Research Group on Combinatorial Optimization, EGRES).

Permission to use, modify and distribute this software is granted
provided that this copyright notice appears in all copies. For
precise terms see the accompanying LICENSE file.

This software is provided "AS IS" with no warranty of any kind,
express or implied, and with no claim as to its suitability for any
purpose.
]]

#[[
This file finds the necessary components for the SOPLEX solver.
]]

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
