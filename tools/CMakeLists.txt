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
This file builds some optional tools.
]]

add_executable(lgf-gen lgf-gen.cc)
target_link_libraries(lgf-gen LEMON::LEMON)

add_executable(dimacs-to-lgf dimacs-to-lgf.cc)
target_link_libraries(dimacs-to-lgf LEMON::LEMON)

add_executable(dimacs-solver dimacs-solver.cc)
target_link_libraries(dimacs-solver LEMON::LEMON)

install(TARGETS lgf-gen dimacs-to-lgf dimacs-solver RUNTIME DESTINATION bin
                                                            COMPONENT bin)

if(NOT WIN32)
  install(
    PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/lemon-0.x-to-1.x.sh
    DESTINATION bin
    COMPONENT bin)
endif()
