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
This file builds an example project consuming LEMON in three different ways:
1. Using locally installed LEMON - find_package.
2. Using add_subdirectory.
3. Using FetchContent.
]]

cmake_minimum_required(VERSION 3.15)
project(
  demo
  VERSION 0.1.0
  LANGUAGES CXX)

# Option 1.
# find_package(LEMON REQUIRED)

# Option 2.
# add_subdirectory(path-to-lemon)

# Option 3.
# FetchContent_Declare(LEMON GIT_REPOSITORY "git@github.com:MultiFlow/LEMON.git")
# FetchContent_MakeAvailable(LEMON)

set(DEMOS arg_parser_demo graph_to_eps_demo lgf_demo)

foreach(DEMO_NAME ${DEMOS})
  add_executable(${DEMO_NAME} ${DEMO_NAME}.cc)
  target_link_libraries(${DEMO_NAME} LEMON::LEMON)
endforeach()
