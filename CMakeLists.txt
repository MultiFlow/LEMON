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
This file defines the build options and orchestrates the build by including
files in the cmake directory.
]]

cmake_minimum_required(VERSION 3.15)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(
  LEMON
  VERSION 1.4.0
  LANGUAGES CXX C)

# OS
if(WIN32)
  set(LEMON_WIN32 ON)
endif()

# Building Options #

# Solvers
option(LEMON_ENABLE_GLPK "Enable GLPK solver backend." OFF)
option(LEMON_ENABLE_ILOG "Enable ILOG (CPLEX) solver backend." OFF)
option(LEMON_ENABLE_COIN "Enable COIN solver backend." OFF)
option(LEMON_ENABLE_SOPLEX "Enable SoPlex solver backend." OFF)

# Maintainer only options
option(LEMON_BUILD_DOCS "" OFF)
option(LEMON_BUILD_TOOLS "" OFF)
option(LEMON_BUILD_DEMO "" OFF)
option(LEMON_BUILD_ARCHIVE "Enable to create archive" OFF)
option(LEMON_BUILD_CPACK "Enable to create cpack" OFF)
option(LEMON_TEST_WITH_VALGRIND "Run the test wit valgrind" OFF)
option(LEMON_ENABLE_DEBUG "Display debug messages" OFF)

if(CMAKE_CONFIGURATION_TYPES)
  list(APPEND CMAKE_CONFIGURATION_TYPES Maintainer)
  list(REMOVE_DUPLICATES CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_CONFIGURATION_TYPES
      "${CMAKE_CONFIGURATION_TYPES}"
      CACHE STRING "Add the configurations that we need" FORCE)
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_BUILD_TYPE
    "${CMAKE_BUILD_TYPE}"
    CACHE
      STRING
      "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel Maintainer."
      FORCE)

enable_testing()

# Add LEMON::LEMON target
include(base)

if(BUILD_TESTING AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  add_subdirectory(test)
endif()

if(LEMON_BUILD_TOOLS AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  add_subdirectory(tools)
endif()

if(LEMON_BUILD_DOCS AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  add_subdirectory(doc)
endif()

if(LEMON_BUILD_DEMO AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  add_subdirectory(demo)
endif()

if(LEMON_BUILD_ARCHIVE AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  include(archive)
endif()

if(LEMON_BUILD_CPACK AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  include(cpack)
endif()
