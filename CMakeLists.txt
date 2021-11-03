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

project(LEMON LANGUAGES CXX C)

# Get version
find_package(Hg)
if(EXISTS ${PROJECT_SOURCE_DIR}/cmake/LEMONVersion.cmake)
  include(${PROJECT_SOURCE_DIR}/cmake/LEMONVersion.cmake)
elseif(DEFINED ENV{LEMON_VERSION})
  set(LEMON_VERSION
      $ENV{LEMON_VERSION}
      CACHE STRING "LEMON version string.")
elseif(HG_FOUND)
  execute_process(
    COMMAND hg log -r. --template "{latesttag}"
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE HG_REVISION_TAG
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND hg log -r. --template "{latesttagdistance}"
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE HG_REVISION_DIST
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(
    COMMAND hg log -r. --template "{node|short}"
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE HG_REVISION_ID
    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(HG_REVISION_TAG STREQUAL "")
    set(HG_REVISION_ID "hg-tip")
  else()
    if(HG_REVISION_TAG STREQUAL "null")
      set(HG_REVISION_TAG "trunk")
    elseif(HG_REVISION_TAG MATCHES "^r")
      string(SUBSTRING ${HG_REVISION_TAG} 1 -1 HG_REVISION_TAG)
    endif()
    if(HG_REVISION_DIST STREQUAL "0")
      set(HG_REVISION ${HG_REVISION_TAG})
    else()
      set(HG_REVISION
          "${HG_REVISION_TAG}+${HG_REVISION_DIST}-${HG_REVISION_ID}")
    endif()
  endif()
  set(LEMON_VERSION ${HG_REVISION})
else() # Hg not found and no environment variable is set
  set(LEMON_VERSION "unknown")
endif()

message(STATUS "LEMON Version: ${LEMON_VERSION}")

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

# Building + Testing
option(LEMON_BUILD_TESTING "Enable unit tests" ON)
option(LEMON_MAINTAINER "Enable Maintainer options." OFF)
option(LEMON_BUILD_DOCS "Enable docs using doxygen and ghostcript" OFF)
option(LEMON_BUILD_TOOLS "Enable project in tools directory" OFF)
option(LEMON_BUILD_DEMO "Enable project in demo directory" OFF)
option(LEMON_BUILD_ARCHIVE "Enable to create archive" OFF)
option(LEMON_BUILD_CPACK "Enable to create cpack" OFF)
option(LEMON_TEST_WITH_VALGRIND "Run the tests with valgrind" OFF)
option(LEMON_ENABLE_DEBUG "Display debug messages" OFF)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_BUILD_TYPE
    "${CMAKE_BUILD_TYPE}"
    CACHE
      STRING
      "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
      FORCE)

enable_testing()

# Add LEMON::LEMON target
include(base)

if(LEMON_MAINTAINER)
  add_custom_target(check ALL COMMAND ${CMAKE_CTEST_COMMAND})
else()
  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND})
endif()

if(LEMON_BUILD_TESTING)
  add_subdirectory(test)
endif()

if(LEMON_BUILD_TOOLS)
  add_subdirectory(tools)
endif()

if(LEMON_BUILD_DOCS)
  add_subdirectory(doc)
endif()

if(LEMON_BUILD_DEMO)
  add_subdirectory(demo)
endif()

if(LEMON_BUILD_ARCHIVE)
  include(archive)
endif()

if(LEMON_BUILD_CPACK)
  include(cpack)
endif()
