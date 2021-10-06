cmake_minimum_required(VERSION 3.14)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(
  LEMON
  VERSION 1.4.0
  LANGUAGES CXX)

option(LEMON_ENABLE_GLPK "Enable GLPK solver backend." OFF)
option(LEMON_ENABLE_ILOG "Enable ILOG (CPLEX) solver backend." OFF)
option(LEMON_ENABLE_COIN "Enable COIN solver backend." OFF)
option(LEMON_ENABLE_SOPLEX "Enable SoPlex solver backend." OFF)

option(BUILD_DOCS "" OFF)
option(BUILD_TOOLS "" OFF)
option(BUILD_ARCHIVE "Enable to create archive" OFF)
option(BUILD_CPACK "Enable to create cpack" OFF)

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

include(base)

if(BUILD_TESTING AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  add_subdirectory(test)
endif()

if(BUILD_TOOLS AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  add_subdirectory(tools)
endif()

if(BUILD_DOCS AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  add_subdirectory(doc)
endif()

if(BUILD_DEMO AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  add_subdirectory(demo)
endif()

if(BUILD_ARCHIVE AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  include(archive)
endif()

if(BUILD_CPACK AND ${CMAKE_BUILD_TYPE} STREQUAL "Maintainer")
  include(cpack)
endif()
