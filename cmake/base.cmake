if(DEFINED ENV{LEMON_CXX_WARNING})
  set(CXX_WARNING $ENV{LEMON_CXX_WARNING})
else()
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(CXX_WARNING
        "-Wall -W -Wunused -Wformat=2 -Wctor-dtor-privacy -Wnon-virtual-dtor -Wno-char-subscripts -Wwrite-strings -Wno-char-subscripts -Wreturn-type -Wcast-qual -Wcast-align -Wsign-promo -Woverloaded-virtual -fno-strict-aliasing -Wold-style-cast -Wno-unknown-pragmas -Wno-unused-local-typedefs"
    )
    set(CMAKE_CXX_FLAGS_DEBUG CACHE STRING "-ggdb")
    set(CMAKE_C_FLAGS_DEBUG CACHE STRING "-ggdb")
  elseif(MSVC)
    # This part is unnecessary 'casue the same is set by the lemon/core.h. Still
    # keep it as an example.
    set(CXX_WARNING "/wd4250 /wd4355 /wd4503 /wd4800 /wd4996")
    # Suppressed warnings: C4250: 'class1' : inherits 'class2::member' via
    # dominance C4355: 'this' : used in base member initializer list C4503:
    # 'function' : decorated name length exceeded, name was truncated C4800:
    # 'type' : forcing value to bool 'true' or 'false' (performance warning)
    # C4996: 'function': was declared deprecated
  else()
    set(CXX_WARNING "-Wall")
  endif()
endif()

set(LEMON_CXX_WARNING_FLAGS
    ${CXX_WARNING}
    CACHE STRING "LEMON warning flags.")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LEMON_CXX_WARNING_FLAGS}")

if(MSVC)
  set(CMAKE_CXX_FLAGS_MAINTAINER
      "/WX ${CMAKE_CXX_FLAGS_DEBUG}"
      CACHE STRING "Flags used by the C++ compiler during maintainer builds.")
  set(CMAKE_C_FLAGS_MAINTAINER
      "/WX ${CMAKE_CXX_FLAGS_DEBUG}"
      CACHE STRING "Flags used by the C compiler during maintainer builds.")
  set(CMAKE_EXE_LINKER_FLAGS_MAINTAINER
      "${CMAKE_EXE_LINKER_FLAGS_DEBUG}"
      CACHE STRING "Flags used for linking binaries during maintainer builds.")
  set(CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
      "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}"
      CACHE
        STRING
        "Flags used by the shared libraries linker during maintainer builds.")
else()
  set(CMAKE_CXX_FLAGS_MAINTAINER
      "-Werror -ggdb -O0"
      CACHE STRING "Flags used by the C++ compiler during maintainer builds.")
  set(CMAKE_C_FLAGS_MAINTAINER
      "-Werror -O0"
      CACHE STRING "Flags used by the C compiler during maintainer builds.")
  set(CMAKE_EXE_LINKER_FLAGS_MAINTAINER
      "${CMAKE_EXE_LINKER_FLAGS_DEBUG}"
      CACHE STRING "Flags used for linking binaries during maintainer builds.")
  set(CMAKE_SHARED_LINKER_FLAGS_MAINTAINER
      "${CMAKE_SHARED_LINKER_FLAGS_DEBUG}"
      CACHE
        STRING
        "Flags used by the shared libraries linker during maintainer builds.")
endif()

mark_as_advanced(
  CMAKE_CXX_FLAGS_MAINTAINER CMAKE_C_FLAGS_MAINTAINER
  CMAKE_EXE_LINKER_FLAGS_MAINTAINER CMAKE_SHARED_LINKER_FLAGS_MAINTAINER)

# include(CheckTypeSize)
# check_type_size("long long" LONG_LONG)
# set(LEMON_HAVE_LONG_LONG ${HAVE_LONG_LONG})

find_package(Threads)

if(NOT LEMON_THREADING)
  if(CMAKE_USE_PTHREADS_INIT)
    set(LEMON_THREADING "Pthread")
  elseif(CMAKE_USE_WIN32_THREADS_INIT)
    set(LEMON_THREADING "Win32")
  else()
    set(LEMON_THREADING "None")
  endif()
endif()

set(LEMON_THREADING
    "${LEMON_THREADING}"
    CACHE STRING
          "Choose the threading library, options are: Pthread Win32 None."
          FORCE)

if(LEMON_THREADING STREQUAL "Pthread")
  set(LEMON_USE_PTHREAD ON)
elseif(LEMON_THREADING STREQUAL "Win32")
  set(LEMON_USE_WIN32_THREADS ON)
endif()

if(LEMON_ENABLE_GLPK)
  find_package(GLPK 4.33)
endif()
if(LEMON_ENABLE_ILOG)
  find_package(ILOG)
endif(LEMON_ENABLE_ILOG)
if()
  find_package(COIN)
endif()
if(LEMON_ENABLE_SOPLEX)
  find_package(SOPLEX)
endif()

if(GLPK_FOUND)
  set(LEMON_HAVE_LP ON)
  set(LEMON_HAVE_MIP ON)
  set(LEMON_HAVE_GLPK ON)
endif()
if(ILOG_FOUND)
  set(LEMON_HAVE_LP ON)
  set(LEMON_HAVE_MIP ON)
  set(LEMON_HAVE_CPLEX ON)
endif()
if(COIN_FOUND)
  set(LEMON_HAVE_LP ON)
  set(LEMON_HAVE_MIP ON)
  set(LEMON_HAVE_CLP ON)
  set(LEMON_HAVE_CBC ON)
endif()
if(SOPLEX_FOUND)
  set(LEMON_HAVE_LP ON)
  set(LEMON_HAVE_SOPLEX ON)
endif()

if(ILOG_FOUND)
  set(DEFAULT_LP "CPLEX")
  set(DEFAULT_MIP "CPLEX")
elseif(COIN_FOUND)
  set(DEFAULT_LP "CLP")
  set(DEFAULT_MIP "CBC")
elseif(GLPK_FOUND)
  set(DEFAULT_LP "GLPK")
  set(DEFAULT_MIP "GLPK")
elseif(SOPLEX_FOUND)
  set(DEFAULT_LP "SOPLEX")
endif()

if(NOT LEMON_DEFAULT_LP
   OR (NOT ILOG_FOUND AND (LEMON_DEFAULT_LP STREQUAL "CPLEX"))
   OR (NOT COIN_FOUND AND (LEMON_DEFAULT_LP STREQUAL "CLP"))
   OR (NOT GLPK_FOUND AND (LEMON_DEFAULT_LP STREQUAL "GLPK"))
   OR (NOT SOPLEX_FOUND AND (LEMON_DEFAULT_LP STREQUAL "SOPLEX")))
  set(LEMON_DEFAULT_LP
      ${DEFAULT_LP}
      CACHE STRING "Default LP solver backend (GLPK, CPLEX, CLP or SOPLEX)"
            FORCE)
else()
  set(LEMON_DEFAULT_LP
      ${DEFAULT_LP}
      CACHE STRING "Default LP solver backend (GLPK, CPLEX, CLP or SOPLEX)")
endif()

if(NOT LEMON_DEFAULT_MIP
   OR (NOT ILOG_FOUND AND (LEMON_DEFAULT_MIP STREQUAL "CPLEX"))
   OR (NOT COIN_FOUND AND (LEMON_DEFAULT_MIP STREQUAL "CBC"))
   OR (NOT GLPK_FOUND AND (LEMON_DEFAULT_MIP STREQUAL "GLPK")))
  set(LEMON_DEFAULT_MIP
      ${DEFAULT_MIP}
      CACHE STRING "Default MIP solver backend (GLPK, CPLEX or CBC)" FORCE)
else()
  set(LEMON_DEFAULT_MIP
      ${DEFAULT_MIP}
      CACHE STRING "Default MIP solver backend (GLPK, CPLEX or CBC)")
endif()

# add_custom_target(check ALL COMMAND ${CMAKE_CTEST_COMMAND})

add_subdirectory(lemon)

if(UNIX)
  install(FILES ${PROJECT_BINARY_DIR}/cmake/Config.cmake
          DESTINATION share/lemon/cmake)
elseif(WIN32)
  install(FILES ${PROJECT_BINARY_DIR}/cmake/LEMONConfig.cmake DESTINATION cmake)
endif()

configure_file(${PROJECT_SOURCE_DIR}/cmake/version.cmake.in
               ${CMAKE_BINARY_DIR}/cmake/version.cmake @ONLY)

include(GNUInstallDirs)
install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE LEMON::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
include(CMakePackageConfigHelpers)
configure_package_config_file(
  cmake/Config.cmake.in "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}")
write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  COMPATIBILITY SameMajorVersion)
install(
  FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
        "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
  COMPONENT Devel)
