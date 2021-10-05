find_path(
  ILOG_ROOT_DIR
  NAMES cplex
  DOC "CPLEX STUDIO root directory"
  PATHS /opt/ibm/ILOG /usr/local/ibm/ILOG /usr/local/ILOG /usr/local/ilog
  PATHS "$ENV{HOME}/ILOG" "$ENV{HOME}/.local/ILOG"
  PATHS "$ENV{HOME}/ibm/ILOG" "$ENV{HOME}/.local/ibm/ILOG"
  PATHS "C:/Program Files/IBM/ILOG"
  PATH_SUFFIXES "CPLEX_Studio126" "CPLEX_Studio125" "CPLEX_Studio124"
                "CPLEX_Studio123" "CPLEX_Studio122" "CPLEX_Studio201"
  NO_DEFAULT_PATH)

if(WIN32)
  if(MSVC_VERSION STREQUAL "1400")
    set(ILOG_WIN_COMPILER "windows_vs2005")
  elseif(MSVC_VERSION STREQUAL "1500")
    set(ILOG_WIN_COMPILER "windows_vs2008")
  elseif(MSVC_VERSION STREQUAL "1600")
    set(ILOG_WIN_COMPILER "windows_vs2010")
  else()
    set(ILOG_WIN_COMPILER "windows_vs2008")
  endif()
  if(CMAKE_CL_64)
    set(ILOG_WIN_COMPILER "x64_${ILOG_WIN_COMPILER}")
    set(ILOG_WIN_PLATFORM "x64_win32")
  else()
    set(ILOG_WIN_COMPILER "x86_${ILOG_WIN_COMPILER}")
    set(ILOG_WIN_PLATFORM "x86_win32")
  endif()
endif()

find_path(
  ILOG_CPLEX_ROOT_DIR
  NAMES include/ilcplex
  HINTS ${ILOG_ROOT_DIR}/cplex ${ILOG_ROOT_DIR}/cplex121
        ${ILOG_ROOT_DIR}/cplex122 ${ILOG_ROOT_DIR}/cplex123
  DOC "CPLEX root directory"
  NO_DEFAULT_PATH)

find_path(
  ILOG_CONCERT_ROOT_DIR
  NAMES include/ilconcert
  HINTS ${ILOG_ROOT_DIR}/concert ${ILOG_ROOT_DIR}/concert29
  DOC "CONCERT root directory"
  NO_DEFAULT_PATH)

find_path(
  ILOG_CPLEX_INCLUDE_DIR ilcplex/cplex.h
  HINTS ${ILOG_CPLEX_ROOT_DIR}/include
  NO_DEFAULT_PATH)

find_path(
  ILOG_CONCERT_INCLUDE_DIR ilconcert/ilobasic.h
  HINTS ${ILOG_CONCERT_ROOT_DIR}/include
  NO_DEFAULT_PATH)

find_library(
  ILOG_CPLEX_LIBRARY cplex cplex121 cplex122 cplex123 cplex124
  HINTS ${ILOG_CPLEX_ROOT_DIR}/lib/x86_sles10_4.1/static_pic
        ${ILOG_CPLEX_ROOT_DIR}/lib/x86-64_sles10_4.1/static_pic
        ${ILOG_CPLEX_ROOT_DIR}/lib/x86_debian4.0_4.1/static_pic
        ${ILOG_CPLEX_ROOT_DIR}/lib/x86-64_debian4.0_4.1/static_pic
        ${ILOG_CPLEX_ROOT_DIR}/lib/x86_linux/static_pic
        ${ILOG_CPLEX_ROOT_DIR}/lib/x86-64_linux/static_pic
        ${ILOG_CPLEX_ROOT_DIR}/lib/${ILOG_WIN_COMPILER}/stat_mda
  NO_DEFAULT_PATH)

find_library(
  ILOG_CONCERT_LIBRARY concert
  HINTS ${ILOG_CONCERT_ROOT_DIR}/lib/x86_sles10_4.1/static_pic
        ${ILOG_CONCERT_ROOT_DIR}/lib/x86-64_sles10_4.1/static_pic
        ${ILOG_CONCERT_ROOT_DIR}/lib/x86_debian4.0_4.1/static_pic
        ${ILOG_CONCERT_ROOT_DIR}/lib/x86-64_debian4.0_4.1/static_pic
        ${ILOG_CONCERT_ROOT_DIR}/lib/x86_linux/static_pic
        ${ILOG_CONCERT_ROOT_DIR}/lib/x86-64_linux/static_pic
        ${ILOG_CONCERT_ROOT_DIR}/lib/${ILOG_WIN_COMPILER}/stat_mda
  NO_DEFAULT_PATH)

find_file(
  ILOG_CPLEX_DLL cplex121.dll cplex122.dll cplex123.dll cplex124.dll
  HINTS ${ILOG_CPLEX_ROOT_DIR}/bin/${ILOG_WIN_PLATFORM}
  NO_DEFAULT_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ILOG DEFAULT_MSG ILOG_CPLEX_LIBRARY
                                  ILOG_CPLEX_INCLUDE_DIR)

if(ILOG_FOUND)
  set(ILOG_INCLUDE_DIRS ${ILOG_CPLEX_INCLUDE_DIR} ${ILOG_CONCERT_INCLUDE_DIR})
  set(ILOG_LIBRARIES ${ILOG_CPLEX_LIBRARY} ${ILOG_CONCERT_LIBRARY})
  if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # SET(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread")
    set(ILOG_LIBRARIES ${ILOG_LIBRARIES} "m" "pthread" "dl")
  endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")
endif(ILOG_FOUND)

mark_as_advanced(ILOG_CPLEX_LIBRARY ILOG_CPLEX_INCLUDE_DIR ILOG_CPLEX_DLL
                 ILOG_CONCERT_LIBRARY ILOG_CONCERT_INCLUDE_DIR ILOG_CONCERT_DLL)
