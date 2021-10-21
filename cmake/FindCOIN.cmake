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
This file finds the necessary components of the COIN-OR libraries.
Cbc, Clp, Cgl, Osi, and Vol
]]

set(COIN_ROOT_DIR
    ""
    CACHE PATH "COIN root directory")

find_path(COIN_INCLUDE_DIR coin/CoinUtilsConfig.h
          HINTS ${COIN_ROOT_DIR}/include)
find_library(
  COIN_CBC_LIBRARY
  NAMES Cbc libCbc
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_CBC_SOLVER_LIBRARY
  NAMES CbcSolver libCbcSolver
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_CGL_LIBRARY
  NAMES Cgl libCgl
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_CLP_LIBRARY
  NAMES Clp libClp
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_COIN_UTILS_LIBRARY
  NAMES CoinUtils libCoinUtils
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_OSI_LIBRARY
  NAMES Osi libOsi
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_OSI_CBC_LIBRARY
  NAMES OsiCbc libOsiCbc
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_OSI_CLP_LIBRARY
  NAMES OsiClp libOsiClp
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_OSI_VOL_LIBRARY
  NAMES OsiVol libOsiVol
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_VOL_LIBRARY
  NAMES Vol libVol
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)

find_library(
  COIN_ZLIB_LIBRARY
  NAMES z libz
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)
find_library(
  COIN_BZ2_LIBRARY
  NAMES bz2 libbz2
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)

find_library(
  COIN_PTHREADS_LIBRARY
  NAMES pthreads libpthreads
  HINTS ${COIN_ROOT_DIR}/lib/coin
  HINTS ${COIN_ROOT_DIR}/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  COIN
  DEFAULT_MSG
  COIN_INCLUDE_DIR
  COIN_CBC_LIBRARY
  COIN_CBC_SOLVER_LIBRARY
  COIN_CGL_LIBRARY
  COIN_CLP_LIBRARY
  COIN_COIN_UTILS_LIBRARY
  COIN_OSI_LIBRARY
  COIN_OSI_CBC_LIBRARY
  COIN_OSI_CLP_LIBRARY
  # COIN_OSI_VOL_LIBRARY
  # COIN_VOL_LIBRARY
)

if(COIN_FOUND)
  set(COIN_INCLUDE_DIRS ${COIN_INCLUDE_DIR})
  set(COIN_CLP_LIBRARIES "${COIN_CLP_LIBRARY};${COIN_COIN_UTILS_LIBRARY}")
  if(COIN_ZLIB_LIBRARY)
    set(COIN_CLP_LIBRARIES "${COIN_CLP_LIBRARIES};${COIN_ZLIB_LIBRARY}")
  endif(COIN_ZLIB_LIBRARY)
  if(COIN_BZ2_LIBRARY)
    set(COIN_CLP_LIBRARIES "${COIN_CLP_LIBRARIES};${COIN_BZ2_LIBRARY}")
  endif(COIN_BZ2_LIBRARY)
  if(COIN_PTHREADS_LIBRARY)
    set(COIN_CLP_LIBRARIES "${COIN_CLP_LIBRARIES};${COIN_PTHREADS_LIBRARY}")
  endif(COIN_PTHREADS_LIBRARY)
  set(COIN_CBC_LIBRARIES
      "${COIN_CBC_LIBRARY};${COIN_CBC_SOLVER_LIBRARY};${COIN_CGL_LIBRARY};${COIN_OSI_LIBRARY};${COIN_OSI_CBC_LIBRARY};${COIN_OSI_CLP_LIBRARY};${COIN_CLP_LIBRARIES}"
  )
  set(COIN_LIBRARIES ${COIN_CBC_LIBRARIES})
endif(COIN_FOUND)

mark_as_advanced(
  COIN_INCLUDE_DIR
  COIN_CBC_LIBRARY
  COIN_CBC_SOLVER_LIBRARY
  COIN_CGL_LIBRARY
  COIN_CLP_LIBRARY
  COIN_COIN_UTILS_LIBRARY
  COIN_OSI_LIBRARY
  COIN_OSI_CBC_LIBRARY
  COIN_OSI_CLP_LIBRARY
  COIN_OSI_VOL_LIBRARY
  COIN_VOL_LIBRARY
  COIN_ZLIB_LIBRARY
  COIN_BZ2_LIBRARY)
