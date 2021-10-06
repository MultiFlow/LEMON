set(NAME lemon)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.in
               ${CMAKE_CURRENT_BINARY_DIR}/config.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/lemon.pc.in
               ${CMAKE_CURRENT_BINARY_DIR}/lemon.pc @ONLY)

set(LEMON_SOURCES
    arg_parser.cc
    base.cc
    color.cc
    lp_base.cc
    lp_skeleton.cc
    random.cc
    bits/windows.cc)

if(LEMON_HAVE_GLPK)
  set(LEMON_SOURCES ${LEMON_SOURCES} glpk.cc)
  include_directories(${GLPK_INCLUDE_DIRS})
  if(WIN32)
    install(FILES ${GLPK_BIN_DIR}/glpk.dll DESTINATION bin)
    install(FILES ${GLPK_BIN_DIR}/libltdl3.dll DESTINATION bin)
    install(FILES ${GLPK_BIN_DIR}/zlib1.dll DESTINATION bin)
  endif()
endif()

if(LEMON_HAVE_CPLEX)
  set(LEMON_SOURCES ${LEMON_SOURCES} cplex.cc)
  include_directories(${ILOG_INCLUDE_DIRS})
endif()

if(LEMON_HAVE_CLP)
  set(LEMON_SOURCES ${LEMON_SOURCES} clp.cc)
  include_directories(${COIN_INCLUDE_DIRS})
endif()

if(LEMON_HAVE_CBC)
  set(LEMON_SOURCES ${LEMON_SOURCES} cbc.cc)
  include_directories(${COIN_INCLUDE_DIRS})
endif()

if(LEMON_HAVE_SOPLEX)
  set(LEMON_SOURCES ${LEMON_SOURCES} soplex.cc)
  include_directories(${SOPLEX_INCLUDE_DIRS})
endif()

file(GLOB _HDRS "*.h" "bits/*h" "concepts/*.h")

add_library(${NAME} ${LEMON_SOURCES} ${_HDRS})

# Compile options (needed for MacOS)
set_target_properties(
  ${NAME}
  PROPERTIES CXX_STANDARD 17
             CXX_STANDARD_REQUIRED ON
             CXX_EXTENSIONS OFF
             POSITION_INDEPENDENT_CODE ON)

target_compile_features(${NAME} PUBLIC cxx_std_17)

target_link_libraries(${NAME} ${GLPK_LIBRARIES} ${COIN_LIBRARIES}
                      ${ILOG_LIBRARIES} ${SOPLEX_LIBRARIES})

set_target_properties(
  ${NAME}
  PROPERTIES PUBLIC_HEADER "${_HDRS};${CMAKE_CURRENT_BINARY_DIR}/config.h"
             OUTPUT_NAME $<IF:$<BOOL:${UNIX}>,emon,${NAME}>
             VERSION ${PROJECT_VERSION}
             SOVERSION ${PROJECT_VERSION_MAJOR})

target_include_directories(
  ${NAME} PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
                 $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>)
# Add alias
add_library(LEMON::LEMON ALIAS ${NAME})

# Install
include(GNUInstallDirs)
install(
  TARGETS ${NAME}
  EXPORT ${PROJECT_NAME}Targets
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lemon
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

# Pkg-config
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/lemon.pc DESTINATION lib/pkgconfig)
