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
This file defines all the test targets for ctest.
Can also be used to run valgrind
]]

if(LEMON_TEST_WITH_VALGRIND)
  set(VALGRIND_FLAGS
      "--leak-check=full"
      CACHE STRING "Valgrind flags used by the tests.")
endif()

set(TESTS
    adaptors_test
    arc_look_up_test
    bellman_ford_test
    bfs_test
    bpgraph_test
    circulation_test
    connectivity_test
    counter_test
    dfs_test
    digraph_test
    dijkstra_test
    dim_test
    edge_set_test
    error_test
    euler_test
    fractional_matching_test
    gomory_hu_test
    graph_copy_test
    graph_test
    graph_utils_test
    hao_orlin_test
    heap_test
    kruskal_test
    lgf_reader_writer_test
    lgf_test
    maps_test
    matching_test
    max_cardinality_search_test
    max_clique_test
    max_flow_test
    min_cost_arborescence_test
    min_cost_flow_test
    min_mean_cycle_test
    multicommodity_flow_test
    nagamochi_ibaraki_test
    path_test
    planarity_test
    radix_sort_test
    random_test
    suurballe_test
    time_measure_test
    tsp_test
    unionfind_test
    vf2_test)

if(LEMON_HAVE_LP)
  if(LEMON_MAINTAINER)
    add_executable(lp_test lp_test.cc)
  else()
    add_executable(lp_test EXCLUDE_FROM_ALL lp_test.cc)
  endif()

  set(LP_TEST_LIBS lemon)

  if(LEMON_HAVE_GLPK)
    set(LP_TEST_LIBS ${LP_TEST_LIBS} ${GLPK_LIBRARIES})
  endif()
  if(LEMON_HAVE_CPLEX)
    set(LP_TEST_LIBS ${LP_TEST_LIBS} ${ILOG_LIBRARIES})
  endif()
  if(LEMON_HAVE_CLP)
    set(LP_TEST_LIBS ${LP_TEST_LIBS} ${COIN_CLP_LIBRARIES})
  endif()
  if(LEMON_HAVE_SOPLEX)
    set(LP_TEST_LIBS ${LP_TEST_LIBS} ${SOPLEX_LIBRARIES})
  endif()

  target_link_libraries(lp_test ${LP_TEST_LIBS})
  add_test(lp_test lp_test)
  # add_dependencies(check lp_test)

  if(WIN32 AND LEMON_HAVE_GLPK)
    get_target_property(TARGET_LOC lp_test LOCATION)
    get_filename_component(TARGET_PATH ${TARGET_LOC} PATH)
    add_custom_command(
      TARGET lp_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${GLPK_BIN_DIR}/glpk.dll ${TARGET_PATH}
      COMMAND ${CMAKE_COMMAND} -E copy ${GLPK_BIN_DIR}/libltdl3.dll
              ${TARGET_PATH}
      COMMAND ${CMAKE_COMMAND} -E copy ${GLPK_BIN_DIR}/zlib1.dll ${TARGET_PATH})
  endif()

  if(WIN32 AND LEMON_HAVE_CPLEX)
    get_target_property(TARGET_LOC lp_test LOCATION)
    get_filename_component(TARGET_PATH ${TARGET_LOC} PATH)
    add_custom_command(
      TARGET lp_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${ILOG_CPLEX_DLL} ${TARGET_PATH})
  endif()
endif()

if(LEMON_HAVE_MIP)
  if(LEMON_MAINTAINER)
    add_executable(mip_test mip_test.cc)
  else()
    add_executable(mip_test EXCLUDE_FROM_ALL mip_test.cc)
  endif()

  set(MIP_TEST_LIBS lemon)

  if(LEMON_HAVE_GLPK)
    set(MIP_TEST_LIBS ${MIP_TEST_LIBS} ${GLPK_LIBRARIES})
  endif()
  if(LEMON_HAVE_CPLEX)
    set(MIP_TEST_LIBS ${MIP_TEST_LIBS} ${ILOG_LIBRARIES})
  endif()
  if(LEMON_HAVE_CBC)
    set(MIP_TEST_LIBS ${MIP_TEST_LIBS} ${COIN_CBC_LIBRARIES})
  endif()

  target_link_libraries(mip_test ${MIP_TEST_LIBS})
  add_test(mip_test mip_test)
  add_dependencies(check mip_test)

  if(WIN32 AND LEMON_HAVE_GLPK)
    get_target_property(TARGET_LOC mip_test LOCATION)
    get_filename_component(TARGET_PATH ${TARGET_LOC} PATH)
    add_custom_command(
      TARGET mip_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${GLPK_BIN_DIR}/glpk.dll ${TARGET_PATH}
      COMMAND ${CMAKE_COMMAND} -E copy ${GLPK_BIN_DIR}/libltdl3.dll
              ${TARGET_PATH}
      COMMAND ${CMAKE_COMMAND} -E copy ${GLPK_BIN_DIR}/zlib1.dll ${TARGET_PATH})
  endif()

  if(WIN32 AND LEMON_HAVE_CPLEX)
    get_target_property(TARGET_LOC mip_test LOCATION)
    get_filename_component(TARGET_PATH ${TARGET_LOC} PATH)
    add_custom_command(
      TARGET mip_test
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy ${ILOG_CPLEX_DLL} ${TARGET_PATH})
  endif()
endif()

foreach(TEST_NAME ${TESTS})
  add_executable(${TEST_NAME} ${TEST_NAME}.cc)
  target_link_libraries(${TEST_NAME} LEMON::LEMON)
  target_include_directories(
    ${TEST_NAME}
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}>
           $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
           $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
           $<INSTALL_INTERFACE:include/lemon>)
  # Inherit lemon target properties (flags, cxx standard)
  if(APPLE)
    set_target_properties(
      ${TEST_NAME} PROPERTIES INSTALL_RPATH
                              "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
  elseif(UNIX AND NOT APPLE)
    set_target_properties(
      ${TEST_NAME} PROPERTIES INSTALL_RPATH
                              "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
  endif()

  if(LEMON_TEST_WITH_VALGRIND)
    add_test(${TEST_NAME} valgrind --error-exitcode=1 ${VALGRIND_FLAGS}
             ${CMAKE_CURRENT_BINARY_DIR}/${TEST_NAME})
  else()
    add_test(${TEST_NAME} ${TEST_NAME})
  endif()
  add_dependencies(check ${TEST_NAME})
endforeach()
