include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})

link_directories(${PROJECT_BINARY_DIR}/lemon)

set(DEMOS arg_parser_demo graph_to_eps_demo lgf_demo)

foreach(DEMO_NAME ${DEMOS})
  add_executable(${DEMO_NAME} ${DEMO_NAME}.cc)
  target_link_libraries(${DEMO_NAME} lemon)
endforeach()
