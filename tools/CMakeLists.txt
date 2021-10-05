include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})

link_directories(${PROJECT_BINARY_DIR}/lemon)

add_executable(lgf-gen lgf-gen.cc)
target_link_libraries(lgf-gen lemon)

add_executable(dimacs-to-lgf dimacs-to-lgf.cc)
target_link_libraries(dimacs-to-lgf lemon)

add_executable(dimacs-solver dimacs-solver.cc)
target_link_libraries(dimacs-solver lemon)

install(TARGETS lgf-gen dimacs-to-lgf dimacs-solver RUNTIME DESTINATION bin
                                                            COMPONENT bin)

if(NOT WIN32)
  install(
    PROGRAMS ${CMAKE_CURRENT_SOURCE_DIR}/lemon-0.x-to-1.x.sh
    DESTINATION bin
    COMPONENT bin)
endif()
