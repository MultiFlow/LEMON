set(ARCHIVE_BASE_NAME ${CMAKE_PROJECT_NAME})
string(TOLOWER ${ARCHIVE_BASE_NAME} ARCHIVE_BASE_NAME)
set(ARCHIVE_NAME ${ARCHIVE_BASE_NAME}-${PROJECT_VERSION})
add_custom_target(
  dist
  COMMAND cmake -E remove_directory ${ARCHIVE_NAME}
  COMMAND hg archive ${ARCHIVE_NAME}
  COMMAND cmake -E copy cmake/version.cmake ${ARCHIVE_NAME}/cmake/version.cmake
  COMMAND tar -czf ${ARCHIVE_BASE_NAME}-nodoc-${PROJECT_VERSION}.tar.gz
          ${ARCHIVE_NAME}
  COMMAND zip -r ${ARCHIVE_BASE_NAME}-nodoc-${PROJECT_VERSION}.zip
          ${ARCHIVE_NAME}
  COMMAND cmake -E copy_directory doc/html ${ARCHIVE_NAME}/doc/html
  COMMAND tar -czf ${ARCHIVE_NAME}.tar.gz ${ARCHIVE_NAME}
  COMMAND zip -r ${ARCHIVE_NAME}.zip ${ARCHIVE_NAME}
  COMMAND cmake -E copy_directory doc/html
          ${ARCHIVE_BASE_NAME}-doc-${PROJECT_VERSION}
  COMMAND tar -czf ${ARCHIVE_BASE_NAME}-doc-${PROJECT_VERSION}.tar.gz
          ${ARCHIVE_BASE_NAME}-doc-${PROJECT_VERSION}
  COMMAND zip -r ${ARCHIVE_BASE_NAME}-doc-${PROJECT_VERSION}.zip
          ${ARCHIVE_BASE_NAME}-doc-${PROJECT_VERSION}
  COMMAND cmake -E remove_directory ${ARCHIVE_NAME}
  COMMAND cmake -E remove_directory ${ARCHIVE_BASE_NAME}-doc-${PROJECT_VERSION}
  DEPENDS html
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR})
