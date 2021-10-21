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
This file creates compressed files of the current version and places them in the
PROJECT_BINARY_DIR as tar.gz and zip files.
]]

set(ARCHIVE_BASE_NAME ${PROJECT_NAME})
string(TOLOWER ${ARCHIVE_BASE_NAME} ARCHIVE_BASE_NAME)
set(ARCHIVE_NAME ${ARCHIVE_BASE_NAME}-${PROJECT_VERSION})
add_custom_target(
  dist
  COMMAND cmake -E remove_directory ${ARCHIVE_NAME}
  COMMAND hg archive ${ARCHIVE_NAME}
  COMMAND cmake -E copy ${PROJECT_NAME}ConfigVersion.cmake
          ${ARCHIVE_NAME}/cmake/${PROJECT_NAME}ConfigVersion.cmake
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
