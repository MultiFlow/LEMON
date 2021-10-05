# CPACK config (Basically for NSIS)
if(${CMAKE_SOURCE_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
  set(CPACK_PACKAGE_NAME ${PROJECT_NAME})
  set(CPACK_PACKAGE_VENDOR "EGRES")
  set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
      "LEMON - Library for Efficient Modeling and Optimization in Networks")
  set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")

  set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})

  set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROJECT_NAME} ${PROJECT_VERSION}")
  set(CPACK_PACKAGE_INSTALL_REGISTRY_KEY "${PROJECT_NAME} ${PROJECT_VERSION}")

  set(CPACK_COMPONENTS_ALL headers library html_documentation bin)

  set(CPACK_COMPONENT_HEADERS_DISPLAY_NAME "C++ headers")
  set(CPACK_COMPONENT_LIBRARY_DISPLAY_NAME "Dynamic-link library")
  set(CPACK_COMPONENT_BIN_DISPLAY_NAME "Command line utilities")
  set(CPACK_COMPONENT_HTML_DOCUMENTATION_DISPLAY_NAME "HTML documentation")

  set(CPACK_COMPONENT_HEADERS_DESCRIPTION "C++ header files")
  set(CPACK_COMPONENT_LIBRARY_DESCRIPTION "DLL and import library")
  set(CPACK_COMPONENT_BIN_DESCRIPTION "Command line utilities")
  set(CPACK_COMPONENT_HTML_DOCUMENTATION_DESCRIPTION
      "Doxygen generated documentation")

  set(CPACK_COMPONENT_HEADERS_DEPENDS library)

  set(CPACK_COMPONENT_HEADERS_GROUP "Development")
  set(CPACK_COMPONENT_LIBRARY_GROUP "Development")
  set(CPACK_COMPONENT_HTML_DOCUMENTATION_GROUP "Documentation")

  set(CPACK_COMPONENT_GROUP_DEVELOPMENT_DESCRIPTION
      "Components needed to develop software using LEMON")
  set(CPACK_COMPONENT_GROUP_DOCUMENTATION_DESCRIPTION "Documentation of LEMON")

  set(CPACK_ALL_INSTALL_TYPES Full Developer)

  set(CPACK_COMPONENT_HEADERS_INSTALL_TYPES Developer Full)
  set(CPACK_COMPONENT_LIBRARY_INSTALL_TYPES Developer Full)
  set(CPACK_COMPONENT_HTML_DOCUMENTATION_INSTALL_TYPES Full)

  set(CPACK_GENERATOR "NSIS")
  set(CPACK_NSIS_MUI_ICON "${PROJECT_SOURCE_DIR}/cmake/nsis/lemon.ico")
  set(CPACK_NSIS_MUI_UNIICON "${PROJECT_SOURCE_DIR}/cmake/nsis/uninstall.ico")
  # SET(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/cmake/nsis\\\\installer.bmp")
  set(CPACK_NSIS_INSTALLED_ICON_NAME "bin\\\\lemon.ico")
  set(CPACK_NSIS_DISPLAY_NAME
      "${CPACK_PACKAGE_INSTALL_DIRECTORY} ${PROJECT_NAME}")
  set(CPACK_NSIS_HELP_LINK "http:\\\\\\\\lemon.cs.elte.hu")
  set(CPACK_NSIS_URL_INFO_ABOUT "http:\\\\\\\\lemon.cs.elte.hu")
  set(CPACK_NSIS_CONTACT "lemon-user@lemon.cs.elte.hu")
  set(CPACK_NSIS_CREATE_ICONS_EXTRA
      "
    CreateShortCut \\\"$SMPROGRAMS\\\\$STARTMENU_FOLDER\\\\Documentation.lnk\\\" \\\"$INSTDIR\\\\share\\\\doc\\\\index.html\\\"
    ")
  set(CPACK_NSIS_DELETE_ICONS_EXTRA
      "
    !insertmacro MUI_STARTMENU_GETFOLDER Application $MUI_TEMP
    Delete \\\"$SMPROGRAMS\\\\$MUI_TEMP\\\\Documentation.lnk\\\"
    ")

  include(CPack)
endif()
