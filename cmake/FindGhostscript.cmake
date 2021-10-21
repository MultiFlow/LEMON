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
This file finds the necessary components for Ghostscript.
]]

include(FindPackageHandleStandardArgs)

find_program(
  GHOSTSCRIPT_EXECUTABLE
  NAMES gs gswin32c
  PATHS "$ENV{ProgramFiles}/gs"
  PATH_SUFFIXES gs8.61/bin gs8.62/bin gs8.63/bin gs8.64/bin gs8.65/bin
  DOC "Ghostscript: PostScript and PDF language interpreter and previewer.")

find_package_handle_standard_args(Ghostscript DEFAULT_MSG
                                  GHOSTSCRIPT_EXECUTABLE)
