###############################################################################
# Find FLANN -- http://www.cs.ubc.ca/research/flann
#
# This sets the following variables:
# FLANN_FOUND        - True if FLANN was found.
# FLANN_INCLUDE_DIRS - Directories containing the FLANN include files.
# FLANN_LIBRARIES    - Libraries needed to use FLANN.
# FLANN_DEFINITIONS  - Compiler flags for FLANN.

find_package(PkgConfig QUIET)
if (PkgConfig_FOUND)
  pkg_check_modules(PC_FLANN flann)
else ()
  set(PC_FLANN_INCLUDEDIR)
  set(PC_FLANN_INCLUDE_DIRS)
  set(PC_FLANN_LIBDIR)
  set(PC_FLANN_LIBRARY_DIRS)
  set(PC_FLANN_CFLAGS_OTHER)
endif ()

find_path(FLANN_INCLUDE_DIR flann/flann.hpp
  HINTS ${PC_FLANN_INCLUDEDIR} ${PC_FLANN_INCLUDE_DIRS}
  PATHS $ENV{FLANN_DIR} ${FLANN_DIR}
  PATH_SUFFIXES include
)

find_library(FLANN_LIBRARY flann
  HINTS ${PC_FLANN_LIBDIR} ${PC_FLANN_LIBRARY_DIRS}
  PATHS $ENV{FLANN_DIR} ${FLANN_DIR}
  PATH_SUFFIXES lib
)

set(FLANN_INCLUDE_DIRS ${FLANN_INCLUDE_DIR})
set(FLANN_LIBRARIES    ${FLANN_LIBRARY})
set(FLANN_DEFINITIONS  ${PC_FLANN_CFLAGS_OTHER})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FLANN DEFAULT_MSG FLANN_LIBRARY FLANN_INCLUDE_DIR)

mark_as_advanced(FLANN_LIBRARY FLANN_INCLUDE_DIR)
