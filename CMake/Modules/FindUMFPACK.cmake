#.rst
# FindUMFPACK
# -----------
#
# Find UMFPACK library which is part of SuiteSparse.
#
# This uses find_package(SuiteSparse COMPONENTS UMFPACK MODULE).

#=============================================================================
# Copyright 2016 Andreas Schuh <andreas.schuh.84@gmail.com>
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# (To distribute this file outside of CMake, substitute the full
#  License text for the above reference.)

# ------------------------------------------------------------------------------
# Look for UMFPACK library which is part of SuiteSparse
if (UMFPACK_ROOT)
  set(SuiteSparse_ROOT ${UMFPACK_ROOT})
elseif (UMFPACK_DIR)
  set(SuiteSparse_ROOT ${UMFPACK_DIR})
endif ()

find_package(SuiteSparse COMPONENTS UMFPACK MODULE QUIET)

set(UMFPACK_VERSION        ${SuiteSparse_VERSION})
set(UMFPACK_VERSION_MAJOR  ${SuiteSparse_VERSION_MAJOR})
set(UMFPACK_VERSION_MINOR  ${SuiteSparse_VERSION_MINOR})
set(UMFPACK_VERSION_PATCH  ${SuiteSparse_VERSION_PATCH})
set(UMFPACK_VERSION_STRING ${SuiteSparse_VERSION_STRING})

set(UMFPACK_INCLUDE_DIR  "${SuiteSparse_UMFPACK_INCLUDE_DIR}")
set(UMFPACK_INCLUDE_DIRS "${SuiteSparse_UMFPACK_INCLUDE_DIRS}")
set(UMFPACK_LIBRARY      "${SuiteSparse_UMFPACK_LIBRARY}")
set(UMFPACK_LIBRARIES    "${SuiteSparse_UMFPACK_LIBRARIES}")

# ------------------------------------------------------------------------------
# Handle QUIET, REQUIRED, and [EXACT] VERSION arguments and set SuiteSparse_FOUND
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  UMFPACK
  VERSION_VAR
    UMFPACK_VERSION
  REQUIRED_VARS
    SuiteSparse_FOUND
    UMFPACK_INCLUDE_DIR
    UMFPACK_LIBRARY
)
