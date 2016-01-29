# ============================================================================
# Copyright (c) 2011-2012 University of Pennsylvania
# Copyright (c) 2013-2014 Andreas Schuh
# All rights reserved.
#
# See COPYING file for license information or visit
# http://opensource.andreasschuh.com/cmake-basis/download.html#license
# ============================================================================

##############################################################################
# @file  FindLIBLINEAR.cmake
# @brief Find LIBLINEAR package.
#
# @par Input varibales:
# <table border="0">
#   <tr>
#     @tp @b LIBLINEAR_DIR @endtp
#     <td>The LIBLINEAR package files are searched primarily under the specified
#         root directory. This variable can be alternatively set as environment
#         variable.</td>
#   </tr>
#   <tr>
#     @tp @b MEX_EXT @endtp
#     <td>The extension of MEX-files. If this variable is not set and the
#         basis_mexext command is available, it is invoked to determine the
#         extension automatically. Otherwise, the MEX extension defaults to
#         "mexa64".</td>
#   </tr>
# </table>
#
# @par Output variables:
# <table border="0">
#   <tr>
#     @tp @b LIBLINEAR_FOUND @endtp
#     <td>Whether the package was found and the following CMake variables are valid.</td>
#   </tr>
#   <tr>
#     @tp @b LIBLINEAR_libsvmwrite_MEX @endtp
#     <td>The libsvmwrite MEX-file.</td>
#   </tr>
#   <tr>
#     @tp @b LIBLINEAR_libsvmread_MEX @endtp
#     <td>The libsvmread MEX-file.</td>
#   </tr>
#   <tr>
#     @tp @b LIBLINEAR_predict_MEX @endtp
#     <td>The predict MEX-file.</td>
#   </tr>
#   <tr>
#     @tp @b LIBLINEAR_train_MEX @endtp
#     <td>The train MEX-file.</td>
#   </tr>
#   <tr>
#     @tp @b LIBLINEAR_MEX_FILES @endtp
#     <td>List of MEX-files (non-cached).</td>
#   </tr>
# </table>
#
# @ingroup CMakeFindModules
##############################################################################

# ----------------------------------------------------------------------------
# initialize search
if (NOT LIBLINEAR_DIR)
  set (LIBLINEAR_DIR "$ENV{LIBLINEAR_DIR}" CACHE PATH "Installation prefix for LIBLINEAR." FORCE)
endif ()

if (NOT MEX_EXT)
  if (COMMAND basis_mexext)
    basis_mexext (MEX_EXT)
  else ()
    set (MEX_EXT "mexa64")
  endif ()
endif ()

#--------------------------------------------------------------
# find paths/files
if (LIBLINEAR_DIR)

  find_file (
    LIBLINEAR_libsvmwrite_MEX
      NAMES         libsvmwrite.${MEX_EXT}
      HINTS         ${LIBLINEAR_DIR}
      PATH_SUFFIXES "matlab"
      DOC           "The libsvmwrite MEX-file of the LIBLINEAR library."
      NO_DEFAULT_PATH
  )

  find_file (
    LIBLINEAR_libsvmread_MEX
      NAMES         libsvmread.${MEX_EXT}
      HINTS         ${LIBLINEAR_DIR}
      PATH_SUFFIXES "matlab"
      DOC           "The libsvmread MEX-file of the LIBLINEAR library."
      NO_DEFAULT_PATH
  )

  find_file (
    LIBLINEAR_train_MEX
      NAMES         train.${MEX_EXT}
      HINTS         ${LIBLINEAR_DIR}
      PATH_SUFFIXES "matlab"
      DOC           "The train MEX-file of the LIBLINEAR library."
      NO_DEFAULT_PATH
  )

  find_file (
    LIBLINEAR_predict_MEX
      NAMES         predict.${MEX_EXT}
      HINTS         ${LIBLINEAR_DIR}
      PATH_SUFFIXES "matlab"
      DOC           "The predict MEX-file of the LIBLINEAR library."
      NO_DEFAULT_PATH
  )

else ()

  find_file (
    LIBLINEAR_libsvmwrite_MEX
      NAMES         libsvmwrite.${MEX_EXT}
      PATH_SUFFIXES "matlab"
      DOC           "The libsvmwrite MEX-file of the LIBLINEAR library."
  )

  find_file (
    LIBLINEAR_libsvmread_MEX
      NAMES         libsvmread.${MEX_EXT}
      PATH_SUFFIXES "matlab"
      DOC           "The libsvmread MEX-file of the LIBLINEAR library."
  )

  find_file (
    LIBLINEAR_train_MEX
      NAMES         train.${MEX_EXT}
      PATH_SUFFIXES "matlab"
      DOC           "The train MEX-file of the LIBLINEAR library."
  )

  find_file (
    LIBLINEAR_predict_MEX
      NAMES         predict.${MEX_EXT}
      PATH_SUFFIXES "matlab"
      DOC           "The predict MEX-file of the LIBLINEAR library."
  )

endif ()

mark_as_advanced (LIBLINEAR_libsvmread_MEX)
mark_as_advanced (LIBLINEAR_libsvmwrite_MEX)
mark_as_advanced (LIBLINEAR_train_MEX)
mark_as_advanced (LIBLINEAR_predict_MEX)

set (LIBLINEAR_MEX_FILES)
if (LIBLINEAR_libsvmread_MEX)
  list (APPEND LIBLINEAR_MEX_FILES "${LIBLINEAR_libsvmread_MEX}")
endif ()
if (LIBLINEAR_libsvmwrite_MEX)
  list (APPEND LIBLINEAR_MEX_FILES "${LIBLINEAR_libsvmwrite_MEX}")
endif ()
if (LIBLINEAR_train_MEX)
  list (APPEND LIBLINEAR_MEX_FILES "${LIBLINEAR_train_MEX}")
endif ()
if (LIBLINEAR_predict_MEX)
  list (APPEND LIBLINEAR_MEX_FILES "${LIBLINEAR_predict_MEX}")
endif ()

# ----------------------------------------------------------------------------
# handle the QUIETLY and REQUIRED arguments and set *_FOUND to TRUE
# if all listed variables are found or TRUE
include (FindPackageHandleStandardArgs)

find_package_handle_standard_args (
  LIBLINEAR
# MESSAGE
    DEFAULT_MSG
# VARIABLES
     LIBLINEAR_libsvmwrite_MEX
     LIBLINEAR_libsvmread_MEX
     LIBLINEAR_predict_MEX
     LIBLINEAR_train_MEX
)

# ----------------------------------------------------------------------------
# set LIBLINEAR_DIR
if (NOT LIBLINEAR_DIR AND LIBLINEAR_FOUND)
  string (REGEX REPLACE "matlab/[^/]+" "" LIBLINEAR_PREFIX "${LIBLINEAR_train_MEX}")
  set (LIBLINEAR_DIR "${LIBLINEAR_PREFIX}" CACHE PATH "Installation prefix for LIBLINEAR." FORCE)
  unset (LIBLINEAR_PREFIX)
endif ()
