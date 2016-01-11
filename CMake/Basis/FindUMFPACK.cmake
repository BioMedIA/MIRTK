# Umfpack lib usually requires linking to a blas library.
# It is up to the user of this module to find a BLAS and link to it.

if (UMFPACK_INCLUDE_DIR AND UMFPACK_LIBRARIES)
  set(UMFPACK_FIND_QUIETLY TRUE)
endif ()

find_path(UMFPACK_INCLUDE_DIR
  NAMES
    umfpack.h
  PATHS
    ${UMFPACK_DIR}/include
    $ENV{UMFPACKDIR}/include
    ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES
    suitesparse
    ufsparse
)

find_library(UMFPACK_LIBRARIES
  NAMES umfpack
  PATHS
    ${UMFPACK_DIR}/lib
    $ENV{UMFPACKDIR}/lib
    ${LIB_INSTALL_DIR}
)

if (UMFPACK_LIBRARIES)

  if (NOT UMFPACK_LIBDIR)
    get_filename_component(UMFPACK_LIBDIR ${UMFPACK_LIBRARIES} PATH)
  endif(NOT UMFPACK_LIBDIR)

  find_library(COLAMD_LIBRARY colamd PATHS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR}/lib ${LIB_INSTALL_DIR})
  if (COLAMD_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${COLAMD_LIBRARY})
  endif (COLAMD_LIBRARY)
  
  find_library(AMD_LIBRARY amd PATHS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR}/lib ${LIB_INSTALL_DIR})
  if (AMD_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${AMD_LIBRARY})
  endif (AMD_LIBRARY)

  find_library(SUITESPARSE_LIBRARY SuiteSparse PATHS ${UMFPACK_LIBDIR} $ENV{UMFPACKDIR}/lib ${LIB_INSTALL_DIR})
  if (SUITESPARSE_LIBRARY)
    set(UMFPACK_LIBRARIES ${UMFPACK_LIBRARIES} ${SUITESPARSE_LIBRARY})
  endif (SUITESPARSE_LIBRARY)

endif ()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UMFPACK DEFAULT_MSG
                                  UMFPACK_INCLUDE_DIR UMFPACK_LIBRARIES)

if (UMFPACK_FOUND)
  string (REGEX REPLACE "/(suitesparse|ufsparse)/?" "" UMFPACK_DIR "${UMFPACK_INCLUDE_DIR}")
  get_filename_component(UMFPACK_DIR "${UMFPACK_DIR}" PATH)
endif ()

mark_as_advanced(UMFPACK_INCLUDE_DIR UMFPACK_LIBRARIES AMD_LIBRARY COLAMD_LIBRARY SUITESPARSE_LIBRARY)
