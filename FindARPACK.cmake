find_library(ARPACK_LIBRARY
  NAMES parpack arpack
  PATHS ${ARPACK_DIR}/lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ARPACK DEFAULT_MSG ARPACK_LIBRARY)

if (ARPACK_FOUND)
  get_filename_component(ARPACK_DIR ${ARPACK_LIBRARY} PATH)
  get_filename_component(ARPACK_DIR ${ARPACK_DIR}     PATH)
  set(ARPACK_INCLUDE_DIR ${ARPACK_DIR}/include CACHE FILEPATH "ARPACK include directory.")
  mark_as_advanced(ARPACK_INCLUDE_DIR)
endif ()
mark_as_advanced(ARPACK_LIBRARY)
