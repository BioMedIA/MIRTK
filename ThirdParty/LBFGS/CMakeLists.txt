basis_include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

basis_add_library(lbfgs STATIC
  src/arithmetic_ansi.h
  src/arithmetic_sse_double.h
  src/arithmetic_sse_float.h
  src/lbfgs.c
  include/lbfgs.h
)

basis_set_target_properties(lbfgs PROPERTIES
  ARCHIVE_OUTPUT_NAME       MIRTKlbfgs
  POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
)
