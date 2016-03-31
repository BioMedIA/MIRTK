cmake_minimum_required(VERSION 2.8.12)

set(VERSION_ARG)
if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
  set(VERSION_ARG VERSION 1.10)
endif ()

project(LBFGS ${VERSION_ARG})
set(PROJECT_SOVERSION ${PROJECT_VERSION_MAJOR})

include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")

basis_add_library(lbfgs STATIC
  src/arithmetic_ansi.h
  src/arithmetic_sse_double.h
  src/arithmetic_sse_float.h
  src/lbfgs.c
  include/lbfgs.h
)

set_target_properties(lbfgs PROPERTIES
  ARCHIVE_OUTPUT_NAME       MIRTKlbfgs
  DEBUG_POSTFIX             D
  POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
)
