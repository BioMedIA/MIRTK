include_directories(
  "${CMAKE_CURRENT_SOURCE_DIR}/znzlib"
  "${CMAKE_CURRENT_SOURCE_DIR}/niftilib"
)

basis_add_library(nifticlib STATIC
  znzlib/znzlib.h
  znzlib/znzlib.c
  niftilib/nifti1.h
  niftilib/nifti1_io.h
  niftilib/nifti1_io.c
)

basis_set_target_properties(nifticlib PROPERTIES
  ARCHIVE_OUTPUT_NAME MIRTKnifticlib
)
