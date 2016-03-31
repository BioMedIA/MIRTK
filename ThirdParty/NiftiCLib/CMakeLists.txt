basis_find_package(ZLIB)
if (ZLIB_FOUND)
  basis_use_package(ZLIB)
endif ()

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
basis_get_target_uid(target_uid nifticlib)

if (MSVC)
  # disable all function not secure warnings
  target_compile_definitions(${target_uid} PRIVATE _CRT_SECURE_NO_WARNINGS)
  # disable all conversion from double to float loss of data warnings
  target_compile_options(${target_uid} PRIVATE /wd4244)
  # disable all conversion from size_t to int loss of data and sign warnings
  target_compile_options(${target_uid} PRIVATE /wd4267)
endif ()

if (ZLIB_FOUND)
  target_link_libraries(${target_uid} ${ZLIB_LIBRARIES})
endif ()

set_target_properties(${target_uid} PROPERTIES
  ARCHIVE_OUTPUT_NAME       MIRTKnifticlib
  DEBUG_POSTFIX             D
  POSITION_INDEPENDENT_CODE ${BUILD_SHARED_LIBS}
)
