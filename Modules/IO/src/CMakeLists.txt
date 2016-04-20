# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2016 Imperial College London
# Copyright 2013-2016 Andreas Schuh
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of MIRTK IO library.
##############################################################################

# ----------------------------------------------------------------------------
# NIfTI I/O library
if (NOT NiftiCLib_FOUND)

  add_library(nifticlib OBJECT
    nifti/nifti1.h
    nifti/nifti2.h
    nifti/nifti2_io.h
    nifti/nifti2_io.cc
    nifti/znzlib.h
    nifti/znzlib.cc
  )
 
  set_target_properties(nifticlib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
 
  if (MSVC)
    # disable all function not secure warnings
    target_compile_definitions(nifticlib PRIVATE _CRT_SECURE_NO_WARNINGS)
    # disable all conversion from double to float loss of data warnings
    target_compile_options(nifticlib PRIVATE /wd4244)
    # disable all conversion from size_t to int loss of data and sign warnings
    target_compile_options(nifticlib PRIVATE /wd4267)
  endif ()

endif ()
# ----------------------------------------------------------------------------
# GIFTI I/O library
if (NOT GiftiCLib_FOUND AND VTK_FOUND AND EXPAT_FOUND)

  add_library(gifticlib OBJECT
    gifti/gifti_io.h
    gifti/gifti_io.cc
    gifti/gifti_xml.h
    gifti/gifti_xml.cc
  )
 
  set_target_properties(gifticlib PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

  if (MSVC)
    # disable all function not secure warnings
    target_compile_definitions(gifticlib PRIVATE _CRT_SECURE_NO_WARNINGS)
    # disable all conversion from double to float loss of data warnings
    target_compile_options(gifticlib PRIVATE /wd4244)
    # disable all conversion from size_t to int loss of data and sign warnings
    target_compile_options(gifticlib PRIVATE /wd4267)
  endif ()
 
endif ()
# ----------------------------------------------------------------------------
# MIRTK I/O library
set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/IOConfig.h
  ${BINARY_INCLUDE_DIR}/mirtk/IOExport.h
  GIPLImageReader.h
  GIPLImageWriter.h
  NiftiImageInfo.h
  NiftiImageReader.h
  NiftiImageWriter.h
  PGMImageReader.h
  PGMImageWriter.h
)

set(SOURCES
  IOConfig.cc
  GIPL.h
  GIPLImageReader.cc
  GIPLImageWriter.cc
  NiftiImage.h
  NiftiImage.cc
  NiftiImageInfo.cc
  NiftiImageReader.cc
  NiftiImageWriter.cc
  PGM.h
  PGMImageReader.cc
  PGMImageWriter.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
)

if (PNG_FOUND)
  add_definitions(${PNG_DEFINITIONS})
  list(APPEND HEADERS PNGImageWriter.h)
  list(APPEND SOURCES PNGImageWriter.cc)
  list(APPEND DEPENDS ${PNG_LIBRARIES})
endif ()

if (NiftiCLib_FOUND)
  list(APPEND DEPENDS ${NiftiCLib_LIBRARIES})
else ()
  list(APPEND SOURCES $<TARGET_OBJECTS:nifticlib>)
endif ()

if (VTK_FOUND)
  list(APPEND HEADERS PointSetIO.h)
  list(APPEND SOURCES PointSetIO.cc brainsuite/dfsurface.h)
  list(APPEND DEPENDS ${VTK_LIBRARIES})
  if (GiftiCLib_FOUND)
    list(APPEND DEPENDS ${GiftiCLib_LIBRARIES})
  elseif (EXPAT_FOUND)
    list(APPEND SOURCES $<TARGET_OBJECTS:gifticlib>)
    list(APPEND DEPENDS ${EXPAT_LIBRARIES})
  endif ()
endif ()

if (ZLIB_FOUND)
  list(APPEND DEPENDS ${ZLIB_LIBRARIES})
endif ()

mirtk_add_library(AUTO_REGISTER)
