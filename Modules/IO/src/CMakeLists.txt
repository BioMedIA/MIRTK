# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
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

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/IOConfig.h
  ${BINARY_INCLUDE_DIR}/mirtk/IOExport.h
  GIPLImageReader.h
  GIPLImageWriter.h
  PGMImageReader.h
  PGMImageWriter.h
)

set(SOURCES
  GIPL.h
  GIPLImageReader.cc
  GIPLImageWriter.cc
  IOConfig.cc
  PGM.h
  PGMImageReader.cc
  PGMImageWriter.cc
)

set(DEPENDS
  LibCommon
  LibNumerics
  LibImage
)

if (NiftiCLib_FOUND)
  list(APPEND HEADERS
    NiftiImageInfo.h
    NiftiImageReader.h
    NiftiImageWriter.h
  )
  list(APPEND SOURCES
    NiftiImage.h
    NiftiImage.cc
    NiftiImageInfo.cc
    NiftiImageReader.cc
    NiftiImageWriter.cc
  )
  list(APPEND DEPENDS ${NiftiCLib_LIBRARIES})
endif ()

if (PNG_FOUND)
  add_definitions(${PNG_DEFINITIONS})
  list(APPEND HEADERS PNGImageWriter.h)
  list(APPEND SOURCES PNGImageWriter.cc)
  list(APPEND DEPENDS ${PNG_LIBRARIES})
endif ()

if (ZLIB_FOUND)
  list(APPEND DEPENDS ${ZLIB_LIBRARIES})
endif ()

if (VTK_FOUND)
  list(APPEND HEADERS PointSetIO.h)
  list(APPEND SOURCES PointSetIO.cc dfsurface.h)
  list(APPEND DEPENDS ${VTK_LIBRARIES})
endif ()

mirtk_add_library(AUTO_REGISTER)
