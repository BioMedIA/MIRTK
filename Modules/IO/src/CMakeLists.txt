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
# GIFTI I/O library
if (NOT ITK_FOUND AND ZLIB_FOUND)

  set(METAIO_SOURCES
    meta/metaUtils.cxx
    meta/metaArrow.cxx
    meta/metaBlob.cxx
    meta/metaCommand.cxx
    meta/metaContour.cxx
    meta/metaDTITube.cxx
    meta/metaEllipse.cxx
    meta/metaGroup.cxx
    meta/metaGaussian.cxx
    meta/metaImage.cxx
    meta/metaImageUtils.cxx
    meta/metaLandmark.cxx
    meta/metaLine.cxx
    meta/metaMesh.cxx
    meta/metaObject.cxx
    meta/metaOutput.cxx
    meta/metaScene.cxx
    meta/metaSurface.cxx
    meta/metaTube.cxx
    meta/metaVesselTube.cxx
    meta/metaTransform.cxx
    meta/metaTubeGraph.cxx
    meta/metaForm.cxx
    meta/metaArray.cxx
    meta/metaFEMObject.cxx
  )

  set(METAIO_HEADERS 
    meta/metaArray.h
    meta/metaArrow.h
    meta/metaBlob.h
    meta/metaCommand.h
    meta/metaContour.h
    meta/metaDTITube.h
    meta/metaEllipse.h
    meta/metaEvent.h
    meta/metaFEMObject.h
    meta/metaForm.h
    meta/metaGaussian.h
    meta/metaGroup.h
    meta/metaImage.h
    meta/metaImageTypes.h
    meta/metaImageUtils.h
    meta/metaIOConfig.h
    meta/metaLandmark.h
    meta/metaLine.h
    meta/metaMesh.h
    meta/metaObject.h
    meta/metaOutput.h
    meta/metaScene.h
    meta/metaSurface.h
    meta/metaTransform.h
    meta/metaTubeGraph.h
    meta/metaTube.h
    meta/metaTypes.h
    meta/metaUtils.h
    meta/metaVesselTube.h
  )

  add_library(metaio OBJECT ${METAIO_HEADERS} ${METAIO_SOURCES})

  set_target_properties(metaio PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

endif ()

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

if (ITK_FOUND OR ZLIB_FOUND)
  list(APPEND HEADERS
    MetaImageReader.h
    MetaImageWriter.h
  )
  list(APPEND SOURCES
    MetaImageReader.cc
    MetaImageWriter.cc
  )
  if (TARGET ITKMetaIO)
    list(APPEND DEPENDS ITKMetaIO)
  else ()
    include_directories("meta")
    list(APPEND SOURCES $<TARGET_OBJECTS:metaio>)
    if (WIN32)
      if (BORLAND)
        list(APPEND DEPENDS import32)
      else ()
        list(APPEND DEPENDS comctl32 wsock32)
      endif ()
    endif ()
  endif ()
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
