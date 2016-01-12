# ============================================================================
# Copyright 2007-2009 Kitware, Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#  notice, this list of conditions and the following disclaimer in the
#  documentation and/or other materials provided with the distribution.
#
# * Neither the names of Kitware, Inc., the Insight Software Consortium,
#  nor the names of their contributors may be used to endorse or promote
#  products derived from this software without specific prior written
#  permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 ============================================================================

##############################################################################
# @file  FindBLAS.cmake
# @brief Find BLAS library.
#
# This module finds an installed fortran library that implements the BLAS
# linear-algebra interface (see http://www.netlib.org/blas/).
# The list of libraries searched for is taken
# from the autoconf macro file, acx_blas.m4 (distributed at
# http://ac-archive.sourceforge.net/ac-archive/acx_blas.html).
#
# Modified by Andreas Schuh to enable the use at SBIA, where an ATLAS C library
# is installed which contains the symbols without trailing _ character, i.e.,
# instead of checking the existence of the cblas_dgemm_ function, the
# existence of the cblas_dgemm function has to be checked. Moreover, added
# code to mark variable as advanced and only show them to the user if
# no required library was found. If the found library is cblas, the corresponding
# header file cblas.h is searched as well. Therefore, added the BLAS_INCLUDE_DIR
# variable which is only defined if required.
#
# This module sets the following variables:
#  BLAS_FOUND - set to true if a library implementing the BLAS interface
#    is found
#  BLAS_LINKER_FLAGS - uncached list of required linker flags (excluding -l
#    and -L).
#  BLAS_LIBRARIES - uncached list of libraries (using full path name) to
#    link against to use BLAS
#  BLAS_INCLUDE_DIR - uncached list of include directories for C libraries
#  BLAS95_LIBRARIES - uncached list of libraries (using full path name)
#    to link against to use BLAS95 interface
#  BLAS95_FOUND - set to true if a library implementing the BLAS f95 interface
#    is found
#  BLA_STATIC  if set on this determines what kind of linkage we do (static)
#  BLA_VENDOR  if set checks only the specified vendor, if not set checks
#     all the possibilities
#  BLA_F95     if set on tries to find the f95 interfaces for BLAS/LAPACK
#
# List of vendors (BLA_VENDOR) valid in this module
# ATLAS, PhiPACK,CXML,DXML,SunPerf,SCSL,SGIMATH,IBMESSL,Intel10_32 (intel mkl v10 32 bit),Intel10_64lp (intel mkl v10 64 bit,lp thread model, lp64 model),
# Intel( older versions of mkl 32 and 64 bit), ACML,ACML_MP,Apple, NAS, Generic
# C/CXX should be enabled to use Intel mkl
#
# @ingroup CMakeFindModules
##############################################################################

set(BLAS_LIBRARIES_VARS) # set by check_fortran_libraries() to a list of
                         # variable names of each searched library such
                         # that these libraries can be made non-advanced
                         # in case no library was found

include(CheckFunctionExists)
include(CheckFortranFunctionExists)

# Check the language being used
get_property( _LANGUAGES_ GLOBAL PROPERTY ENABLED_LANGUAGES )
if( _LANGUAGES_ MATCHES Fortran )
  set( _CHECK_FORTRAN TRUE )
elseif( (_LANGUAGES_ MATCHES C) OR (_LANGUAGES_ MATCHES CXX) )
  set( _CHECK_FORTRAN FALSE )
else()
  if(BLAS_FIND_REQUIRED)
    message(FATAL_ERROR "FindBLAS requires Fortran, C, or C++ to be enabled.")
  else(BLAS_FIND_REQUIRED)
    message(STATUS "Looking for BLAS... - NOT found (Unsupported languages)")
    return()
  endif(BLAS_FIND_REQUIRED)
endif( )

macro(Check_Fortran_Libraries LIBRARIES _prefix _name _flags _list _threads)
# This macro checks for the existence of the combination of fortran libraries
# given by _list.  If the combination is found, this macro checks (using the
# Check_Fortran_Function_Exists macro) whether can link against that library
# combination using the name of a routine given by _name using the linker
# flags given by _flags.  If the combination of libraries is found and passes
# the link test, LIBRARIES is set to the list of complete library paths that
# have been found.  Otherwise, LIBRARIES is set to FALSE.

# N.B. _prefix is the prefix applied to the names of all cached variables that
# are generated internally and marked advanced by this macro.

set(_libraries_work TRUE)
set(${LIBRARIES})
set(_combined_name)
foreach(_library ${_list})
  set(_combined_name ${_combined_name}_${_library})

  if(_libraries_work)
   if ( WIN32 )
    if(BLA_STATIC)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".lib;.dll")
    endif(BLA_STATIC)
    find_library(${_prefix}_${_library}_LIBRARY
    NAMES ${_library}
    PATHS ENV LIB
    )
   endif ( WIN32 )

   if ( APPLE )
    if(BLA_STATIC)
     set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.dylib")
    endif(BLA_STATIC)
    find_library(${_prefix}_${_library}_LIBRARY
    NAMES ${_library}
    PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV DYLD_LIBRARY_PATH
    )

   else ( APPLE )
    if(BLA_STATIC)
      set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;.so")
    endif(BLA_STATIC)
    find_library(${_prefix}_${_library}_LIBRARY
    NAMES ${_library}
    PATHS /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 ENV LD_LIBRARY_PATH
    )
   endif( APPLE )
    mark_as_advanced(${_prefix}_${_library}_LIBRARY)
    list (APPEND BLAS_LIBRARIES_VARS ${_prefix}_${_library}_LIBRARY)
    set(${LIBRARIES} ${${LIBRARIES}} ${${_prefix}_${_library}_LIBRARY})
    set(_libraries_work ${${_prefix}_${_library}_LIBRARY})
  endif(_libraries_work)
endforeach(_library ${_list})
if(_libraries_work)
  # Test this combination of libraries.
  set(CMAKE_REQUIRED_LIBRARIES ${_flags} ${${LIBRARIES}} ${_threads})
#  message("DEBUG: CMAKE_REQUIRED_LIBRARIES = ${CMAKE_REQUIRED_LIBRARIES}")
  if (_CHECK_FORTRAN)
    check_fortran_function_exists("${_name}" ${_prefix}${_combined_name}_WORKS)
  else()
    check_function_exists("${_name}" ${_prefix}${_combined_name}_WORKS)
    if (NOT ${_prefix}${_combined_name}_WORKS)
      check_function_exists("${_name}_" ${_prefix}${_combined_name}_WORKS)
    endif ()
  endif()
  set(CMAKE_REQUIRED_LIBRARIES)
  mark_as_advanced(${_prefix}${_combined_name}_WORKS)
  set(_libraries_work ${${_prefix}${_combined_name}_WORKS})
endif(_libraries_work)
if(NOT _libraries_work)
  set(${LIBRARIES} FALSE)
endif(NOT _libraries_work)
#message("DEBUG: ${LIBRARIES} = ${${LIBRARIES}}")
endmacro(Check_Fortran_Libraries)

set(BLAS_LINKER_FLAGS)
set(BLAS_LIBRARIES)
set(BLAS95_LIBRARIES)
if ($ENV{BLA_VENDOR} MATCHES ".+")
  set(BLA_VENDOR $ENV{BLA_VENDOR})
else ($ENV{BLA_VENDOR} MATCHES ".+")
  if(NOT BLA_VENDOR)
    set(BLA_VENDOR "All")
  endif(NOT BLA_VENDOR)
endif ($ENV{BLA_VENDOR} MATCHES ".+")

if (BLA_VENDOR STREQUAL "ATLAS" OR BLA_VENDOR STREQUAL "All")
 if(NOT BLAS_LIBRARIES)
  # BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  cblas_dgemm
  ""
  "cblas;f77blas;atlas"
  ""
  )
  if (BLAS_LIBRARIES_VARS MATCHES "cblas")
    find_path (BLAS_INCLUDE_DIR NAMES cblas.h DOC "Include directory of the cblas.h header file.")
  endif ()
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "ATLAS" OR BLA_VENDOR STREQUAL "All")

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if (BLA_VENDOR STREQUAL "PhiPACK" OR BLA_VENDOR STREQUAL "All")
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  ""
  "sgemm;dgemm;blas"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "PhiPACK" OR BLA_VENDOR STREQUAL "All")

# BLAS in Alpha CXML library?
if (BLA_VENDOR STREQUAL "CXML" OR BLA_VENDOR STREQUAL "All")
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  ""
  "cxml"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "CXML" OR BLA_VENDOR STREQUAL "All")

# BLAS in Alpha DXML library? (now called CXML, see above)
if (BLA_VENDOR STREQUAL "DXML" OR BLA_VENDOR STREQUAL "All")
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  ""
  "dxml"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "DXML" OR BLA_VENDOR STREQUAL "All")

# BLAS in Sun Performance library?
if (BLA_VENDOR STREQUAL "SunPerf" OR BLA_VENDOR STREQUAL "All")
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  "-xlic_lib=sunperf"
  "sunperf;sunmath"
  ""
  )
  if(BLAS_LIBRARIES)
    set(BLAS_LINKER_FLAGS "-xlic_lib=sunperf")
  endif(BLAS_LIBRARIES)
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "SunPerf" OR BLA_VENDOR STREQUAL "All")

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if (BLA_VENDOR STREQUAL "SCSL" OR BLA_VENDOR STREQUAL "All")
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  ""
  "scsl"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "SCSL" OR BLA_VENDOR STREQUAL "All")

# BLAS in SGIMATH library?
if (BLA_VENDOR STREQUAL "SGIMATH" OR BLA_VENDOR STREQUAL "All")
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  ""
  "complib.sgimath"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "SGIMATH" OR BLA_VENDOR STREQUAL "All")

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if (BLA_VENDOR STREQUAL "IBMESSL" OR BLA_VENDOR STREQUAL "All")
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  ""
  "essl;blas"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "IBMESSL" OR BLA_VENDOR STREQUAL "All")

#BLAS in acml library?
if (BLA_VENDOR STREQUAL "ACML" OR BLA_VENDOR STREQUAL "ACML_MP" OR BLA_VENDOR STREQUAL "All")
# the patch from Chuck Atkins:
 if( ((_BLAS_VENDOR STREQUAL "ACML") AND (NOT BLAS_ACML_LIB_DIRS)) OR
    ((_BLAS_VENDOR STREQUAL "ACML_MP") AND (NOT BLAS_ACML_MP_LIB_DIRS)) )
   if( WIN32 )
    file( GLOB _ACML_ROOT "C:/AMD/acml*/ACML-EULA.txt" )
   else()
    file( GLOB _ACML_ROOT "/opt/acml*/ACML-EULA.txt" )
   endif()
   if( _ACML_ROOT )
    get_filename_component( _ACML_ROOT ${_ACML_ROOT} PATH )
    if( SIZEOF_INTEGER EQUAL 8 )
     set( _ACML_PATH_SUFFIX "_int64" )
    else()
    set( _ACML_PATH_SUFFIX "" )
   endif()
   if( CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" )
    set( _ACML_COMPILER32 "ifort32" )
    set( _ACML_COMPILER64 "ifort64" )
   elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro" )
    set( _ACML_COMPILER32 "sun32" )
    set( _ACML_COMPILER64 "sun64" )
   elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "PGI" )
    set( _ACML_COMPILER32 "pgi32" )
    if( WIN32 )
     set( _ACML_COMPILER64 "win64" )
    else()
     set( _ACML_COMPILER64 "pgi64" )
    endif()
   elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "Open64" )
    # 32 bit builds not supported on Open64 but for code simplicity
    # We'll just use the same directory twice
    set( _ACML_COMPILER32 "open64_64" )
    set( _ACML_COMPILER64 "open64_64" )
   elseif( CMAKE_Fortran_COMPILER_ID STREQUAL "NAG" )
    set( _ACML_COMPILER32 "nag32" )
    set( _ACML_COMPILER64 "nag64" )
   else() #if( CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" )
    set( _ACML_COMPILER32 "gfortran32" )
    set( _ACML_COMPILER64 "gfortran64" )
   endif()

   if( _BLAS_VENDOR STREQUAL "ACML_MP" )
    set(_ACML_MP_LIB_DIRS
     "${_ACML_ROOT}/${_ACML_COMPILER32}_mp${_ACML_PATH_SUFFIX}/lib"
     "${_ACML_ROOT}/${_ACML_COMPILER64}_mp${_ACML_PATH_SUFFIX}/lib" )
   else() #if( _BLAS_VENDOR STREQUAL "ACML" )
    set(_ACML_LIB_DIRS
     "${_ACML_ROOT}/${_ACML_COMPILER32}${_ACML_PATH_SUFFIX}/lib"
     "${_ACML_ROOT}/${_ACML_COMPILER64}${_ACML_PATH_SUFFIX}/lib" )
   endif()
  endif()
 endif()

 if( _BLAS_VENDOR STREQUAL "ACML_MP" )
  foreach( BLAS_ACML_MP_LIB_DIRS ${_ACML_MP_LIB_DIRS} )
   _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "acml_mp;acml_mv" "" )
   if( BLAS_${_BLAS_VENDOR}_FOUND )
    break()
   endif()
  endforeach()
 else() #if( _BLAS_VENDOR STREQUAL "ACML" )
  foreach( BLAS_ACML_LIB_DIRS ${_ACML_LIB_DIRS} )
   _BLAS_LOCATE_AND_TEST( ${_BLAS_VENDOR} "acml;acml_mv" "" )
   if( BLAS_${_BLAS_VENDOR}_FOUND )
    break()
   endif()
  endforeach()
 endif()

 # Either acml or acml_mp should be in LD_LIBRARY_PATH but not both
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  ""
  "acml;acml_mv"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  ""
  "acml_mp;acml_mv"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
endif () # ACML

# Apple BLAS library?
if (BLA_VENDOR STREQUAL "Apple" OR BLA_VENDOR STREQUAL "All")
if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  cblas_dgemm
  ""
  "Accelerate"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "Apple" OR BLA_VENDOR STREQUAL "All")

if (BLA_VENDOR STREQUAL "NAS" OR BLA_VENDOR STREQUAL "All")
 if ( NOT BLAS_LIBRARIES )
    check_fortran_libraries(
    BLAS_LIBRARIES
    BLAS
    cblas_dgemm
    ""
    "vecLib"
    ""
    )
 endif ( NOT BLAS_LIBRARIES )
endif (BLA_VENDOR STREQUAL "NAS" OR BLA_VENDOR STREQUAL "All")
# Generic BLAS library?
if (BLA_VENDOR STREQUAL "Generic" OR BLA_VENDOR STREQUAL "All")
 if(NOT BLAS_LIBRARIES)
  check_fortran_libraries(
  BLAS_LIBRARIES
  BLAS
  sgemm
  ""
  "blas"
  ""
  )
 endif(NOT BLAS_LIBRARIES)
endif (BLA_VENDOR STREQUAL "Generic" OR BLA_VENDOR STREQUAL "All")

#BLAS in intel mkl 10 library? (em64t 64bit)
if (BLA_VENDOR MATCHES "Intel*" OR BLA_VENDOR STREQUAL "All")
 if (NOT WIN32)
  set(LM "-lm")
 endif ()
 if (_LANGUAGES_ MATCHES C OR _LANGUAGES_ MATCHES CXX)
  if(BLAS_FIND_QUIETLY OR NOT BLAS_FIND_REQUIRED)
    find_package(Threads)
  else(BLAS_FIND_QUIETLY OR NOT BLAS_FIND_REQUIRED)
    find_package(Threads REQUIRED)
  endif(BLAS_FIND_QUIETLY OR NOT BLAS_FIND_REQUIRED)
  if (WIN32)
  if(BLA_F95)
    if(NOT BLAS95_LIBRARIES)
    check_fortran_libraries(
    BLAS95_LIBRARIES
    BLAS
    sgemm
    ""
    "mkl_blas95;mkl_intel_c;mkl_intel_thread;mkl_core;libguide40"
    ""
    )
    endif(NOT BLAS95_LIBRARIES)
  else(BLA_F95)
    if(NOT BLAS_LIBRARIES)
    check_fortran_libraries(
    BLAS_LIBRARIES
    BLAS
    SGEMM
    ""
    "mkl_c_dll;mkl_intel_thread_dll;mkl_core_dll;libguide40"
    ""
    )
    endif(NOT BLAS_LIBRARIES)
  endif(BLA_F95)
  else(WIN32)
  if (BLA_VENDOR STREQUAL "Intel10_32" OR BLA_VENDOR STREQUAL "All")
    if(BLA_F95)
      if(NOT BLAS95_LIBRARIES)
      check_fortran_libraries(
      BLAS95_LIBRARIES
      BLAS
      sgemm
      ""
      "mkl_blas95;mkl_intel;mkl_intel_thread;mkl_core;guide"
      "${CMAKE_THREAD_LIBS_INIT};${LM}"
      )
      endif(NOT BLAS95_LIBRARIES)
    else(BLA_F95)
    if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
      BLAS_LIBRARIES
      BLAS
      sgemm
      ""
      "mkl_intel;mkl_intel_thread;mkl_core;guide"
      "${CMAKE_THREAD_LIBS_INIT}"
      "${LM}"
      )
      endif(NOT BLAS_LIBRARIES)
    endif(BLA_F95)
  endif (BLA_VENDOR STREQUAL "Intel10_32" OR BLA_VENDOR STREQUAL "All")
  if (BLA_VENDOR STREQUAL "Intel10_64lp" OR BLA_VENDOR STREQUAL "All")
   if(BLA_F95)
    if(NOT BLAS95_LIBRARIES)
      check_fortran_libraries(
      BLAS95_LIBRARIES
      BLAS
      sgemm
      ""
      "mkl_blas95;mkl_intel_lp64;mkl_intel_thread;mkl_core;guide"
      "${CMAKE_THREAD_LIBS_INIT};${LM}"
      )
    endif(NOT BLAS95_LIBRARIES)
   else(BLA_F95)
     if(NOT BLAS_LIBRARIES)
      check_fortran_libraries(
      BLAS_LIBRARIES
      BLAS
      sgemm
      ""
      "mkl_intel_lp64;mkl_intel_thread;mkl_core;guide"
      "${CMAKE_THREAD_LIBS_INIT};${LM}"
      )
     endif(NOT BLAS_LIBRARIES)
   endif(BLA_F95)
  endif (BLA_VENDOR STREQUAL "Intel10_64lp" OR BLA_VENDOR STREQUAL "All")
  endif (WIN32)
  #older vesions of intel mkl libs
  # BLAS in intel mkl library? (shared)
  if(NOT BLAS_LIBRARIES)
    check_fortran_libraries(
    BLAS_LIBRARIES
    BLAS
    sgemm
    ""
    "mkl;guide"
    "${CMAKE_THREAD_LIBS_INIT};${LM}"
    )
  endif(NOT BLAS_LIBRARIES)
  #BLAS in intel mkl library? (static, 32bit)
  if(NOT BLAS_LIBRARIES)
    check_fortran_libraries(
    BLAS_LIBRARIES
    BLAS
    sgemm
    ""
    "mkl_ia32;guide"
    "${CMAKE_THREAD_LIBS_INIT};${LM}"
    )
  endif(NOT BLAS_LIBRARIES)
  #BLAS in intel mkl library? (static, em64t 64bit)
  if(NOT BLAS_LIBRARIES)
    check_fortran_libraries(
    BLAS_LIBRARIES
    BLAS
    sgemm
    ""
    "mkl_em64t;guide"
    "${CMAKE_THREAD_LIBS_INIT};${LM}"
    )
  endif(NOT BLAS_LIBRARIES)
 endif (_LANGUAGES_ MATCHES C OR _LANGUAGES_ MATCHES CXX)
endif (BLA_VENDOR MATCHES "Intel*" OR BLA_VENDOR STREQUAL "All")

if (BLAS_LIBRARIES_VARS)
  mark_as_advanced (FORCE ${BLAS_LIBRARIES_VARS})
endif ()
if (BLAS95_LIBRARIES)
  mark_as_advanced (FORCE ${BLAS95_LIBRARIES})
endif ()
if (BLAS_INCLUDE_DIR)
  mark_as_advanced (FORCE BLAS_INCLUDE_DIR)
endif ()

if(BLA_F95)
 if(BLAS95_LIBRARIES)
    set(BLAS95_FOUND TRUE)
  else(BLAS95_LIBRARIES)
    set(BLAS95_FOUND FALSE)
  endif(BLAS95_LIBRARIES)

  if(NOT BLAS_FIND_QUIETLY)
    if(BLAS95_FOUND)
      mark_as_advanced (${BLAS95_LIBRARIES})
      message(STATUS "A library with BLAS95 API found.")
    else(BLAS95_FOUND)
      if(BLAS_FIND_REQUIRED)
        mark_as_advanced (CLEAR ${BLAS95_LIBRARIES})
        message(FATAL_ERROR
        "A required library with BLAS95 API not found. Please specify library location.")
      else(BLAS_FIND_REQUIRED)
        message(STATUS
        "A library with BLAS95 API not found. Please specify library location.")
      endif(BLAS_FIND_REQUIRED)
    endif(BLAS95_FOUND)
  endif(NOT BLAS_FIND_QUIETLY)
  set(BLAS_FOUND TRUE)
  set(BLAS_LIBRARIES "${BLAS95_LIBRARIES}")
else(BLA_F95)
  if(BLAS_LIBRARIES)
    if(BLAS_LIBRARIES_VARS MATCHES "cblas" AND NOT BLAS_INCLUDE_DIR)
      set(BLAS_FOUND FALSE)
    else()
      set(BLAS_FOUND TRUE)
    endif()
  else(BLAS_LIBRARIES)
    set(BLAS_FOUND FALSE)
  endif(BLAS_LIBRARIES)

  if(NOT BLAS_FIND_QUIETLY)
    if(BLAS_FOUND)
      message(STATUS "A library with BLAS API found.")
    else(BLAS_FOUND)
      if(BLAS_LIBRARIES AND BLAS_cblas_LIBRARY AND NOT BLAS_INCLUDE_DIR)
        if(BLAS_FIND_REQUIRED)
          mark_as_advanced (CLEAR BLAS_INCLUDE_DIR)
          message(FATAL_ERROR
            "Location of cblas.h header file for C BLAS library not found!"
            " Please specify the directory containing this file for library ${BLAS_cblas_LIBRARY}"
            " using the BLAS_INCLUDE_DIR variable."
          )
        else(BLAS_FIND_REQUIRED)
          message(STATUS
          "Location of cblas.h header file for C BLAS library not found! Please specify header location."
          )
        endif(BLAS_FIND_REQUIRED)
      else()
        if(BLAS_FIND_REQUIRED)
          mark_as_advanced (CLEAR ${BLAS_LIBRARIES_VARS})
          message(FATAL_ERROR
            "A required library with BLAS API not found. Please specify library location"
            " by setting the following variables: [${BLAS_LIBRARIES_VARS}]."
          )
        else(BLAS_FIND_REQUIRED)
          message(STATUS
          "A library with BLAS API not found. Please specify library location."
          )
        endif(BLAS_FIND_REQUIRED)
      endif()
    endif(BLAS_FOUND)
  endif(NOT BLAS_FIND_QUIETLY)
endif(BLA_F95)

unset(BLAS_LIBRARIES_VARS)
