# ============================================================================
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2013-2015 Imperial College London
# Copyright 2013-2015 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ============================================================================

##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of MIRTK Common library.
##############################################################################

set(HEADERS
  ${BINARY_INCLUDE_DIR}/mirtk/CommonConfig.h
  ${BINARY_INCLUDE_DIR}/mirtk/CommonExport.h
  ${BINARY_INCLUDE_DIR}/mirtk/VersionInfo.h
  Algorithm.h
  Allocate.h
  Array.h
  ArrayHeap.h
  Assert.h
  Cfstream.h
  Cifstream.h
  Cofstream.h
  Common.h
  Config.h
  Configurable.h
  Cuda.h
  CudaRuntime.h
  CutilMath.h
  Deallocate.h
  EnergyMeasure.h
  Event.h
  EventDelegate.h
  Exception.h
  FastDelegate.h
  Indent.h
  List.h
  Math.h
  Memory.h
  Numeric.h
  Object.h
  ObjectFactory.h
  Observable.h
  Observer.h
  Options.h
  OrderedMap.h
  OrderedSet.h
  Pair.h
  Parallel.h
  Path.h
  Profiling.h
  PriorityQueue.h
  Queue.h
  Random.h
  Stack.h
  Status.h
  Stream.h
  String.h
  System.h
  Terminal.h
  TestProd.h
  UnorderedMap.h
  UnorderedSet.h
  Utils.h
  Version.h
)

set(SOURCES
  Cifstream.cc
  Cofstream.cc
  Configurable.cc
  Math.cc
  Memory.cc
  Observer.cc
  Options.cc
  Parallel.cc
  Path.cc
  Profiling.cc
  String.cc
  System.cc
  Terminal.cc
  Version.cc
)

if (VTK_FOUND)
  list(APPEND HEADERS Vtk.h VtkMath.h)
  list(APPEND SOURCES Vtk.cc)
  list(APPEND DEPENDS ${VTK_LIBRARIES})
endif ()

if (MATLAB_FOUND)
  list(APPEND HEADERS Matlab.h)
  list(APPEND SOURCES Matlab.cc)
  list(APPEND DEPENDS ${MATLAB_mwmclmcrrt_LIBRARY})
endif ()

if (TARGET TBB::tbb)
  list(APPEND DEPENDS TBB::tbb)
endif ()

if (ZLIB_FOUND)
  list(APPEND DEPENDS ${ZLIB_LIBRARIES})
endif ()

mirtk_add_library(HEADERS ${HEADERS} SOURCES ${SOURCES} DEPENDS ${DEPENDS})
