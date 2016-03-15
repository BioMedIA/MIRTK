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
  mirtkAlgorithm.h
  mirtkAllocate.h
  mirtkArray.h
  mirtkArrayHeap.h
  mirtkAssert.h
  mirtkCfstream.h
  mirtkCifstream.h
  mirtkCofstream.h
  mirtkCommon.h
  mirtkConfig.h
  mirtkConfigurable.h
  mirtkCuda.h
  mirtkCudaRuntime.h
  mirtkCutilMath.h
  mirtkDeallocate.h
  mirtkEnergyMeasure.h
  mirtkEvent.h
  mirtkEventDelegate.h
  mirtkFastDelegate.h
  mirtkIndent.h
  mirtkList.h
  mirtkMath.h
  mirtkMemory.h
  mirtkObject.h
  mirtkObjectFactory.h
  mirtkObservable.h
  mirtkObserver.h
  mirtkOptions.h
  mirtkOrderedMap.h
  mirtkOrderedSet.h
  mirtkPair.h
  mirtkParallel.h
  mirtkPath.h
  mirtkProfiling.h
  mirtkPriorityQueue.h
  mirtkQueue.h
  mirtkRandom.h
  mirtkStack.h
  mirtkStatus.h
  mirtkStream.h
  mirtkString.h
  mirtkTerminal.h
  mirtkTestProd.h
  mirtkUnorderedMap.h
  mirtkUnorderedSet.h
  mirtkUtils.h
  mirtkVersion.h
  ${BINARY_INCLUDE_DIR}/mirtkVersionInfo.h
)

set(SOURCES
  mirtkCifstream.cc
  mirtkCofstream.cc
  mirtkConfigurable.cc
  mirtkMath.cc
  mirtkMemory.cc
  mirtkObserver.cc
  mirtkOptions.cc
  mirtkParallel.cc
  mirtkPath.cc
  mirtkProfiling.cc
  mirtkString.cc
  mirtkTerminal.cc
  mirtkVersion.cc
)

if (VTK_FOUND)
  list(APPEND HEADERS mirtkVtk.h mirtkVtkMath.h)
  list(APPEND SOURCES mirtkVtk.cc)
  list(APPEND DEPENDS ${VTK_LIBRARIES})
endif ()

if (MATLAB_FOUND)
  list(APPEND HEADERS mirtkMatlab.h)
  list(APPEND SOURCES mirtkMatlab.cc)
  list(APPEND DEPENDS ${MATLAB_mwmclmcrrt_LIBRARY})
endif ()

if (TBB_FOUND)
  list(APPEND DEPENDS ${TBB_LIBRARIES})
endif ()

if (ZLIB_FOUND)
  list(APPEND DEPENDS ${ZLIB_LIBRARIES})
endif ()

mirtk_add_library(HEADERS ${HEADERS} SOURCES ${SOURCES} DEPENDS ${DEPENDS})
