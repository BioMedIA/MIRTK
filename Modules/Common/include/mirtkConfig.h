/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Imperial College London
 * Copyright 2013-2015 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef MIRTK_Config_H
#define MIRTK_Config_H


// ===========================================================================
// General
// ===========================================================================

/// Whether to build for execution on Microsoft Windows
#ifndef WINDOWS
#  if defined(_WIN32) || defined(_WIN64) || defined(_WINDOWS) || defined(WIN32) 
#    define WINDOWS
#  endif
#endif

/// Precision of floating point types to use by default
/// 0: single-precision 1: double-precision
#define MIRTK_USE_FLOAT_BY_DEFAULT 0

// ===========================================================================
// CUDA
// ===========================================================================

// ---------------------------------------------------------------------------
#ifndef MIRTKCU_API
#  if __CUDACC__
#    define MIRTKCU_API __device__ __host__
#  else
#    define MIRTKCU_API
#  endif
#endif

// ---------------------------------------------------------------------------
#ifndef MIRTKCU_HOST_API
#  if __CUDACC__
#    define MIRTKCU_HOST_API __host__
#  else
#    define MIRTKCU_HOST_API
#  endif
#endif

// ---------------------------------------------------------------------------
#ifndef MIRTKCU_DEVICE_API
#  if __CUDACC__
#    define MIRTKCU_DEVICE_API __device__
#  else
#    define MIRTKCU_DEVICE_API
#  endif
#endif


#endif // MIRTK_Config_H
