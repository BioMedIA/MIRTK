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

#ifndef MIRTK_CudaRuntime_H
#define MIRTK_CudaRuntime_H


#ifdef HAVE_CUDA
#  include <cuda_runtime.h>
#endif


namespace mirtk {


// =============================================================================
// vector_types.h substitutes in mirtk namespace
// =============================================================================

typedef unsigned char  uchar;
typedef unsigned int   uint;
typedef unsigned short ushort;

#if defined(__VECTOR_TYPES_H__)

#define __vector1_type__(tag, type) using tag
#define __vector2_type__(tag, type) using tag
#define __vector3_type__(tag, type) using tag
#define __vector4_type__(tag, type) using tag

#else // __VECTOR_TYPES_H__

#define __vector1_type__(tag, type) \
  struct tag { type x; }; \
  typedef struct tag tag
#define __vector2_type__(tag, type) \
  struct tag { type x; type y; }; \
  typedef struct tag tag
#define __vector3_type__(tag, type) \
  struct tag { type x; type y; type z; }; \
  typedef struct tag tag
#define __vector4_type__(tag, type) \
  struct tag { type x; type y; type z; type w; }; \
  typedef struct tag tag

#endif // __VECTOR_TYPES_H__


#define __vector_type__(type) \
  __vector1_type__(type##1, type); \
  __vector2_type__(type##2, type); \
  __vector3_type__(type##3, type); \
  __vector4_type__(type##4, type)

__vector_type__(char);
__vector_type__(uchar);
__vector_type__(short);
__vector_type__(ushort);
__vector_type__(int);
__vector_type__(uint);
__vector_type__(float);
__vector_type__(double);

#undef __vector_type__
#undef __vector1_type__
#undef __vector2_type__
#undef __vector3_type__
#undef __vector4_type__


#if defined(__VECTOR_TYPES_H__)
using dim3;
#else
struct dim3
{
  uint x, y, z;
  dim3(uint vx = 1, uint vy = 1, uint vz = 1) : x(vx), y(vy), z(vz) {}
  dim3(uint3 v) : x(v.x), y(v.y), z(v.z) {}
  operator uint3() { uint3 t; t.x = x; t.y = y; t.z = z; return t; }
};
typedef struct dim3 dim3;
#endif // !defined(__VECTOR_TYPES_H__)

// =============================================================================
// vector_functions.h substitutes in mirtk namespace
// =============================================================================

#define __vector1_func__(type) \
  inline type##1 make_##type##1(type x) \
  { \
    type##1 t; t.x = x; return t; \
  }

#define __vector2_func__(type) \
  inline type##2 make_##type##2(type x, type y) \
  { \
    type##2 t; t.x = x, t.y = y; return t; \
  }

#define __vector3_func__(type) \
  inline type##3 make_##type##3(type x, type y, type z) \
  { \
    type##3 t; t.x = x, t.y = y, t.z = z; return t; \
  }

#define __vector4_func__(type) \
  inline type##4 make_##type##4(type x, type y, type z, type w) \
  { \
    type##4 t; t.x = x, t.y = y, t.z = z, t.w = w; return t; \
  }

#define __vector_func__(type) \
  __vector1_func__(type); \
  __vector2_func__(type); \
  __vector3_func__(type); \
  __vector4_func__(type)

__vector_func__(char);
__vector_func__(uchar);
__vector_func__(short);
__vector_func__(ushort);
__vector_func__(int);
__vector_func__(uint);
__vector_func__(float);
__vector_func__(double);

#undef __vector_func__
#undef __vector1_func__
#undef __vector2_func__
#undef __vector3_func__
#undef __vector4_func__


} // namespace mirtk

#endif // MIRTK_CudaRuntime_H
