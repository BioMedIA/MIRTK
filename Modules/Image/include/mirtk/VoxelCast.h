/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#ifndef MIRTK_VoxelCast_H
#define MIRTK_VoxelCast_H

#include "mirtk/Voxel.h"
#include "mirtk/Stream.h"
#include "mirtk/Math.h"


// Overloading the C++ conversion operators for vector types would be possible
// as well, however, this can lead to ambiguitis in mathematical expressions
// when a vector can be constructed from a scalar and converted to a scalar at
// the same time. Then it is no longer clear whether the arithmetic operation
// should be performed between scalars, vectors, or a mix of both. Therefore,
// only support the needed voxel type conversions as template specializations
// of our own voxel_cast template function.

namespace mirtk {


// -----------------------------------------------------------------------------
/// Auxiliary template class for partial voxel_cast specialization
template <class TIn, class TOut>
struct VoxelCaster
{
  /// By default, use static_cast to convert voxel value from one type to another.
  /// If the value is out of the range which can be represented by the output
  /// type, the minium/maximum value of the output type is returned instead.
  static TOut Convert(const TIn &value)
  {
    if (static_cast<double>(value) < voxel_limits<TOut>::min()) {
      return voxel_limits<TOut>::min_value();
    } else if (static_cast<double>(value) > voxel_limits<TOut>::max()) {
      return voxel_limits<TOut>::max_value();
    }
    return static_cast <TOut>(value);
  }
};

// -----------------------------------------------------------------------------
template <class T>
struct VoxelCaster<T, T>
{
  static T Convert(const T &value)
  {
    return value;
  }
};


// =============================================================================
// Scalar types
// =============================================================================

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float, float>
{
  static float Convert(const float &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float, double>
{
  static double Convert(const float &value)
  {
    return static_cast<double>(value);
  }
};

// -----------------------------------------------------------------------------
#define _MIRTK_VOXELCAST_FLOAT_TO_INT(TOut) \
    template <> \
    struct VoxelCaster<float, TOut> \
    { \
      static TOut Convert(const float &value) \
      { \
        if (static_cast<double>(value) < voxel_limits<TOut>::min()) { \
          return voxel_limits<TOut>::min_value(); \
        } else if (static_cast<double>(value) > voxel_limits<TOut>::max()) { \
          return voxel_limits<TOut>::max_value(); \
        } \
        return static_cast<TOut>(round(value)); \
      } \
    }
_MIRTK_VOXELCAST_FLOAT_TO_INT(unsigned char);
_MIRTK_VOXELCAST_FLOAT_TO_INT(char);
_MIRTK_VOXELCAST_FLOAT_TO_INT(unsigned short);
_MIRTK_VOXELCAST_FLOAT_TO_INT(short);
_MIRTK_VOXELCAST_FLOAT_TO_INT(int);
_MIRTK_VOXELCAST_FLOAT_TO_INT(unsigned int);
_MIRTK_VOXELCAST_FLOAT_TO_INT(long);
_MIRTK_VOXELCAST_FLOAT_TO_INT(unsigned long);
#undef _MIRTK_VOXELCAST_FLOAT_TO_INT

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double, float>
{
  static float Convert(const double &value)
  {
    return static_cast<float>(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double, double>
{
  static double Convert(const double &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
#define _MIRTK_VOXELCAST_DOUBLE_TO_INT(TOut) \
    template <> \
    struct VoxelCaster<double, TOut> \
    { \
      static TOut Convert(const double &value) \
      { \
        if (value < voxel_limits<TOut>::min()) { \
          return voxel_limits<TOut>::min_value(); \
        } else if (value > voxel_limits<TOut>::max()) { \
          return voxel_limits<TOut>::max_value(); \
        } \
        return static_cast<TOut>(round(value)); \
      } \
    }
_MIRTK_VOXELCAST_DOUBLE_TO_INT(unsigned char);
_MIRTK_VOXELCAST_DOUBLE_TO_INT(char);
_MIRTK_VOXELCAST_DOUBLE_TO_INT(unsigned short);
_MIRTK_VOXELCAST_DOUBLE_TO_INT(short);
_MIRTK_VOXELCAST_DOUBLE_TO_INT(int);
_MIRTK_VOXELCAST_DOUBLE_TO_INT(unsigned int);
_MIRTK_VOXELCAST_DOUBLE_TO_INT(long);
_MIRTK_VOXELCAST_DOUBLE_TO_INT(unsigned long);
#undef _MIRTK_VOXELCAST_DOUBLE_TO_INT

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float1>
{
  static float1 Convert(const TIn &value)
  {
    return make_float1(value);
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double1>
{
  static double1 Convert(const TIn &value)
  {
    return make_double1(static_cast<double>(value));
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float2>
{
  static float2 Convert(const TIn &value)
  {
    return make_float2(static_cast<float>(value));
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double2>
{
  static double2 Convert(const TIn &value)
  {
    return make_double2(value);
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float3>
{
  static float3 Convert(const TIn &value)
  {
    return make_float3(static_cast<float>(value));
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double3>
{
  static double3 Convert(const TIn &value)
  {
    return make_double3(static_cast<double>(value));
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, Vector>
{
  static Vector Convert(const TIn &value)
  {
    return Vector(1, VoxelCaster<TIn, double>::Convert(value));
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<int, Vector3D<TOut> >
{
  static Vector3D<TOut> Convert(const int &value)
  {
    return Vector3D<TOut>(VoxelCaster<int, TOut>::Convert(value));
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float4>
{
  static float4 Convert(const TIn &value)
  {
    return make_float4(static_cast<float>(value));
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double4>
{
  static double4 Convert(const TIn &value)
  {
    return make_double4(static_cast<double>(value));
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, float3x3>
{
  static float3x3 Convert(const TIn &value)
  {
    return make_float3x3(static_cast<float>(value));
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<TIn, double3x3>
{
  static double3x3 Convert(const TIn &value)
  {
    return make_double3x3(static_cast<double>(value));
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<int, VectorND<9, TOut> >
{
  static VectorND<9, TOut> Convert(const int &value)
  {
    return VectorND<9, TOut>(VoxelCaster<int, TOut>::Convert(value));
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double, VectorND<9, TOut> >
{
  static VectorND<9, TOut> Convert(const double &value)
  {
    return VectorND<9, TOut>(VoxelCaster<double, TOut>::Convert(value));
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float, Vector>
{
  static Vector Convert(const float &value)
  {
    return Vector(1, static_cast<double>(value));
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double, Vector>
{
  static Vector Convert(const double &value)
  {
    return Vector(1, value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double, Vector3D<TOut> >
{
  static Vector3D<TOut> Convert(const double &value)
  {
    return Vector3D<TOut>(VoxelCaster<double, TOut>::Convert(value));
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<int, Vector4D<TOut> >
{
  static Vector4D<TOut> Convert(const int &value)
  {
    return Vector4D<TOut>(VoxelCaster<int, TOut>::Convert(value));
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double, Vector4D<TOut> >
{
  static Vector4D<TOut> Convert(const double &value)
  {
    return Vector4D<TOut>(VoxelCaster<double, TOut>::Convert(value));
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float1, TOut>
{
  static TOut Convert(const float1 &value)
  {
    return VoxelCaster<float, TOut>::Convert(value.x);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float1, float1>
{
  static float1 Convert(const float1 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float1, double1>
{
  static double1 Convert(const float1 &value)
  {
    return make_double1(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float1, Vector>
{
  static Vector Convert(const float1 &value)
  {
    Vector v(1);
    v(0) = value.x;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double1, TOut>
{
  static TOut Convert(const double1 &value)
  {
    return VoxelCaster<double, TOut>::Convert(value.x);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double1, float1>
{
  static float1 Convert(const double1 &value)
  {
    return make_float1(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double1, double1>
{
  static double1 Convert(const double1 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double1, Vector>
{
  static Vector Convert(const double1 &value)
  {
    Vector v(1);
    v(0) = value.x;
    return v;
  }
};


// =============================================================================
// Vector for generic BaseImage interface
// =============================================================================

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<Vector, TOut>
{
  static TOut Convert(const Vector &value)
  {
    if (value.Rows() == 1) return VoxelCaster<double, TOut>::Convert(value(0));
    cerr << "Can only cast vector with exactly one element to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, Vector>
{
  static Vector Convert(const Vector &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, float1>
{
  static float1 Convert(const Vector &value)
  {
    if (value.Rows() == 1) return make_float1(value(0));
    cerr << "Can only cast vector with exactly one element to a 1D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, double1>
{
  static double1 Convert(const Vector &value)
  {
    if (value.Rows() == 1) return make_double1(value(0));
    cerr << "Can only cast vector with exactly one element to a 1D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, float2>
{
  static float2 Convert(const Vector &value)
  {
    if (value.Rows() == 2) return make_float2(value(0), value(1));
    cerr << "Can only cast vector with exactly two elements to a 2D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, double2>
{
  static double2 Convert(const Vector &value)
  {
    if (value.Rows() == 2) return make_double2(value(0), value(1));
    cerr << "Can only cast vector with exactly two elements to a 2D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, float3>
{
  static float3 Convert(const Vector &value)
  {
    if (value.Rows() == 3) return make_float3(value(0), value(1), value(2));
    cerr << "Can only cast vector with exactly three elements to a 3D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, double3>
{
  static double3 Convert(const Vector &value)
  {
    if (value.Rows() == 3) return make_double3(value(0), value(1), value(2));
    cerr << "Can only cast vector with exactly three elements to a 3D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<Vector, Vector3D<TOut> >
{
  static Vector3D<TOut> Convert(const Vector &value)
  {
    if (value.Rows() == 3) {
      return Vector3D<TOut>(VoxelCaster<double, TOut>::Convert(value(0)),
                            VoxelCaster<double, TOut>::Convert(value(1)),
                            VoxelCaster<double, TOut>::Convert(value(2)));
    }
    cerr << "Can only cast vector with exactly three elements to a 3D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, float4>
{
  static float4 Convert(const Vector &value)
  {
    if (value.Rows() == 4) return make_float4(value(0), value(1), value(2), value(3));
    cerr << "Can only cast vector with exactly four elements to a 4D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, double4>
{
  static double4 Convert(const Vector &value)
  {
    if (value.Rows() == 4) return make_double4(value(0), value(1), value(2), value(3));
    cerr << "Can only cast vector with exactly four elements to a 4D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<Vector, Vector4D<TOut> >
{
  static Vector4D<TOut> Convert(const Vector &value)
  {
    if (value.Rows() == 4) {
      return Vector4D<TOut>(VoxelCaster<double, TOut>::Convert(value(0)),
                            VoxelCaster<double, TOut>::Convert(value(1)),
                            VoxelCaster<double, TOut>::Convert(value(2)),
                            VoxelCaster<double, TOut>::Convert(value(3)));
    }
    cerr << "Can only cast vector with exactly four elements to a 4D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<Vector, VectorND<9, TOut> >
{
  static VectorND<9, TOut> Convert(const Vector &value)
  {
    if (value.Rows() == 9) {
      VectorND<9, TOut> v;
      for (int i = 0; i < 9; ++i) {
        v(i) = VoxelCaster<double, TOut>::Convert(value(i));
      }
      return v;
    }
    cerr << "Can only cast vector with exactly nine elements to a 9D vector!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, float3x3>
{
  static float3x3 Convert(const Vector &v)
  {
    float3x3 m;
    if (v.Rows() == 9) {
      m.a.x = VoxelCaster<double, float>::Convert(v(0));
      m.a.y = VoxelCaster<double, float>::Convert(v(1));
      m.a.z = VoxelCaster<double, float>::Convert(v(2));
      m.b.x = VoxelCaster<double, float>::Convert(v(3));
      m.b.y = VoxelCaster<double, float>::Convert(v(4));
      m.b.z = VoxelCaster<double, float>::Convert(v(5));
      m.c.x = VoxelCaster<double, float>::Convert(v(6));
      m.c.y = VoxelCaster<double, float>::Convert(v(7));
      m.c.z = VoxelCaster<double, float>::Convert(v(8));
    } else if (v.Rows() == 6) {
      m.a.x = VoxelCaster<double, float>::Convert(v(0));
      m.a.y = VoxelCaster<double, float>::Convert(v(1));
      m.a.z = VoxelCaster<double, float>::Convert(v(2));
      m.b.x = m.a.y;
      m.b.y = VoxelCaster<double, float>::Convert(v(3));
      m.b.z = VoxelCaster<double, float>::Convert(v(4));
      m.c.x = m.a.z;
      m.c.y = m.b.z;
      m.c.z = VoxelCaster<double, float>::Convert(v(5));
    } else {
      cerr << "Can only cast vector of size 6 or 9 to a 3x3 matrix!" << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
    }
    return m;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<Vector, double3x3>
{
  static double3x3 Convert(const Vector &v)
  {
    double3x3 m;
    if (v.Rows() == 9) {
      m.a.x = v(0);
      m.a.y = v(1);
      m.a.z = v(2);
      m.b.x = v(3);
      m.b.y = v(4);
      m.b.z = v(5);
      m.c.x = v(6);
      m.c.y = v(7);
      m.c.z = v(8);
    } else if (v.Rows() == 6) {
      m.a.x = v(0);
      m.a.y = v(1);
      m.a.z = v(2);
      m.b.x = m.a.y;
      m.b.y = v(3);
      m.b.z = v(4);
      m.c.x = m.a.z;
      m.c.y = m.b.z;
      m.c.z = v(5);
    } else {
      cerr << "Can only cast vector of size 6 or 9 to a 3x3 matrix!" << endl;
      cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
      exit(1);
    }
    return m;
  }
};


// =============================================================================
// 2D vector types
// =============================================================================

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float2, TOut>
{
  static TOut Convert(const float2 &)
  {
    cerr << "Cannot cast 2D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float2, float2>
{
  static float2 Convert(const float2 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float2, double2>
{
  static double2 Convert(const float2 &value)
  {
    return make_double2(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float2, Vector>
{
  static Vector Convert(const float2 &value)
  {
    Vector v(2);
    v(0) = value.x;
    v(1) = value.y;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double2, TOut>
{
  static TOut Convert(const double2 &)
  {
    cerr << "Cannot cast 2D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double2, float2>
{
  static float2 Convert(const double2 &value)
  {
    return make_float2(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double2, double2>
{
  static double2 Convert(const double2 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double2, Vector>
{
  static Vector Convert(const double2 &value)
  {
    Vector v(2);
    v(0) = value.x;
    v(1) = value.y;
    return v;
  }
};


// =============================================================================
// 3D vector types
// =============================================================================

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float3, TOut>
{
  static TOut Convert(const float3 &)
  {
    cerr << "Cannot cast 3D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3, float3>
{
  static float3 Convert(const float3 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3, double3>
{
  static double3 Convert(const float3 &value)
  {
    return make_double3(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3, Vector>
{
  static Vector Convert(const float3 &value)
  {
    Vector v(3);
    v(0) = value.x;
    v(1) = value.y;
    v(2) = value.z;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double3, TOut>
{
  static TOut Convert(const double3 &)
  {
    cerr << "Cannot cast 3D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3, float3>
{
  static float3 Convert(const double3 &value)
  {
    return make_float3(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3, double3>
{
  static double3 Convert(const double3 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3, Vector>
{
  static Vector Convert(const double3 &value)
  {
    Vector v(3);
    v(0) = value.x;
    v(1) = value.y;
    v(2) = value.z;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<Vector3D<TIn>, Vector>
{
  static Vector Convert(const Vector3D<TIn> &value)
  {
    Vector v(3);
    v.Put(value);
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TIn, class TOut>
struct VoxelCaster<Vector3D<TIn>, TOut>
{
  static TOut Convert(const Vector3D<TIn> &)
  {
    cerr << "Cannot cast 3D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TIn, class TOut>
struct VoxelCaster<Vector3D<TIn>, Vector3D<TOut> >
{
  static Vector3D<TOut> Convert(const Vector3D<TIn> &value)
  {
    return Vector3D<TOut>(VoxelCaster<TIn, TOut>::Convert(value._x),
                          VoxelCaster<TIn, TOut>::Convert(value._y),
                          VoxelCaster<TIn, TOut>::Convert(value._z));
  }
};

// -----------------------------------------------------------------------------
template <class T>
struct VoxelCaster<Vector3D<T>, Vector3D<T> >
{
  static Vector3D<T> Convert(const Vector3D<T> &value)
  {
    return value;
  }
};


// =============================================================================
// 4D vector types
// =============================================================================

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float4, TOut>
{
  static TOut Convert(const float4 &)
  {
    cerr << "Cannot cast 4D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float4, float4>
{
  static float4 Convert(const float4 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float4, double4>
{
  static double4 Convert(const float4 &value)
  {
    return make_double4(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float4, Vector>
{
  static Vector Convert(const float4 &value)
  {
    Vector v(4);
    v(0) = value.x;
    v(1) = value.y;
    v(2) = value.z;
    v(3) = value.w;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double4, TOut>
{
  static TOut Convert(const double4 &)
  {
    cerr << "Cannot cast 4D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double4, float4>
{
  static float4 Convert(const double4 &value)
  {
    return make_float4(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double4, double4>
{
  static double4 Convert(const double4 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double4, Vector>
{
  static Vector Convert(const double4 &value)
  {
    Vector v(4);
    v(0) = value.x;
    v(1) = value.y;
    v(2) = value.z;
    v(3) = value.w;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TIn, class TOut>
struct VoxelCaster<Vector4D<TIn>, TOut>
{
  static TOut Convert(const Vector4D<TIn> &)
  {
    cerr << "Cannot cast 4D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TIn, class TOut>
struct VoxelCaster<Vector4D<TIn>, Vector4D<TOut> >
{
  static Vector4D<TOut> Convert(const Vector4D<TIn> &value)
  {
    return Vector4D<TOut>(VoxelCaster<TIn, TOut>::Convert(value._x),
                          VoxelCaster<TIn, TOut>::Convert(value._y),
                          VoxelCaster<TIn, TOut>::Convert(value._z),
                          VoxelCaster<TIn, TOut>::Convert(value._t));
  }
};

// -----------------------------------------------------------------------------
template <class T>
struct VoxelCaster<Vector4D<T>, Vector4D<T> >
{
  static Vector4D<T> Convert(const Vector4D<T> &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<Vector4D<TIn>, Vector>
{
  static Vector Convert(const Vector4D<TIn> &value)
  {
    Vector v(4);
    v.Put(value);
    return v;
  }
};


// =============================================================================
// 9D vector types
// =============================================================================

// -----------------------------------------------------------------------------
template <class TIn, class TOut>
struct VoxelCaster<VectorND<9, TIn>, TOut>
{
  static TOut Convert(const VectorND<9, TIn> &)
  {
    cerr << "Cannot cast 9D vector to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <class TIn, class TOut>
struct VoxelCaster<VectorND<9, TIn>, VectorND<9, TOut> >
{
  static VectorND<9, TOut> Convert(const VectorND<9, TIn> &value)
  {
    VectorND<9, TOut> v;
    for (int i = 0; i < 9; ++i) {
      v(i) = VoxelCaster<TIn, TOut>::Convert(value(i));
    }
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class T>
struct VoxelCaster<VectorND<9, T>, VectorND<9, T> >
{
  static VectorND<9, T> Convert(const VectorND<9, T> &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<VectorND<9, TIn>, Vector>
{
  static Vector Convert(const VectorND<9, TIn> &value)
  {
    Vector v(9);
    for (int i = 0; i < 9; ++i) {
      v(i) = VoxelCaster<TIn, double>::Convert(value(i));
    }
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<VectorND<9, TIn>, float3x3>
{
  static float3x3 Convert(const VectorND<9, TIn> &v)
  {
    float3x3 m;
    m.a.x = VoxelCaster<TIn, float>::Convert(v(0));
    m.a.y = VoxelCaster<TIn, float>::Convert(v(1));
    m.a.z = VoxelCaster<TIn, float>::Convert(v(2));
    m.b.x = VoxelCaster<TIn, float>::Convert(v(3));
    m.b.y = VoxelCaster<TIn, float>::Convert(v(4));
    m.b.z = VoxelCaster<TIn, float>::Convert(v(5));
    m.c.x = VoxelCaster<TIn, float>::Convert(v(6));
    m.c.y = VoxelCaster<TIn, float>::Convert(v(7));
    m.c.z = VoxelCaster<TIn, float>::Convert(v(8));
    return m;
  }
};

// -----------------------------------------------------------------------------
template <class TIn>
struct VoxelCaster<VectorND<9, TIn>, double3x3>
{
  static double3x3 Convert(const VectorND<9, TIn> &v)
  {
    double3x3 m;
    m.a.x = VoxelCaster<TIn, double>::Convert(v(0));
    m.a.y = VoxelCaster<TIn, double>::Convert(v(1));
    m.a.z = VoxelCaster<TIn, double>::Convert(v(2));
    m.b.x = VoxelCaster<TIn, double>::Convert(v(3));
    m.b.y = VoxelCaster<TIn, double>::Convert(v(4));
    m.b.z = VoxelCaster<TIn, double>::Convert(v(5));
    m.c.x = VoxelCaster<TIn, double>::Convert(v(6));
    m.c.y = VoxelCaster<TIn, double>::Convert(v(7));
    m.c.z = VoxelCaster<TIn, double>::Convert(v(8));
    return m;
  }
};


// =============================================================================
// Matrix types
// =============================================================================

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float3x3, TOut>
{
  static TOut Convert(const float3x3 &)
  {
    cerr << "Cannot cast 3x3 matrix to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3x3, float3x3>
{
  static float3x3 Convert(const float3x3 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3x3, double3x3>
{
  static double3x3 Convert(const float3x3 &value)
  {
    return make_double3x3(value);
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<float3x3, VectorND<9, TOut> >
{
  static VectorND<9, TOut> Convert(const float3x3 &value)
  {
    VectorND<9, TOut> v;
    v(0) = VoxelCaster<float, TOut>::Convert(value.a.x);
    v(1) = VoxelCaster<float, TOut>::Convert(value.a.y);
    v(2) = VoxelCaster<float, TOut>::Convert(value.a.z);
    v(3) = VoxelCaster<float, TOut>::Convert(value.b.x);
    v(4) = VoxelCaster<float, TOut>::Convert(value.b.y);
    v(5) = VoxelCaster<float, TOut>::Convert(value.b.z);
    v(6) = VoxelCaster<float, TOut>::Convert(value.c.x);
    v(7) = VoxelCaster<float, TOut>::Convert(value.c.y);
    v(8) = VoxelCaster<float, TOut>::Convert(value.c.z);
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<float3x3, Vector>
{
  static Vector Convert(const float3x3 &value)
  {
    Vector v(9);
    v(0) = value.a.x;
    v(1) = value.a.y;
    v(2) = value.a.z;
    v(3) = value.b.x;
    v(4) = value.b.y;
    v(5) = value.b.z;
    v(6) = value.c.x;
    v(7) = value.c.y;
    v(8) = value.c.z;
    return v;
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double3x3, TOut>
{
  static TOut Convert(const double3x3 &)
  {
    cerr << "Cannot cast 3x3 matrix to a scalar!" << endl;
    cerr << "Set breakpoint in " << __FILE__ << ":" << __LINE__ << " to debug." << endl;
    exit(1);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3x3, float3x3>
{
  static float3x3 Convert(const double3x3 &value)
  {
    return make_float3x3(value);
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3x3, double3x3>
{
  static double3x3 Convert(const double3x3 &value)
  {
    return value;
  }
};

// -----------------------------------------------------------------------------
template <class TOut>
struct VoxelCaster<double3x3, VectorND<9, TOut> >
{
  static VectorND<9, TOut> Convert(const double3x3 &value)
  {
    VectorND<9, TOut> v;
    v(0) = VoxelCaster<double, TOut>::Convert(value.a.x);
    v(1) = VoxelCaster<double, TOut>::Convert(value.a.y);
    v(2) = VoxelCaster<double, TOut>::Convert(value.a.z);
    v(3) = VoxelCaster<double, TOut>::Convert(value.b.x);
    v(4) = VoxelCaster<double, TOut>::Convert(value.b.y);
    v(5) = VoxelCaster<double, TOut>::Convert(value.b.z);
    v(6) = VoxelCaster<double, TOut>::Convert(value.c.x);
    v(7) = VoxelCaster<double, TOut>::Convert(value.c.y);
    v(8) = VoxelCaster<double, TOut>::Convert(value.c.z);
    return v;
  }
};

// -----------------------------------------------------------------------------
template <>
struct VoxelCaster<double3x3, Vector>
{
  static Vector Convert(const double3x3 &value)
  {
    Vector v(9);
    v(0) = value.a.x;
    v(1) = value.a.y;
    v(2) = value.a.z;
    v(3) = value.b.x;
    v(4) = value.b.y;
    v(5) = value.b.z;
    v(6) = value.c.x;
    v(7) = value.c.y;
    v(8) = value.c.z;
    return v;
  }
};


// =============================================================================
// voxel_cast function
// =============================================================================

// -----------------------------------------------------------------------------
template <class TOut, class TIn>
TOut voxel_cast(const TIn &value)
{
  return VoxelCaster<TIn, TOut>::Convert(value);
}


} // namespace mirtk

#endif // MIRTK_VoxelCast_H
