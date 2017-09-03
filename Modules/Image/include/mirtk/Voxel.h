/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_Voxel_H
#define MIRTK_Voxel_H

#include "mirtk/Math.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Vector4D.h"
#include "mirtk/VectorND.h"
#include "mirtk/String.h"

#ifdef HAVE_VTK 
  #include "vtkType.h"
#endif


namespace mirtk {

// =============================================================================
// Voxel type enumeration
// =============================================================================

/// Enumeration of voxel data types
enum ImageDataType
{
  MIRTK_VOXEL_UNKNOWN,
  MIRTK_VOXEL_CHAR,
  MIRTK_VOXEL_UNSIGNED_CHAR,
  MIRTK_VOXEL_SHORT,
  MIRTK_VOXEL_UNSIGNED_SHORT,
  MIRTK_VOXEL_INT,
  MIRTK_VOXEL_UNSIGNED_INT,
  MIRTK_VOXEL_FLOAT,
  MIRTK_VOXEL_DOUBLE,
  MIRTK_VOXEL_RGB,
  MIRTK_VOXEL_FLOAT1, // unused
  MIRTK_VOXEL_FLOAT2,
  MIRTK_VOXEL_FLOAT3,
  MIRTK_VOXEL_FLOAT4,
  MIRTK_VOXEL_FLOAT9,
  MIRTK_VOXEL_DOUBLE1, // unused
  MIRTK_VOXEL_DOUBLE2,
  MIRTK_VOXEL_DOUBLE3,
  MIRTK_VOXEL_DOUBLE4,
  MIRTK_VOXEL_DOUBLE9,
  MIRTK_VOXEL_FLOAT1x1, // unused
  MIRTK_VOXEL_FLOAT2x2,
  MIRTK_VOXEL_FLOAT3x3,
  MIRTK_VOXEL_FLOAT3x4,
  MIRTK_VOXEL_FLOAT4x4,
  MIRTK_VOXEL_DOUBLE1x1, // unused
  MIRTK_VOXEL_DOUBLE2x2,
  MIRTK_VOXEL_DOUBLE3x3,
  MIRTK_VOXEL_DOUBLE3x4,
  MIRTK_VOXEL_DOUBLE4x4,
  // Last entry of unique enumeration values
  MIRKT_VOXEL_LAST,
  // Aliases for common voxel types
  MIRTK_VOXEL_BINARY = MIRTK_VOXEL_UNSIGNED_CHAR,
  MIRTK_VOXEL_BYTE   = MIRTK_VOXEL_UNSIGNED_CHAR,
  MIRTK_VOXEL_GREY   = MIRTK_VOXEL_SHORT,
#if MIRTK_USE_FLOAT_BY_DEFAULT
  MIRTK_VOXEL_REAL    = MIRTK_VOXEL_FLOAT,
  MIRTK_VOXEL_REAL1   = MIRTK_VOXEL_FLOAT1,
  MIRTK_VOXEL_REAL2   = MIRTK_VOXEL_FLOAT2,
  MIRTK_VOXEL_REAL3   = MIRTK_VOXEL_FLOAT3,
  MIRTK_VOXEL_REAL4   = MIRTK_VOXEL_FLOAT4,
  MIRTK_VOXEL_REAL9   = MIRTK_VOXEL_FLOAT9,
  MIRTK_VOXEL_REAL1x1 = MIRTK_VOXEL_FLOAT1x1,
  MIRTK_VOXEL_REAL2x2 = MIRTK_VOXEL_FLOAT2x2,
  MIRTK_VOXEL_REAL3x3 = MIRTK_VOXEL_FLOAT3x3,
  MIRTK_VOXEL_REAL3x4 = MIRTK_VOXEL_FLOAT3x4,
  MIRTK_VOXEL_REAL4x4 = MIRTK_VOXEL_FLOAT4x4
#else // MIRTK_USE_FLOAT_BY_DEFAULT
  MIRTK_VOXEL_REAL    = MIRTK_VOXEL_DOUBLE,
  MIRTK_VOXEL_REAL1   = MIRTK_VOXEL_DOUBLE1,
  MIRTK_VOXEL_REAL2   = MIRTK_VOXEL_DOUBLE2,
  MIRTK_VOXEL_REAL3   = MIRTK_VOXEL_DOUBLE3,
  MIRTK_VOXEL_REAL4   = MIRTK_VOXEL_DOUBLE4,
  MIRTK_VOXEL_REAL9   = MIRTK_VOXEL_DOUBLE9,
  MIRTK_VOXEL_REAL1x1 = MIRTK_VOXEL_DOUBLE1x1,
  MIRTK_VOXEL_REAL2x2 = MIRTK_VOXEL_DOUBLE2x2,
  MIRTK_VOXEL_REAL3x3 = MIRTK_VOXEL_DOUBLE3x3,
  MIRTK_VOXEL_REAL3x4 = MIRTK_VOXEL_DOUBLE3x4,
  MIRTK_VOXEL_REAL4x4 = MIRTK_VOXEL_DOUBLE4x4
#endif // MIRTK_USE_FLOAT_BY_DEFAULT
};

/// Convert image data type enumeration value to string
template <> string ToString(const ImageDataType &, int, char, bool);

/// Convert string to image data type enumeration value
template <> bool FromString(const char *, ImageDataType &);

// =============================================================================
// Conversion from/to VTK data types
// =============================================================================
#ifdef HAVE_VTK

// -----------------------------------------------------------------------------
/// Get VTK data type from IRTK voxel type
inline int ToVTKDataType(int type)
{
  switch (type) {
    case MIRTK_VOXEL_CHAR:           return VTK_CHAR;
    case MIRTK_VOXEL_UNSIGNED_CHAR:  return VTK_UNSIGNED_CHAR;
    case MIRTK_VOXEL_SHORT:          return VTK_SHORT;
    case MIRTK_VOXEL_UNSIGNED_SHORT: return VTK_UNSIGNED_SHORT;
    case MIRTK_VOXEL_INT:            return VTK_INT;
    case MIRTK_VOXEL_UNSIGNED_INT:   return VTK_UNSIGNED_INT;
    case MIRTK_VOXEL_FLOAT:          return VTK_FLOAT;
    case MIRTK_VOXEL_DOUBLE:         return VTK_DOUBLE;
    default:                         return VTK_VOID;
  }
}

// -----------------------------------------------------------------------------
/// Get IRTK voxel type from VTK data type
inline int FromVTKDataType(int type)
{
  switch (type) {
    case VTK_CHAR:           return MIRTK_VOXEL_CHAR;
    case VTK_UNSIGNED_CHAR:  return MIRTK_VOXEL_UNSIGNED_CHAR;
    case VTK_SHORT:          return MIRTK_VOXEL_SHORT;
    case VTK_UNSIGNED_SHORT: return MIRTK_VOXEL_UNSIGNED_SHORT;
    case VTK_INT:            return MIRTK_VOXEL_INT;
    case VTK_UNSIGNED_INT:   return MIRTK_VOXEL_UNSIGNED_INT;
    case VTK_FLOAT:          return MIRTK_VOXEL_FLOAT;
    case VTK_DOUBLE:         return MIRTK_VOXEL_DOUBLE;
    default:                 return MIRTK_VOXEL_UNKNOWN;
  }
}

#endif // HAVE_VTK
// =============================================================================
// Voxel types
// =============================================================================

typedef unsigned char BinaryPixel;
typedef unsigned char BytePixel;
typedef short         GreyPixel;
typedef realt         RealPixel;

typedef Vector3D<float>     Float3;
typedef Vector4D<float>     Float4;
typedef VectorND<9, float>  Float9;
typedef Vector3D<double>    Double3;
typedef Vector4D<double>    Double4;
typedef VectorND<9, double> Double9;

// =============================================================================
// Voxel type limits
// =============================================================================

const GreyPixel MIN_GREY = numeric_limits<GreyPixel>::min();
const GreyPixel MAX_GREY = numeric_limits<GreyPixel>::max();

// -----------------------------------------------------------------------------
template <class T>
struct voxel_limits
{
  /// Minimum value that can be represented by this voxel type as double
  /// \note For vector types, this corresponds to the minium value that can
  ///       be represented by each component of the vector.
  static double min() throw();
  /// Maximum value that can be represented by this voxel type as double
  /// \note For vector types, this corresponds to the maxium value that can
  ///       be represented by each component of the vector.
  static double max() throw();
  /// Minimum value that can be represented by this voxel type
  static T min_value() throw();
  /// Maximum value that can be represented by this voxel type
  static T max_value() throw();
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<char>
{
  static char   min_value() throw() { return numeric_limits<char>::lowest(); }
  static char   max_value() throw() { return numeric_limits<char>::max(); }
  static double min()       throw() { return static_cast<double>(min_value()); }
  static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<unsigned char>
{
  static unsigned char min_value() throw() { return numeric_limits<unsigned char>::lowest(); }
  static unsigned char max_value() throw() { return numeric_limits<unsigned char>::max(); }
  static double        min()       throw() { return static_cast<double>(min_value()); }
  static double        max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<short>
{
  static short  min_value() throw() { return numeric_limits<short>::lowest(); }
  static short  max_value() throw() { return numeric_limits<short>::max(); }
  static double min()       throw() { return static_cast<double>(min_value()); }
  static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<unsigned short>
{
  static unsigned short min_value() throw() { return numeric_limits<unsigned short>::lowest(); }
  static unsigned short max_value() throw() { return numeric_limits<unsigned short>::max(); }
  static double         min()       throw() { return static_cast<double>(min_value()); }
  static double         max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<int> {
  static int    min_value() throw() { return numeric_limits<int>::lowest(); }
  static int    max_value() throw() { return numeric_limits<int>::max();  }
  static double min()       throw() { return static_cast<double>(min_value()); }
  static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<unsigned int>
{
  static unsigned int min_value() throw() { return numeric_limits<unsigned int>::lowest(); }
  static unsigned int max_value() throw() { return numeric_limits<unsigned int>::max(); }
  static double       min()       throw() { return static_cast<double      >(min_value()); }
  static double       max()       throw() { return static_cast<double      >(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float>
{
  static float  min_value() throw() { return numeric_limits<float>::lowest(); }
  static float  max_value() throw() { return numeric_limits<float>::max(); }
  static double min()       throw() { return static_cast<double>(min_value()); }
  static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double>
{
  static double min_value() throw() { return numeric_limits<double>::lowest(); }
  static double max_value() throw() { return numeric_limits<double>::max(); }
  static double min()       throw() { return static_cast<double>(min_value()); }
  static double max()       throw() { return static_cast<double>(max_value()); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float1>
{
  static float1 min_value() throw() { return make_float1(voxel_limits<float>::min_value()); }
  static float1 max_value() throw() { return make_float1(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
  static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float2>
{
  static float2 min_value() throw() { return make_float2(voxel_limits<float>::min_value()); }
  static float2 max_value() throw() { return make_float2(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
  static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float3>
{
  static float3 min_value() throw() { return make_float3(voxel_limits<float>::min_value()); }
  static float3 max_value() throw() { return make_float3(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
  static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float4>
{
  static float4 min_value() throw() { return make_float4(voxel_limits<float>::min_value()); }
  static float4 max_value() throw() { return make_float4(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
  static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<float3x3>
{
  static float3x3 min_value() throw() { return make_float3x3(voxel_limits<float>::min_value()); }
  static float3x3 max_value() throw() { return make_float3x3(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
  static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double1>
{
  static double1 min_value() throw() { return make_double1(voxel_limits<double>::min_value()); }
  static double1 max_value() throw() { return make_double1(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
  static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double2>
{
  static double2 min_value() throw() { return make_double2(voxel_limits<double>::min_value()); }
  static double2 max_value() throw() { return make_double2(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
  static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double3>
{
  static double3 min_value() throw() { return make_double3(voxel_limits<double>::min_value()); }
  static double3 max_value() throw() { return make_double3(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
  static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double4>
{
  static double4 min_value() throw() { return make_double4(voxel_limits<double>::min_value()); }
  static double4 max_value() throw() { return make_double4(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
  static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<double3x3>
{
  static double3x3 min_value() throw() { return make_double3x3(voxel_limits<double>::min_value()); }
  static double3x3 max_value() throw() { return make_double3x3(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
  static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<Float3>
{
  static Float3 min_value() throw()
  { return Float3(voxel_limits<float>::min_value()); }
  static Float3 max_value() throw()
  { return Float3(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
  static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<Double3>
{
  static Double3 min_value() throw()
  { return Double3(voxel_limits<double>::min_value()); }
  static Double3 max_value() throw()
  { return Double3(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
  static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<Float4>
{
  static Float4 min_value() throw()
  { return Float4(voxel_limits<float>::min_value()); }
  static Float4 max_value() throw()
  { return Float4(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
  static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<Double4>
{
  static Double4 min_value() throw()
  { return Double4(voxel_limits<double>::min_value()); }
  static Double4 max_value() throw()
  { return Double4(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
  static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<Float9>
{
  static Float9 min_value() throw()
  { return Float9(voxel_limits<float>::min_value()); }
  static Float9 max_value() throw()
  { return Float9(voxel_limits<float>::max_value()); }
  static double min() throw() { return voxel_limits<float>::min(); }
  static double max() throw() { return voxel_limits<float>::max(); }
};

// -----------------------------------------------------------------------------
template <> struct voxel_limits<Double9>
{
  static Double9 min_value() throw()
  { return Double9(voxel_limits<double>::min_value()); }
  static Double9 max_value() throw()
  { return Double9(voxel_limits<double>::max_value()); }
  static double min() throw() { return voxel_limits<double>::min(); }
  static double max() throw() { return voxel_limits<double>::max(); }
};

// -----------------------------------------------------------------------------
// Variable length vector type not allowed as actual voxel type of an image
// instance. Only used as voxel type by base class methods and general
// interpolators. Treat voxel type as if it was a scalar type here.
template <> struct voxel_limits<Vector>
{
  static Vector min_value() throw() { return Vector(1, min()); }
  static Vector max_value() throw() { return Vector(1, max()); }
  static double min()       throw() { return numeric_limits<double>::lowest(); }
  static double max()       throw() { return numeric_limits<double>::max(); }
};

// =============================================================================
// Voxel type information
// =============================================================================

// -----------------------------------------------------------------------------
int DataTypeSize(int);

// -----------------------------------------------------------------------------
inline string DataTypeName(int type)
{
  return ToString(static_cast<ImageDataType>(type));
}

// -----------------------------------------------------------------------------
inline int ToDataType(const char *str)
{
  ImageDataType type;
  if (FromString(str, type)) return type;
  return MIRTK_VOXEL_UNKNOWN;
}

// -----------------------------------------------------------------------------
inline int ToDataType(const string &str)
{
  return ToDataType(str.c_str());
}

// -----------------------------------------------------------------------------
template <class T>
struct voxel_info : public voxel_limits<T>
{
  /// Scalar type compatible with this voxel type
  typedef T ScalarType;
  /// Floating point type compatible with this voxel type
  typedef RealPixel RealType;
  /// Number of (vector) elements stored by this voxel
  static int vector_size() throw();
  /// Enumeration value corresponding to voxel type of (vector) elements
  static int element_type() throw();
  /// Enumeration value corresponding to this voxel type
  static int type() throw();
  /// Minimum value that can be represented by this voxel type as double
  /// \note For vector types, this corresponds to the minium value that can
  ///       be represented by each component of the vector.
  using voxel_limits<T>::min;
  /// Maximum value that can be represented by this voxel type as double
  /// \note For vector types, this corresponds to the maxium value that can
  ///       be represented by each component of the vector.
  using voxel_limits<T>::max;
  /// Minimum value that can be represented by this voxel type
  using voxel_limits<T>::min_value;
  /// Maximum value that can be represented by this voxel type
  using voxel_limits<T>::max_value;
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<char>
{
  typedef char        ScalarType;
  typedef RealPixel   RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_CHAR; }
  static int type()         throw() { return MIRTK_VOXEL_CHAR; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<unsigned char>
{
  typedef unsigned char ScalarType;
  typedef RealPixel     RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_UNSIGNED_CHAR; }
  static int type()         throw() { return MIRTK_VOXEL_UNSIGNED_CHAR; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<short>
{
  typedef short     ScalarType;
  typedef RealPixel RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_SHORT; }
  static int type()         throw() { return MIRTK_VOXEL_SHORT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<unsigned short>
{
  typedef unsigned short ScalarType;
  typedef RealPixel      RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_UNSIGNED_SHORT; }
  static int type()         throw() { return MIRTK_VOXEL_UNSIGNED_SHORT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<int>
{
  typedef int       ScalarType;
  typedef RealPixel RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_INT; }
  static int type()         throw() { return MIRTK_VOXEL_INT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<unsigned int>
{
  typedef unsigned int ScalarType;
  typedef RealPixel    RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_UNSIGNED_INT; }
  static int type()         throw() { return MIRTK_VOXEL_UNSIGNED_INT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float>
{
  typedef float1   Type1;
  typedef float2   Type2;
  typedef float3   Type3;
  typedef float4   Type4;
  typedef float2x2 Type2x2;
  typedef float3x3 Type3x3;
  typedef float3x4 Type3x4;
  typedef float4x4 Type4x4;
  typedef float    ScalarType;
  typedef float    RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_FLOAT; }
  static int type()         throw() { return MIRTK_VOXEL_FLOAT; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float1>
{
  typedef float  ScalarType;
  typedef float1 RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_FLOAT; }
  static int type()         throw() { return MIRTK_VOXEL_FLOAT1; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float2>
{
  typedef float  ScalarType;
  typedef float2 RealType;
  static int vector_size()  throw() { return 2; }
  static int element_type() throw() { return MIRTK_VOXEL_FLOAT; }
  static int type()         throw() { return MIRTK_VOXEL_FLOAT2; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float3>
{
  typedef float  ScalarType;
  typedef float3 RealType;
  static int vector_size()  throw() { return 3; }
  static int element_type() throw() { return MIRTK_VOXEL_FLOAT; }
  static int type()         throw() { return MIRTK_VOXEL_FLOAT3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float4>
{
  typedef float  ScalarType;
  typedef float4 RealType;
  static int vector_size()  throw() { return 4; }
  static int element_type() throw() { return MIRTK_VOXEL_FLOAT; }
  static int type()         throw() { return MIRTK_VOXEL_FLOAT4; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<float3x3>
{
  typedef float    ScalarType;
  typedef float3x3 RealType;
  static int vector_size()  throw() { return 9; }
  static int element_type() throw() { return MIRTK_VOXEL_FLOAT; }
  static int type()         throw() { return MIRTK_VOXEL_FLOAT3x3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double>
{
  typedef double1   Type1;
  typedef double2   Type2;
  typedef double3   Type3;
  typedef double4   Type4;
  typedef double2x2 Type2x2;
  typedef double3x3 Type3x3;
  typedef double3x4 Type3x4;
  typedef double4x4 Type4x4;
  typedef double    ScalarType;
  typedef double    RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double1>
{
  typedef double  ScalarType;
  typedef double1 RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE1; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double2>
{
  typedef double  ScalarType;
  typedef double2 RealType;
  static int vector_size()  throw() { return 2; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE2; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double3>
{
  typedef double  ScalarType;
  typedef double3 RealType;
  static int vector_size()  throw() { return 3; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double4>
{
  typedef double  ScalarType;
  typedef double4 RealType;
  static int vector_size()  throw() { return 4; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE4; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<double3x3>
{
  typedef double    ScalarType;
  typedef double3x3 RealType;
  static int vector_size()  throw() { return 9; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE3x3; }
};

// -----------------------------------------------------------------------------
// Variable length vector type not allowed as actual voxel type of an image
// instance. Only used as voxel type by base class methods and general
// interpolators. Treat voxel type as if it was a scalar type here.
template <> struct voxel_info<Vector>
{
  typedef double ScalarType;
  typedef double RealType;
  static int vector_size()  throw() { return 1; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<Float3>
{
  typedef float  ScalarType;
  typedef Float3 RealType;
  static int vector_size()  throw() { return 3; }
  static int element_type() throw() { return MIRTK_VOXEL_FLOAT; }
  static int type()         throw() { return MIRTK_VOXEL_FLOAT3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<Double3>
{
  typedef double  ScalarType;
  typedef Double3 RealType;
  static int vector_size()  throw() { return 3; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE3; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<Float4>
{
  typedef float  ScalarType;
  typedef Float4 RealType;
  static int vector_size()  throw() { return 4; }
  static int element_type() throw() { return MIRTK_VOXEL_FLOAT; }
  static int type()         throw() { return MIRTK_VOXEL_FLOAT4; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<Double4>
{
  typedef double  ScalarType;
  typedef Double4 RealType;
  static int vector_size()  throw() { return 4; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE4; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<Float9>
{
  typedef float  ScalarType;
  typedef Float9 RealType;
  static int vector_size()  throw() { return 9; }
  static int element_type() throw() { return MIRTK_VOXEL_FLOAT; }
  static int type()         throw() { return MIRTK_VOXEL_FLOAT4; }
};

// -----------------------------------------------------------------------------
template <> struct voxel_info<Double9>
{
  typedef double  ScalarType;
  typedef Double9 RealType;
  static int vector_size()  throw() { return 9; }
  static int element_type() throw() { return MIRTK_VOXEL_DOUBLE; }
  static int type()         throw() { return MIRTK_VOXEL_DOUBLE4; }
};


} // namespace mirtk

#endif // MIRTK_Voxel_H
