/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#include <mirtkVoxel.h>


namespace mirtk {


// -----------------------------------------------------------------------------
template <> string ToString(const ImageDataType &value, int w, char c, bool left)
{
  const char *str = "unknown";
  switch (value) {
    case MIRTK_VOXEL_CHAR:           str = "char"; break;
    case MIRTK_VOXEL_UNSIGNED_CHAR:  str = "uchar"; break;
    case MIRTK_VOXEL_SHORT:          str = "short"; break;
    case MIRTK_VOXEL_UNSIGNED_SHORT: str = "ushort"; break;
    case MIRTK_VOXEL_INT:            str = "int"; break;
    case MIRTK_VOXEL_UNSIGNED_INT:   str = "uint"; break;
    case MIRTK_VOXEL_FLOAT:          str = "float"; break;
    case MIRTK_VOXEL_DOUBLE:         str = "double"; break;
    case MIRTK_VOXEL_RGB:            str = "RGB"; break;
    case MIRTK_VOXEL_FLOAT1:         str = "float1"; break;
    case MIRTK_VOXEL_FLOAT2:         str = "float2"; break;
    case MIRTK_VOXEL_FLOAT3:         str = "float3"; break;
    case MIRTK_VOXEL_FLOAT4:         str = "float4"; break;
    case MIRTK_VOXEL_FLOAT1x1:       str = "float1x1"; break;
    case MIRTK_VOXEL_FLOAT2x2:       str = "float2x2"; break;
    case MIRTK_VOXEL_FLOAT3x3:       str = "float3x3"; break;
    case MIRTK_VOXEL_FLOAT3x4:       str = "float3x4"; break;
    case MIRTK_VOXEL_FLOAT4x4:       str = "float4x4"; break;
    case MIRTK_VOXEL_DOUBLE1:        str = "double1"; break;
    case MIRTK_VOXEL_DOUBLE2:        str = "double2"; break;
    case MIRTK_VOXEL_DOUBLE3:        str = "double3"; break;
    case MIRTK_VOXEL_DOUBLE4:        str = "double4"; break;
    case MIRTK_VOXEL_DOUBLE1x1:      str = "double1x1"; break;
    case MIRTK_VOXEL_DOUBLE2x2:      str = "double2x2"; break;
    case MIRTK_VOXEL_DOUBLE3x3:      str = "double3x3"; break;
    case MIRTK_VOXEL_DOUBLE3x4:      str = "double3x4"; break;
    case MIRTK_VOXEL_DOUBLE4x4:      str = "double4x4"; break;
    default:                         str = "unknown"; break;
  }
  return ToString(str, c, w, left);
}

// -----------------------------------------------------------------------------
template <> bool FromString(const char *str, ImageDataType &value)
{
  value = static_cast<ImageDataType>(MIRKT_VOXEL_LAST - 1);
  while (value != MIRTK_VOXEL_UNKNOWN) {
    if (ToString(value) == str) break;
    value = static_cast<ImageDataType>(value - 1);
  }
  return value != MIRTK_VOXEL_UNKNOWN;
}

// -----------------------------------------------------------------------------
int DataTypeSize(int type)
{
  switch (type){
    case MIRTK_VOXEL_CHAR:           return sizeof(char);
    case MIRTK_VOXEL_UNSIGNED_CHAR:  return sizeof(unsigned char);
    case MIRTK_VOXEL_SHORT:          return sizeof(short);
    case MIRTK_VOXEL_UNSIGNED_SHORT: return sizeof(unsigned short);
    case MIRTK_VOXEL_INT:            return sizeof(int);
    case MIRTK_VOXEL_UNSIGNED_INT:   return sizeof(unsigned int);
    case MIRTK_VOXEL_FLOAT:          return sizeof(float);
    case MIRTK_VOXEL_DOUBLE:         return sizeof(double);
    case MIRTK_VOXEL_FLOAT1:         return sizeof(float1);
    case MIRTK_VOXEL_FLOAT2:         return sizeof(float2);
    case MIRTK_VOXEL_FLOAT3:         return sizeof(float3);
    case MIRTK_VOXEL_FLOAT4:         return sizeof(float4);
    case MIRTK_VOXEL_FLOAT2x2:       return sizeof(float2x2);
    case MIRTK_VOXEL_FLOAT3x3:       return sizeof(float3x3);
    case MIRTK_VOXEL_FLOAT3x4:       return sizeof(float3x4);
    case MIRTK_VOXEL_FLOAT4x4:       return sizeof(float4x4);
    case MIRTK_VOXEL_DOUBLE1:        return sizeof(double1);
    case MIRTK_VOXEL_DOUBLE2:        return sizeof(double2);
    case MIRTK_VOXEL_DOUBLE3:        return sizeof(double3);
    case MIRTK_VOXEL_DOUBLE4:        return sizeof(double4);
    case MIRTK_VOXEL_DOUBLE2x2:      return sizeof(double2x2);
    case MIRTK_VOXEL_DOUBLE3x3:      return sizeof(double3x3);
    case MIRTK_VOXEL_DOUBLE3x4:      return sizeof(double3x4);
    case MIRTK_VOXEL_DOUBLE4x4:      return sizeof(double4x4);
    default:                         return 0;
  }
}


} // namespace mirtk
