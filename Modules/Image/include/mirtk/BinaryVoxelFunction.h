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

#ifndef MIRTK_BinaryVoxelFunction_H
#define MIRTK_BinaryVoxelFunction_H

#include "mirtk/VoxelFunction.h"

#include "mirtk/NeighborhoodOffsets.h"
#include "mirtk/InterpolateImageFunction.h"


/**
 * These basic binary voxel functions can be used as VoxelFunc template parameter
 * of the binary ForEachVoxel function templates as follows:
 *
 * \code
 * // Add one image to another in-place
 * GreyImage input (attr);
 * GreyImage output(attr);
 * // Set image in-place to maximum of both images
 * BinaryVoxelFunction::Max max;
 * ForEachVoxel(input, output, max);
 * BinaryVoxelFunction::Add add;
 * ParallelForEachVoxel(input, output, add);
 * // Compute sum-of-squared differences (SSD)
 * GreyImage target(attr);
 * GreyImage source(attr);
 * BinaryVoxelFunction::SSD ssd;
 * ForEachVoxel(target, source, ssd);
 * printf("SSD=%f\n", ssd.value);
 * \endcode
 *
 */

namespace mirtk { namespace BinaryVoxelFunction {


// =============================================================================
// Copy
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * Copies the voxel value of one image to the corresponding voxel of another
 */
struct Copy : public VoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out) { *out = *in; }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// =============================================================================
// Basic mathematical voxel-wise in-place operations
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * Add image to another
 */
struct Add : public VoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out)
  {
    *out = static_cast<T2>(static_cast<double>(*out) + static_cast<double>(*in));
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// -----------------------------------------------------------------------------
/**
 * Subtract image from another
 */
struct Sub : public VoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out)
  {
    *out = static_cast<T2>(static_cast<double>(*out) - static_cast<double>(*in));
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// -----------------------------------------------------------------------------
/**
 * Multiplies the voxel value of one image by the value of another
 */
struct Mul : public VoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out)
  {
    *out = static_cast<T2>(static_cast<double>(*out) * static_cast<double>(*in));
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// -----------------------------------------------------------------------------
/**
 * Divides the voxel value of one image by the value of another
 */
struct Div : public VoxelFunction
{
  template <class T1, class T2>
  void operator ()(const T1 *in, T2 *out)
  {
    double divisor = static_cast<double>(*in);
    if (divisor == .0) {
      *out = static_cast<T2>(0);
    } else {
      *out = static_cast<T2>(static_cast<double>(*out) / divisor);
    }
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// =============================================================================
// Morphological operators
// =============================================================================

/**
 * Replace voxel value by minimum of neighboring voxels
 */
class Erode : public VoxelFunction
{
  const NeighborhoodOffsets &_Offsets;

public:

  Erode(const NeighborhoodOffsets &offsets) : _Offsets(offsets) {}

  template <class T>
  void operator ()(int i, int j, int k, int, const T *in, T *out)
  {
    const T *ptr2offset;
    T &value = *out;
    if (_Domain->IsBoundary(i, j, k)) {
      value = *in;
    } else {
      value = *in;
      for (int n = 0; n < _Offsets.Size(); ++n) {
        ptr2offset = in + _Offsets(n);
        if (*ptr2offset < value) value = *ptr2offset;
      }
    }
  }
};

/**
 * Replace voxel value by maximum of neighboring voxels
 */
class Dilate : public VoxelFunction
{
  const NeighborhoodOffsets &_Offsets;

public:

  Dilate(const NeighborhoodOffsets &offsets) : _Offsets(offsets) {}

  template <class T>
  void operator ()(int i, int j, int k, int, const T *in, T *out)
  {
    const T *ptr2offset;
    T &value = *out;
    if (_Domain->IsBoundary(i, j, k)) {
      value = *in;
    } else {
      value = *in;
      for (int n = 0; n < _Offsets.Size(); ++n) {
        ptr2offset = in + _Offsets(n);
        if (*ptr2offset > value) value = *ptr2offset;
      }
    }
  }
};

// =============================================================================
// Image similarity measures
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * Compute sum-of-squared differences (SSD)
 */
struct SSD : public VoxelReduction
{
  double value;

  SSD()       : value(.0)      {}
  SSD(SSD &o) : value(o.value) {}

  void split(const SSD &)    { value = .0; }
  void join (const SSD &rhs) { value += rhs.value; }

  template <class T1, class T2>
  void operator ()(const T1 *in1, const T2 *in2)
  {
    double diff = static_cast<double>(*in1) - static_cast<double>(*in2);
    value += diff * diff;
  }

  template <class TImage, class T1, class T2>
  void operator ()(const TImage&, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }

  template <class T1, class T2>
  void operator ()(int, int, int, int, const T1 *in, T2 *out)
  {
    this->operator ()(in, out);
  }
};

// =============================================================================
// Composition
// =============================================================================
  
// -----------------------------------------------------------------------------
/// Compose two 2D displacement fields: D = D1 ° D2
template <class TReal, class TInterpolator = InterpolateImageFunction>
struct ComposeDisplacementFields2D : VoxelFunction
{
  // ---------------------------------------------------------------------------
  ComposeDisplacementFields2D(TInterpolator *d1, const GenericImage<TReal> *d2)
  :
    _D1(d1->Input()),
    _D1Interpolator(d1),
    _D2(d2),
    _y(d2->X() * d2->Y())
  {
  }
  
  // ---------------------------------------------------------------------------
  void operator()(int i, int j, int k, int, const TReal *d2, TReal *dout)
  {
    double d[3]; // 2D displacement field can have constant third component!
    double x1 = i, y1 = j, z1 = k;
    _D2->ImageToWorld(x1, y1, z1);
    double x2 = x1 + d2[_x]; double x = x2;
    double y2 = y1 + d2[_y]; double y = y2;
    _D1            ->WorldToImage(x, y, z1);
    _D1Interpolator->Evaluate (d, x, y, z1);
    x2 = x2 + d[0];
    y2 = y2 + d[1];
    dout[_x] = static_cast<TReal>(x2 - x1);
    dout[_y] = static_cast<TReal>(y2 - y1);
  }
  
private:
  const BaseImage           *_D1;             ///< First displacement field
  TInterpolator             *_D1Interpolator; ///< Interpolator of first displacement field
  const GenericImage<TReal> *_D2;             ///< Second displacement field
  
  static const int _x = 0; ///< Offset of x component
  int              _y;     ///< Offset of y component
};

// -----------------------------------------------------------------------------
/// Compose two 3D displacement fields: D = D1 ° D2
template <class TReal, class TInterpolator = InterpolateImageFunction>
struct ComposeDisplacementFields3D : VoxelFunction
{
  // ---------------------------------------------------------------------------
  ComposeDisplacementFields3D(TInterpolator *d1, const GenericImage<TReal> *d2)
  :
    _D1(d1->Input()),
    _D1Interpolator(d1),
    _D2(d2),
    _y(d2->X() * d2->Y() * d2->Z()),
    _z(2 * _y)
  {
  }
  
  // ---------------------------------------------------------------------------
  void operator()(int i, int j, int k, int, const TReal *d2, TReal *dout)
  {
    double d[3];
    double x1 = i, y1 = j, z1 = k;
    _D2->ImageToWorld(x1, y1, z1);
    double x2 = x1 + d2[_x]; double x = x2;
    double y2 = y1 + d2[_y]; double y = y2;
    double z2 = z1 + d2[_z]; double z = z2;
    _D1            ->WorldToImage(x, y, z);
    _D1Interpolator->Evaluate (d, x, y, z);
    x2 = x2 + d[0];
    y2 = y2 + d[1];
    z2 = z2 + d[2];
    dout[_x] = static_cast<TReal>(x2 - x1);
    dout[_y] = static_cast<TReal>(y2 - y1);
    dout[_z] = static_cast<TReal>(z2 - z1);
  }
  
private:
  const BaseImage           *_D1;             ///< First displacement field
  TInterpolator             *_D1Interpolator; ///< Interpolator of first displacement field
  const GenericImage<TReal> *_D2;             ///< Second displacement field
  
  static const int _x = 0; ///< Offset of x component
  int              _y;     ///< Offset of y component
  int              _z;     ///< Offset of z component
};


} } // namespace mirtk::BinaryVoxelFunction

#endif
