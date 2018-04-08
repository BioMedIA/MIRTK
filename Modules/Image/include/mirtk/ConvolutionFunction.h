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

#ifndef MIRTK_ConvolutionFunction_H
#define MIRTK_ConvolutionFunction_H

#include "mirtk/Voxel.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/BaseImage.h"
#include "mirtk/GenericImage.h"


namespace mirtk {

/// Voxel functions implementing convolution of an image with a discrete kernel
namespace ConvolutionFunction {


/// Tolerance used for floating point comparisons
inline double Epsilon()
{
  return 1e-6;
}


// =============================================================================
// Boundary conditions for image domain
// =============================================================================

// -----------------------------------------------------------------------------
/// Mirror image at boundary
struct MirrorBoundaryCondition
{
  /// Apply mirror boundary condition to voxel index
  ///
  /// \param[in] i Voxel index along image dimension
  /// \param[in] N Number of voxels in image dimension
  int operator ()(int i, int N) const
  {
    N -= 1;
    int m, n;
    if (i < 0) {
      i = -i;
      n = i / N;
      m = i - n * N;
      if (n % 2) i = N - 1 - m;
      else       i = m;
    } else if (i > N) {
      i = i - N;
      n = i / N;
      m = i - n * N;
      if (n % 2) i = m;
      else       i = N - m;
    }
    return i;
  }

  /// Apply mirror boundary condition to voxel pointer
  template <class T>
  const T *operator ()(int i, int N, const T *p, int stride)
  {
    return p + (this->operator ()(i, N) - i) * stride;
  }
};

// =============================================================================
// Convolve at boundary truncated image with 1D kernel
// =============================================================================

// -----------------------------------------------------------------------------
/// Perform convolution along x with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveInX : public VoxelFunction
{
  const TKernel *_Kernel;    ///< Convolution kernel
  int            _Radius;    ///< Radius of kernel
  bool           _Normalize; ///< Wether to divide by sum of kernel weights
  int            _Offset;    ///< Vector component offset
  int            _X;         ///< Number of voxels in x

  ConvolveInX(const BaseImage *image, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    _Kernel(kernel), _Radius((size - 1) / 2), _Normalize(norm),
    _Offset(l * image->NumberOfSpatialVoxels()), _X(image->X())
  {}

  template <class T1, class T2>
  void operator ()(int i, int, int, int, const T1 *in, T2 *out) const
  {
    typedef typename voxel_info<T2>::RealType         AccType;
    typedef typename voxel_info<AccType>::ScalarType  RealType;
    // Convolve l-th vector component
    in += _Offset, out += _Offset;
    // Go to start of input and end of kernel
    i    -= _Radius;
    in   -= _Radius;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (i < 0) n += i, in -= i, i = 0;
    // Inside image domain
    AccType  v, acc(0.);
    RealType w, sum = .0;
    while (i < _X && n >= 0) {
      w = static_cast<RealType>(_Kernel[n]);
      v = voxel_cast<AccType>(*in);
      acc += (v *= w);
      sum += w;
      ++i, ++in, --n;
    }
    // Output result of convolution
    if (_Normalize && !AreEqual(sum, 0., Epsilon())) acc /= sum;
    *out = voxel_cast<T2>(acc);
  }
};

// -----------------------------------------------------------------------------
/// Perform convolution along y with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveInY : public VoxelFunction
{
  const TKernel *_Kernel;    ///< Convolution kernel
  int            _Radius;    ///< Radius of kernel
  bool           _Normalize; ///< Wether to divide by sum of kernel weights
  int            _Offset;    ///< Vector component offset
  int            _X;         ///< Number of voxels in x
  int            _Y;         ///< Number of voxels in y

  ConvolveInY(const BaseImage *image, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    _Kernel(kernel), _Radius((size - 1) / 2), _Normalize(norm),
    _Offset(l * image->NumberOfSpatialVoxels()), _X(image->X()), _Y(image->Y())
  {}

  template <class T1, class T2>
  void operator ()(int, int j, int, int, const T1 *in, T2 *out) const
  {
    typedef typename voxel_info<T2>::RealType         AccType;
    typedef typename voxel_info<AccType>::ScalarType  RealType;
    // Convolve l-th vector component
    in += _Offset, out += _Offset;
    // Go to start of input and end of kernel
    j    -= _Radius;
    in   -= _Radius * _X;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (j < 0) n += j, in -= j * _X, j = 0;
    // Inside image domain
    AccType  v, acc = 0.;
    RealType w, sum = 0.;
    while (j < _Y && n >= 0) {
      w = static_cast<RealType>(_Kernel[n]);
      v = voxel_cast<AccType>(*in);
      acc += (v *= w);
      sum += w;
      ++j, in += _X, --n;
    }
    // Output result of convolution
    if (_Normalize && !AreEqual(sum, 0., Epsilon())) acc /= sum;
    *out = voxel_cast<T2>(acc);
  }
};

// -----------------------------------------------------------------------------
/// Perform convolution along z with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveInZ : public VoxelFunction
{
  const TKernel *_Kernel;    ///< Convolution kernel
  int            _Radius;    ///< Radius of kernel
  bool           _Normalize; ///< Wether to divide by sum of kernel weights
  int            _Offset;    ///< Vector component offset
  int            _XY;        ///< Number of voxels in slice (nx * ny)
  int            _Z;         ///< Number of voxels in z

  ConvolveInZ(const BaseImage *image, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    _Kernel(kernel), _Radius((size - 1) / 2), _Normalize(norm),
    _Offset(l * image->NumberOfSpatialVoxels()), _XY(image->X() * image->Y()), _Z(image->Z())
  {}

  template <class T1, class T2>
  void operator ()(int, int, int k, int, const T1 *in, T2 *out) const
  {
    typedef typename voxel_info<T2>::RealType         AccType;
    typedef typename voxel_info<AccType>::ScalarType  RealType;
    // Convolve l-th vector component
    in += _Offset, out += _Offset;
    // Go to start of input and end of kernel
    k    -= _Radius;
    in   -= _Radius * _XY;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (k < 0) n += k, in -= k * _XY, k = 0;
    // Inside image domain
    AccType  v, acc = 0.;
    RealType w, sum = 0.;
    while (k < _Z && n >= 0) {
      w = static_cast<RealType>(_Kernel[n]);
      v = voxel_cast<AccType>(*in);
      acc += (v *= w);
      sum += w;
      ++k, in += _XY, --n;
    }
    // Output result of convolution
    if (_Normalize && !AreEqual(sum, 0., Epsilon())) acc /= sum;
    *out = voxel_cast<T2>(acc);
  }
};

// -----------------------------------------------------------------------------
/// Perform convolution along t with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveInT : public VoxelFunction
{
  const TKernel *_Kernel;    ///< Convolution kernel
  int            _Radius;    ///< Radius of kernel
  bool           _Normalize; ///< Wether to divide by sum of kernel weights
  int            _XYZ;       ///< Number of voxels in volume (nx * ny * nz)
  int            _T;         ///< Number of voxels in t

  ConvolveInT(const BaseImage *image, const TKernel *kernel, int size, bool norm = true)
  :
    _Kernel(kernel), _Radius((size - 1) / 2), _Normalize(norm),
    _XYZ(image->X() * image->Y() * image->Z()), _T(image->T())
  {}

  template <class T1, class T2>
  void operator ()(int, int, int, int l, const T1 *in, T2 *out) const
  {
    typedef typename voxel_info<T2>::RealType         AccType;
    typedef typename voxel_info<AccType>::ScalarType  RealType;
    // Go to start of input and end of kernel
    l    -= _Radius;
    in   -= _Radius * _XYZ;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (l < 0) n += l, in -= l * _XYZ, l = 0;
    // Inside image domain
    AccType  v, acc = 0.;
    RealType w, sum = 0.;
    while (l < _T && n >= 0) {
      w = static_cast<RealType>(_Kernel[n]);
      v = voxel_cast<AccType>(*in);
      acc += (v *= w);
      sum += w;
      ++l, in += _XYZ, --n;
    }
    // Output result of convolution
    if (_Normalize && !AreEqual(sum, 0., Epsilon())) acc /= sum;
    *out = voxel_cast<T2>(acc);
  }
};

// =============================================================================
// Convolve at boundary truncated image foreground with 1D kernel
// =============================================================================

// -----------------------------------------------------------------------------
/// Perform convolution of image foreground along x with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveForegroundInX : public VoxelFunction
{
  const BaseImage *_Image;     ///< Image defining foreground region
  const TKernel   *_Kernel;    ///< Convolution kernel
  int              _Radius;    ///< Radius of kernel
  bool             _Normalize; ///< Wether to divide by sum of kernel weights
  int              _Offset;    ///< Vector component offset
  int              _X;         ///< Number of voxels in x

  ConvolveForegroundInX(const BaseImage *image, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    _Image(image), _Kernel(kernel), _Radius((size - 1) / 2), _Normalize(norm),
    _Offset(l * image->NumberOfSpatialVoxels()), _X(image->X())
  {}

  template <class T1, class T2>
  void operator ()(int i, int j, int k, int l, const T1 *in, T2 *out) const
  {
    typedef typename voxel_info<T2>::RealType         AccType;
    typedef typename voxel_info<AccType>::ScalarType  RealType;
    // Convolve l-th vector component
    in += _Offset, out += _Offset;
    // Go to start of input and end of kernel
    i    -= _Radius;
    in   -= _Radius;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (i < 0) n += i, in -= i, i = 0;
    // Inside image domain
    AccType  v, acc = 0.;
    RealType w, sum = 0.;
    while (i < _X && n >= 0) {
      if (_Image->IsForeground(i, j, k, l)) {
        w = static_cast<RealType>(_Kernel[n]);
        v = voxel_cast<AccType>(*in);
        acc += (v += w);
        sum += w;
      }
      ++i, ++in, --n;
    }
    // Output result of convolution
    if (_Normalize && !AreEqual(sum, 0., Epsilon())) acc /= sum;
    *out = voxel_cast<T2>(acc);
  }
};

// -----------------------------------------------------------------------------
/// Perform convolution of image foreground along y with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveForegroundInY : public VoxelFunction
{
  const BaseImage *_Image;     ///< Image defining foreground region
  const TKernel   *_Kernel;    ///< Convolution kernel
  int              _Radius;    ///< Radius of kernel
  bool             _Normalize; ///< Wether to divide by sum of kernel weights
  int              _Offset;    ///< Vector component offset
  int              _X;         ///< Number of voxels in x
  int              _Y;         ///< Number of voxels in y

  ConvolveForegroundInY(const BaseImage *image, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    _Image(image), _Kernel(kernel), _Radius((size - 1) / 2), _Normalize(norm),
    _Offset(l * image->NumberOfSpatialVoxels()), _X(image->X()), _Y(image->Y())
  {}

  template <class T1, class T2>
  void operator ()(int i, int j, int k, int l, const T1 *in, T2 *out) const
  {
    typedef typename voxel_info<T2>::RealType         AccType;
    typedef typename voxel_info<AccType>::ScalarType  RealType;
    // Convolve l-th vector component
    in += _Offset, out += _Offset;
    // Go to start of input and end of kernel
    j    -= _Radius;
    in   -= _Radius * _X;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (j < 0) n += j, in -= j * _X, j = 0;
    // Inside image domain
    AccType  v, acc = 0.;
    RealType w, sum = 0.;
    while (j < _Y && n >= 0) {
      if (_Image->IsForeground(i, j, k, l)) {
        w = static_cast<RealType>(_Kernel[n]);
        v = voxel_cast<AccType>(*in);
        acc += (v *= w);
        sum += w;
      }
      ++j, in += _X, --n;
    }
    // Output result of convolution
    if (_Normalize && !AreEqual(sum, 0., Epsilon())) acc /= sum;
    *out = voxel_cast<T2>(acc);
  }
};

// -----------------------------------------------------------------------------
/// Perform convolution of image foreground along z with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveForegroundInZ : public VoxelFunction
{
  const BaseImage *_Image;     ///< Image defining foreground region
  const TKernel   *_Kernel;    ///< Convolution kernel
  int              _Radius;    ///< Radius of kernel
  bool             _Normalize; ///< Wether to divide by sum of kernel weights
  int              _Offset;    ///< Vector component offset
  int              _XY;        ///< Number of voxels in slice (nx * ny)
  int              _Z;         ///< Number of voxels in z

  ConvolveForegroundInZ(const BaseImage *image, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    _Image(image), _Kernel(kernel), _Radius((size - 1) / 2), _Normalize(norm),
    _Offset(l * image->NumberOfSpatialVoxels()), _XY(image->X() * image->Y()), _Z(image->Z())
  {}

  template <class T1, class T2>
  void operator ()(int i, int j, int k, int l, const T1 *in, T2 *out) const
  {
    typedef typename voxel_info<T2>::RealType         AccType;
    typedef typename voxel_info<AccType>::ScalarType  RealType;
    // Convolve l-th vector component
    in += _Offset, out += _Offset;
    // Go to start of input and end of kernel
    k    -= _Radius;
    in   -= _Radius * _XY;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (k < 0) n += k, in -= k * _XY, k = 0;
    // Inside image domain
    AccType  v, acc = 0.;
    RealType w, sum = 0.;
    while (k < _Z && n >= 0) {
      if (_Image->IsForeground(i, j, k, l)) {
        w = static_cast<RealType>(_Kernel[n]);
        v = voxel_cast<AccType>(*in);
        acc += (v *= w);
        sum += w;
      }
      ++k, in += _XY, --n;
    }
    // Output result of convolution
    if (_Normalize && !AreEqual(sum, 0., Epsilon())) acc /= sum;
    *out = voxel_cast<T2>(acc);
  }
};

// -----------------------------------------------------------------------------
/// Perform convolution of image foreground along t with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveForegroundInT : public VoxelFunction
{
  const BaseImage *_Image;     ///< Image defining foreground region
  const TKernel   *_Kernel;    ///< Convolution kernel
  int              _Radius;    ///< Radius of kernel
  bool             _Normalize; ///< Wether to divide by sum of kernel weights
  int              _XYZ;       ///< Number of voxels in volume (nx * ny * nz)
  int              _T;         ///< Number of voxels in t

  ConvolveForegroundInT(const BaseImage *image, const TKernel *kernel, int size, bool norm = true)
  :
    _Image(image), _Kernel(kernel), _Radius((size - 1) / 2), _Normalize(norm),
    _XYZ(image->X() * image->Y() * image->Z()), _T(image->T())
  {}

  template <class T1, class T2>
  void operator ()(int i, int j, int k, int l, const T1 *in, T2 *out) const
  {
    typedef typename voxel_info<T2>::RealType         AccType;
    typedef typename voxel_info<AccType>::ScalarType  RealType;
    // Go to start of input and end of kernel
    l    -= _Radius;
    in   -= _Radius * _XYZ;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (l < 0) n += l, in -= l * _XYZ, l = 0;
    // Inside image domain
    AccType  v, acc = 0.;
    RealType w, sum = 0.;
    while (l < _T && n >= 0) {
      if (_Image->IsForeground(i, j, k, l)) {
        w = static_cast<RealType>(_Kernel[n]);
        v = voxel_cast<AccType>(*in);
        acc += (v *= w);
        sum += w;
      }
      ++l, in += _XYZ, --n;
    }
    // Output result of convolution
    if (_Normalize && !AreEqual(sum, 0., Epsilon())) acc /= sum;
    *out = voxel_cast<T2>(acc);
  }
};

// =============================================================================
// Convolve at boundary truncated voxel-wise weighted image with 1D kernel
// =============================================================================

// -----------------------------------------------------------------------------
/// Perform convolution of weighted image along x with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveWeightedImageInX : public VoxelFunction
{
  const TKernel *_Kernel;  ///< Convolution kernel
  int            _Radius;  ///< Radius of kernel
  int            _Offset1; ///< Image  component offset
  int            _Offset2; ///< Weight component offset
  int            _X;       ///< Number of voxels in x

  ConvolveWeightedImageInX(const BaseImage *image, const TKernel *kernel,
                           int size, int l1 = 0, int l2 = 0)
  :
    _Kernel(kernel), _Radius((size - 1) / 2),
    _Offset1(l1 * image->NumberOfSpatialVoxels()),
    _Offset2(l2 * image->NumberOfSpatialVoxels()),
    _X(image->X())
  {}

  template <class T1, class T2, class T3>
  void operator ()(int i, int, int, int, const T1 *in, const T2 *win, T3 *out) const
  {
    // Convolve l1-th vector component using l2-th weight
    in += _Offset1, win += _Offset2, out += _Offset1;
    // Go to start of input and end of kernel
    i    -= _Radius;
    in   -= _Radius;
    win  -= _Radius;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (i < 0) n += i, in -= i, win -= i, i = 0;
    // Inside image domain
    double w, acc = .0, sum = .0;
    while (i < _X && n >= 0) {
      w = static_cast<double>(*win) * static_cast<double>(_Kernel[n]);
      acc += w * static_cast<double>(*in);
      sum += w;
      ++i, ++in, ++win, --n;
    }
    // Output result of convolution
    acc = (sum > Epsilon() ? acc / sum : 0.);
    *out = static_cast<T3>(acc);
  }
};

// -----------------------------------------------------------------------------
/// Perform convolution of weighted image along y with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveWeightedImageInY : public VoxelFunction
{
  const TKernel *_Kernel;  ///< Convolution kernel
  int            _Radius;  ///< Radius of kernel
  int            _Offset1; ///< Image  component offset
  int            _Offset2; ///< Weight component offset
  int            _X;       ///< Number of voxels in x
  int            _Y;       ///< Number of voxels in y

  ConvolveWeightedImageInY(const BaseImage *image, const TKernel *kernel, int size, int l1 = 0, int l2 = 0)
  :
    _Kernel(kernel), _Radius((size - 1) / 2),
    _Offset1(l1 * image->NumberOfSpatialVoxels()),
    _Offset2(l2 * image->NumberOfSpatialVoxels()),
    _X(image->X()), _Y(image->Y())
  {}

  template <class T1, class T2, class T3>
  void operator ()(int, int j, int, int, const T1 *in, const T2 *win, T3 *out) const
  {
    // Convolve l1-th vector component using l2-th weight
    in += _Offset1, win += _Offset2, out += _Offset1;
    // Go to start of input and end of kernel
    j    -= _Radius;
    in   -= _Radius * _X;
    win  -= _Radius * _X;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (j < 0) n += j, in -= j * _X, win -= j * _X, j = 0;
    // Inside image domain
    double w, acc = .0, sum = .0;
    while (j < _Y && n >= 0) {
      w = static_cast<double>(_Kernel[n]) * static_cast<double>(*win);
      acc += w * static_cast<double>(*in);
      sum += w;
      ++j, in += _X, win += _X, --n;
    }
    // Output result of convolution
    acc = (sum > Epsilon() ? acc / sum : 0.);
    *out = static_cast<T3>(acc);
  }
};

// -----------------------------------------------------------------------------
/// Perform convolution of weighted image along z with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveWeightedImageInZ : public VoxelFunction
{
  const TKernel *_Kernel;  ///< Convolution kernel
  int            _Radius;  ///< Radius of kernel
  int            _Offset1; ///< Image  component offset
  int            _Offset2; ///< Weight component offset
  int            _XY;      ///< Number of voxels in slice (nx * ny)
  int            _Z;       ///< Number of voxels in z

  ConvolveWeightedImageInZ(const BaseImage *image, const TKernel *kernel, int size, int l1 = 0, int l2 = 0)
  :
    _Kernel(kernel), _Radius((size - 1) / 2),
    _Offset1(l1 * image->NumberOfSpatialVoxels()),
    _Offset2(l2 * image->NumberOfSpatialVoxels()),
    _XY(image->X() * image->Y()), _Z(image->Z())
  {}

  template <class T1, class T2, class T3>
  void operator ()(int, int, int k, int, const T1 *in, const T2 *win, T3 *out) const
  {
    // Convolve l1-th vector component using l2-th weight
    in += _Offset1, win += _Offset2, out += _Offset1;
    // Go to start of input and end of kernel
    k    -= _Radius;
    in   -= _Radius * _XY;
    win  -= _Radius * _XY;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (k < 0) n += k, in -= k * _XY, win -= k * _XY, k = 0;
    // Inside image domain
    double w, acc = .0, sum = .0;
    while (k < _Z && n >= 0) {
      w = static_cast<double>(_Kernel[n]) * static_cast<double>(*win);
      acc += w * static_cast<double>(*in);
      sum += w;
      ++k, in += _XY, win += _XY, --n;
    }
    // Output result of convolution
    acc = (sum > Epsilon() ? acc / sum : 0.);
    *out = static_cast<T3>(acc);
  }
};

// -----------------------------------------------------------------------------
/// Perform convolution of weighted image along t with truncation of kernel at boundary
template <class TKernel = double>
struct ConvolveWeightedImageInT : public VoxelFunction
{
  const TKernel *_Kernel; ///< Convolution kernel
  int            _Radius; ///< Radius of kernel
  int            _XYZ;    ///< Number of voxels in volume (nx * ny * nz)
  int            _T;      ///< Number of voxels in t

  ConvolveWeightedImageInT(const BaseImage *image, const TKernel *kernel, int size)
  :
    _Kernel(kernel), _Radius((size - 1) / 2),
    _XYZ(image->X() * image->Y() * image->Z()), _T(image->T())
  {}

  template <class T1, class T2, class T3>
  void operator ()(int, int, int, int l, const T1 *in, const T2 *win, T3 *out) const
  {
    // Go to start of input and end of kernel
    l    -= _Radius;
    in   -= _Radius * _XYZ;
    win  -= _Radius * _XYZ;
    int n = _Radius + _Radius/* + 1 - 1*/;
    // Outside left boundary
    if (l < 0) n += l, in -= l * _XYZ, win -= l * _XYZ, l = 0;
    // Inside image domain
    double w, acc = .0, sum = .0;
    while (l < _T && n >= 0) {
      w = static_cast<double>(_Kernel[n]) * static_cast<double>(*win);
      acc += w * static_cast<double>(*in);
      sum += w;
      ++l, in += _XYZ, win += _XYZ, --n;
    }
    // Output result of convolution
    acc = (sum > Epsilon() ? acc / sum : 0.);
    *out = static_cast<T3>(acc);
  }
};

// =============================================================================
// Convolve at truncated image foreground with 1D kernel
// =============================================================================

// -----------------------------------------------------------------------------
/// Base class of 1D convolution functions which truncate the kernel at the foreground boundary
template <class TKernel = double>
class TruncatedForegroundConvolution1D : public VoxelFunction
{
protected:

  // ---------------------------------------------------------------------------
  /// Constructor
  TruncatedForegroundConvolution1D(const BaseImage *image, const TKernel *kernel, int size, bool norm = true)
  :
    _Kernel(kernel), _Size(size), _Radius((size - 1) / 2), _Normalize(norm)
  {
    if (image->HasBackgroundValue()) {
      _Background = image->GetBackgroundValueAsDouble();
    } else {
      _Background = NaN;
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply kernel initially at center voxel
  template <class T>
  bool ConvolveCenterVoxel(const T *in, double &acc, double &sum) const
  {
    if (AreEqualOrNaN(static_cast<double>(*in), _Background, Epsilon())) {
      acc = _Background;
      sum = .0;
      return false;
    } else {
      acc = static_cast<double>(_Kernel[_Radius]) * static_cast<double>(*in);
      sum = static_cast<double>(_Kernel[_Radius]);
      return true;
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply kernel to left neighbors of given voxel
  template <class T>
  void ConvolveLeftNeighbors(int i, int n, const T *in, int s, double &acc, double &sum) const
  {
    int di = -1;
    for (int k = _Radius + 1; k < _Size; ++k) {
      // Move voxel index/pointer
      i  += di;
      in -= s;
      // Stop if outside image/foreground
      if (i < 0 || i >= n || AreEqualOrNaN(static_cast<double>(*in), _Background, Epsilon())) break;
      // Multiply with kernel weight
      acc += static_cast<double>(_Kernel[k]) * static_cast<double>(*in);
      sum += static_cast<double>(_Kernel[k]);
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply kernel to right neighbors of given voxel
  template <class T>
  void ConvolveRightNeighbors(int i, int n, const T *in, int s, double &acc, double &sum) const
  {
    int di = +1;
    for (int k = _Radius - 1; k >= 0; --k) {
      // Move voxel index/pointer
      i  += di;
      in += s;
      // Stop if outside image/foreground
      if (i < 0 || i >= n || AreEqualOrNaN(static_cast<double>(*in), _Background, Epsilon())) break;
      // Multiply with kernel weight
      acc += static_cast<double>(_Kernel[k]) * static_cast<double>(*in);
      sum += static_cast<double>(_Kernel[k]);
    }
  }

  // ---------------------------------------------------------------------------
  template <class T>
  void Put(T *out, double acc, double sum) const
  {
    if (_Normalize && !AreEqual(sum, 0., Epsilon())) acc /= sum;
    (*out) = static_cast<T>(acc);
  }

  // ---------------------------------------------------------------------------
  // Data members
protected:

  const TKernel *_Kernel;     ///< Convolution kernel
  int            _Size;       ///< Size of kernel = 2 * _Radius + 1
  int            _Radius;     ///< Radius of kernel
  bool           _Normalize;  ///< Whether to normalize by sum of overlapping kernel weights
  double         _Background; ///< Background value (padding)
};


// -----------------------------------------------------------------------------
/// Compute convolution of voxel with kernel in x dimension
template <class TKernel = double>
struct ConvolveTruncatedForegroundInX : public TruncatedForegroundConvolution1D<TKernel>
{
  ConvolveTruncatedForegroundInX(const BaseImage *image, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _Offset(l * image->NumberOfSpatialVoxels()), _X(image->X())
  {}

  ConvolveTruncatedForegroundInX(const BaseImage *image, double bg, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _Offset(l * image->NumberOfSpatialVoxels()), _X(image->X())
  {
    this->_Background = bg;
  }

  ConvolveTruncatedForegroundInX(const BaseImage *image, const GenericImage<TKernel> *kernel, bool norm = true, int l = 0)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _Offset(l * image->NumberOfSpatialVoxels()), _X(image->X())
  {}

  template <class T1, class T2>
  void operator ()(int i, int, int, int, const T1 *in, T2 *out) const
  {
    in += _Offset, out += _Offset;
    double acc, sum;
    if (this->ConvolveCenterVoxel(in, acc, sum)) {
      this->ConvolveLeftNeighbors (i, _X, in, 1, acc, sum);
      this->ConvolveRightNeighbors(i, _X, in, 1, acc, sum);
    }
    this->Put(out, acc, sum);
  }

  int _Offset; ///< Vector component offset
  int _X;      ///< Number of voxels in x dimension (nx)
};

// -----------------------------------------------------------------------------
/// Compute convolution for given voxel along y dimension
template <class TKernel = double>
struct ConvolveTruncatedForegroundInY : public TruncatedForegroundConvolution1D<TKernel>
{
  ConvolveTruncatedForegroundInY(const BaseImage *image, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _Offset(l * image->NumberOfSpatialVoxels()),
    _X(image->X()),
    _Y(image->Y())
  {}

  ConvolveTruncatedForegroundInY(const BaseImage *image, double bg, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _Offset(l * image->NumberOfSpatialVoxels()),
    _X(image->X()),
    _Y(image->Y())
  {
    this->_Background = bg;
  }

  ConvolveTruncatedForegroundInY(const BaseImage *image, const GenericImage<TKernel> *kernel, bool norm = true, int l = 0)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _Offset(l * image->NumberOfSpatialVoxels()),
    _X(image->X()),
    _Y(image->Y())
  {}

  template <class T1, class T2>
  void operator ()(int, int j, int, int, const T1 *in, T2 *out) const
  {
    in += _Offset, out += _Offset;
    double acc, sum;
    if (this->ConvolveCenterVoxel(in, acc, sum)) {
      this->ConvolveLeftNeighbors (j, _Y, in, _X, acc, sum);
      this->ConvolveRightNeighbors(j, _Y, in, _X, acc, sum);
    }
    this->Put(out, acc, sum);
  }

  int _Offset; ///< Vector component offset
  int _X;      ///< Number of voxels in x dimension (nx)
  int _Y;      ///< Number of voxels in y dimension (ny)
};

// -----------------------------------------------------------------------------
/// Compute convolution for given voxel along z dimension
template <class TKernel = double>
struct ConvolveTruncatedForegroundInZ : public TruncatedForegroundConvolution1D<TKernel>
{
  ConvolveTruncatedForegroundInZ(const BaseImage *image, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _Offset(l * image->NumberOfSpatialVoxels()),
    _XY(image->X() * image->Y()),
    _Z (image->Z())
  {}

  ConvolveTruncatedForegroundInZ(const BaseImage *image, double bg, const TKernel *kernel, int size, bool norm = true, int l = 0)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _Offset(l * image->NumberOfSpatialVoxels()),
    _XY(image->X() * image->Y()),
    _Z (image->Z())
  {
    this->_Background = bg;
  }

  ConvolveTruncatedForegroundInZ(const BaseImage *image, const GenericImage<TKernel> *kernel, bool norm = true, int l = 0)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _Offset(l * image->NumberOfSpatialVoxels()),
    _XY(image->X() * image->Y()),
    _Z (image->Z())
  {}

  template <class T1, class T2>
  void operator ()(int, int, int k, int, const T1 *in, T2 *out) const
  {
    in += _Offset, out += _Offset;
    double acc, sum;
    if (this->ConvolveCenterVoxel(in, acc, sum)) {
      this->ConvolveLeftNeighbors (k, _Z, in, _XY, acc, sum);
      this->ConvolveRightNeighbors(k, _Z, in, _XY, acc, sum);
    }
    this->Put(out, acc, sum);
  }

  int _Offset; ///< Vector component offset
  int _XY;     ///< Number of voxels in slice (nx * ny)
  int _Z;      ///< Number of voxels in z dimension (nz)
};

// -----------------------------------------------------------------------------
/// Compute convolution for given voxel along t dimension
template <class TKernel = double>
struct ConvolveTruncatedForegroundInT : public TruncatedForegroundConvolution1D<TKernel>
{
  ConvolveTruncatedForegroundInT(const BaseImage *image, const TKernel *kernel, int size, bool norm = true)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _XYZ(image->X() * image->Y() * image->Z()),
    _T  (image->T())
  {}

  ConvolveTruncatedForegroundInT(const BaseImage *image, double bg, const TKernel *kernel, int size, bool norm = true)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _XYZ(image->X() * image->Y() * image->Z()),
    _T  (image->T())
  {
    this->_Background = bg;
  }

  ConvolveTruncatedForegroundInT(const BaseImage *image, const GenericImage<TKernel> *kernel, bool norm = true)
  :
    TruncatedForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _XYZ(image->X() * image->Y() * image->Z()),
    _T  (image->T())
  {}

  template <class T1, class T2>
  void operator ()(int, int, int, int t, const T1 *in, T2 *out) const
  {
    double acc, sum;
    if (this->ConvolveCenterVoxel(in, acc, sum)) {
      this->ConvolveLeftNeighbors (t, _T, in, _XYZ, acc, sum);
      this->ConvolveRightNeighbors(t, _T, in, _XYZ, acc, sum);
    }
    this->Put(out, acc, sum);
  }

  int _XYZ; ///< Number of voxels in frame (nx * ny * nz)
  int _T;   ///< Number of voxels in t dimension (nt)
};

// =============================================================================
// Convolve at boundary mirrored image foreground with 1D kernel
// =============================================================================

// -----------------------------------------------------------------------------
/// Base class of 1D convolution functions which mirror the foreground into background
template <class TKernel = double>
class MirroredForegroundConvolution1D : public VoxelFunction
{
protected:

  // ---------------------------------------------------------------------------
  /// Constructor
  MirroredForegroundConvolution1D(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
  :
    _Kernel(kernel), _Size(size), _Radius((size - 1) / 2), _Norm(norm)
  {
    if (image->HasBackgroundValue()) {
      _Background = image->GetBackgroundValueAsDouble();
    } else {
      _Background = NaN;
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply kernel initially at center voxel
  template <class T>
  bool ConvolveCenterVoxel(const T *in, double &acc) const
  {
    if (AreEqualOrNaN(static_cast<double>(*in), _Background, Epsilon())) {
      acc = _Background / _Norm;
      return false;
    } else {
      acc = static_cast<double>(*in) * static_cast<double>(_Kernel[_Radius]);
      return true;
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply kernel to left neighbors of given voxel
  template <class T>
  void ConvolveLeftNeighbors(int i, int n, const T *in, int s, double &acc) const
  {
    int di = -1;
    for (int k = _Radius + 1; k < _Size; ++k) {
      // Move voxel index/pointer
      i  += di;
      in -= s;
      // Reverse if outside image/foreground
      if (i < 0 || i >= n || AreEqualOrNaN(static_cast<double>(*in), _Background, Epsilon())) {
        di  = -di;
        s   = -s;
        i  += di;
        in -= s;
      }
      // Multiply with kernel weight
      acc += static_cast<double>(*in) * static_cast<double>(_Kernel[k]);
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply kernel to right neighbors of given voxel
  template <class T>
  void ConvolveRightNeighbors(int i, int n, const T *in, int s, double &acc) const
  {
    int di = +1;
    for (int k = _Radius - 1; k >= 0; --k) {
      // Move voxel index/pointer
      i  += di;
      in += s;
      // Reverse if outside image/foreground
      if (i < 0 || i >= n || AreEqualOrNaN(static_cast<double>(*in), _Background, Epsilon())) {
        di  = -di;
        s   = -s;
        i  += di;
        in +=  s;
      }
      // Multiply with kernel weight
      acc += static_cast<double>(*in) * static_cast<double>(_Kernel[k]);
    }
  }

  // ---------------------------------------------------------------------------
  template <class T>
  void Put(T *out, double acc) const
  {
    (*out) = static_cast<T>(_Norm * acc);
  }

  // ---------------------------------------------------------------------------
  // Data members
private:

  const TKernel *_Kernel;     ///< Convolution kernel
  int            _Size;       ///< Size of kernel = 2 * _Radius + 1
  int            _Radius;     ///< Radius of kernel
  double         _Norm;       ///< Normalization factor
  double         _Background; ///< Background value (padding)
};


// -----------------------------------------------------------------------------
/// Compute convolution of voxel with kernel in x dimension
template <class TKernel = double>
struct ConvolveMirroredForegroundInX : public MirroredForegroundConvolution1D<TKernel>
{
  ConvolveMirroredForegroundInX(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
    : MirroredForegroundConvolution1D<TKernel>(image, kernel, size, norm), _X(image->X()) {}

  ConvolveMirroredForegroundInX(const BaseImage *image, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    MirroredForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _X(image->X())
  {}

  template <class T1, class T2>
  void operator ()(int i, int, int, int, const T1 *in, T2 *out) const
  {
    double acc;
    if (this->ConvolveCenterVoxel(in, acc)) {
      this->ConvolveLeftNeighbors (i, _X, in, 1, acc);
      this->ConvolveRightNeighbors(i, _X, in, 1, acc);
    }
    this->Put(out, acc);
  }

  int _X; ///< Number of voxels in x dimension (nx)
};

// -----------------------------------------------------------------------------
/// Compute convolution for given voxel along y dimension
template <class TKernel = double>
struct ConvolveMirroredForegroundInY : public MirroredForegroundConvolution1D<TKernel>
{
  ConvolveMirroredForegroundInY(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
  :
    MirroredForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _X(image->X()),
    _Y(image->Y())
  {}

  ConvolveMirroredForegroundInY(const BaseImage *image, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    MirroredForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _X(image->X()),
    _Y(image->Y())
  {}

  template <class T1, class T2>
  void operator ()(int, int j, int, int, const T1 *in, T2 *out) const
  {
    double acc;
    if (this->ConvolveCenterVoxel(in, acc)) {
      this->ConvolveLeftNeighbors (j, _Y, in, _X, acc);
      this->ConvolveRightNeighbors(j, _Y, in, _X, acc);
    }
    this->Put(out, acc);
  }

  int _X; ///< Number of voxels in x dimension (nx)
  int _Y; ///< Number of voxels in y dimension (ny)
};

// -----------------------------------------------------------------------------
/// Compute convolution for given voxel along z dimension
template <class TKernel = double>
struct ConvolveMirroredForegroundInZ : public MirroredForegroundConvolution1D<TKernel>
{
  ConvolveMirroredForegroundInZ(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
  :
    MirroredForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _XY(image->X() * image->Y()),
    _Z (image->Z())
  {}

  ConvolveMirroredForegroundInZ(const BaseImage *image, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    MirroredForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _XY(image->X() * image->Y()),
    _Z (image->Z())
  {}

  template <class T1, class T2>
  void operator ()(int, int, int k, int, const T1 *in, T2 *out) const
  {
    double acc;
    if (this->ConvolveCenterVoxel(in, acc)) {
      this->ConvolveLeftNeighbors (k, _Z, in, _XY, acc);
      this->ConvolveRightNeighbors(k, _Z, in, _XY, acc);
    }
    this->Put(out, acc);
  }

  int _XY; ///< Number of voxels in slice (nx * ny)
  int _Z;  ///< Number of voxels in z dimension (nz)
};

// -----------------------------------------------------------------------------
/// Compute convolution for given voxel along t dimension
template <class TKernel = double>
struct ConvolveMirroredForegroundInT : public MirroredForegroundConvolution1D<TKernel>
{
  ConvolveMirroredForegroundInT(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
  :
    MirroredForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _XYZ(image->X() * image->Y() * image->Z()),
    _T  (image->T())
  {}

  ConvolveMirroredForegroundInT(const BaseImage *image, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    MirroredForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _XYZ(image->X() * image->Y() * image->Z()),
    _T  (image->T())
  {}

  template <class T1, class T2>
  void operator ()(int, int, int, int t, const T1 *in, T2 *out) const
  {
    double acc;
    if (this->ConvolveCenterVoxel(in, acc)) {
      this->ConvolveLeftNeighbors (t, _T, in, _XYZ, acc);
      this->ConvolveRightNeighbors(t, _T, in, _XYZ, acc);
    }
    this->Put(out, acc);
  }

  int _XYZ; ///< Number of voxels in frame (nx * ny * nz)
  int _T;   ///< Number of voxels in t dimension (nt)
};

// =============================================================================
// Convolve at boundary extended image foreground with 1D kernel
// =============================================================================

// -----------------------------------------------------------------------------
/// Base class of 1D convolution functions which extends the foreground into background
template <class TKernel = double>
class ExtendedForegroundConvolution1D : public VoxelFunction
{
protected:

  // ---------------------------------------------------------------------------
  /// Constructor
  ExtendedForegroundConvolution1D(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
  :
    _Kernel(kernel), _Size(size), _Radius((size - 1) / 2), _Norm(norm)
  {
    if (image->HasBackgroundValue()) {
      _Background = image->GetBackgroundValueAsDouble();
    } else {
      _Background = NaN;
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply kernel initially at center voxel
  template <class T>
  bool ConvolveCenterVoxel(const T *in, double &acc) const
  {
    if (AreEqualOrNaN(static_cast<double>(*in), _Background, Epsilon())) {
      acc = _Background / _Norm;
      return false;
    } else {
      acc = (*in) * _Kernel[_Radius];
      return true;
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply kernel to left neighbors of given voxel
  template <class T>
  void ConvolveLeftNeighbors(int i, int n, const T *in, int s, double &acc) const
  {
    int di = 1;
    for (int k = _Radius + 1; k < _Size; ++k) {
      i  -= di;
      in -= s;
      if (i < 0 || AreEqualOrNaN(static_cast<double>(*in), _Background, Epsilon())) {
        i  += di;
        in += s;
        di = s = 0;
      }
      // Multiply with kernel weight
      acc += static_cast<double>(*in) * static_cast<double>(_Kernel[k]);
    }
  }

  // ---------------------------------------------------------------------------
  /// Apply kernel to right neighbors of given voxel
  template <class T>
  void ConvolveRightNeighbors(int i, int n, const T *in, int s, double &acc) const
  {
    int di = 1;
    for (int k = _Radius - 1; k >= 0; --k) {
      i  += di;
      in += s;
      if (n <= i || AreEqualOrNaN(static_cast<double>(*in), _Background, Epsilon())) {
        i  -= di;
        in -= s;
        di = s = 0;
      }
      // Multiply with kernel weight
      acc += static_cast<double>(*in) * static_cast<double>(_Kernel[k]);
    }
  }

  // ---------------------------------------------------------------------------
  template <class T>
  void Put(T *out, double acc) const
  {
    (*out) = static_cast<T>(_Norm * acc);
  }

  // ---------------------------------------------------------------------------
  // Data members
private:

  const TKernel *_Kernel;     ///< Convolution kernel
  int            _Size;       ///< Size of kernel = 2 * _Radius + 1
  int            _Radius;     ///< Radius of kernel
  double         _Norm;       ///< Normalization factor
  double         _Background; ///< Background value (padding)
};


// -----------------------------------------------------------------------------
/// Compute convolution of voxel with kernel in x dimension
template <class TKernel = double>
struct ConvolveExtendedForegroundInX : public ExtendedForegroundConvolution1D<TKernel>
{
  ConvolveExtendedForegroundInX(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
    : ExtendedForegroundConvolution1D<TKernel>(image, kernel, size, norm), _X(image->X()) {}

  ConvolveExtendedForegroundInX(const BaseImage *image, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ExtendedForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _X(image->X())
  {}

  template <class T1, class T2>
  void operator ()(int i, int, int, int, const T1 *in, T2 *out) const
  {
    double acc;
    if (this->ConvolveCenterVoxel(in, acc)) {
      this->ConvolveLeftNeighbors (i, _X, in, 1, acc);
      this->ConvolveRightNeighbors(i, _X, in, 1, acc);
    }
    this->Put(out, acc);
  }

  int _X; ///< Number of voxels in x dimension (nx)
};

// -----------------------------------------------------------------------------
/// Compute convolution for given voxel along y dimension
template <class TKernel = double>
struct ConvolveExtendedForegroundInY : public ExtendedForegroundConvolution1D<TKernel>
{
  ConvolveExtendedForegroundInY(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
  :
    ExtendedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _X(image->X()),
    _Y(image->Y())
  {}

  ConvolveExtendedForegroundInY(const BaseImage *image, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ExtendedForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _X(image->X()),
    _Y(image->Y())
  {}

  template <class T1, class T2>
  void operator ()(int, int j, int, int, const T1 *in, T2 *out) const
  {
    double acc;
    if (this->ConvolveCenterVoxel(in, acc)) {
      this->ConvolveLeftNeighbors (j, _Y, in, _X, acc);
      this->ConvolveRightNeighbors(j, _Y, in, _X, acc);
    }
    this->Put(out, acc);
  }

  int _X; ///< Number of voxels in x dimension (nx)
  int _Y; ///< Number of voxels in y dimension (ny)
};

// -----------------------------------------------------------------------------
/// Compute convolution for given voxel along z dimension
template <class TKernel = double>
struct ConvolveExtendedForegroundInZ : public ExtendedForegroundConvolution1D<TKernel>
{
  ConvolveExtendedForegroundInZ(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
  :
    ExtendedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _XY(image->X() * image->Y()),
    _Z (image->Z())
  {}

  ConvolveExtendedForegroundInZ(const BaseImage *image, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ExtendedForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _XY(image->X() * image->Y()),
    _Z (image->Z())
  {}

  template <class T1, class T2>
  void operator ()(int, int, int k, int, const T1 *in, T2 *out) const
  {
    double acc;
    if (this->ConvolveCenterVoxel(in, acc)) {
      this->ConvolveLeftNeighbors (k, _Z, in, _XY, acc);
      this->ConvolveRightNeighbors(k, _Z, in, _XY, acc);
    }
    this->Put(out, acc);
  }

  int _XY; ///< Number of voxels in slice (nx * ny)
  int _Z;  ///< Number of voxels in z dimension (nz)
};

// -----------------------------------------------------------------------------
/// Compute convolution for given voxel along t dimension
template <class TKernel = double>
struct ConvolveExtendedForegroundInT : public ExtendedForegroundConvolution1D<TKernel>
{
  ConvolveExtendedForegroundInT(const BaseImage *image, const TKernel *kernel, int size, double norm = 1.0)
  :
    ExtendedForegroundConvolution1D<TKernel>(image, kernel, size, norm),
    _XYZ(image->X() * image->Y() * image->Z()),
    _T  (image->T())
  {}

  ConvolveExtendedForegroundInT(const BaseImage *image, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ExtendedForegroundConvolution1D<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _XYZ(image->X() * image->Y() * image->Z()),
    _T  (image->T())
  {}

  template <class T1, class T2>
  void operator ()(int, int, int, int t, const T1 *in, T2 *out) const
  {
    double acc;
    if (this->ConvolveCenterVoxel(in, acc)) {
      this->ConvolveLeftNeighbors (t, _T, in, _XYZ, acc);
      this->ConvolveRightNeighbors(t, _T, in, _XYZ, acc);
    }
    this->Put(out, acc);
  }

  int _XYZ; ///< Number of voxels in frame (nx * ny * nz)
  int _T;   ///< Number of voxels in t dimension (nt)
};

// =============================================================================
// Smooth and downsample mirrored foreground in one step
// =============================================================================

// -----------------------------------------------------------------------------
/// Downsample image using convolution of original voxels with kernel in x dimension
template <class TVoxel, class TKernel = RealPixel>
struct DownsampleConvolvedMirroredForegroundInX : public ConvolveMirroredForegroundInX<TKernel>
{
  DownsampleConvolvedMirroredForegroundInX(const GenericImage<TVoxel> *image, int m, const TKernel *kernel, int size, double norm = 1.0)
  :
    ConvolveMirroredForegroundInX<TKernel>(image, kernel, size, norm),
    _Input(image), _Offset((image->X() % m + 1) / 2), _Factor(m)
  {}

  DownsampleConvolvedMirroredForegroundInX(const GenericImage<TVoxel> *image, int m, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ConvolveMirroredForegroundInX<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _Input(image), _Offset((image->X() % m + 1) / 2), _Factor(m)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int l, T *out) const
  {
    i = _Offset + i * _Factor;
    ConvolveMirroredForegroundInX<TKernel>::operator ()(i, j, k, l, _Input->Data(i, j, k, l), out);
  }

  const GenericImage<TVoxel> *_Input;  ///< Image to downsample
  int                         _Offset; ///< Offset of first voxel
  int                         _Factor; ///< Downsampling factor
};

// -----------------------------------------------------------------------------
/// Downsample image using convolution of original voxels with kernel in y dimension
template <class TVoxel, class TKernel = RealPixel>
struct DownsampleConvolvedMirroredForegroundInY : public ConvolveMirroredForegroundInY<TKernel>
{
  DownsampleConvolvedMirroredForegroundInY(const GenericImage<TVoxel> *image, int m, const TKernel *kernel, int size, double norm = 1.0)
  :
    ConvolveMirroredForegroundInY<TKernel>(image, kernel, size, norm),
    _Input(image), _Offset((image->Y() % m + 1) / 2), _Factor(m)
  {}

  DownsampleConvolvedMirroredForegroundInY(const GenericImage<TVoxel> *image, int m, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ConvolveMirroredForegroundInY<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _Input(image), _Offset((image->Y() % m + 1) / 2), _Factor(m)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int l, T *out) const
  {
    j = _Offset + j * _Factor;
    ConvolveMirroredForegroundInY<TKernel>::operator ()(i, j, k, l, _Input->Data(i, j, k, l), out);
  }

  const GenericImage<TVoxel> *_Input;  ///< Image to downsample
  int                         _Offset; ///< Offset of first voxel
  int                         _Factor; ///< Downsampling factor
};

// -----------------------------------------------------------------------------
/// Downsample image using convolution of original voxels with kernel in z dimension
template <class TVoxel, class TKernel = RealPixel>
struct DownsampleConvolvedMirroredForegroundInZ : public ConvolveMirroredForegroundInZ<TKernel>
{
  DownsampleConvolvedMirroredForegroundInZ(const GenericImage<TVoxel> *image, int m, const TKernel *kernel, int size, double norm = 1.0)
  :
    ConvolveMirroredForegroundInZ<TKernel>(image, kernel, size, norm),
    _Input(image), _Offset((image->Z() % m + 1) / 2), _Factor(m)
  {}

  DownsampleConvolvedMirroredForegroundInZ(const GenericImage<TVoxel> *image, int m, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ConvolveMirroredForegroundInZ<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _Input(image), _Offset((image->Z() % m + 1) / 2), _Factor(m)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int l, T *out) const
  {
    k = _Offset + k * _Factor;
    ConvolveMirroredForegroundInZ<TKernel>::operator ()(i, j, k, l, _Input->Data(i, j, k, l), out);
  }

  const GenericImage<TVoxel> *_Input;  ///< Image to downsample
  int                         _Offset; ///< Offset of first voxel
  int                         _Factor; ///< Downsampling factor
};

// =============================================================================
// Smooth and downsample extended foreground in one step
// =============================================================================

// -----------------------------------------------------------------------------
/// Downsample image using convolution of original voxels with kernel in x dimension
template <class TVoxel, class TKernel = RealPixel>
struct DownsampleConvolvedExtendedForegroundInX : public ConvolveExtendedForegroundInX<TKernel>
{
  DownsampleConvolvedExtendedForegroundInX(const GenericImage<TVoxel> *image, int m, const TKernel *kernel, int size, double norm = 1.0)
  :
    ConvolveExtendedForegroundInX<TKernel>(image, kernel, size, norm),
    _Input(image), _Offset((image->X() % m + 1) / 2), _Factor(m)
  {}

  DownsampleConvolvedExtendedForegroundInX(const GenericImage<TVoxel> *image, int m, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ConvolveExtendedForegroundInX<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _Input(image), _Offset((image->X() % m + 1) / 2), _Factor(m)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int l, T *out) const
  {
    i = _Offset + i * _Factor;
    ConvolveExtendedForegroundInX<TKernel>::operator ()(i, j, k, l, _Input->Data(i, j, k, l), out);
  }

  const GenericImage<TVoxel> *_Input;  ///< Image to downsample
  int                         _Offset; ///< Offset of first voxel
  int                         _Factor; ///< Downsampling factor
};

// -----------------------------------------------------------------------------
/// Downsample image using convolution of original voxels with kernel in y dimension
template <class TVoxel, class TKernel = RealPixel>
struct DownsampleConvolvedExtendedForegroundInY : public ConvolveExtendedForegroundInY<TKernel>
{
  DownsampleConvolvedExtendedForegroundInY(const GenericImage<TVoxel> *image, int m, const TKernel *kernel, int size, double norm = 1.0)
  :
    ConvolveExtendedForegroundInY<TKernel>(image, kernel, size, norm),
    _Input(image), _Offset((image->Y() % m + 1) / 2), _Factor(m)
  {}

  DownsampleConvolvedExtendedForegroundInY(const GenericImage<TVoxel> *image, int m, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ConvolveExtendedForegroundInY<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _Input(image), _Offset((image->Y() % m + 1) / 2), _Factor(m)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int l, T *out) const
  {
    j = _Offset + j * _Factor;
    ConvolveExtendedForegroundInY<TKernel>::operator ()(i, j, k, l, _Input->Data(i, j, k, l), out);
  }

  const GenericImage<TVoxel> *_Input;  ///< Image to downsample
  int                         _Offset; ///< Offset of first voxel
  int                         _Factor; ///< Downsampling factor
};

// -----------------------------------------------------------------------------
/// Downsample image using convolution of original voxels with kernel in z dimension
template <class TVoxel, class TKernel = RealPixel>
struct DownsampleConvolvedExtendedForegroundInZ : public ConvolveExtendedForegroundInZ<TKernel>
{
  DownsampleConvolvedExtendedForegroundInZ(const GenericImage<TVoxel> *image, int m, const TKernel *kernel, int size, double norm = 1.0)
  :
    ConvolveExtendedForegroundInZ<TKernel>(image, kernel, size, norm),
    _Input(image), _Offset((image->Z() % m + 1) / 2), _Factor(m)
  {}

  DownsampleConvolvedExtendedForegroundInZ(const GenericImage<TVoxel> *image, int m, const GenericImage<TKernel> *kernel, double norm = 1.0)
  :
    ConvolveExtendedForegroundInZ<TKernel>(image, kernel->Data(), kernel->NumberOfVoxels(), norm),
    _Input(image), _Offset((image->Z() % m) / 2), _Factor(m)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int l, T *out) const
  {
    k = _Offset + k * _Factor;
    ConvolveExtendedForegroundInZ<TKernel>::operator ()(i, j, k, l, _Input->Data(i, j, k, l), out);
  }

  const GenericImage<TVoxel> *_Input;  ///< Image to downsample
  int                         _Offset; ///< Offset of first voxel
  int                         _Factor; ///< Downsampling factor
};


} } // namespace mirtk::VoxelConvolution

#endif // MIRTK_ConvolutionFunction_H
