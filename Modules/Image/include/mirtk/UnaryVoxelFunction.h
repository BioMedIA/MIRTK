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

#ifndef MIRTK_UnaryVoxelFunction_H
#define MIRTK_UnaryVoxelFunction_H

#include "mirtk/VoxelFunction.h"

#include "mirtk/Assert.h"
#include "mirtk/Voxel.h"
#include "mirtk/VoxelCast.h"
#include "mirtk/BaseImage.h"
#include "mirtk/GenericImage.h"
#include "mirtk/Math.h"


namespace mirtk {


// Forward declaration
class InterpolateImageFunction;


/**
 * These basic unary voxel functions can be used as VoxelFunc template parameter
 * of the unary ForEachVoxel function templates as follows:
 *
 * \code
 * GreyImage image(attr);
 * // Clamp intensities such that they are in the range [0, 255]
 * UnaryVoxelFunction::Clamp clamp(0, 255);
 * ForEachVoxel(image, clamp);
 * // Determine intensity range of grayscale image
 * UnaryVoxelFunction::GetMinMax minmax;
 * ParallelForEachVoxel(image, minmax);
 * double min = minmax.GetMinAsDouble();
 * double max = minmax.GetMaxAsDouble();
 * \endcode
 */
namespace UnaryVoxelFunction {


// =============================================================================
// Minimum/maximum
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * Finds the minimum intensity in an image
 */
struct GetMin : public VoxelReduction
{
  void Reset() { _Min = voxel_limits<double>::max(); }

  GetMin() { Reset(); }
  GetMin(const GetMin &o) : _Min(o._Min) {}

  void split(const GetMin &rhs) { _Min = rhs._Min; }
  void join (const GetMin &rhs) { if (rhs._Min < _Min) _Min = rhs._Min; }

  template <class TImage, class T>
  void operator ()(const TImage&, int, const T *p)
  {
    if (static_cast<double>(*p) < _Min) _Min = static_cast<double>(*p);
  }

  template <class T>
  void operator ()(int, int, int, int, const T *p)
  {
    if (static_cast<double>(*p) < _Min) _Min = static_cast<double>(*p);
  }

  double GetMinAsDouble() const
  {
    if (_Min == voxel_limits<double>::max()) {
      return numeric_limits<double>::quiet_NaN();
    }
    return _Min;
  }

protected:
  double _Min;
};

// -----------------------------------------------------------------------------
/**
 * Finds the maximum intensity in an image
 */
struct GetMax : public VoxelReduction
{
  void Reset() { _Max = voxel_limits<double>::min(); }

  GetMax() { Reset(); }
  GetMax(const GetMax &o) : _Max(o._Max) {}

  void split(const GetMax &rhs) { _Max = rhs._Max; }
  void join (const GetMax &rhs) { if (rhs._Max < _Max) _Max = rhs._Max; }

  template <class TImage, class T>
  void operator ()(const TImage&, int, const T *p)
  {
    if (static_cast<double>(*p) > _Max) _Max = static_cast<double>(*p);
  }

  template <class T>
  void operator ()(int, int, int, int, const T *p)
  {
    if (static_cast<double>(*p) > _Max) _Max = static_cast<double>(*p);
  }

  double GetMaxAsDouble() const
  {
    if (_Max == voxel_limits<double>::min()) {
      return numeric_limits<double>::quiet_NaN();
    }
    return _Max;
  }

protected:
  double _Max;
};

// -----------------------------------------------------------------------------
/**
 * Finds the minimum and maximum intensities in an image
 */
struct GetMinMax : public VoxelReduction
{
  void Reset()
  {
    _Min = voxel_limits<double>::max();
    _Max = voxel_limits<double>::min();
  }

  GetMinMax() { Reset(); }
  GetMinMax(const GetMinMax &o) : _Min(o._Min), _Max(o._Max) {}

  void split(const GetMinMax &rhs)
  {
    _Min = rhs._Min;
    _Max = rhs._Max;
  }

  void join(const GetMinMax &rhs)
  {
    if (rhs._Min < _Min) _Min = rhs._Min;
    if (rhs._Max > _Max) _Max = rhs._Max;
  }

  template <class TImage, class T>
  void operator ()(const TImage&, int, const T *p)
  {
    if (static_cast<double>(*p) < _Min) _Min = static_cast<double>(*p);
    if (static_cast<double>(*p) > _Max) _Max = static_cast<double>(*p);
  }

  template <class T>
  void operator ()(int, int, int, int, const T *p)
  {
    if (static_cast<double>(*p) < _Min) _Min = static_cast<double>(*p);
    if (static_cast<double>(*p) > _Max) _Max = static_cast<double>(*p);
  }

  double GetMinAsDouble() const
  {
    if (_Min == voxel_limits<double>::max()) {
      return numeric_limits<double>::quiet_NaN();
    }
    return _Min;
  }

  double GetMaxAsDouble() const
  {
    if (_Max == voxel_limits<double>::min()) {
      return numeric_limits<double>::quiet_NaN();
    }
    return _Max;
  }

protected:
  double _Min, _Max;
};

// =============================================================================
// Thresholding
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * Sets all voxel values below a given threshold to this threshold value
 */
template <class T>
struct LowerThreshold : public VoxelFunction
{
  LowerThreshold(T threshold) : _LowerThreshold(threshold) {}

  void operator ()(const GenericImage<T>&, int, T *p)
  {
    if (*p < _LowerThreshold) *p = _LowerThreshold;
  }

  void operator ()(int, int, int, int, T *p)
  {
    if (*p < _LowerThreshold) *p = _LowerThreshold;
  }

  T _LowerThreshold; ///< Lower threshold value
};

// -----------------------------------------------------------------------------
/**
 * Sets all voxel values exceeding a given threshold to this threshold value
 */
template <class T>
struct UpperThreshold : public VoxelFunction
{
  UpperThreshold(T threshold) : _UpperThreshold(threshold) {}

  void operator ()(const GenericImage<T>&, int, T *p)
  {
    if (*p > _UpperThreshold) *p = _UpperThreshold;
  }

  void operator ()(int, int, int, int, T *p)
  {
    if (*p > _UpperThreshold) *p = _UpperThreshold;
  }

  T _UpperThreshold; ///< Upper threshold value
};

// -----------------------------------------------------------------------------
/**
 * Sets voxel values to a minimum/maximum value if below/above minimum/maximum
 */
template <class T>
struct Clamp : public VoxelFunction
{
  Clamp(T min, T max) : _LowerThreshold(min), _UpperThreshold(max) {}

  void operator ()(const GenericImage<T>&, int, T *p)
  {
    if      (*p < _LowerThreshold) *p = _LowerThreshold;
    else if (*p > _UpperThreshold) *p = _UpperThreshold;
  }

  void operator ()(int, int, int, int, T *p)
  {
    if      (*p < _LowerThreshold) *p = _LowerThreshold;
    else if (*p > _UpperThreshold) *p = _UpperThreshold;
  }

  T _LowerThreshold; ///< Lower threshold value
  T _UpperThreshold; ///< Upper threshold value
};

// =============================================================================
// Mathematical operations
// =============================================================================

// -----------------------------------------------------------------------------
/**
 * Takes the square root of the intensities
 */
struct Sqrt : public VoxelFunction
{
  template <class TImage, class T>
  void operator ()(const TImage&, int, T *p)
  {
    *p = static_cast<T>(sqrt(static_cast<double>(*p)));
  }

  template <class T>
  void operator ()(int, int, int, int, const T *p)
  {
    *p = static_cast<T>(sqrt(static_cast<double>(*p)));
  }
};

// =============================================================================
// Casts
// =============================================================================

// -----------------------------------------------------------------------------
/// Casts intensities to grey values
struct CastToGreyValue : public VoxelFunction
{
  template <class TImage, class T>
  void operator ()(const TImage &, int, T *v)
  {
    *v = static_cast<GreyPixel>(*v);
  }
};

// =============================================================================
// Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
/// Interpolate scalar or vector image
template <class TInterpolator = InterpolateImageFunction>
struct InterpolateImage : public VoxelFunction
{
  const BaseImage     *_Input;
  const TInterpolator *_Interpolator;
  BaseImage           *_Output;

  InterpolateImage(const TInterpolator *input, BaseImage *output)
  :
    _Input       (input->Input()),
    _Interpolator(input),
    _Output      (output)
  {
    mirtkAssert(input->Input()->T() == 1 || input->Input()->GetTSize() != .0,
                "Input image is not a multi-channel image");
  }

  InterpolateImage(const InterpolateImage &other)
  :
    _Input       (other._Input),
    _Interpolator(other._Interpolator),
    _Output      (other._Output)
  {}

  /// Interpolate scalar/vector image
  template <class T>
  void operator()(int i, int j, int k, int l, T *o)
  {
    // Convert output voxel to world coordinates
    double x = i, y = j, z = k, t = l;
        _Output->ImageToWorld(x, y, z);
    t = _Output->ImageToTime (t);
    // Convert world coordinates to input voxel coordinates
        _Input->WorldToImage(x, y, z);
    t = _Input->TimeToImage (t);
    // Interpolate input scalar or vector
    Vector v;
    _Interpolator->Evaluate(v, x, y, z, t);
    (*o) = voxel_cast<T>(v); // Convert to output type
  }
};

// -----------------------------------------------------------------------------
/// Interpolate scalar image
template <class TInterpolator = InterpolateImageFunction>
struct InterpolateScalarImage : public VoxelFunction
{
  const BaseImage     *_Input;
  const TInterpolator *_Interpolator;
  BaseImage           *_Output;

  InterpolateScalarImage(const TInterpolator *input, BaseImage *output)
  :
    _Input       (input->Input()),
    _Interpolator(input),
    _Output      (output)
  {
    mirtkAssert(input->Input()->N() == 1, "Input image has scalar data type");
  }

  InterpolateScalarImage(const InterpolateScalarImage &other)
  :
    _Input       (other._Input),
    _Interpolator(other._Interpolator),
    _Output      (other._Output)
  {}

  /// Interpolate scalar at (i, j, k, l)
  template <class T>
  void operator()(int i, int j, int k, int l, T *o)
  {
    // Convert output voxel to world coordinates
    double x = i, y = j, z = k, t = l;
    _Output->ImageToWorld(x, y, z);
    t = _Output->ImageToTime (t);
    // Convert world coordinates to input voxel coordinates
    _Input->WorldToImage(x, y, z);
    t = _Input->TimeToImage (t);
    // Interpolate input scalar and convert to output type
    (*o) = voxel_cast<T>(_Interpolator->Evaluate(x, y, z, t));
  }
};

// -----------------------------------------------------------------------------
/// Interpolate multi-channel image (3D+c)
template <class TInterpolator = InterpolateImageFunction>
struct InterpolateMultiChannelImage : public VoxelFunction
{
  const BaseImage     *_Input;
  const TInterpolator *_Interpolator;
  BaseImage           *_Output;
  int                  _NumberOfVoxels;

  InterpolateMultiChannelImage(const TInterpolator *input, BaseImage *output)
  :
    _Input         (input->Input()),
    _Interpolator  (input),
    _Output        (output),
    _NumberOfVoxels(output->NumberOfSpatialVoxels())
  {
    mirtkAssert(input->Input() != NULL, "Interpolator is initialized");
    mirtkAssert(output->T() == input->Input()->T(),
                "Input and output images have same number of channels");
  }

  InterpolateMultiChannelImage(const InterpolateMultiChannelImage &other)
  :
    _Input         (other._Input),
    _Interpolator  (other._Interpolator),
    _Output        (other._Output),
    _NumberOfVoxels(other._NumberOfVoxels)
  {}

  /// Interpolate all channels at (i, j, k) at once
  template <class T>
  void operator()(int i, int j, int k, int, T *o)
  {
    // Convert output voxel to world coordinates
    double x = i, y = j, z = k;
    _Output->ImageToWorld(x, y, z);
    // Convert world coordinates to input voxel coordinates
    _Input->WorldToImage(x, y, z);
    // Interpolate input channels
    double *v = new double[_Output->T()];
    _Interpolator->Evaluate(v, x, y, z);
    for (int l = 0; l < _Output->T(); ++l, o += _NumberOfVoxels) {
      (*o) = voxel_cast<T>(v[l]); // Convert to output type
    }
    delete[] v;
  }

  /// Interpolate all channels at (i, j, k) at once
  void operator()(int i, int j, int k, int, double *o)
  {
    // Convert output voxel to world coordinates
    double x = i, y = j, z = k;
    _Output->ImageToWorld(x, y, z);
    // Convert world coordinates to input voxel coordinates
    _Input->WorldToImage(x, y, z);
    // Interpolate input scalar or vector
    _Interpolator->Evaluate(o, x, y, z, _NumberOfVoxels);
  }
};

// =============================================================================
// Downsampling
// =============================================================================

// -----------------------------------------------------------------------------
template <class TVoxel>
struct DownsampleX : public VoxelFunction
{
  const GenericImage<TVoxel> *_Input;
  int                         _Offset;
  int                         _Factor;

  DownsampleX(const GenericImage<TVoxel> *input, int m)
    : _Input(input), _Offset((input->GetX() % m) / 2), _Factor(m)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int l, T *out)
  {
    *out = _Input->Get(_Offset + i * _Factor, j, k, l);
  }
};

// -----------------------------------------------------------------------------
template <class TVoxel>
struct DownsampleY : public VoxelFunction
{
  const GenericImage<TVoxel> *_Input;
  int                         _Offset;
  int                         _Factor;

  DownsampleY(const GenericImage<TVoxel> *input, int m)
    : _Input(input), _Offset((input->GetY() % m) / 2), _Factor(m)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int l, T *out)
  {
    *out = _Input->Get(i, _Offset + j * _Factor, k, l);
  }
};

// -----------------------------------------------------------------------------
template <class TVoxel>
struct DownsampleZ : public VoxelFunction
{
  const GenericImage<TVoxel> *_Input;
  int                         _Offset;
  int                         _Factor;

  DownsampleZ(const GenericImage<TVoxel> *input, int m)
    : _Input(input), _Offset((input->GetZ() % m) / 2), _Factor(m)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int l, T *out)
  {
    *out = _Input->Get(i, j, _Offset + k * _Factor, l);
  }
};


} } // namespace mirtk::UnaryVoxelFunction

#endif // MIRTK_UnaryVoxelFunction_H
