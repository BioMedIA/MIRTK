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

#ifndef MIRTK_VoxelDomain_H
#define MIRTK_VoxelDomain_H


namespace mirtk {


/**
 * Types which implement inside/outside image domain checks
 *
 * This namespace defines auxiliary types which implement inside/outside image
 * domain checks for voxels as static methods. These types can be used as
 * Domain template argument as used by the ForEachVoxelIf and ParallelForEachVoxelIf
 * template functions to restrict the evaluation of the voxel function to a specific
 * subdomain of the image(s). Moreover, if it is known how the subdomain is
 * defined, e.g., by a binary mask, a background value, a lower foreground threshold,
 * or an interpolation domain, a specialized inside/outside domain check is more
 * efficient as it only needs to test one particular condition. The more generic
 * default domain check has to evaluate multiple conditions to first decide
 * how the domain is being defined and then evaluate the actual inside/outside
 * condition for this domain. This requires multiple if-statements and can be
 * expensive when all these if-statements have to be evaluated for each voxel.
 * Therefore, whenever possible use a specialized domain check to focus on the
 * condition that matters.
 *
 * Example:
 * \code
 * #include "mirtk/GenericImage.h"
 * #include "mirtk/VoxelFunction.h"
 * #include "mirtk/VoxelDomain.h"
 *
 * namespace mirtk {
 *
 *
 * void ProcessImageVoxelByVoxel(GreyImage &image)
 * {
 *   using namespace ForEachVoxelDomain;
 *   ForEachVoxelIf<AboveBackgroundLevel>(image, func);
 * }
 *
 *
 * } // namespace mirtk
 * \endcode
 *
 * \note The last image passed to ForEachVoxelIf/ParallelForEachVoxelIf is
 *       the one which defines the image domain over which to iterate as it
 *       is often an output image. Therefore, this last image is passed on
 *       to the domain checks defined within the ImageDomain namespace.
 *
 */
namespace ForEachVoxelDomain {


// -----------------------------------------------------------------------------
/**
 * Checks if voxel is in foreground using BaseImage::IsForeground
 */
struct Foreground
{
  static inline bool IsInside(const BaseImage &image, int idx, const void *)
  {
    return image.IsForeground(idx);
  }

  static inline bool IsInside(const BaseImage &image, int i, int j, int k, int l, const void *)
  {
    return image.IsForeground(i, j, k, l);
  }
};

// -----------------------------------------------------------------------------
/**
 * Checks if voxel is in foreground using background value of image
 */
struct NotBackgroundValue
{
  template <class T>
  static inline bool IsInside(const BaseImage &image, int, const T *p)
  {
    const double bg = image.GetBackgroundValueAsDouble();
    return (*p != bg) && (bg == bg || *p == *p /*, i.e., not both NaN */);
  }

  template <class T>
  static inline bool IsInside(const BaseImage &image, int, int, int, int, const T *p)
  {
    const double bg = image.GetBackgroundValueAsDouble();
    return (*p != bg) && (bg == bg || *p == *p /*, i.e., not both NaN */);
  }
};

// -----------------------------------------------------------------------------
/**
 * Checks if voxel is in foreground using background value of image as threshold
 */
struct AboveBackgroundLevel
{
  template <class T>
  static inline bool IsInside(const BaseImage &image, int, const T *p)
  {
    const double bg = image.GetBackgroundValueAsDouble();
    return (*p > bg) && (bg == bg || *p == *p /*, i.e., not both NaN */);
  }

  template <class T>
  static inline bool IsInside(const BaseImage &image, int, int, int, int, const T *p)
  {
    const double bg = image.GetBackgroundValueAsDouble();
    return (*p > bg) && (bg == bg || *p == *p /*, i.e., not both NaN */);
  }
};

// -----------------------------------------------------------------------------
/**
 * Checks if voxel is in background using BaseImage::IsBackground
 */
struct Background
{
  static inline bool IsInside(const BaseImage &image, int idx, const void *)
  {
    return image.IsBackground(idx);
  }

  static inline bool IsInside(const BaseImage &image, int i, int j, int k, int l, const void *)
  {
    return image.IsBackground(i, j, k, l);
  }
};

// -----------------------------------------------------------------------------
/**
 * Checks if voxel is in background using background value of image
 */
struct BackgroundValue
{
  template <class T>
  static inline bool IsInside(const BaseImage &image, int, const T *p)
  {
    const double bg = image.GetBackgroundValueAsDouble();
    return (*p == bg) || (bg != bg && *p != *p /*, i.e., both NaN */);
  }

  template <class T>
  static inline bool IsInside(const BaseImage &image, int, int, int, int, const T *p)
  {
    const double bg = image.GetBackgroundValueAsDouble();
    return (*p == bg) || (bg != bg && *p != *p /*, i.e., both NaN */);
  }
};

// -----------------------------------------------------------------------------
/**
 * Checks if voxel is in image domain using mask of image
 */
struct InMask
{
  static inline bool IsInside(const BaseImage &image, int idx, const void *)
  {
    const BinaryImage * const mask = image.GetMask();
    if (mask->GetT() > 1) {
      return (mask->Get(idx) != BinaryPixel(0));
    } else {
      return (mask->Get(idx % (image.GetX() * image.GetY() * image.GetZ())) != BinaryPixel(0));
    }
  }

  static inline bool IsInside(const BaseImage &image, int i, int j, int k, int l, const void *)
  {
    const BinaryImage * const mask = image.GetMask();
    return (mask->Get(i, j, k, mask->GetT() > 1 ? l : 0) != BinaryPixel(0));
  }
};

// -----------------------------------------------------------------------------
/**
 * Checks if voxel is in image domain using mask of image (BaseImage::GetMask)
 *
 * Use this domain check if the foreground mask is never four-dimensional even
 * when the image itself has multiple frames, channels, or vector components.
 */
struct InSpatialMask
{
  static inline bool IsInside(const BaseImage &image, int idx, const void *)
  {
    return (image.GetMask()->Get(idx % (image.GetX() * image.GetY() * image.GetZ())) != BinaryPixel(0));
  }

  static inline bool IsInside(const BaseImage &image, int i, int j, int k, int, const void *)
  {
    return (image.GetMask()->Get(i, j, k) != BinaryPixel(0));
  }
};


} // namespace ForEachVoxelDomain


/**
 * Types which implement inside/outside image interpolation domain checks
 *
 * This namespace defines auxiliary types which implement inside/outside
 * interpolation domain checks for (transformed) voxel coordinates as
 * static methods. These types can be used, in particular, as InputDomain
 * template argument of the voxel transformation functions derived from
 * MultipleVoxelTransformation::BaseTransform.
 */
namespace InterpolationDomain {


// -----------------------------------------------------------------------------
/**
 * Wraps call to InterpolateImageFunction::IsForeground
 *
 * Checks if voxel coordinates are inside the domain for which the interpolation
 * is defined and also if all values required for interpolation are part of the
 * foreground region of the image.
 */
struct Foreground
{
  template <class InterpolateImageFunction>
  static inline bool IsInside(InterpolateImageFunction *interpolator, const double &x, const double &y, const double &z)
  {
    return interpolator->IsForeground(x, y, z);
  }
};

// -----------------------------------------------------------------------------
/**
 * Wraps call to InterpolateImageFunction::IsInside
 *
 * Use as argument for the InputDomain template parameter of the voxel
 * transformation functors to interpolate the input at a given transformed
 * voxel when its coordinates are inside the domain for which the interpolation
 * is well defined.
 *
 * More efficient then the default Foreground domain check, but only ensures that
 * the transformed coordinates are within the interpolation domain of the input
 * image(s). Can be used when background value of input images is very dominant,
 * e.g., minimum negative value supported by scalar type, and thus interpolated
 * values are below a certain threshold, or if influence of background values
 * is considered negligible.
 */
struct Inside
{
  template <class InterpolateImageFunction>
  static inline bool IsInside(InterpolateImageFunction *interpolator, const double &x, const double &y, const double &z)
  {
    return interpolator->IsInside(x, y, z);
  }
};


} // namespace InterpolationDomain


} // namespace mirtk

#endif // MIRTK_VoxelDomain_H
