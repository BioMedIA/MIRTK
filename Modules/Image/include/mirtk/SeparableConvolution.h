/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2017 Imperial College London
 * Copyright 2017 Andreas Schuh
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

#ifndef MIRTK_SeparableConvolution_H
#define MIRTK_SeparableConvolution_H

#include "mirtk/GenericImage.h"
#include "mirtk/ImageToImage.h"


namespace mirtk {


/**
 * Generic filter for convolution of an image with a separable kernel
 */
template <class TVoxel, class TKernel = RealPixel>
class SeparableConvolution : public ImageToImage<TVoxel>
{
  mirtkInPlaceImageFilterMacro(SeparableConvolution, TVoxel);

public:

  /// Type of convolution kernel
  typedef GenericImage<TKernel> KernelType;

protected:

  /// Convolution kernel along x axis
  mirtkPublicAggregateMacro(const KernelType, KernelX);

  /// Convolution kernel along y axis
  mirtkPublicAggregateMacro(const KernelType, KernelY);

  /// Convolution kernel along z axis
  mirtkPublicAggregateMacro(const KernelType, KernelZ);

  /// Convolution kernel along t axis
  mirtkPublicAggregateMacro(const KernelType, KernelT);

  /// Whether to normalize kernels
  mirtkPublicAttributeMacro(bool, Normalize);

  /// Whether to truncate convolution at background specified by image mask
  ///
  /// When this option is set and the image has a background mask, unmasked
  /// values of the convolution sum (i.e., background values) are ignored.
  mirtkPublicAttributeMacro(bool, UseBackgroundMask);

  /// Whether to truncate convolution at input background value if set
  ///
  /// This option is ignored when UseBackgroundMask is set and the image
  /// has a background mask image set.
  mirtkPublicAttributeMacro(bool, UseBackgroundValue);

  /// Whether to use padding value
  ///
  /// This option is ignored when UseBackgroundMask is set and the image
  /// has a background mask image set.
  ///
  /// When also UseBackground is set, the padding value is used instead
  /// when the image has a background value set. Otherwise the padding
  /// value is used when this option is turned on.
  mirtkPublicAttributeMacro(bool, UsePaddingValue);

  /// Padding value
  mirtkPublicAttributeMacro(double, PaddingValue);

  /// Number of vector components of input image to convolve
  ///
  /// By default, when the input image as temporal voxel size 0, each
  /// vector component stored in the temporal dimension is convolved
  /// separately. When this attribute is set to a positive value,
  /// only the first vector components up to the specified number
  /// of components are being convolved, and the output image may
  /// have less vector components stored in the temporal dimension
  /// than the input image.
  mirtkPublicAttributeMacro(int, Components);

protected:

  /// Initialize filter
  virtual void Initialize();

  /// Finalize filter
  virtual void Finalize();

public:

  /// Constructor
  ///
  /// \param[in] k Isotropic 1D convolution kernel.
  SeparableConvolution(const KernelType *k = nullptr);

  /// Constructor
  ///
  /// \param[in] kx Convolution kernel along x axis.
  /// \param[in] ky Convolution kernel along y axis.
  /// \param[in] kz Convolution kernel along z axis.
  /// \param[in] kt Convolution kernel along t axis.
  SeparableConvolution(const KernelType *kx,
                       const KernelType *ky,
                       const KernelType *kz = nullptr,
                       const KernelType *kt = nullptr);

  /// Destructor
  ~SeparableConvolution();

  /// Set isotropic convolution kernel for all dimensions
  void Kernel(const KernelType *);

  /// Check if given kernel is valid
  bool CheckKernel(const KernelType *) const;

  /// Convolve image
  virtual void Run();

  /// Convolve image along x only
  virtual void RunX();

  /// Convolve image along y only
  virtual void RunY();

  /// Convolve image along z only
  virtual void RunZ();

  /// Convolve image along t only
  virtual void RunT();

};


} // namespace mirtk

#endif // MIRTK_SeparableConvolution_H
