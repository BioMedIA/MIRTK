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

#ifndef MIRTK_GaussianBlurring_H
#define MIRTK_GaussianBlurring_H

#include "mirtk/SeparableConvolution.h"


namespace mirtk {


/**
 * Class for Gaussian blurring of images
 *
 * This class defines and implements the Gaussian blurring of images. The
 * blurring is implemented by three successive 1D convolutions with a 1D
 * Gaussian kernel.
 *
 * By default, if one isotropic Gaussian standard deviation is specified,
 * the first three dimensions of the image are blurred only. Otherwise,
 * the 1D convolution with a 1D Gaussian kernel is performed only for
 * dimensions of more than one voxel size and for which a non-zero
 * standard deviation for the Gaussian kernel has been set.
 */
template <class TVoxel>
class GaussianBlurring : public SeparableConvolution<TVoxel>
{
  mirtkInPlaceImageFilterMacro(GaussianBlurring, TVoxel);

protected:

  /// Type of convolution kernels
  typedef typename SeparableConvolution<TVoxel>::KernelType KernelType;

  /// Standard deviation of Gaussian kernel in x
  mirtkAttributeMacro(double, SigmaX);

  /// Standard deviation of Gaussian kernel in y
  mirtkAttributeMacro(double, SigmaY);

  /// Standard deviation of Gaussian kernel in z
  mirtkAttributeMacro(double, SigmaZ);

  /// Standard deviation of Gaussian kernel in t
  mirtkAttributeMacro(double, SigmaT);

protected:

  // Base class setters unused, should not be called by user
  virtual void KernelX(const KernelType *) {}
  virtual void KernelY(const KernelType *) {}
  virtual void KernelZ(const KernelType *) {}
  virtual void KernelT(const KernelType *) {}

  // Instantiated Gaussian kernels
  UniquePtr<KernelType> _GaussianKernel[4];

  /// Initialize 1D Gaussian kernel with sigma given in voxel units
  UniquePtr<KernelType> InitializeKernel(double);

  /// Initialize filter
  virtual void Initialize();

public:

  /// Constructor
  GaussianBlurring(double = 1.);

  /// Constructor
  GaussianBlurring(double, double, double = 0., double = 0.);

  /// Destructor
  ~GaussianBlurring();

  /// Set sigma
  virtual void SetSigma(double);

  /// Set sigma
  virtual void SetSigma(double, double, double = 0., double = 0.);

  /// Kernel size used for a given sigma (divided by voxel size)
  static int KernelSize(double);

};


} // namespace mirtk

#endif // MIRTK_GaussianBlurring_H
