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

#ifndef MIRTK_GaussianBlurring_H
#define MIRTK_GaussianBlurring_H

#include <mirtkGenericImage.h>
#include <mirtkImageToImage.h>


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
template <class VoxelType>
class GaussianBlurring : public ImageToImage<VoxelType>
{
  mirtkInPlaceImageFilterMacro(GaussianBlurring);

  /// Standard deviation of Gaussian kernel in x
  mirtkAttributeMacro(double, SigmaX);

  /// Standard deviation of Gaussian kernel in y
  mirtkAttributeMacro(double, SigmaY);

  /// Standard deviation of Gaussian kernel in z
  mirtkAttributeMacro(double, SigmaZ);

  /// Standard deviation of Gaussian kernel in t
  mirtkAttributeMacro(double, SigmaT);

protected:

  /// Gaussian convolution kernel
  GenericImage<RealPixel> *_Kernel;

  /// Initialize filter
  virtual void Initialize();

  /// Initialize 1D Gaussian kernel with sigma given in voxel units
  virtual void InitializeKernel(double);

  /// Finalize filter
  virtual void Finalize();

public:

  /// Constructor
  GaussianBlurring(double = 1.0);

  /// Constructor
  GaussianBlurring(double, double, double = .0, double = .0);

  /// Destructor
  ~GaussianBlurring();

  /// Run Gaussian blurring
  virtual void Run();

  /// Run Gaussian blurring along x only
  virtual void RunX();

  /// Run Gaussian blurring along y only
  virtual void RunY();

  /// Run Gaussian blurring along z only
  virtual void RunZ();

  /// Run Gaussian blurring along t only
  virtual void RunT();

  /// Set sigma
  virtual void SetSigma(double);

  /// Set sigma
  virtual void SetSigma(double, double, double = .0, double = .0);

  /// Kernel size used for a given sigma (divided by voxel size)
  static int KernelSize(double);

};


} // namespace mirtk

#endif // MIRTK_GaussianBlurring_H
