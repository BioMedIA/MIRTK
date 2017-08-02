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

#ifndef MIRTK_CubicBSplineConvolution_H
#define MIRTK_CubicBSplineConvolution_H

#include "mirtk/SeparableConvolution.h"


namespace mirtk {


/**
 * Convolves image with separable cubic B-spline kernel
 */
template <class TVoxel>
class CubicBSplineConvolution : public SeparableConvolution<TVoxel>
{
  mirtkInPlaceImageFilterMacro(CubicBSplineConvolution, TVoxel);

protected:

  /// Type of convolution kernels
  typedef typename SeparableConvolution<TVoxel>::KernelType KernelType;

  /// Radius of cubic B-spline support region along x axis
  ///
  /// - r<0: voxel units
  /// - r=0: no smoothing in x
  /// - r>0: mm units
  mirtkAttributeMacro(double, RadiusX);

  /// Radius of cubic B-spline support region along y axis
  ///
  /// - r<0: voxel units
  /// - r=0: no smoothing in y
  /// - r>0: mm units
  mirtkAttributeMacro(double, RadiusY);

  /// Radius of cubic B-spline support region along z axis
  ///
  /// - r<0: voxel units
  /// - r=0: no smoothing in z
  /// - r>0: mm units
  mirtkAttributeMacro(double, RadiusZ);

  /// Radius of cubic B-spline support region along t axis
  ///
  /// - r<0: voxel units
  /// - r=0: no smoothing in t
  /// - r>0: mm units
  mirtkAttributeMacro(double, RadiusT);

protected:

  // Base class setters unused, should not be called by user
  virtual void KernelX(const KernelType *) {}
  virtual void KernelY(const KernelType *) {}
  virtual void KernelZ(const KernelType *) {}
  virtual void KernelT(const KernelType *) {}

  // Instantiated cubic B-spline kernels
  UniquePtr<KernelType> _BSplineKernel[4];

  /// Initialize 1D cubic B-spline kernel with radius given in voxel units
  UniquePtr<KernelType> InitializeKernel(double);

  /// Initialize filter
  virtual void Initialize();

public:

  /// Constructor
  CubicBSplineConvolution(double = 2.);

  /// Constructor
  CubicBSplineConvolution(double, double, double = 0., double = 0.);

  /// Destructor
  ~CubicBSplineConvolution();

  /// Set isotropic radius of spatial cubic B-spline support region
  virtual void Radius(double);

  /// Set radius of cubic B-spline support region
  virtual void Radius(double, double, double = 0., double = 0.);

  /// Kernel size used for a given radius (divided by voxel size)
  static int KernelSize(double);

};


} // namespace mirtk

#endif // MIRTK_CubicBSplineConvolution_H
