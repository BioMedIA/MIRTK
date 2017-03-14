/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
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

#ifndef MIRTK_NormalizedIntensityCrossCorrelation_H
#define MIRTK_NormalizedIntensityCrossCorrelation_H

#include "mirtk/ImageSimilarity.h"

#include "mirtk/Voxel.h"
#include "mirtk/GenericImage.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Parallel.h"


namespace mirtk {


/**
 * (Local) Normalized cross correlation image similarity measure
 */
class NormalizedIntensityCrossCorrelation : public ImageSimilarity
{
  mirtkEnergyTermMacro(NormalizedIntensityCrossCorrelation, EM_NCC);

  // ---------------------------------------------------------------------------
  // Types
public:

  typedef GradientImageType::VoxelType   RealType;
  typedef GenericImage<RealType>         RealImage;
  typedef GenericImage<RealType>         KernelImage;

  /// Enumeration of local window size units
  enum Units { UNITS_Default, UNITS_MM, UNITS_Voxel };

private:

  /// Enumeration of local window types
  enum KernelType { CustomKernel /**< reserved */, BoxWindow, GaussianKernel };

  // ---------------------------------------------------------------------------
  // Attributes

  /// Type of kernel for output in parameter file
  mirtkAttributeMacro(enum KernelType, KernelType);

  /// Local window kernel in x dimension
  /// \note The kernel image is deleted by the destructor of this class.
  mirtkComponentMacro(KernelImage, KernelX);

  /// Local window kernel in y dimension
  /// \note The kernel image is deleted by the destructor of this class.
  mirtkComponentMacro(KernelImage, KernelY);

  /// Local window kernel in z dimension
  /// \note The kernel image is deleted by the destructor of this class.
  mirtkComponentMacro(KernelImage, KernelZ);

  /// Term A of normalized cross-correlation
  mirtkComponentMacro(RealImage, A);

  /// Term B of normalized cross-correlation
  mirtkComponentMacro(RealImage, B);

  /// Term C of normalized cross-correlation
  mirtkComponentMacro(RealImage, C);

  /// Mean normalized source image
  mirtkComponentMacro(RealImage, S);

  /// Mean normalized target image
  mirtkComponentMacro(RealImage, T);

  /// Global inner products computed over whole image domain
  double _GlobalA, _GlobalB, _GlobalC, _GlobalS, _GlobalT;

  /// Temporary image used by ComputeWeightedAverage and WriteDataSets
  mutable RealImage _Temp;

  /// Size of box window or Full Width at Tenth Maximum (FWTM) of Gaussian kernel.
  /// A negative value indicates voxel units, while a positive value corresponds
  /// to world units (i.e., mm). When 0, the global normalized cross correlation
  /// is computed instead of local cross correlation within a neighborhood.
  mirtkPublicAttributeMacro(Vector3D<double>, NeighborhoodSize);

  /// Radius of local neighborhood window in number of voxels
  mirtkReadOnlyAttributeMacro(Vector3D<int>, NeighborhoodRadius);

  /// Sum of local correlation coefficients / global normalized cross correlation
  double _Sum;

  /// Size of overlap domain
  int _N;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Create 1D Gaussian kernel with given standard deviation
  static KernelImage *CreateGaussianKernel(double);

  /// Reset local window kernel
  virtual void ClearKernel();

public:

  /// Constructor
  NormalizedIntensityCrossCorrelation(const char * = "");

  /// Copy constructor
  NormalizedIntensityCrossCorrelation(const NormalizedIntensityCrossCorrelation &);

  /// Assignment operator
  NormalizedIntensityCrossCorrelation &operator =(const NormalizedIntensityCrossCorrelation &);

  /// Destructor
  ~NormalizedIntensityCrossCorrelation();

  // ---------------------------------------------------------------------------
  // Parameters

  /// Set kernel to a box window with given radius
  virtual void SetKernelToBoxWindow(double, double = -1, double = -1, Units = UNITS_MM);

  /// Set kernel to a Gaussian with given standard deviation
  virtual void SetKernelToGaussian(double, double = -1, double = -1, Units = UNITS_MM);

protected:

  /// Set parameter value from string
  virtual bool SetWithPrefix(const char *, const char *);

public:

  // Import other overloads
  using ImageSimilarity::Parameter;

  /// Get parameter key/value as string map
  virtual ParameterList Parameter() const;

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize similarity measure once input and parameters have been set
  virtual void Initialize();

  // ---------------------------------------------------------------------------
  // Evaluation

  /// Update moving image and internal state of similarity measure
  virtual void Update(bool = true);

  /// Exclude region from similarity evaluation
  virtual void Exclude(const blocked_range3d<int> &);

  /// Include region in similarity evaluation
  virtual void Include(const blocked_range3d<int> &);

protected:

  /// Convolve image with separable window kernel
  virtual void ComputeWeightedAverage(const blocked_range3d<int> &, RealImage *);

  /// Compute statistics of local intensity distribution
  virtual void ComputeStatistics(const blocked_range3d<int> &,
                                 const RegisteredImage *,
                                 RealImage *, RealImage *);

  /// Evaluate similarity of images
  virtual double Evaluate();

  /// Evaluate non-parametric similarity gradient w.r.t the given image
  virtual bool NonParametricGradient(const RegisteredImage *, GradientImageType *);

  // ---------------------------------------------------------------------------
  // Debugging
public:

  /// Return unweighted and unnormalized raw energy term value
  /// \remarks Use for progress reporting only.
  virtual double RawValue(double) const;

  /// Print debug information
  virtual void Print(Indent = 0) const;

  /// Write input of data fidelity term
  virtual void WriteDataSets(const char *, const char *, bool = true) const;

};


} // namespace mirtk

#endif // MIRTK_NormalizedIntensityCrossCorrelation_H
