/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2015 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_ScalarGaussian_H
#define MIRTK_ScalarGaussian_H

#include "mirtk/ScalarFunction.h"
#include "mirtk/Vector4D.h"


namespace mirtk {


/**
 * Scalar Gaussian function.
 */
class ScalarGaussian : public ScalarFunction
{
  mirtkObjectMacro(ScalarGaussian);

protected:

  /// Variance of the Gaussian
  Vector4D<double> _Variance;

  /// Center of Gaussian
  Vector4D<double> _Center;

  /// Normalization of ND Gaussian function
  double _Norm[4];

  /// Pre-compute normalization factor(s)
  void ComputeNorm();

public:

  /// Construct Gaussian with isotropic sigma of 1 and center at origin
  ScalarGaussian();

  /// Construct Gaussian with isotropic sigma and center at origin
  ScalarGaussian(double);

  /// Construct Gaussian with isotropic sigma and specific center
  ScalarGaussian(double, double, double, double = .0, double = .0);

  /// Construct 3D Gaussian with anisotropic sigma and specific center
  ScalarGaussian(double, double, double,
                 double, double, double);

  /// Construct 4D Gaussian with anisotropic sigma and specific center
  ScalarGaussian(double, double, double, double,
                 double, double, double, double);

  /// Destructor
  virtual ~ScalarGaussian();

  /// Evaluate 1D Gaussian
  virtual double Evaluate(double);

  /// Evaluate 2D Gaussian
  virtual double Evaluate(double, double);

  /// Evaluate 3D Gaussian
  virtual double Evaluate(double, double, double);

  /// Evaluate 4D Gaussian
  virtual double Evaluate(double, double, double, double);
};


} // namespace mirtk

#endif // MIRTK_ScalarGaussian_H
