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

#include "mirtk/ScalarGaussian.h"
#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void ScalarGaussian::ComputeNorm()
{
  _Norm[0] = (_Variance._x > .0 ? (1.0      / sqrt(two_pi * _Variance._x)) : 1.0);
  _Norm[1] = (_Variance._y > .0 ? (_Norm[0] / sqrt(two_pi * _Variance._y)) : 1.0);
  _Norm[2] = (_Variance._z > .0 ? (_Norm[1] / sqrt(two_pi * _Variance._z)) : 1.0);
  _Norm[3] = (_Variance._t > .0 ? (_Norm[2] / sqrt(two_pi * _Variance._t)) : 1.0);
  if (_Variance._x <= .0) _Variance._x = 1.0;
  if (_Variance._y <= .0) _Variance._y = 1.0;
  if (_Variance._z <= .0) _Variance._z = 1.0;
  if (_Variance._t <= .0) _Variance._t = 1.0;
}

// -----------------------------------------------------------------------------
ScalarGaussian::ScalarGaussian()
:
  _Variance(1.0),
  _Center(.0)
{
  ComputeNorm();
}

// -----------------------------------------------------------------------------
ScalarGaussian::ScalarGaussian(double sigma)
:
  _Variance(sigma * sigma),
  _Center(.0)
{
  ComputeNorm();
}

// -----------------------------------------------------------------------------
ScalarGaussian::ScalarGaussian(double sigma, double x_0, double y_0, double z_0, double t_0)
:
  _Variance(sigma * sigma),
  _Center(x_0, y_0, z_0, t_0)
{
  ComputeNorm();
}

// -----------------------------------------------------------------------------
ScalarGaussian::ScalarGaussian(double sigma_x, double sigma_y, double sigma_z,
                               double x_0, double y_0, double z_0)
:
  _Variance(sigma_x*sigma_x, sigma_y*sigma_y, sigma_z*sigma_z, .0),
  _Center(x_0, y_0, z_0, .0)
{
  ComputeNorm();
}

// -----------------------------------------------------------------------------
ScalarGaussian::ScalarGaussian(double sigma_x, double sigma_y, double sigma_z, double sigma_t,
                               double x_0, double y_0, double z_0, double t_0)
:
  _Variance(sigma_x*sigma_x, sigma_y*sigma_y, sigma_z*sigma_z, sigma_t*sigma_t),
  _Center(x_0, y_0, z_0, t_0)
{
  ComputeNorm();
}

// -----------------------------------------------------------------------------
ScalarGaussian::~ScalarGaussian()
{
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
double ScalarGaussian::Evaluate(double x)
{
  x -= _Center._x;
  return _Norm[0] * exp(- 0.5 * (x * x) / _Variance._x);
}

// -----------------------------------------------------------------------------
double ScalarGaussian::Evaluate(double x, double y)
{
  x -= _Center._x, y -= _Center._y;
  return _Norm[1] * exp(- 0.5 * (  (x * x) / _Variance._x
                                 + (y * y) / _Variance._y));
}

// -----------------------------------------------------------------------------
double ScalarGaussian::Evaluate(double x, double y, double z)
{
  x -= _Center._x, y -= _Center._y, z -= _Center._z;
  return _Norm[2] * exp(- 0.5 * (  (x * x) / _Variance._x
                                 + (y * y) / _Variance._y
                                 + (z * z) / _Variance._z));
}

// -----------------------------------------------------------------------------
double ScalarGaussian::Evaluate(double x, double y, double z, double t)
{
  x -= _Center._x, y -= _Center._y, z -= _Center._z, t -= _Center._t;
  return _Norm[3] * exp(- 0.5 * (  (x * x) / _Variance._x
                                 + (y * y) / _Variance._y
                                 + (z * z) / _Variance._z
                                 + (t * t) / _Variance._t));
}


} // namespace mirtk
