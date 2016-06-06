/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/PolynomialSolvers.h"

#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Root finding
// =============================================================================

// -----------------------------------------------------------------------------
int SolveCubicEquation(double x[3], double a, double b, double c)
{
	double a2 = a*a;
  double q  = (a2 - 3.*b)/9.;
	double r  = (a*(2.*a2 - 9.*b) + 27.*c) / 54.;
  double r2 = r*r;
	double q3 = q*q*q;
	double A, B;
  if(r2<q3) {
    double t = acos(clamp(r / sqrt(q3), -1., 1.));
    a /= 3.;
    q = -2. * sqrt(q);
    x[0] = q * cos( t           / 3.) - a;
    x[1] = q * cos((t + two_pi) / 3.) - a;
    x[2] = q * cos((t - two_pi) / 3.) - a;
    return 3;
  } else {
    A = -pow(abs(r) + sqrt(r2 - q3), 1./3.);
    if (r < 0.) A = -A;
    B = (A == 0. ? 0 : B = q / A);
    a /= 3.;
    x[0] = (A + B) - a;
    x[1] = -.5 * (A + B) -a;
    x[2] =  .5 * sqrt(3.) * (A - B);
    if (abs(x[2]) < 1e-14) {
      x[2] = x[1];
      return 2;
    }
    return 1;
  }
}

// =============================================================================
// Minimization
// =============================================================================

// -----------------------------------------------------------------------------
double MinimumOf4thDegreePolynomial(double a, double b, double c)
{
  int    i, j, n;
  double x[3], x2;

  // Find roots of derivative
  n = SolveCubicEquation(x, 4., 3. * a, 2. * b, c);

  // Determine global minimum
  double value, min_value = numeric_limits<double>::max();
  for (i = j = 0; i < n; ++i) {
    x2 = x[i] * x[i];
    value = x2 * x2 + a * x2 * x[i] + b * x2 + c * x[i];
    if (value < min_value) {
      min_value = value, j = i;
    }
  }

  return x[j];
}


} // namespace mirtk
