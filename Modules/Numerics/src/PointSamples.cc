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

// Disable MSVC warning caused by <boost/random.hpp>
// (see https://svn.boost.org/trac/boost/ticket/11426)
#ifdef _MSC_VER
#  pragma warning(disable:4996)
#endif

#include "mirtk/PointSamples.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Array.h"

#include "boost/random.hpp"


namespace mirtk {


// =============================================================================
// Auxiliaries
// =============================================================================

namespace PointSamplesUtils {


typedef boost::mt19937 RandomNumberGeneratorType;


} // namespace PointSamplesUtils
using namespace PointSamplesUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
PointSamples::PointSamples(int n, int seed)
:
  PointSet(n), _RandomNumberGenerator(new RandomNumberGeneratorType())
{
  RandomNumberGeneratorType &rng = *(RandomNumberGeneratorType*)_RandomNumberGenerator;
  if (seed < 0) rng.seed(static_cast<unsigned int>(std::time(0)));
  else          rng.seed(seed);
}

// -----------------------------------------------------------------------------
PointSamples::~PointSamples()
{
  delete (RandomNumberGeneratorType*)_RandomNumberGenerator;
}

// =============================================================================
// Uniform grid sampling
// =============================================================================

// -----------------------------------------------------------------------------
void PointSamples::SampleGrid(const Point &p1, const Point &p2,
                              int nx, int ny, int nz)
{
  SampleGrid(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z, nx, ny, nz);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleGrid(const Point &p1, const Point &p2,
                              double dx, double dy, double dz)
{
  SampleGrid(p1._x, p1._y, p1._z, p2._x, p2._y, p2._z, dx, dy, dz);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleGrid(double x1, double y1, double z1,
                              double x2, double y2, double z2,
                              int    nx, int    ny, int    nz)
{
  if (x2 < x1) swap(x1, x2);
  if (y2 < y1) swap(y1, y2);
  if (z2 < z1) swap(z1, z2);
  if (nx <= 0) nx = 1;
  if (ny <= 0) ny = 1;
  if (nz <= 0) nz = 1;
  const double sx = x2 - x1;
  const double sy = x2 - x1;
  const double sz = x2 - x1;
  const double dx = sx / nx;
  const double dy = sy / ny;
  const double dz = sz / nz;
  Size(nx * ny * nz);
  int n = 0;
  for (int k = 0; k < nz; ++k)
  for (int j = 0; j < ny; ++j)
  for (int i = 0; i < nx; ++i, ++n) {
    _data[n] = Point(x1 + i * dx, y1 + j * dy, z1 + k * dz);
  }
}

// -----------------------------------------------------------------------------
void PointSamples::SampleGrid(double x1, double y1, double z1,
                              double x2, double y2, double z2,
                              double dx, double dy, double dz)
{
  if (x2 < x1) swap(x1, x2);
  if (y2 < y1) swap(y1, y2);
  if (z2 < z1) swap(z1, z2);
  if (dx < 0) dx = .0;
  if (dy < 0) dy = .0;
  if (dz < 0) dz = .0;
  const double sx = x2 - x1;
  const double sy = x2 - x1;
  const double sz = x2 - x1;
  const int    nx = (dx > .0 ? iround(sx / dx) : 1);
  const int    ny = (dy > .0 ? iround(sy / dy) : 1);
  const int    nz = (dz > .0 ? iround(sz / dz) : 1);
  Size(nx * ny * nz);
  int n = 0;
  for (int k = 0; k < nz; ++k)
  for (int j = 0; j < ny; ++j)
  for (int i = 0; i < nx; ++i, ++n) {
    _data[n] = Point(x1 + i * dx, y1 + j * dy, z1 + k * dz);
  }
}

// =============================================================================
// Uniform spherical distribution
// =============================================================================

// -----------------------------------------------------------------------------
void PointSamples::SampleSphere(double r)
{
  SampleSphere(.0, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleSphere(double c, double r)
{
  SampleSphere(c, c, c, r, r, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleSphere(const Point & c, double r)
{
  SampleSphere(c._x, c._y, c._z, r, r, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleSphere(const Point & c, double rx, double ry, double rz)
{
  SampleSphere(c._x, c._y, c._z, rx, ry, rz);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleSphere(double cx, double cy, double cz, double r)
{
  SampleSphere(cx, cy, cz, r, r, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleSphere(double cx, double cy, double cz,
                                double rx, double ry, double rz)
{
  typedef boost::uniform_on_sphere<double>                                   Distribution;
  typedef boost::variate_generator<RandomNumberGeneratorType&, Distribution> Generator;

  RandomNumberGeneratorType &rng = *(RandomNumberGeneratorType*)_RandomNumberGenerator;

  Distribution dist(3);
  Generator next(rng, dist);

  Array<double> p(3);
  for (int i = 0; i < _n; ++i) {
    p = next();
    _data[i]._x = cx + rx * p[0];
    _data[i]._y = cy + ry * p[1];
    _data[i]._z = cz + rz * p[2];
  }
}

// =============================================================================
// Regular spherical samples
// =============================================================================

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularSphere(double r)
{
  SampleRegularSphere(.0, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularSphere(double c, double r)
{
  SampleRegularSphere(c, c, c, r, r, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularSphere(const Point & c, double r)
{
  SampleRegularSphere(c._x, c._y, c._z, r, r, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularSphere(const Point & c, double rx, double ry, double rz)
{
  SampleRegularSphere(c._x, c._y, c._z, rx, ry, rz);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularSphere(double cx, double cy, double cz, double r)
{
  SampleRegularSphere(cx, cy, cz, r, r, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularSphere(double cx, double cy, double cz,
                                       double rx, double ry, double rz)
{
  cerr << "PointSamples::SampleRegularSphere: Not implemented" << endl;
  exit(1);

  const int n = _n;

  // Sample half sphere
  _n = (n + 1) / 2;
  SampleRegularHalfSphere(cx, cy, cz, rx, ry, rz);

  // TODO: Duplicate samples with opposite sign
  _n = n;
}

// =============================================================================
// Regular spherical sampling of half sphere
// =============================================================================

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularHalfSphere(double r)
{
  SampleRegularHalfSphere(.0, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularHalfSphere(double c, double r)
{
  SampleRegularHalfSphere(c, c, c, r, r, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularHalfSphere(const Point & c, double r)
{
  SampleRegularHalfSphere(c._x, c._y, c._z, r, r, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularHalfSphere(const Point & c, double rx, double ry, double rz)
{
  SampleRegularHalfSphere(c._x, c._y, c._z, rx, ry, rz);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularHalfSphere(double cx, double cy, double cz, double r)
{
  SampleRegularHalfSphere(cx, cy, cz, r, r, r);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleRegularHalfSphere(double cx, double cy, double cz,
                                           double rx, double ry, double rz)
{
  const double epsilon = 1e-6;

  // Determine angular sampling rate from desired number of samples (_n)
  const double delta    = 1.23 * sqrt(4.0 / _n);
  const double dTheta   = two_pi / round(two_pi / clamp(delta, .0, two_pi));
  const int    maxTheta = iround(two_pi / dTheta);

  Array<int>    maxPhi(maxTheta + 1);
  Array<double> dPhi  (maxTheta + 1);

  int    n     = 0; // actual number of samples
  double theta = .0;
  for (int i = 0; i <= maxTheta; ++i, theta += dTheta) {
    // On the tip, the circle is just a single point
    if (theta < epsilon) {
      maxPhi[i] = 0;
      dPhi  [i] = .0;
    }
    // On the equator, we only need half of the directions
    else if (theta >= (two_pi - epsilon)) {
      dPhi  [i] = min(two_pi, dTheta / sin(theta));
      maxPhi[i] = max(1, static_cast<int>(pi / dPhi[i]) - 1);
      dPhi  [i] = pi / (maxPhi[i] + 1);
    }
    // Otherwise, let's sample the entire circle parallel to the equator
    else {
      dPhi  [i] = min(two_pi, dTheta / sin(theta));
      maxPhi[i] = max(3, 2 * static_cast<int>(pi / dPhi[i]) - 1);
      dPhi  [i] = two_pi / (maxPhi[i] + 1);
    }
    n += maxPhi[i] + 1;
  }

  Size(n);

  theta = .0;
  for (int i = 0, k = 0; i <= maxTheta; ++i, theta += dTheta) {
    double phi = .0;
    for (int j = 0; j <= maxPhi[i]; ++j, phi += dPhi[i]) {
      double sin_theta = sin(theta);
      _data[k++] = Point(cx + rx * sin_theta * cos(phi),
                         cy + ry * sin_theta * sin(phi),
                         cz + rz * cos(theta));
    }
  }
}

// =============================================================================
// Normal distribution
// =============================================================================

// -----------------------------------------------------------------------------
void PointSamples::SampleGaussian(double s)
{
  SampleGaussian(.0, s);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleGaussian(double m, double s)
{
  SampleGaussian(m, m, m, s, s, s);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleGaussian(const Point & m, double s)
{
  SampleGaussian(m._x, m._y, m._z, s, s, s);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleGaussian(const Point & m, double sx, double sy, double sz)
{
  SampleGaussian(m._x, m._y, m._z, sx, sy, sz);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleGaussian(double mx, double my, double mz, double s)
{
  SampleGaussian(mx, my, mz, s, s, s);
}

// -----------------------------------------------------------------------------
void PointSamples::SampleGaussian(double mx, double my, double mz,
                                  double sx, double sy, double sz)
{
  typedef boost::normal_distribution<double>                                 Distribution;
  typedef boost::variate_generator<RandomNumberGeneratorType&, Distribution> Generator;

  RandomNumberGeneratorType &rng = *(RandomNumberGeneratorType*)_RandomNumberGenerator;

  Distribution distx(mx, sx);
  Distribution disty(my, sy);
  Distribution distz(mz, sz);

  Generator x(rng, distx);
  Generator y(rng, disty);
  Generator z(rng, distz);

  for (int i = 0; i < _n; ++i) {
    _data[i]._x = x();
    _data[i]._y = y();
    _data[i]._z = z();
  }
}


} // namespace mirtk
