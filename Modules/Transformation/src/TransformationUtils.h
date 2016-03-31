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

#ifndef MIRTK_TransformationUtils_H
#define MIRTK_TransformationUtils_H

#include "mirtk/Transformation.h"
#include "mirtk/GenericImage.h"
#include "mirtk/Parallel.h"


namespace mirtk {
namespace TransformationUtils {


// -----------------------------------------------------------------------------
/// Parallel helper for Transformation::Approximate implementations
class SubDisplacements
{
  const Transformation *_Transformation;
  const double         *_x,  *_y,  *_z, *_t;
  double               *_dx, *_dy, *_dz;

public:

  SubDisplacements(const Transformation *transformation,
                   const double *x, const double *y, const double *z, const double *t,
                   double *dx, double *dy, double *dz)
  :
    _Transformation(transformation), _x(x), _y(y), _z(z), _t(t), _dx(dx), _dy(dy), _dz(dz)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    const double *x  = _x  + re.begin();
    const double *y  = _y  + re.begin();
    const double *z  = _z  + re.begin();
    const double *t  = _t  + re.begin();
    double       *dx = _dx + re.begin();
    double       *dy = _dy + re.begin();
    double       *dz = _dz + re.begin();

    double tx, ty, tz;
    for (int idx = re.begin(); idx < re.end(); ++idx, ++x, ++y, ++z, ++t, ++dx, ++dy, ++dz) {
      tx = *x, ty = *y, tz = *z;
      _Transformation->Displacement(tx, ty, tz, *t);
      *dx -= tx, *dy -= ty, *dz -= tz;
    }
  }
};

// -----------------------------------------------------------------------------
/// Body of Transformation::Transform(int, double *, double *, double *, ...)
class TransformPoints
{
  const Transformation *_Transformation;
  double               *_x, *_y, *_z, _t0, _t1;
  const double         *_t;

public:

  TransformPoints(const Transformation *transformation, double *x, double *y, double *z, double t, double t0)
  :
    _Transformation(transformation),
    _x(x), _y(y), _z(z), _t0(t0), _t1(t), _t(NULL)
  {}

  TransformPoints(const Transformation *transformation, double *x, double *y, double *z, const double *t, double t0)
  :
    _Transformation(transformation),
    _x(x), _y(y), _z(z), _t0(t0), _t1(.0), _t(t)
  {}

  void operator ()(const blocked_range<int> &idx) const
  {
    double *x = _x + idx.begin();
    double *y = _y + idx.begin();
    double *z = _z + idx.begin();
    for (int i = idx.begin(); i != idx.end(); ++i, ++x, ++y, ++z) {
      _Transformation->Transform(*x, *y, *z, (_t ? _t[i] : _t1), _t0);
    }
  }
};

// -----------------------------------------------------------------------------
/// Body of Transformation::Transform(WorldCoordsImage &)
class TransformWorldCoords
{
  const Transformation *_Transformation;
  WorldCoordsImage     *_Coords;
  const int             _NumberOfPoints;
  double                _t, _t0;

public:

  TransformWorldCoords(const Transformation *t, WorldCoordsImage &coords, double t0)
  :
    _Transformation(t),
    _Coords(&coords),
    _NumberOfPoints(coords.NumberOfSpatialVoxels()),
    _t(coords.GetTOrigin()), _t0(t0)
  {}

  void operator ()(const blocked_range<int> &idx) const
  {
    WorldCoordsImage::VoxelType *wx, *wy, *wz;
    wx = _Coords->Data() + idx.begin();
    wy = wx + _NumberOfPoints;
    wz = wy + _NumberOfPoints;
    for (int i = idx.begin(); i != idx.end(); ++i, ++wx, ++wy, ++wz) {
      _Transformation->Transform(*wx, *wy, *wz, _t, _t0);
    }
  }
};


} } // namespace mirtk::TransformationUtils

#endif // MIRTK_TransformationUtils_H
