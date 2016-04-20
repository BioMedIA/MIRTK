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

#ifndef MIRTK_NaryVoxelFunction_H
#define MIRTK_NaryVoxelFunction_H

#include "mirtk/VoxelFunction.h"
#include "mirtk/BaseImage.h"
#include "mirtk/InterpolateImageFunction.h"


namespace mirtk {


/**
 * Namespace containing voxel functions with different number of arguments
 *
 * The voxel functions defined in this namespace may provide different overloads
 * of the operator() for different number of image arguments to the ForEachVoxel
 * template functions.
 *
 */
namespace NaryVoxelFunction {


// =============================================================================
// Stationary velocity fields
// =============================================================================

// -----------------------------------------------------------------------------
/// Voxel function for exponentiation of 2D velocity field
template <class TInterpolator = InterpolateImageFunction>
struct ExpVelocityFieldEuler2D : public VoxelFunction
{
  const BaseImage *_VelocityField;        ///< Velocity field
  TInterpolator   *_VelocityInterpolator; ///< Velocity field interpolator
  int              _NumberOfSteps;        ///< Number of integration steps
 
  double           _dt;    ///< Temporal step length
  static const int _x = 0; ///< Offset of x component
  int              _y;     ///< Offset of y component
  
  
  ExpVelocityFieldEuler2D(TInterpolator *v, int n = 100, double t = 1.0)
  :
    _VelocityField(v->Input()),
    _VelocityInterpolator (v),
    _NumberOfSteps(n),
    _dt(t / _NumberOfSteps),
    _y(_VelocityField->GetX() * _VelocityField->GetY())
  {}

  template <class T>
  void operator ()(int i, int j, int k, int, T *d)
  {
    double x, y, z, v[3]; // 2D vector field could have constant third component!
    double x1 = i, y1 = j, z1 = k;
    _VelocityField->ImageToWorld(x1, y1, z1);
    double x2 = x1, y2 = y1;
    for (int s = 0; s < _NumberOfSteps; ++s) {
      x = x2, y = y2, z = z1;
      _VelocityField      ->WorldToImage(x, y, z);
      _VelocityInterpolator->Evaluate(v, x, y, z);
      x2 += v[0] * _dt;
      y2 += v[1] * _dt;
    }
    d[_x] = static_cast<T>(x2 - x1);
    d[_y] = static_cast<T>(y2 - y1);
  }
  
  template <class T>
  void operator ()(int i, int j, int k, int, const T *din, T *dout)
  {
    double x, y, z, v[3]; // 2D vector field could have constant third component!
    double x1 = i, y1 = j, z1 = k;
    _VelocityField->ImageToWorld(x1, y1, z1);
    double x2 = x1 + din[_x];
    double y2 = y1 + din[_y];
    for (int s = 0; s < _NumberOfSteps; ++s) {
      x = x2, y = y2, z = z1;
      _VelocityField      ->WorldToImage(x, y, z);
      _VelocityInterpolator->Evaluate(v, x, y, z);
      x2 += v[0] * _dt;
      y2 += v[1] * _dt;
    }
    dout[_x] = static_cast<T>(x2 - x1);
    dout[_y] = static_cast<T>(y2 - y1);
  }
};

// -----------------------------------------------------------------------------
/// Voxel function for exponentiation of 3D velocity field
template <class TInterpolator = InterpolateImageFunction>
struct ExpVelocityFieldEuler3D : public VoxelFunction
{
  const BaseImage *_VelocityField;        ///< Velocity field
  TInterpolator   *_VelocityInterpolator; ///< Velocity field interpolator
  int              _NumberOfSteps;        ///< Number of integration steps
 
  double           _dt;    ///< Temporal step length
  static const int _x = 0; ///< Offset of x component
  int              _y, _z; ///< Offset of y and z components
  
  
  ExpVelocityFieldEuler3D(TInterpolator *v, int n = 100, double t = 1.0)
  :
    _VelocityField(v->Input()),
    _VelocityInterpolator (v),
    _NumberOfSteps(n),
    _dt(t / _NumberOfSteps),
    _y(_VelocityField->X() * _VelocityField->Y() * _VelocityField->Z()),
    _z(2 * _y)
  {}

  template <class T>
  void operator ()(int i, int j, int k, int, T *d)
  {
    double x, y, z, v[3];
    double x1 = i, y1 = j, z1 = k;
    _VelocityField->ImageToWorld(x1, y1, z1);
    double x2 = x1, y2 = y1, z2 = z1;
    for (int s = 0; s < _NumberOfSteps; ++s) {
      x = x2, y = y2, z = z2;
      _VelocityField      ->WorldToImage(x, y, z);
      _VelocityInterpolator->Evaluate(v, x, y, z);
      x2 += v[0] * _dt;
      y2 += v[1] * _dt;
      z2 += v[2] * _dt;
    }
    d[_x] = static_cast<T>(x2 - x1);
    d[_y] = static_cast<T>(y2 - y1);
    d[_z] = static_cast<T>(z2 - z1);
  }
  
  template <class T>
  void operator ()(int i, int j, int k, int, const T *din, T *dout)
  {
    double x, y, z, v[3];
    double x1 = i, y1 = j, z1 = k;
    _VelocityField->ImageToWorld(x1, y1, z1);
    double x2 = x1 + din[_x];
    double y2 = y1 + din[_y];
    double z2 = z1 + din[_z];
    for (int s = 0; s < _NumberOfSteps; ++s) {
      x = x2, y = y2, z = z2;
      _VelocityField      ->WorldToImage(x, y, z);
      _VelocityInterpolator->Evaluate(v, x, y, z);
      x2 += v[0] * _dt;
      y2 += v[1] * _dt;
      z2 += v[2] * _dt;
    }
    dout[_x] = static_cast<T>(x2 - x1);
    dout[_y] = static_cast<T>(y2 - y1);
    dout[_z] = static_cast<T>(z2 - z1);
  }
};

// -----------------------------------------------------------------------------
/// Compute voxel-wise weighted sum
struct VoxelWiseWeightedSum : public VoxelFunction
{
  const double *_Weight;

  template <class TImage, class T>
  void operator()(const TImage &, int, const T *a, const T *b, T *out)
  {
    (*out) = static_cast<T>(_Weight[0] * (*a) + _Weight[1] * (*b));
  }

  template <class TImage, class T>
  void operator()(const TImage &, int, const T *a, const T *b, const T *c, T *out)
  {
    (*out) = static_cast<T>(_Weight[0] * (*a) + _Weight[1] * (*b) + _Weight[2] * (*c));
  }

  template <class TImage, class T>
  void operator()(const TImage &, int, const T *a, const T *b, const T *c, const T *d, T *out)
  {
    (*out) = static_cast<T>(_Weight[0] * (*a) + _Weight[1] * (*b) + _Weight[2] * (*c) + _Weight[3] * (*d));
  }

  template <class TImage, class T>
  void operator()(const TImage &, int, const T *a, const T *b, const T *c, const T *d, const T *e, T *out)
  {
    (*out) = static_cast<T>(_Weight[0] * (*a) + _Weight[1] * (*b) + _Weight[2] * (*c) + _Weight[3] * (*d) + _Weight[4] * (*e));
  }

  template <class TImage, class T>
  void operator()(const TImage &, int, const T *a, const T *b, const T *c, const T *d, const T *e, const T *f, T *out)
  {
    (*out) = static_cast<T>(_Weight[0] * (*a) + _Weight[1] * (*b) + _Weight[2] * (*c) + _Weight[3] * (*d) + _Weight[4] * (*e) + _Weight[5] * (*f));
  }

  template <class TImage, class T>
  void operator()(const TImage &, int, const T *a, const T *b, const T *c, const T *d, const T *e, const T *f, const T *g, T *out)
  {
    (*out) = static_cast<T>(_Weight[0] * (*a) + _Weight[1] * (*b) + _Weight[2] * (*c) + _Weight[3] * (*d) + _Weight[4] * (*e) + _Weight[5] * (*f) + _Weight[6] * (*g));
  }
};

// -----------------------------------------------------------------------------
/// Evaluate Baker-Campbell-Hausdorff (BCH) formula
struct EvaluateBCHFormula : public VoxelFunction
{
  template <class TImage, class T>
  void operator()(const TImage &, int, const T *v, const T *d, T *out)
  {
    (*out) = (*v) + (*d);
  }
  
  template <class TImage, class T>
  void operator()(const TImage &, int, const T *v, const T *d, const T *l1, T *out)
  {
    (*out) = static_cast<T>((*v) + (*d) + (*l1) / 2.0);
  }
  
  template <class TImage, class T>
  void operator()(const TImage &, int, const T *v, const T *d, const T *l1, const T *l2, T *out)
  {
    (*out) = static_cast<T>((*v) + (*d) + (*l1) / 2.0 + (*l2) / 12.0);
  }
  
  template <class TImage, class T>
  void operator()(const TImage &, int, const T *v, const T *d, const T *l1, const T *l2, const T *l3, T *out)
  {
    (*out) = static_cast<T>((*v) + (*d) + (*l1) / 2.0 + (*l2) / 12.0 - (*l3) / 12.0);
  }
  
  template <class TImage, class T>
  void operator()(const TImage &, int, const T *v, const T *d, const T *l1, const T *l2, const T *l3, const T *l4, T *out)
  {
    (*out) = static_cast<T>((*v) + (*d) + (*l1) / 2.0 + (*l2) / 12.0 - (*l3) / 12.0 - (*l4) / 24.0);
  }
};

// -----------------------------------------------------------------------------
/// Evaluate velocity update using Baker-Campbell-Hausdorff (BCH) formula
struct EvaluateBCHUpdate : public VoxelFunction
{
  template <class TImage, class T>
  void operator()(const TImage &, int, const T *d, const T *l1, T *out)
  {
    (*out) = static_cast<T>((*d) + (*l1) / 2.0);
  }

  template <class TImage, class T>
  void operator()(const TImage &, int, const T *d, const T *l1, const T *l2, T *out)
  {
    (*out) = static_cast<T>((*d) + (*l1) / 2.0 + (*l2) / 12.0);
  }

  template <class TImage, class T>
  void operator()(const TImage &, int, const T *d, const T *l1, const T *l2, const T *l3, T *out)
  {
    (*out) = static_cast<T>((*d) + (*l1) / 2.0 + (*l2) / 12.0 - (*l3) / 12.0);
  }

  template <class TImage, class T>
  void operator()(const TImage &, int, const T *d, const T *l1, const T *l2, const T *l3, const T *l4, T *out)
  {
    (*out) = static_cast<T>((*d) + (*l1) / 2.0 + (*l2) / 12.0 - (*l3) / 12.0 - (*l4) / 24.0);
  }
};


} } // namespace mirtk::NaryVoxelFunction

#endif // MIRTK_NaryVoxelFunction_H
