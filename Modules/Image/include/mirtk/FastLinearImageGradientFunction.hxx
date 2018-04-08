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

#ifndef MIRTK_FastLinearImageGradientFunction_HXX
#define MIRTK_FastLinearImageGradientFunction_HXX

#include "mirtk/FastLinearImageGradientFunction.h"
#include "mirtk/ImageGradientFunction.hxx"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
GenericFastLinearImageGradientFunction<TImage>
::GenericFastLinearImageGradientFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
GenericFastLinearImageGradientFunction<TImage>
::~GenericFastLinearImageGradientFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastLinearImageGradientFunction<TImage>::Initialize(bool coeff)
{
  // Initialize base class
  Superclass::Initialize(coeff);

  // Domain on which linear interpolation is defined: [0, x)
  switch (this->NumberOfDimensions()) {
    case 4:
      this->_t1 = 0.;
      this->_t2 = fdec(this->Input()->T() - 1);
    case 3:
      this->_z1 = 0.;
      this->_z2 = fdec(this->Input()->Z() - 1);
    default:
      this->_y1 = 0.;
      this->_y2 = fdec(this->Input()->Y() - 1);
      this->_x1 = 0.;
      this->_x2 = fdec(this->Input()->X() - 1);
  }
}

// =============================================================================
// Domain checks
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastLinearImageGradientFunction<TImage>
::BoundingInterval(double x, int &i, int &I) const
{
  i = ifloor(x), I = i + 1;
}

// =============================================================================
// Interpolation weights
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline int GenericFastLinearImageGradientFunction<TImage>
::ComputeWeights(double x, Real w[2])
{
  const int i = ifloor(x);
  w[1] = static_cast<Real>(x) - static_cast<Real>(i);
  w[0] = static_cast<Real>(1) - w[1];
  return i;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::Get2D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = iround(z);
  const int l = iround(t);

  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return this->DefaultValue();
  }

  GradientType val(.0);
  Real         nrm(0), coeff;

  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    if (0 <= jb && jb < this->Input()->Y()) {
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        if (0 <= ia && ia < this->Input()->X()) {
          coeff = voxel_cast<Real>(this->Input()->Get(ia, jb, k, l));
          val._x += wd[a] * wy[b] * coeff;
          val._y += wx[a] * wd[b] * coeff;
          nrm    += wx[a] * wy[b];
        }
      }
    }
  }

  if (nrm > 1e-3) val /= nrm;
  else val = this->DefaultValue();

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPadding2D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = iround(z);
  const int l = iround(t);

  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return this->DefaultValue();
  }

  GradientType val(.0);
  Real         coeff;

  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    for (int a = 0; a <= 1; ++a) {
      ia = i + a;
      if (this->Input()->IsInsideForeground(ia, jb, k, l)) {
        coeff = voxel_cast<Real>(this->Input()->Get(ia, jb, k, l));
        val._x += wd[a] * wy[b] * coeff;
        val._y += wx[a] * wd[b] * coeff;
      } else {
        return this->DefaultValue();
      }
    }
  }

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::Get2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = iround(z);
  const int l = iround(t);

  GradientType val(.0);
  Real         coeff;

  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    for (int a = 0; a <= 1; ++a) {
      ia = i + a;
      coeff = voxel_cast<Real>(input->Get(ia, jb, k, l));
      val._x += wd[a] * wy[b] * coeff;
      val._y += wx[a] * wd[b] * coeff;
    }
  }

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPadding2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = iround(z);
  const int l = iround(t);

  if (k < 0 || k >= input->Z() ||
      l < 0 || l >= input->T()) {
    return this->DefaultValue();
  }

  GradientType val(.0);
  Real         coeff;

  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    for (int a = 0; a <= 1; ++a) {
      ia = i + a;
      if (input->IsForeground(ia, jb, k, l)) {
        coeff = voxel_cast<Real>(input->Get(ia, jb, k, l));
        val._x += wd[a] * wy[b] * coeff;
        val._y += wx[a] * wd[b] * coeff;
      } else {
        return this->DefaultValue();
      }
    }
  }

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::Get3D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = iround(t);

  if (l < 0 || l >= this->Input()->T()) {
    return this->DefaultValue();
  }

  GradientType val(.0);
  Real         nrm(0), coeff;

  int ia, jb, kc;
  for (int c = 0; c <= 1; ++c) {
    kc = k + c;
    if (0 <= kc && kc < this->Input()->Z()) {
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        if (0 <= jb && jb < this->Input()->Y()) {
          for (int a = 0; a <= 1; ++a) {
            ia = i + a;
            if (0 <= ia && ia < this->Input()->X()) {
              coeff = voxel_cast<Real>(this->Input()->Get(ia, jb, kc, l));
              val._x += wd[a] * wy[b] * wz[c] * coeff;
              val._y += wx[a] * wd[b] * wz[c] * coeff;
              val._z += wx[a] * wy[b] * wd[c] * coeff;
              nrm    += wx[a] * wy[b] * wz[c];
            }
          }
        }
      }
    }
  }

  if (nrm > 1e-3) val /= nrm;
  else val = this->DefaultValue();

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPadding3D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = iround(t);

  if (l < 0 || l >= this->Input()->T()) {
    return this->DefaultValue();
  }

  GradientType val(.0);
  Real         coeff;

  int ia, jb, kc;
  for (int c = 0; c <= 1; ++c) {
    kc = k + c;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        if (this->Input()->IsInsideForeground(ia, jb, kc, l)) {
          coeff = voxel_cast<Real>(this->Input()->Get(ia, jb, kc, l));
          val._x += wd[a] * wy[b] * wz[c] * coeff;
          val._y += wx[a] * wd[b] * wz[c] * coeff;
          val._z += wx[a] * wy[b] * wd[c] * coeff;
        } else {
          return this->DefaultValue();
        }
      }
    }
  }

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::Get3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = iround(t);

  GradientType val(.0);
  Real         coeff;

  int ia, jb, kc;
  for (int c = 0; c <= 1; ++c) {
    kc = k + c;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        coeff = voxel_cast<Real>(input->Get(ia, jb, kc, l));
        val._x += wd[a] * wy[b] * wz[c] * coeff;
        val._y += wx[a] * wd[b] * wz[c] * coeff;
        val._z += wx[a] * wy[b] * wd[c] * coeff;
      }
    }
  }

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPadding3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = iround(t);

  if (l < 0 || l >= input->T()) {
    return this->DefaultValue();
  }

  GradientType val(.0);
  Real         coeff;

  int ia, jb, kc;
  for (int c = 0; c <= 1; ++c) {
    kc = k + c;
    for (int b = 0; b <= 1; ++b) {
      jb  = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        if (input->IsForeground(ia, jb, kc, l)) {
          coeff = voxel_cast<Real>(input->Get(ia, jb, kc, l));
          val._x += wd[a] * wy[b] * wz[c] * coeff;
          val._y += wx[a] * wd[b] * wz[c] * coeff;
          val._z += wx[a] * wy[b] * wd[c] * coeff;
        } else {
          return this->DefaultValue();
        }
      }
    }
  }

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::Get4D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wt[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(t, wt);

  GradientType val(.0);
  Real         nrm(0), coeff;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    if (0 <= ld && ld < this->Input()->T()) {
      for (int c = 0; c <= 1; ++c) {
        kc = k + c;
        if (0 <= kc && kc < this->Input()->Z()) {
          for (int b = 0; b <= 1; ++b) {
            jb = j + b;
            if (0 <= jb && jb < this->Input()->Y()) {
              for (int a = 0; a <= 1; ++a) {
                ia = i + a;
                if (0 <= ia && ia < this->Input()->X()) {
                  coeff = voxel_cast<Real>(this->Input()->Get(ia, jb, kc, ld));
                  val._x += wd[a] * wy[b] * wz[c] * wt[d] * coeff;
                  val._y += wx[a] * wd[b] * wz[c] * wt[d] * coeff;
                  val._z += wx[a] * wy[b] * wd[c] * wt[d] * coeff;
                  nrm += wx[a] * wy[b] * wz[c] * wt[d];
                }
              }
            }
          }
        }
      }
    }
  }

  if (nrm > 1e-3) val /= nrm;
  else val = this->DefaultValue();

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPadding4D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wt[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(t, wt);

  GradientType val(.0);
  Real         coeff;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          if (this->Input()->IsInsideForeground(ia, jb, kc, ld)) {
            coeff = voxel_cast<Real>(this->Input()->Get(ia, jb, kc, ld));
            val._x += wd[a] * wy[b] * wz[c] * wt[d] * coeff;
            val._y += wx[a] * wd[b] * wz[c] * wt[d] * coeff;
            val._z += wx[a] * wy[b] * wd[c] * wt[d] * coeff;
          } else {
            return this->DefaultValue();
          }
        }
      }
    }
  }

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::Get4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wt[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(t, wt);

  GradientType val(.0);
  Real         coeff;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          coeff = voxel_cast<Real>(input->Get(ia, jb, kc, ld));
          val._x += wd[a] * wy[b] * wz[c] * wt[d] * coeff;
          val._y += wx[a] * wd[b] * wz[c] * wt[d] * coeff;
          val._z += wx[a] * wy[b] * wd[c] * wt[d] * coeff;
        }
      }
    }
  }

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPadding4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wt[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(t, wt);

  GradientType val(.0);
  Real         coeff;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          if (input->IsForeground(ia, jb, kc, ld)) {
            coeff = voxel_cast<Real>(input->Get(ia, jb, kc, ld));
            val._x += wd[a] * wy[b] * wz[c] * wt[d] * coeff;
            val._y += wx[a] * wd[b] * wz[c] * wt[d] * coeff;
            val._z += wx[a] * wy[b] * wd[c] * wt[d] * coeff;
          } else {
            return this->DefaultValue();
          }
        }
      }
    }
  }

  return val;
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::Get(double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return Get3D(x, y, z, t);
    case 2:  return Get2D(x, y, z, t);
    default: return Get4D(x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPadding(double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return GetWithPadding3D(x, y, z, t);
    case 2:  return GetWithPadding2D(x, y, z, t);
    default: return GetWithPadding4D(x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::Get(const TOtherImage *input, double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return Get3D(input, x, y, z, t);
    case 2:  return Get2D(input, x, y, z, t);
    default: return Get4D(input, x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPadding(const TOtherImage *input, double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return GetWithPadding3D(input, x, y, z, t);
    case 2:  return GetWithPadding2D(input, x, y, z, t);
    default: return GetWithPadding4D(input, x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetInside(double x, double y, double z, double t) const
{
  return Get(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return Get(this->Extrapolator(), x, y, z, t);
  } else {
    return Get(x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return GetWithPadding(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastLinearImageGradientFunction<TImage>::GradientType
GenericFastLinearImageGradientFunction<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return GetWithPadding(this->Extrapolator(), x, y, z, t);
  } else {
    return GetWithPadding(x, y, z, t);
  }
}


} // namespace mirtk

#endif // MIRTK_FastLinearImageGradientFunction_HXX
