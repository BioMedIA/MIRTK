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

#ifndef MIRTK_CSplineInterpolateImageFunction_HXX
#define MIRTK_CSplineInterpolateImageFunction_HXX

#include "mirtk/CSplineInterpolateImageFunction.h"

#include "mirtk/Math.h"
#include "mirtk/VoxelCast.h"
#include "mirtk/InterpolateImageFunction.hxx"


namespace mirtk {


// =============================================================================
// Cubic spline function
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::Real
GenericCSplineInterpolateImageFunction<TImage>::CSpline(Real x)
{
  const Real xi   = abs(x);
  const Real xii  = xi  * xi;
  const Real xiii = xii * xi;
  if      (xi < 1) return  xiii - 2 * xii + 1;
  else if (xi < 2) return -xiii + 5 * xii - 8 * xi + 4;
  return 0;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
GenericCSplineInterpolateImageFunction<TImage>
::GenericCSplineInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
GenericCSplineInterpolateImageFunction<TImage>
::~GenericCSplineInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericCSplineInterpolateImageFunction<TImage>::Initialize(bool coeff)
{
  // Initialize base class
  Superclass::Initialize(coeff);

  // Domain on which the C-spline interpolation is defined
  switch (this->NumberOfDimensions()) {
    case 4:
      this->_t1 = 1;
      this->_t2 = this->Input()->T() - 3;
    case 3:
      this->_z1 = 1;
      this->_z2 = this->Input()->Z() - 3;
    default:
      this->_y1 = 1;
      this->_y2 = this->Input()->Y() - 3;
      this->_x1 = 1;
      this->_x2 = this->Input()->X() - 3;
  }
}

// =============================================================================
// Domain checks
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
void GenericCSplineInterpolateImageFunction<TImage>
::BoundingInterval(double x, int &i, int &I) const
{
  i = ifloor(x) - 1, I = i + 3;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::Get2D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), wy, w;

  for (j = j1; j <= j2; ++j) {
    if (0 <= j && j < this->Input()->Y()) {
      wy = CSpline(Real(y - j));
      for (i = i1; i <= i2; ++i) {
        if (0 <= i && i < this->Input()->X()) {
          w    = CSpline(Real(x - i)) * wy;
          val += w * voxel_cast<RealType>(this->Input()->Get(i, j, k, l));
          nrm += w;
        }
      }
    }
  }

  if (nrm > 1e-3) val /= nrm;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::GetWithPadding2D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wy, w;

  for (j = j1; j <= j2; ++j) {
    wy = CSpline(Real(y - j));
    for (i = i1; i <= i2; ++i) {
      w = CSpline(Real(x - i)) * wy;
      if (this->Input()->IsInsideForeground(i, j, k, l)) {
        val += w * voxel_cast<RealType>(this->Input()->Get(i, j, k, l));
        fgw += w;
      } else {
        bgw += w;
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val / fgw);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::Get2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), wy, w;

  for (j = j1; j <= j2; ++j) {
    wy = CSpline(Real(y - j));
    for (i = i1; i <= i2; ++i) {
      w    = CSpline(Real(x - i)) * wy;
      val += w * voxel_cast<RealType>(input->Get(i, j, k, l));
      nrm += w;
    }
  }

  if (nrm > 1e-3) val /= nrm;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::GetWithPadding2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wy, w;

  for (j = j1; j <= j2; ++j) {
    wy = CSpline(Real(y - j));
    for (i = i1; i <= i2; ++i) {
      w    = CSpline(Real(x - i)) * wy;
      if (input->IsForeground(i, j, k, l)) {
        val += w * voxel_cast<RealType>(input->Get(i, j, k, l));
        fgw += w;
      } else {
        bgw += w;
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val / fgw);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::Get3D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  if (l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int k1 = k - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;
  const int k2 = k + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), wz, wyz, w;

  for (k = k1; k <= k2; ++k) {
    if (0 <= k && k < this->Input()->Z()) {
      wz = CSpline(Real(z - k));
      for (j = j1; j <= j2; ++j) {
        if (0 <= j && j < this->Input()->Y()) {
          wyz = CSpline(Real(y - j)) * wz;
          for (i = i1; i <= i2; ++i) {
            if (0 <= i && i < this->Input()->X()) {
              w    = CSpline(Real(x - i)) * wyz;
              val += w * voxel_cast<RealType>(this->Input()->Get(i, j, k, l));
              nrm += w;
            }
          }
        }
      }
    }
  }

  if (nrm > 1e-3) val /= nrm;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::GetWithPadding3D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  if (l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int k1 = k - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;
  const int k2 = k + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wz, wyz, w;

  for (k = k1; k <= k2; ++k) {
    wz = CSpline(Real(z - k));
    for (j = j1; j <= j2; ++j) {
      wyz = CSpline(Real(y - j)) * wz;
      for (i = i1; i <= i2; ++i) {
        w = CSpline(Real(x - i)) * wyz;
        if (this->Input()->IsInsideForeground(i, j, k, l)) {
          val += w * voxel_cast<RealType>(this->Input()->Get(i, j, k, l));
          fgw += w;
        } else {
          bgw += w;
        }
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val / fgw);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::Get3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int k1 = k - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;
  const int k2 = k + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), wz, wyz, w;

  for (k = k1; k <= k2; ++k) {
    wz = CSpline(Real(z - k));
    for (j = j1; j <= j2; ++j) {
      wyz = CSpline(Real(y - j)) * wz;
      for (i = i1; i < i2; ++i) {
        w    = CSpline(Real(x - i)) * wyz;
        val += w * voxel_cast<RealType>(input->Get(i, j, k, l));
        nrm += w;
      }
    }
  }

  if (nrm > 1e-3) val /= nrm;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::GetWithPadding3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int k1 = k - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;
  const int k2 = k + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wz, wyz, w;

  for (k = k1; k <= k2; ++k) {
    wz = CSpline(Real(z - k));
    for (j = j1; j <= j2; ++j) {
      wyz = CSpline(Real(y - j)) * wz;
      for (i = i1; i < i2; ++i) {
        w = CSpline(Real(x - i)) * wyz;
        if (input->IsForeground(i, j, k, l)) {
          val += w * voxel_cast<RealType>(input->Get(i, j, k, l));
          fgw += w;
        } else {
          bgw += w;
        }
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val / fgw);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::Get4D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int k1 = k - 1;
  const int l1 = l - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;
  const int k2 = k + 2;
  const int l2 = l + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), wt, wzt, wyzt, w;

  for (l = l1; l <= l2; ++l) {
    if (0 <= l && l < this->Input()->T()) {
      wt = CSpline(Real(t - l));
      for (k = k1; k <= k2; ++k) {
        if (0 <= k && k < this->Input()->Z()) {
          wzt = CSpline(Real(z - k)) * wt;
          for (j = j1; j <= j2; ++j) {
            if (0 <= j && j < this->Input()->Y()) {
              wyzt = CSpline(Real(y - j)) * wzt;
              for (i = i1; i <= i2; ++i) {
                if (0 <= i && i < this->Input()->X()) {
                  w    = CSpline(Real(x - i)) * wyzt;
                  val += w * voxel_cast<RealType>(this->Input()->Get(i, j, k, l));
                  nrm += w;
                }
              }
            }
          }
        }
      }
    }
  }

  if (nrm > 1e-3) val /= nrm;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::GetWithPadding4D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int k1 = k - 1;
  const int l1 = l - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;
  const int k2 = k + 2;
  const int l2 = l + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wt, wzt, wyzt, w;

  for (l = l1; l <= l2; ++l) {
    wt = CSpline(Real(t - l));
    for (k = k1; k <= k2; ++k) {
      wzt = CSpline(Real(z - k)) * wt;
      for (j = j1; j <= j2; ++j) {
        wyzt = CSpline(Real(y - j)) * wzt;
        for (i = i1; i <= i2; ++i) {
          w = CSpline(Real(x - i)) * wyzt;
          if (this->Input()->IsInsideForeground(i, j, k, l)) {
            val += w * voxel_cast<RealType>(this->Input()->Get(i, j, k, l));
            fgw += w;
          } else {
            bgw += w;
          }
        }
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val / fgw);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::Get4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int k1 = k - 1;
  const int l1 = l - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;
  const int k2 = k + 2;
  const int l2 = l + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), wt, wzt, wyzt, w;

  for (l = l1; l <= l2; ++l) {
    wt = CSpline(Real(t - l));
    for (k = k1; k <= k2; ++k) {
      wzt = CSpline(Real(z - k)) * wt;
      for (j = j1; j <= j2; ++j) {
        wyzt = CSpline(Real(y - j)) * wzt;
        for (i = i1; i <= i2; ++i) {
          w    = CSpline(Real(x - i)) * wyzt;
          val += w * voxel_cast<RealType>(input->Get(i, j, k, l));
          nrm += w;
        }
      }
    }
  }

  if (nrm > 1e-3) val /= nrm;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::GetWithPadding4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int i1 = i - 1;
  const int j1 = j - 1;
  const int k1 = k - 1;
  const int l1 = l - 1;
  const int i2 = i + 2;
  const int j2 = j + 2;
  const int k2 = k + 2;
  const int l2 = l + 2;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wt, wzt, wyzt, w;

  for (l = l1; l <= l2; ++l) {
    wt = CSpline(Real(t - l));
    for (k = k1; k <= k2; ++k) {
      wzt = CSpline(Real(z - k)) * wt;
      for (j = j1; j <= j2; ++j) {
        wyzt = CSpline(Real(y - j)) * wzt;
        for (i = i1; i <= i2; ++i) {
          w = CSpline(Real(x - i)) * wyzt;
          if (input->IsForeground(i, j, k, l)) {
            val += w * voxel_cast<RealType>(input->Get(i, j, k, l));
            fgw += w;
          } else {
            bgw += w;
          }
        }
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val / fgw);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
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
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
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
inline typename TOtherImage::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
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
inline typename TOtherImage::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
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
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::GetInside(double x, double y, double z, double t) const
{
  return Get(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
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
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return GetWithPadding(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericCSplineInterpolateImageFunction<TImage>::VoxelType
GenericCSplineInterpolateImageFunction<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return GetWithPadding(this->Extrapolator(), x, y, z, t);
  } else {
    return GetWithPadding(x, y, z, t);
  }
}


} // namespace mirtk

#endif // MIRTK_CSplineInterpolateImageFunction_HXX
