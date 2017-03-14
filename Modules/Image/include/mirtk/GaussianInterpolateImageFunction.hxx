/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#ifndef MIRTK_GaussianInterpolateImageFunction_HXX
#define MIRTK_GaussianInterpolateImageFunction_HXX

#include "mirtk/GaussianInterpolateImageFunction.h"
#include "mirtk/InterpolateImageFunction.hxx"

#include "mirtk/Math.h"
#include "mirtk/VoxelCast.h"
#include "mirtk/ScalarGaussian.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
GenericGaussianInterpolateImageFunction<TImage>
::GenericGaussianInterpolateImageFunction(double sigma)
:
  _Sigma(sigma),
  _RadiusX(.0), _RadiusY(.0), _RadiusZ(.0), _RadiusT(.0),
  _dx     (.0), _dy     (.0), _dz     (.0), _dt     (.0)
{
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericGaussianInterpolateImageFunction<TImage>::Initialize(bool coeff)
{
  // Initialize base class
  Superclass::Initialize(coeff);

  // Get input voxel size
  this->Input()->GetPixelSize(_dx, _dy, _dz, _dt);

  // Compute filter extent
  _RadiusX = _RadiusY = _RadiusZ = _RadiusT = .0;
  switch (this->NumberOfDimensions()) {
    case 4:  _RadiusT = round(ceil(3.0 * _Sigma/_dt));
    case 3:  _RadiusZ = round(ceil(3.0 * _Sigma/_dz));
    default: _RadiusY = round(ceil(3.0 * _Sigma/_dy));
             _RadiusX = round(ceil(3.0 * _Sigma/_dx));
  }

  // Domain on which the C-spline interpolation is defined
  this->_x1 = _RadiusX;
  this->_x2 = this->Input()->X() - _RadiusX - 1;
  this->_y1 = _RadiusY;
  this->_y2 = this->Input()->Y() - _RadiusY - 1;
  this->_z1 = _RadiusZ;
  this->_z2 = this->Input()->Z() - _RadiusZ - 1;
  this->_t1 = _RadiusT;
  this->_t2 = this->Input()->T() - _RadiusT - 1;
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericGaussianInterpolateImageFunction<TImage>
::BoundingInterval(double x, int &i, int &I) const
{
  double r = max(_RadiusX, max(_RadiusY, max(_RadiusZ, _RadiusT)));
  i = ifloor(x - r), I = ifloor(x + r);
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericGaussianInterpolateImageFunction<TImage>
::BoundingBox(double x, double y, int &i1, int &i2,
                                  int &j1, int &j2) const
{
  i1 = ifloor(x - _RadiusX), i2 = ifloor(x + _RadiusX);
  j1 = ifloor(y - _RadiusY), j2 = ifloor(y + _RadiusY);
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericGaussianInterpolateImageFunction<TImage>
::BoundingBox(double x, double y, double z, int &i1, int &i2, int &k1,
                                            int &j1, int &j2, int &k2) const
{
  i1 = ifloor(x - _RadiusX), i2 = ifloor(x + _RadiusX);
  j1 = ifloor(y - _RadiusY), j2 = ifloor(y + _RadiusY);
  k1 = ifloor(z - _RadiusZ), k2 = ifloor(z + _RadiusZ);
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericGaussianInterpolateImageFunction<TImage>
::BoundingBox(double x, double y, double z, double t,
              int &i1, int &i2, int &k1, int &l1,
              int &j1, int &j2, int &k2, int &l2) const
{
  i1 = ifloor(x - _RadiusX), i2 = ifloor(x + _RadiusX);
  j1 = ifloor(y - _RadiusY), j2 = ifloor(y + _RadiusY);
  k1 = ifloor(z - _RadiusZ), k2 = ifloor(z + _RadiusZ);
  l1 = ifloor(t - _RadiusT), l2 = ifloor(t + _RadiusT);
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
::Get2D(double x, double y, double z, double t) const
{
  const int k = iround(z);
  const int l = iround(t);

  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, .0, x, y, .0);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), w;

  for (int j = j1; j <= j2; ++j) {
    if (0 <= j && j < this->Input()->Y()) {
      for (int i = i1; i <= i2; ++i) {
        if (0 <= i && i < this->Input()->X()) {
          w   = static_cast<Real>(kernel.Evaluate(i, j));
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
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
::GetWithPadding2D(double x, double y, double z, double t) const
{
  const int k = iround(z);
  const int l = iround(t);

  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, .0, x, y, .0);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  for (int j = j1; j <= j2; ++j) {
    for (int i = i1; i <= i2; ++i) {
      w = static_cast<Real>(kernel.Evaluate(i, j));
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
GenericGaussianInterpolateImageFunction<TImage>
::Get2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int k = iround(z);
  const int l = iround(t);

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, .0, x, y, .0);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), w;

  for (int j = j1; j <= j2; ++j) {
    for (int i = i1; i <= i2; ++i) {
      w   = static_cast<Real>(kernel.Evaluate(i, j));
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
GenericGaussianInterpolateImageFunction<TImage>
::GetWithPadding2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int k = iround(z);
  const int l = iround(t);

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, .0, x, y, .0);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  for (int j = j1; j <= j2; ++j) {
    for (int i = i1; i <= i2; ++i) {
      w = static_cast<Real>(kernel.Evaluate(i, j));
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
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
::Get3D(double x, double y, double z, double t) const
{
  const int l = iround(t);

  if (l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int k1 = ifloor(z - _RadiusZ);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);
  const int k2 = ifloor(z + _RadiusZ);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, _Sigma/_dz, x, y, z);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), w;

  for (int k = k1; k <= k2; ++k) {
    if (0 <= k && k < this->Input()->Z()) {
      for (int j = j1; j <= j2; ++j) {
        if (0 <= j && j < this->Input()->Y()) {
          for (int i = i1; i <= i2; ++i) {
            if (0 <= i && i < this->Input()->X()) {
              w   = static_cast<Real>(kernel.Evaluate(i, j, k));
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
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
::GetWithPadding3D(double x, double y, double z, double t) const
{
  const int l = iround(t);

  if (l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int k1 = ifloor(z - _RadiusZ);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);
  const int k2 = ifloor(z + _RadiusZ);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, _Sigma/_dz, x, y, z);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  for (int k = k1; k <= k2; ++k) {
    for (int j = j1; j <= j2; ++j) {
      for (int i = i1; i <= i2; ++i) {
        w = static_cast<Real>(kernel.Evaluate(i, j, k));
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
GenericGaussianInterpolateImageFunction<TImage>
::Get3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int l = iround(t);

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int k1 = ifloor(z - _RadiusZ);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);
  const int k2 = ifloor(z + _RadiusZ);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, _Sigma/_dz, x, y, z);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), w;

  for (int k = k1; k <= k2; ++k) {
    for (int j = j1; j <= j2; ++j) {
      for (int i = i1; i <= i2; ++i) {
        w   = static_cast<Real>(kernel.Evaluate(i, j, k));
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
GenericGaussianInterpolateImageFunction<TImage>
::GetWithPadding3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int l = iround(t);

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int k1 = ifloor(z - _RadiusZ);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);
  const int k2 = ifloor(z + _RadiusZ);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, _Sigma/_dz, x, y, z);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  for (int k = k1; k <= k2; ++k) {
    for (int j = j1; j <= j2; ++j) {
      for (int i = i1; i <= i2; ++i) {
        w = static_cast<Real>(kernel.Evaluate(i, j, k));
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
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
::Get4D(double x, double y, double z, double t) const
{
  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int k1 = ifloor(z - _RadiusZ);
  const int l1 = ifloor(t - _RadiusT);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);
  const int k2 = ifloor(z + _RadiusZ);
  const int l2 = ifloor(t + _RadiusT);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, _Sigma/_dz, _Sigma/_dt, x, y, z, t);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), w;

  for (int l = l1; l <= l2; ++l) {
    if (0 <= l && l < this->Input()->T()) {
      for (int k = k1; k <= k2; ++k) {
        if (0 <= k && k < this->Input()->Z()) {
          for (int j = j1; j <= j2; ++j) {
            if (0 <= j && j < this->Input()->Y()) {
              for (int i = i1; i <= i2; ++i) {
                if (0 <= i && i < this->Input()->X()) {
                  w   = static_cast<Real>(kernel.Evaluate(i, j, k, l));
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
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
::GetWithPadding4D(double x, double y, double z, double t) const
{
  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int k1 = ifloor(z - _RadiusZ);
  const int l1 = ifloor(t - _RadiusT);
  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);
  const int k2 = ifloor(z + _RadiusZ);
  const int l2 = ifloor(t + _RadiusT);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, _Sigma/_dz, _Sigma/_dt, x, y, z, t);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  for (int l = l1; l <= l2; ++l) {
    for (int k = k1; k <= k2; ++k) {
      for (int j = j1; j <= j2; ++j) {
        for (int i = i1; i <= i2; ++i) {
          w = static_cast<Real>(kernel.Evaluate(i, j, k, l));
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
GenericGaussianInterpolateImageFunction<TImage>
::Get4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int k1 = ifloor(z - _RadiusZ);
  const int l1 = ifloor(t - _RadiusT);

  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);
  const int k2 = ifloor(z + _RadiusZ);
  const int l2 = ifloor(t + _RadiusT);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, _Sigma/_dz, _Sigma/_dt, x, y, z, t);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), w;

  for (int l = l1; l <= l2; ++l) {
    for (int k = k1; k <= k2; ++k) {
      for (int j = j1; j <= j2; ++j) {
        for (int i = i1; i <= i2; ++i) {
          w   = static_cast<Real>(kernel.Evaluate(i, j, k, l));
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
GenericGaussianInterpolateImageFunction<TImage>
::GetWithPadding4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int i1 = ifloor(x - _RadiusX);
  const int j1 = ifloor(y - _RadiusY);
  const int k1 = ifloor(z - _RadiusZ);
  const int l1 = ifloor(t - _RadiusT);

  const int i2 = ifloor(x + _RadiusX);
  const int j2 = ifloor(y + _RadiusY);
  const int k2 = ifloor(z + _RadiusZ);
  const int l2 = ifloor(t + _RadiusT);

  // Create Gaussian interpolation kernel
  ScalarGaussian kernel(_Sigma/_dx, _Sigma/_dy, _Sigma/_dz, _Sigma/_dt, x, y, z, t);

  // Perform Gaussian interpolation
  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  for (int l = l1; l <= l2; ++l) {
    for (int k = k1; k <= k2; ++k) {
      for (int j = j1; j <= j2; ++j) {
        for (int i = i1; i <= i2; ++i) {
          w = static_cast<Real>(kernel.Evaluate(i, j, k, l));
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
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
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
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
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
GenericGaussianInterpolateImageFunction<TImage>
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
GenericGaussianInterpolateImageFunction<TImage>
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
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
::GetInside(double x, double y, double z, double t) const
{
  return Get(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
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
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return GetWithPadding(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericGaussianInterpolateImageFunction<TImage>::VoxelType
GenericGaussianInterpolateImageFunction<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return GetWithPadding(this->Extrapolator(), x, y, z, t);
  } else {
    return GetWithPadding(x, y, z, t);
  }
}


} // namespace mirtk

#endif // MIRTK_GaussianInterpolateImageFunction_HXX
