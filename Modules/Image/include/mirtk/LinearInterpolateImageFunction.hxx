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

#ifndef MIRTK_LinearInterpolateImageFunction_HXX
#define MIRTK_LinearInterpolateImageFunction_HXX

#include "mirtk/LinearInterpolateImageFunction.h"
#include "mirtk/InterpolateImageFunction.hxx"

#include "mirtk/Math.h"
#include "mirtk/VoxelCast.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
GenericLinearInterpolateImageFunction<TImage>
::GenericLinearInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
GenericLinearInterpolateImageFunction<TImage>
::~GenericLinearInterpolateImageFunction()
{
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>::Initialize(bool coeff)
{
  /// Initialize base class
  Superclass::Initialize(coeff);

  // Domain on which linear interpolation is defined: [0, x)
  switch (this->NumberOfDimensions()) {
    case 4:
      this->_t1 = .0;
      this->_t2 = fdec(this->Input()->T() - 1);
    case 3:
      this->_z1 = .0;
      this->_z2 = fdec(this->Input()->Z() - 1);
    default:
      this->_y1 = .0;
      this->_y2 = fdec(this->Input()->Y() - 1);
      this->_x1 = .0;
      this->_x2 = fdec(this->Input()->X() - 1);
  }

  // Calculate offsets for fast pixel access
  const int xoffset[2] = {0, 1};
  const int yoffset[2] = {0, this->Input()->X()};
  const int zoffset[2] = {0, this->Input()->X() * this->Input()->Y()};
  const int toffset[2] = {0, this->Input()->X() * this->Input()->Y() * this->Input()->Z()};

  int idx = 0;
  for (int d = 0; d <= 1; ++d) {
    for (int c = 0; c <= 1; ++c) {
      for (int b = 0; b <= 1; ++b) {
        for (int a = 0; a <= 1; ++a) {
          _Offset[idx++] = xoffset[a] + yoffset[b] + zoffset[c] + toffset[d];
        }
      }
    }
  }
}

// =============================================================================
// Domain checks
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::BoundingInterval(double x, int &i, int &I) const
{
  i = static_cast<int>(floor(x)), I = i + 1;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get2D(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(round(z));
  const int l = static_cast<int>(round(t));

  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2];
  wx[1] = Real(x - i); wx[0] = Real(1) - wx[1];
  wy[1] = Real(y - j); wy[0] = Real(1) - wy[1];

  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), w;
  
  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    if (0 <= jb && jb < this->Input()->Y()) {
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        if (0 <= ia && ia < this->Input()->X()) {
          w    = wx[a] * wy[b];
          val += w * this->Input()->Get(ia, jb, k, l);
          nrm += w;
        }
      }
    }
  }

  if (nrm) val /= nrm;
  else     val  = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding2D(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(round(z));
  const int l = static_cast<int>(round(t));

  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2];
  wx[1] = Real(x - i); wx[0] = Real(1) - wx[1];
  wy[1] = Real(y - j); wy[0] = Real(1) - wy[1];

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    for (int a = 0; a <= 1; ++a) {
      ia = i + a;
      w  = wx[a] * wy[b];
      if (this->Input()->IsInsideForeground(ia, jb, k, l)) {
        val += w * this->Input()->Get(ia, jb, k, l);
        fgw += w;
      } else {
        bgw += w;
      }
    }
  }

  if (fgw > bgw) val /= fgw;
  else           val  = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(round(z));
  const int l = static_cast<int>(round(t));
  const int I = i + 1;
  const int J = j + 1;

  const Real A = Real(x - i);
  const Real B = Real(y - j);
  const Real a = Real(1) - A;
  const Real b = Real(1) - B;

  typename TOtherImage::RealType val;
  val = (b * (a * input->Get(i, j, k, l) + A * input->Get(I, j, k, l))  +
         B * (a * input->Get(i, J, k, l) + A * input->Get(I, J, k, l)));
  return voxel_cast<typename TOtherImage::VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(round(z));
  const int l = static_cast<int>(round(t));

  if (k < 0 || k >= input->Z() ||
      l < 0 || l >= input->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2];
  wx[1] = Real(x - i); wx[0] = Real(1) - wx[1];
  wy[1] = Real(y - j); wy[0] = Real(1) - wy[1];

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    for (int a = 0; a <= 1; ++a) {
      ia = i + a;
      w  = wx[a] * wy[b];
      if (input->IsForeground(ia, jb, k, l)) {
        val += w * input->Get(ia, jb, k, l);
        fgw += w;
      } else {
        bgw += w;
      }
    }
  }

  if (fgw > bgw) val /= fgw;
  else           val  = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get3D(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(floor(z));
  const int l = static_cast<int>(round(t));

  if (l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2], wz[2];
  wx[1] = Real(x - i); wx[0] = Real(1) - wx[1];
  wy[1] = Real(y - j); wy[0] = Real(1) - wy[1];
  wz[1] = Real(z - k); wz[0] = Real(1) - wz[1];

  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), w;

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
              w    = wx[a] * wy[b] * wz[c];
              val += w * this->Input()->Get(ia, jb, kc, l);
              nrm += w;
            }
          }
        }
      }
    }
  }

  if (nrm) val /= nrm;
  else     val  = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding3D(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(floor(z));
  const int l = static_cast<int>(round(t));

  if (l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2], wz[2];
  wx[1] = Real(x - i); wx[0] = Real(1) - wx[1];
  wy[1] = Real(y - j); wy[0] = Real(1) - wy[1];
  wz[1] = Real(z - k); wz[0] = Real(1) - wz[1];

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  int ia, jb, kc;
  for (int c = 0; c <= 1; ++c) {
    kc = k + c;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        w  = wx[a] * wy[b] * wz[c];
        if (this->Input()->IsInsideForeground(ia, jb, kc, l)) {
          val += w * this->Input()->Get(ia, jb, kc, l);
          fgw += w;
        } else {
          bgw += w;
        }
      }
    }
  }

  if (fgw > bgw) val /= fgw;
  else           val  = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(floor(z));
  const int l = static_cast<int>(round(t));
  const int I = i + 1;
  const int J = j + 1;
  const int K = k + 1;

  const Real A = Real(x - i);
  const Real B = Real(y - j);
  const Real C = Real(z - k);
  const Real a = Real(1) - A;
  const Real b = Real(1) - B;
  const Real c = Real(1) - C;

  typename TOtherImage::RealType val;
  val = (c * (b * (a * input->Get(i, j, k, l) + A * input->Get(I, j, k, l))  +
              B * (a * input->Get(i, J, k, l) + A * input->Get(I, J, k, l))) +
         C * (b * (a * input->Get(i, j, K, l) + A * input->Get(I, j, K, l))  +
              B * (a * input->Get(i, J, K, l) + A * input->Get(I, J, K, l))));
  return voxel_cast<typename TOtherImage::VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(floor(z));
  const int l = static_cast<int>(round(t));

  if (l < 0 || l >= input->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2], wz[2];
  wx[1] = Real(x - i); wx[0] = Real(1) - wx[1];
  wy[1] = Real(y - j); wy[0] = Real(1) - wy[1];
  wz[1] = Real(z - k); wz[0] = Real(1) - wz[1];

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wyz, w;

  int ia, jb, kc;
  for (int c = 0; c <= 1; ++c) {
    kc = k + c;
    for (int b = 0; b <= 1; ++b) {
      jb  = j + b;
      wyz = wy[b] * wz[c];
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        w  = wx[a] * wyz;
        if (input->IsForeground(ia, jb, kc, l)) {
          val += w * input->Get(ia, jb, kc, l);
          fgw += w;
        } else {
          bgw += w;
        }
      }
    }
  }

  if (fgw > bgw) val /= fgw;
  else           val  = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get4D(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(floor(z));
  const int l = static_cast<int>(floor(t));

  Real wx[2], wy[2], wz[2], wt[2];
  wx[1] = Real(x - i); wx[0] = Real(1) - wx[1];
  wy[1] = Real(y - j); wy[0] = Real(1) - wy[1];
  wz[1] = Real(z - k); wz[0] = Real(1) - wz[1];
  wt[1] = Real(t - l); wt[0] = Real(1) - wt[1];

  RealType val = voxel_cast<RealType>(0);
  Real     nrm(0), w;

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
                  w    = wx[a] * wy[b] * wz[c] * wt[d];
                  val += w * this->Input()->Get(ia, jb, kc, ld);
                  nrm += w;
                }
              }
            }
          }
        }
      }
    }
  }

  if (nrm) val /= nrm;
  else     val  = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding4D(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(floor(z));
  const int l = static_cast<int>(floor(t));

  Real wx[2], wy[2], wz[2], wt[2];
  wx[1] = Real(x - i); wx[0] = Real(1) - wx[1];
  wy[1] = Real(y - j); wy[0] = Real(1) - wy[1];
  wz[1] = Real(z - k); wz[0] = Real(1) - wz[1];
  wt[1] = Real(t - l); wt[0] = Real(1) - wt[1];

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          w  = wx[a] * wy[b] * wz[c] * wt[d];
          if (this->Input()->IsInsideForeground(ia, jb, kc, ld)) {
            val += w * this->Input()->Get(ia, jb, kc, ld);
            fgw += w;
          } else {
            bgw += w;
          }
        }
      }
    }
  }

  if (fgw > bgw) val /= fgw;
  else           val  = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(floor(z));
  const int l = static_cast<int>(floor(t));
  const int I = i + 1;
  const int J = j + 1;
  const int K = k + 1;
  const int L = l + 1;

  const Real A = Real(x - i);
  const Real B = Real(y - j);
  const Real C = Real(z - k);
  const Real D = Real(t - l);
  const Real a = Real(1) - A;
  const Real b = Real(1) - B;
  const Real c = Real(1) - C;
  const Real d = Real(1) - D;

  typename TOtherImage::RealType val;
  val = (d * (c * (b * (a * input->Get(i, j, k, l) + A * input->Get(I, j, k, l))   +
                   B * (a * input->Get(i, J, k, l) + A * input->Get(I, J, k, l)))  +
              C * (b * (a * input->Get(i, j, K, l) + A * input->Get(I, j, K, l))   +
                   B * (a * input->Get(i, J, K, l) + A * input->Get(I, J, K, l)))) +
         D * (c * (b * (a * input->Get(i, j, k, L) + A * input->Get(I, j, k, L))   +
                   B * (a * input->Get(i, J, k, L) + A * input->Get(I, J, k, L)))  +
              C * (b * (a * input->Get(i, j, K, L) + A * input->Get(I, j, K, L))   +
                   B * (a * input->Get(i, J, K, L) + A * input->Get(I, J, K, L)))));
  return voxel_cast<typename TOtherImage::VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int i = static_cast<int>(floor(x));
  const int j = static_cast<int>(floor(y));
  const int k = static_cast<int>(floor(z));
  const int l = static_cast<int>(floor(t));

  Real wx[2], wy[2], wz[2], wt[2];
  wx[1] = Real(x - i); wx[0] = Real(1) - wx[1];
  wy[1] = Real(y - j); wy[0] = Real(1) - wy[1];
  wz[1] = Real(z - k); wz[0] = Real(1) - wz[1];
  wt[1] = Real(t - l); wt[0] = Real(1) - wt[1];

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wzt, wyzt, w;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    for (int c = 0; c <= 1; ++c) {
      kc  = k + c;
      wzt = wz[c] * wt[d];
      for (int b = 0; b <= 1; ++b) {
        jb  = j + b;
        wyzt = wy[b] * wzt;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          w  = wx[a] * wyzt;
          if (input->IsForeground(ia, jb, kc, ld)) {
            val += w * input->Get(ia, jb, kc, ld);
            fgw += w;
          } else {
            bgw += w;
          }
        }
      }
    }
  }

  if (fgw > bgw) val /= fgw;
  else           val  = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
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
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
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
GenericLinearInterpolateImageFunction<TImage>
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
GenericLinearInterpolateImageFunction<TImage>
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
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetInside2D(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(x);
  const int j = static_cast<int>(y);
  const int k = static_cast<int>(round(z));
  const int l = static_cast<int>(round(t));

  const Real A = Real(x - i);
  const Real B = Real(y - j);
  const Real a = Real(1) - A;
  const Real b = Real(1) - B;

  const VoxelType *img;
  img = reinterpret_cast<const VoxelType *>(this->Input()->GetDataPointer(i, j, k, l));
  return voxel_cast<VoxelType>(b * (a * img[_Offset[0]] + A * img[_Offset[1]]) +
                               B * (a * img[_Offset[2]] + A * img[_Offset[3]]));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetInside3D(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(x);
  const int j = static_cast<int>(y);
  const int k = static_cast<int>(z);
  const int l = static_cast<int>(round(t));

  const Real A = Real(x - i);
  const Real B = Real(y - j);
  const Real C = Real(z - k);
  const Real a = Real(1) - A;
  const Real b = Real(1) - B;
  const Real c = Real(1) - C;

  const VoxelType *img;
  img = reinterpret_cast<const VoxelType *>(this->Input()->GetDataPointer(i, j, k, l));
  return voxel_cast<VoxelType>(c * (b * (a * img[_Offset[0]] + A * img[_Offset[1]])  +
                                    B * (a * img[_Offset[2]] + A * img[_Offset[3]])) +
                               C * (b * (a * img[_Offset[4]] + A * img[_Offset[5]])  +
                                    B * (a * img[_Offset[6]] + A * img[_Offset[7]])));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetInside4D(double x, double y, double z, double t) const
{
  const int i = static_cast<int>(x);
  const int j = static_cast<int>(y);
  const int k = static_cast<int>(z);
  const int l = static_cast<int>(t);

  const Real A = Real(x - i);
  const Real B = Real(y - j);
  const Real C = Real(z - k);
  const Real D = Real(t - l);
  const Real a = Real(1) - A;
  const Real b = Real(1) - B;
  const Real c = Real(1) - C;
  const Real d = Real(1) - D;

  const VoxelType *img;
  img = reinterpret_cast<const VoxelType *>(this->Input()->GetDataPointer(i, j, k, l));
  return voxel_cast<VoxelType>(d * (c * (b * (a * img[_Offset[ 0]] + A * img[_Offset[ 1]])   +
                                         B * (a * img[_Offset[ 2]] + A * img[_Offset[ 3]]))  +
                                    C * (b * (a * img[_Offset[ 4]] + A * img[_Offset[ 5]])   +
                                         B * (a * img[_Offset[ 6]] + A * img[_Offset[ 7]]))) +
                               D * (c * (b * (a * img[_Offset[ 8]] + A * img[_Offset[ 9]])   +
                                         B * (a * img[_Offset[10]] + A * img[_Offset[11]]))  +
                                    C * (b * (a * img[_Offset[12]] + A * img[_Offset[13]])   +
                                         B * (a * img[_Offset[14]] + A * img[_Offset[15]]))));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetInside(double x, double y, double z, double t) const
{
  // Use faster pixel access than Get(Input(), x, y, z, t), which requires
  // a 4D array lookup with three indirections instead of a direct 1D lookup
  switch (this->NumberOfDimensions()) {
    case 3:  return GetInside3D(x, y, z, t);
    case 2:  return GetInside2D(x, y, z, t);
    default: return GetInside4D(x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <>
inline typename GenericLinearInterpolateImageFunction<BaseImage>::VoxelType
GenericLinearInterpolateImageFunction<BaseImage>
::GetInside(double x, double y, double z, double t) const
{
  return Get(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
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
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  return GetWithPadding(this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    return GetWithPadding(this->Extrapolator(), x, y, z, t);
  } else {
    return GetWithPadding(x, y, z, t);
  }
}


} // namespace mirtk

#endif // MIRTK_LinearInterpolateImageFunction_HXX
