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

  // Calculate offsets for fast pixel access
  const int xoffset[2] = {0, 1};
  const int yoffset[2] = {0, this->Input()->X()};
  const int zoffset[2] = {0, this->Input()->X() * this->Input()->Y()};
  const int toffset[2] = {0, this->Input()->X() * this->Input()->Y() * this->Input()->Z()};

  int idx = 0;
  for (int d = 0; d <= 1; ++d)
  for (int c = 0; c <= 1; ++c)
  for (int b = 0; b <= 1; ++b)
  for (int a = 0; a <= 1; ++a, ++idx) {
    _Offset[idx] = xoffset[a] + yoffset[b] + zoffset[c] + toffset[d];
  }
}

// =============================================================================
// Interpolation weights
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline int GenericLinearInterpolateImageFunction<TImage>
::ComputeWeights(double x, Real w[2])
{
  const int i = ifloor(x);
  w[1] = static_cast<Real>(x) - static_cast<Real>(i);
  w[0] = static_cast<Real>(1) - w[1];
  return i;
}

// =============================================================================
// Domain checks
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::BoundingInterval(double x, int &i, int &I) const
{
  i = ifloor(x), I = i + 1;
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
  const int k = iround(z);
  const int l = iround(t);
  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);

  RealType val = voxel_cast<RealType>(0);
  Real nrm(0), w;

  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    if (0 <= jb && jb < this->Input()->Y()) {
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        if (0 <= ia && ia < this->Input()->X()) {
          w = wx[a] * wy[b];
          val += w * voxel_cast<RealType>(this->Input()->Get(ia, jb, k, l));
          nrm += w;
        }
      }
    }
  }

  if (nrm > Real(1e-3)) val /= nrm;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding2D(double x, double y, double z, double t) const
{
  const int k = iround(z);
  const int l = iround(t);
  if (k < 0 || k >= this->Input()->Z() ||
      l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);

  RealType val = voxel_cast<RealType>(0);
  Real fgw(0), w;

  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    for (int a = 0; a <= 1; ++a) {
      ia = i + a;
      w  = wx[a] * wy[b];
      if (this->Input()->IsInsideForeground(ia, jb, k, l)) {
        val += w * voxel_cast<RealType>(this->Input()->Get(ia, jb, k, l));
        fgw += w;
      }
    }
  }

  if (ifloor(fgw * 1.e3) > 500) val /= fgw;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  Real wx[2], wy[2];
  const int i = ComputeWeights(x, wx), I = i + 1;
  const int j = ComputeWeights(y, wy), J = j + 1;
  const int k = iround(z);
  const int l = iround(t);

  const auto v00 = voxel_cast<RealType>(input->Get(i, j, k, l));
  const auto v01 = voxel_cast<RealType>(input->Get(I, j, k, l));
  const auto v10 = voxel_cast<RealType>(input->Get(i, J, k, l));
  const auto v11 = voxel_cast<RealType>(input->Get(I, J, k, l));

  return voxel_cast<VoxelType>(wy[0] * (wx[0] * v00 + wx[1] * v01)  +
                               wy[1] * (wx[0] * v10 + wx[1] * v11));
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding2D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int k = iround(z);
  const int l = iround(t);
  if (k < 0 || k >= input->Z() ||
      l < 0 || l >= input->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);

  RealType val = voxel_cast<RealType>(0);
  Real fgw(0), w;

  int ia, jb;
  for (int b = 0; b <= 1; ++b) {
    jb = j + b;
    for (int a = 0; a <= 1; ++a) {
      ia = i + a;
      w  = wx[a] * wy[b];
      if (input->IsForeground(ia, jb, k, l)) {
        val += w * voxel_cast<RealType>(input->Get(ia, jb, k, l));
        fgw += w;
      }
    }
  }

  if (ifloor(fgw * 1.e3) > 500) val /= fgw;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get3D(double x, double y, double z, double t) const
{
  const int l = iround(t);
  if (l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2], wz[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);

  RealType val = voxel_cast<RealType>(0);
  Real nrm(0), w;

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
              val += w * voxel_cast<RealType>(this->Input()->Get(ia, jb, kc, l));
              nrm += w;
            }
          }
        }
      }
    }
  }

  if (nrm > Real(1e-3)) val /= nrm;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding3D(double x, double y, double z, double t) const
{
  const int l = iround(t);
  if (l < 0 || l >= this->Input()->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2], wz[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);

  RealType val = voxel_cast<RealType>(0);
  Real fgw(0), w;

  int ia, jb, kc;
  for (int c = 0; c <= 1; ++c) {
    kc = k + c;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        w  = wx[a] * wy[b] * wz[c];
        if (this->Input()->IsInsideForeground(ia, jb, kc, l)) {
          val += w * voxel_cast<RealType>(this->Input()->Get(ia, jb, kc, l));
          fgw += w;
        }
      }
    }
  }

  if (ifloor(fgw * 1.e3) > 500) val /= fgw;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  Real wx[2], wy[2], wz[2];
  const int i = ComputeWeights(x, wx), I = i + 1;
  const int j = ComputeWeights(y, wy), J = j + 1;
  const int k = ComputeWeights(z, wz), K = k + 1;
  const int l = iround(t);

  const auto v000 = voxel_cast<RealType>(input->Get(i, j, k, l));
  const auto v001 = voxel_cast<RealType>(input->Get(I, j, k, l));
  const auto v010 = voxel_cast<RealType>(input->Get(i, J, k, l));
  const auto v011 = voxel_cast<RealType>(input->Get(I, J, k, l));
  const auto v100 = voxel_cast<RealType>(input->Get(i, j, K, l));
  const auto v101 = voxel_cast<RealType>(input->Get(I, j, K, l));
  const auto v110 = voxel_cast<RealType>(input->Get(i, J, K, l));
  const auto v111 = voxel_cast<RealType>(input->Get(I, J, K, l));

  return voxel_cast<VoxelType>(wz[0] * (wy[0] * (wx[0] * v000 + wx[1] * v001)  +
                                        wy[1] * (wx[0] * v010 + wx[1] * v011)) +
                               wz[1] * (wy[0] * (wx[0] * v100 + wx[1] * v101)  +
                                        wy[1] * (wx[0] * v110 + wx[1] * v111)));
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding3D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  const int l = iround(t);
  if (l < 0 || l >= input->T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  Real wx[2], wy[2], wz[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);

  RealType val = voxel_cast<RealType>(0);
  Real fgw(0), wyz, w;

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
          val += w * static_cast<RealType>(input->Get(ia, jb, kc, l));
          fgw += w;
        }
      }
    }
  }

  if (ifloor(fgw * 1.e3) > 500) val /= fgw;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get4D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wt[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(t, wt);

  RealType val = voxel_cast<RealType>(0);
  Real nrm(0), w;

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
                  w = wx[a] * wy[b] * wz[c] * wt[d];
                  val += w * voxel_cast<RealType>(this->Input()->Get(ia, jb, kc, ld));
                  nrm += w;
                }
              }
            }
          }
        }
      }
    }
  }

  if (nrm > Real(1e-3)) val /= nrm;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding4D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wt[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(t, wt);

  RealType val = voxel_cast<RealType>(0);
  Real fgw(0), w;

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
            val += w * voxel_cast<RealType>(this->Input()->Get(ia, jb, kc, ld));
            fgw += w;
          }
        }
      }
    }
  }

  if (ifloor(fgw * 1.e3) > 500) val /= fgw;
  else val = voxel_cast<RealType>(this->DefaultValue());

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::Get4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  Real wx[2], wy[2], wz[2], wt[2];
  const int i = ComputeWeights(x, wx), I = i + 1;
  const int j = ComputeWeights(y, wy), J = j + 1;
  const int k = ComputeWeights(z, wz), K = k + 1;
  const int l = ComputeWeights(t, wt), L = l + 1;

  const auto v0000 = voxel_cast<RealType>(input->Get(i, j, k, l));
  const auto v0001 = voxel_cast<RealType>(input->Get(I, j, k, l));
  const auto v0010 = voxel_cast<RealType>(input->Get(i, J, k, l));
  const auto v0011 = voxel_cast<RealType>(input->Get(I, J, k, l));
  const auto v0100 = voxel_cast<RealType>(input->Get(i, j, K, l));
  const auto v0101 = voxel_cast<RealType>(input->Get(I, j, K, l));
  const auto v0110 = voxel_cast<RealType>(input->Get(i, J, K, l));
  const auto v0111 = voxel_cast<RealType>(input->Get(I, J, K, l));
  const auto v1000 = voxel_cast<RealType>(input->Get(i, j, k, L));
  const auto v1001 = voxel_cast<RealType>(input->Get(I, j, k, L));
  const auto v1010 = voxel_cast<RealType>(input->Get(i, J, k, L));
  const auto v1011 = voxel_cast<RealType>(input->Get(I, J, k, L));
  const auto v1100 = voxel_cast<RealType>(input->Get(i, j, K, L));
  const auto v1101 = voxel_cast<RealType>(input->Get(I, j, K, L));
  const auto v1110 = voxel_cast<RealType>(input->Get(i, J, K, L));
  const auto v1111 = voxel_cast<RealType>(input->Get(I, J, K, L));

  return voxel_cast<VoxelType>(wt[0] * (wz[0] * (wy[0] * (wx[0] * v0000 + wx[1] * v0001)   +
                                                 wy[1] * (wx[0] * v0010 + wx[1] * v0011))  +
                                        wz[1] * (wy[0] * (wx[0] * v0100 + wx[1] * v0101)   +
                                                 wy[1] * (wx[0] * v0110 + wx[1] * v0111))) +
                               wt[1] * (wz[0] * (wy[0] * (wx[0] * v1000 + wx[1] * v1001)   +
                                                 wy[1] * (wx[0] * v1010 + wx[1] * v1011))  +
                                        wz[1] * (wy[0] * (wx[0] * v1100 + wx[1] * v1101)   +
                                                 wy[1] * (wx[0] * v1110 + wx[1] * v1111))));
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
inline typename TOtherImage::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetWithPadding4D(const TOtherImage *input, double x, double y, double z, double t) const
{
  typedef typename TOtherImage::VoxelType VoxelType;
  typedef typename TOtherImage::RealType  RealType;

  Real wx[2], wy[2], wz[2], wt[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(t, wt);

  RealType val = voxel_cast<RealType>(0);
  Real fgw(0), wzt, wyzt, w;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    for (int c = 0; c <= 1; ++c) {
      kc  = k + c;
      wzt = wz[c] * wt[d];
      for (int b = 0; b <= 1; ++b) {
        jb   = j + b;
        wyzt = wy[b] * wzt;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          w  = wx[a] * wyzt;
          if (input->IsForeground(ia, jb, kc, ld)) {
            val += w * voxel_cast<RealType>(input->Get(ia, jb, kc, ld));
            fgw += w;
          }
        }
      }
    }
  }

  if (ifloor(fgw * 1.e3) > 500) val /= fgw;
  else val = voxel_cast<RealType>(this->DefaultValue());

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
  Real wx[2], wy[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = iround(z);
  const int l = iround(t);

  auto img = reinterpret_cast<const VoxelType *>(this->Input()->GetDataPointer(i, j, k, l));

  const auto v00 = voxel_cast<RealType>(img[_Offset[0]]);
  const auto v01 = voxel_cast<RealType>(img[_Offset[1]]);
  const auto v10 = voxel_cast<RealType>(img[_Offset[2]]);
  const auto v11 = voxel_cast<RealType>(img[_Offset[3]]);

  return voxel_cast<VoxelType>(wy[0] * (wx[0] * v00 + wx[1] * v01)  +
                               wy[1] * (wx[0] * v10 + wx[1] * v11));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetInside3D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = iround(t);

  auto img = reinterpret_cast<const VoxelType *>(this->Input()->GetDataPointer(i, j, k, l));

  const auto v000 = voxel_cast<RealType>(img[_Offset[0]]);
  const auto v001 = voxel_cast<RealType>(img[_Offset[1]]);
  const auto v010 = voxel_cast<RealType>(img[_Offset[2]]);
  const auto v011 = voxel_cast<RealType>(img[_Offset[3]]);
  const auto v100 = voxel_cast<RealType>(img[_Offset[4]]);
  const auto v101 = voxel_cast<RealType>(img[_Offset[5]]);
  const auto v110 = voxel_cast<RealType>(img[_Offset[6]]);
  const auto v111 = voxel_cast<RealType>(img[_Offset[7]]);

  return voxel_cast<VoxelType>(wz[0] * (wy[0] * (wx[0] * v000 + wx[1] * v001)  +
                                        wy[1] * (wx[0] * v010 + wx[1] * v011)) +
                               wz[1] * (wy[0] * (wx[0] * v100 + wx[1] * v101)  +
                                        wy[1] * (wx[0] * v110 + wx[1] * v111)));
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericLinearInterpolateImageFunction<TImage>::VoxelType
GenericLinearInterpolateImageFunction<TImage>
::GetInside4D(double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wt[2];
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(t, wt);

  auto img = reinterpret_cast<const VoxelType *>(this->Input()->GetDataPointer(i, j, k, l));

  const auto v0000 = voxel_cast<RealType>(img[_Offset[ 0]]);
  const auto v0001 = voxel_cast<RealType>(img[_Offset[ 1]]);
  const auto v0010 = voxel_cast<RealType>(img[_Offset[ 2]]);
  const auto v0011 = voxel_cast<RealType>(img[_Offset[ 3]]);
  const auto v0100 = voxel_cast<RealType>(img[_Offset[ 4]]);
  const auto v0101 = voxel_cast<RealType>(img[_Offset[ 5]]);
  const auto v0110 = voxel_cast<RealType>(img[_Offset[ 6]]);
  const auto v0111 = voxel_cast<RealType>(img[_Offset[ 7]]);
  const auto v1000 = voxel_cast<RealType>(img[_Offset[ 8]]);
  const auto v1001 = voxel_cast<RealType>(img[_Offset[ 9]]);
  const auto v1010 = voxel_cast<RealType>(img[_Offset[10]]);
  const auto v1011 = voxel_cast<RealType>(img[_Offset[11]]);
  const auto v1100 = voxel_cast<RealType>(img[_Offset[12]]);
  const auto v1101 = voxel_cast<RealType>(img[_Offset[13]]);
  const auto v1110 = voxel_cast<RealType>(img[_Offset[14]]);
  const auto v1111 = voxel_cast<RealType>(img[_Offset[15]]);

  return voxel_cast<VoxelType>(wt[0] * (wz[0] * (wy[0] * (wx[0] * v0000 + wx[1] * v0001)   +
                                                 wy[1] * (wx[0] * v0010 + wx[1] * v0011))  +
                                        wz[1] * (wy[0] * (wx[0] * v0100 + wx[1] * v0101)   +
                                                 wy[1] * (wx[0] * v0110 + wx[1] * v0111))) +
                               wt[1] * (wz[0] * (wy[0] * (wx[0] * v1000 + wx[1] * v1001)   +
                                                 wy[1] * (wx[0] * v1010 + wx[1] * v1011))  +
                                        wz[1] * (wy[0] * (wx[0] * v1100 + wx[1] * v1101)   +
                                                 wy[1] * (wx[0] * v1110 + wx[1] * v1111))));
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

// =============================================================================
// Evaluation of derivatives
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::EvaluateJacobianInside(Matrix &jac, double x, double y, double z, double t) const
{
  Jacobian(jac, this->Input(), x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::EvaluateJacobianOutside(Matrix &jac, double x, double y, double z, double t) const
{
  if (this->Extrapolator()) {
    Jacobian(jac, this->Extrapolator(), x, y, z, t);
  } else {
    Jacobian(jac, x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::EvaluateJacobianWithPaddingInside(Matrix &jac, double x, double y, double z, double t) const
{
  JacobianWithPadding(jac, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::EvaluateJacobianWithPaddingOutside(Matrix &jac, double x, double y, double z, double t) const
{
  JacobianWithPadding(jac, x, y, z, t);
}

// -----------------------------------------------------------------------------
// without use of extrapolation
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::Jacobian2D(Matrix &jac, double x, double y, double z, double t) const
{
  const TImage * const input = this->Input();

  Real wx[2], wy[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = iround(z);

  if (IsNaN(t) && IsZero(input->TSize()) && input->N() == 1) {

    jac.Initialize(input->T(), 2);

    if (i < 0 || i >= input->X() - 1 ||
        j < 0 || j >= input->Y() - 1 ||
        k < 0 || k >= input->Z()) {
      return; // zero
    }

    Real dv, w[2];
    int ia, jb;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        w[0] = wd[a] * wy[b];
        w[1] = wx[a] * wd[b];
        for (int l = 0; l < input->T(); ++l) {
          dv = voxel_cast<Real>(input->Get(ia, jb, k, l));
          jac(l, 0) += w[0] * dv;
          jac(l, 1) += w[1] * dv;
        }
      }
    }

  } else {

    jac.Initialize(input->N(), 2);

    const int l = (IsNaN(t) ? 0 : iround(t));
    if (i < 0 || i >= input->X() - 1 ||
        j < 0 || j >= input->Y() - 1 ||
        k < 0 || k >= input->Z() ||
        l < 0 || l >= input->T()) {
      return; // zero
    }

    RealType dv; // input value at discrete point
    RealType dx = voxel_cast<RealType>(0); // derivative(s) w.r.t. x
    RealType dy = voxel_cast<RealType>(0); // derivative(s) w.r.t. y

    int ia, jb;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        dv = voxel_cast<RealType>(input->Get(ia, jb, k, l));
        dx += (wd[a] * wy[b]) * dv;
        dy += (wx[a] * wd[b]) * dv;
      }
    }

    for (int r = 0; r < input->N(); ++r) {
      jac(r, 0) = get(dx, r);
      jac(r, 1) = get(dy, r);
    }

  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::Jacobian3D(Matrix &jac, double x, double y, double z, double t) const
{
  const TImage * const input = this->Input();

  Real wx[2], wy[2], wz[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);

  if (IsNaN(t) && IsZero(input->TSize()) && input->N() == 1) {

    jac.Initialize(input->T(), 3);

    if (i < 0 || i >= input->X() - 1 ||
        j < 0 || j >= input->Y() - 1 ||
        k < 0 || k >= input->Z() - 1) {
      return; // zero
    }

    Real w[3], dv;
    int ia, jb, kc;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          w[0] = wd[a] * wy[b] * wz[c];
          w[1] = wx[a] * wd[b] * wz[c];
          w[2] = wx[a] * wy[b] * wd[c];
          for (int l = 0; l < input->T(); ++l) {
            dv = voxel_cast<Real>(input->Get(ia, jb, kc, l));
            jac(l, 0) += w[0] * dv;
            jac(l, 1) += w[1] * dv;
            jac(l, 2) += w[2] * dv;
          }
        }
      }
    }

  } else {

    jac.Initialize(input->N(), 3);

    const int l = iround(t);
    if (i < 0 || i >= input->X() - 1 ||
        j < 0 || j >= input->Y() - 1 ||
        k < 0 || k >= input->Z() - 1 ||
        l < 0 || l >= input->T()) {
      return; // zero
    }

    RealType dv; // input value at discrete point
    RealType dx = voxel_cast<RealType>(0); // derivative(s) w.r.t. x
    RealType dy = voxel_cast<RealType>(0); // derivative(s) w.r.t. y
    RealType dz = voxel_cast<RealType>(0); // derivative(s) w.r.t. z

    int ia, jb, kc;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          dv = voxel_cast<RealType>(input->Get(ia, jb, kc, l));
          dx += (wd[a] * wy[b] * wz[c]) * dv;
          dy += (wx[a] * wd[b] * wz[c]) * dv;
          dz += (wx[a] * wy[b] * wd[c]) * dv;
        }
      }
    }

    for (int r = 0; r < input->N(); ++r) {
      jac(r, 0) = get(dx, r);
      jac(r, 1) = get(dy, r);
      jac(r, 2) = get(dz, r);
    }

  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::Jacobian4D(Matrix &jac, double x, double y, double z, double t) const
{
  const TImage * const input = this->Input();

  jac.Initialize(input->N(), 4);

  Real wx[2], wy[2], wz[2], wt[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(IsNaN(t) ? 0. : t, wt);

  if (i < 0 || i >= input->X() - 1 ||
      j < 0 || j >= input->Y() - 1 ||
      k < 0 || k >= input->Z() - 1 ||
      l < 0 || l >= input->T() - 1) {
    return; // zero
  }

  RealType dv; // input value at discrete point
  RealType dx = voxel_cast<RealType>(0); // derivative(s) w.r.t. x
  RealType dy = voxel_cast<RealType>(0); // derivative(s) w.r.t. y
  RealType dz = voxel_cast<RealType>(0); // derivative(s) w.r.t. z
  RealType dt = voxel_cast<RealType>(0); // derivative(s) w.r.t. t

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          dv = voxel_cast<RealType>(input->Get(ia, jb, kc, ld));
          dx += (wd[a] * wy[b] * wz[c] * wt[d]) * dv;
          dy += (wx[a] * wd[b] * wz[c] * wt[d]) * dv;
          dz += (wx[a] * wy[b] * wd[c] * wt[d]) * dv;
          dt += (wx[a] * wy[b] * wz[c] * wd[d]) * dv;
        }
      }
    }
  }

  for (int r = 0; r < input->N(); ++r) {
    jac(r, 0) = get(dx, r);
    jac(r, 1) = get(dy, r);
    jac(r, 2) = get(dz, r);
    jac(r, 3) = get(dt, r);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::Jacobian(Matrix &jac, double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return Jacobian3D(jac, x, y, z, t);
    case 2:  return Jacobian2D(jac, x, y, z, t);
    default: return Jacobian4D(jac, x, y, z, t);
  }
}


// -----------------------------------------------------------------------------
// without use of extrapolation, inside image foreground
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::JacobianWithPadding2D(Matrix &jac, double x, double y, double z, double t) const
{
  const TImage * const input = this->Input();

  Real wx[2], wy[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = iround(z);

  if (IsNaN(t) && IsZero(input->TSize()) && input->N() == 1) {

    jac.Initialize(input->T(), 2);

    if (i < 0 || i >= input->X() - 1 ||
        j < 0 || j >= input->Y() - 1 ||
        k < 0 || k >= input->Z()) {
      return; // zero
    }

    Real w[2], dv;
    int ia, jb;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        if (!input->IsForeground(ia, jb, k)) {
          jac.Initialize(input->T(), 2);
          return; // zero
        }
        w[0] = wd[a] * wy[b];
        w[1] = wx[a] * wd[b];
        for (int l = 0; l < input->T(); ++l) {
          dv = voxel_cast<Real>(input->Get(ia, jb, k, l));
          jac(l, 0) += w[0] * dv;
          jac(l, 1) += w[1] * dv;
        }
      }
    }

  } else {

    jac.Initialize(input->N(), 2);

    const int l = (IsNaN(t) ? 0 : iround(t));
    if (i < 0 || i >= input->X() - 1 ||
        j < 0 || j >= input->Y() - 1 ||
        k < 0 || k >= input->Z() ||
        l < 0 || l >= input->T()) {
      return; // zero
    }

    RealType dv; // input value at discrete point
    RealType dx = voxel_cast<RealType>(0); // derivative(s) w.r.t. x
    RealType dy = voxel_cast<RealType>(0); // derivative(s) w.r.t. y

    int ia, jb;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        if (!input->IsForeground(ia, jb, k, l)) {
          return; // zero
        }
        dv = voxel_cast<RealType>(input->Get(ia, jb, k, l));
        dx += (wd[a] * wy[b]) * dv;
        dy += (wx[a] * wd[b]) * dv;
      }
    }

    for (int r = 0; r < input->N(); ++r) {
      jac(r, 0) = get(dx, r);
      jac(r, 1) = get(dy, r);
    }

  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::JacobianWithPadding3D(Matrix &jac, double x, double y, double z, double t) const
{
  const TImage * const input = this->Input();

  Real wx[2], wy[2], wz[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);

  if (IsNaN(t) && IsZero(input->TSize()) && input->N() == 1) {

    jac.Initialize(input->T(), 3);

    if (i < 0 || i >= input->X() - 1 ||
        j < 0 || j >= input->Y() - 1 ||
        k < 0 || k >= input->Z() - 1) {
      return; // zero
    }

    Real w[3], dv;
    int ia, jb, kc;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          if (!input->IsForeground(ia, jb, kc)) {
            jac.Initialize(input->T(), 3);
            return; // zero
          }
          w[0] = wd[a] * wy[b] * wz[c];
          w[1] = wx[a] * wd[b] * wz[c];
          w[2] = wx[a] * wy[b] * wd[c];
          for (int l = 0; l < input->T(); ++l) {
            dv = voxel_cast<Real>(input->Get(ia, jb, kc, l));
            jac(l, 0) += w[0] * dv;
            jac(l, 1) += w[1] * dv;
            jac(l, 2) += w[2] * dv;
          }
        }
      }
    }

  } else {

    jac.Initialize(input->N(), 3);

    const int l = (IsNaN(t) ? 0 : iround(t));
    if (i < 0 || i >= input->X() - 1 ||
        j < 0 || j >= input->Y() - 1 ||
        k < 0 || k >= input->Z() - 1 ||
        l < 0 || l >= input->T()) {
      return; // zero
    }

    RealType dv; // input value at discrete point
    RealType dx = voxel_cast<RealType>(0); // derivative(s) w.r.t. x
    RealType dy = voxel_cast<RealType>(0); // derivative(s) w.r.t. y
    RealType dz = voxel_cast<RealType>(0); // derivative(s) w.r.t. z

    int ia, jb, kc;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          if (!input->IsForeground(ia, jb, kc, l)) {
            return; // zero
          }
          dv = voxel_cast<RealType>(input->Get(ia, jb, kc, l));
          dx += (wd[a] * wy[b] * wz[c]) * dv;
          dy += (wx[a] * wd[b] * wz[c]) * dv;
          dz += (wx[a] * wy[b] * wd[c]) * dv;
        }
      }
    }

    for (int r = 0; r < input->N(); ++r) {
      jac(r, 0) = get(dx, r);
      jac(r, 1) = get(dy, r);
      jac(r, 2) = get(dz, r);
    }

  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::JacobianWithPadding4D(Matrix &jac, double x, double y, double z, double t) const
{
  const TImage * const input = this->Input();

  jac.Initialize(input->N(), 4);

  Real wx[2], wy[2], wz[2], wt[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(IsNaN(t) ? 0. : t, wt);

  if (i < 0 || i >= input->X() - 1 ||
      j < 0 || j >= input->Y() - 1 ||
      k < 0 || k >= input->Z() - 1 ||
      l < 0 || l >= input->T() - 1) {
    return; // zero
  }

  RealType dv; // input value at discrete point
  RealType dx = voxel_cast<RealType>(0); // derivative(s) w.r.t. x
  RealType dy = voxel_cast<RealType>(0); // derivative(s) w.r.t. y
  RealType dz = voxel_cast<RealType>(0); // derivative(s) w.r.t. z
  RealType dt = voxel_cast<RealType>(0); // derivative(s) w.r.t. t

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          if (!input->IsForeground(ia, jb, kc, ld)) {
            return; // zero
          }
          dv = voxel_cast<RealType>(input->Get(ia, jb, kc, ld));
          dx += (wd[a] * wy[b] * wz[c] * wt[d]) * dv;
          dy += (wx[a] * wd[b] * wz[c] * wt[d]) * dv;
          dz += (wx[a] * wy[b] * wd[c] * wt[d]) * dv;
          dt += (wx[a] * wy[b] * wz[c] * wd[d]) * dv;
        }
      }
    }
  }

  for (int r = 0; r < input->N(); ++r) {
    jac(r, 0) = get(dx, r);
    jac(r, 1) = get(dy, r);
    jac(r, 2) = get(dz, r);
    jac(r, 3) = get(dt, r);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericLinearInterpolateImageFunction<TImage>
::JacobianWithPadding(Matrix &jac, double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return JacobianWithPadding3D(jac, x, y, z, t);
    case 2:  return JacobianWithPadding2D(jac, x, y, z, t);
    default: return JacobianWithPadding4D(jac, x, y, z, t);
  }
}


// -----------------------------------------------------------------------------
// inside image domain or with extrapolation function
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
void GenericLinearInterpolateImageFunction<TImage>
::Jacobian2D(Matrix &jac, const TOtherImage *input, double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = iround(z);

  if (IsNaN(t) && IsZero(input->TSize()) && input->N() == 1) {

    jac.Initialize(input->T(), 2);

    if (k < 0 || k >= input->Z()) {
      return; // zero
    }

    Real w[2], dv;
    int ia, jb;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        w[0] = wd[a] * wy[b];
        w[1] = wx[a] * wd[b];
        for (int l = 0; l < input->T(); ++l) {
          dv = voxel_cast<Real>(input->Get(ia, jb, k, l));
          jac(l, 0) += w[0] * dv;
          jac(l, 1) += w[1] * dv;
        }
      }
    }

  } else {

    jac.Initialize(input->N(), 2);

    const int l = (IsNaN(t) ? 0 : iround(t));
    if (k < 0 || k >= input->Z() ||
        l < 0 || l >= input->T()) {
      return; // zero
    }

    RealType dv; // input value at discrete point
    RealType dx = voxel_cast<RealType>(0); // derivative(s) w.r.t. x
    RealType dy = voxel_cast<RealType>(0); // derivative(s) w.r.t. y

    int ia, jb;
    for (int b = 0; b <= 1; ++b) {
      jb = j + b;
      for (int a = 0; a <= 1; ++a) {
        ia = i + a;
        dv = voxel_cast<RealType>(input->Get(ia, jb, k, l));
        dx += (wd[a] * wy[b]) * dv;
        dy += (wx[a] * wd[b]) * dv;
      }
    }

    for (int r = 0; r < input->N(); ++r) {
      jac(r, 0) = get(dx, r);
      jac(r, 1) = get(dy, r);
    }

  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
void GenericLinearInterpolateImageFunction<TImage>
::Jacobian3D(Matrix &jac, const TOtherImage *input,
             double x, double y, double z, double t) const
{
  Real wx[2], wy[2], wz[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);

  if (IsNaN(t) && IsZero(input->TSize()) && input->N() == 1) {

    jac.Initialize(input->T(), 3);

    Real w[3], dv;
    int ia, jb, kc;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          w[0] = wd[a] * wy[b] * wz[c];
          w[1] = wx[a] * wd[b] * wz[c];
          w[2] = wx[a] * wy[b] * wd[c];
          for (int l = 0; l < input->T(); ++l) {
            dv = voxel_cast<Real>(input->Get(ia, jb, kc, l));
            jac(l, 0) += w[0] * dv;
            jac(l, 1) += w[1] * dv;
            jac(l, 2) += w[2] * dv;
          }
        }
      }
    }

  } else {

    jac.Initialize(input->N(), 3);

    const int l = (IsNaN(t) ? 0 : iround(t));
    if (l < 0 || l >= input->T()) {
      return; // zero
    }

    RealType dv; // input value at discrete point
    RealType dx = voxel_cast<RealType>(0); // derivative(s) w.r.t. x
    RealType dy = voxel_cast<RealType>(0); // derivative(s) w.r.t. y
    RealType dz = voxel_cast<RealType>(0); // derivative(s) w.r.t. z

    int ia, jb, kc;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          dv = voxel_cast<RealType>(input->Get(ia, jb, kc, l));
          dx += (wd[a] * wy[b] * wz[c]) * dv;
          dy += (wx[a] * wd[b] * wz[c]) * dv;
          dz += (wx[a] * wy[b] * wd[c]) * dv;
        }
      }
    }

    for (int r = 0; r < input->N(); ++r) {
      jac(r, 0) = get(dx, r);
      jac(r, 1) = get(dy, r);
      jac(r, 2) = get(dz, r);
    }

  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
void GenericLinearInterpolateImageFunction<TImage>
::Jacobian4D(Matrix &jac, const TOtherImage *input,
             double x, double y, double z, double t) const
{
  jac.Initialize(input->N(), 4);

  Real wx[2], wy[2], wz[2], wt[2], wd[2] = {Real(-1), Real(1)};
  const int i = ComputeWeights(x, wx);
  const int j = ComputeWeights(y, wy);
  const int k = ComputeWeights(z, wz);
  const int l = ComputeWeights(IsNaN(t) ? 0. : t, wt);

  RealType dv; // input value at discrete point
  RealType dx = voxel_cast<RealType>(0); // derivative(s) w.r.t. x
  RealType dy = voxel_cast<RealType>(0); // derivative(s) w.r.t. y
  RealType dz = voxel_cast<RealType>(0); // derivative(s) w.r.t. z
  RealType dt = voxel_cast<RealType>(0); // derivative(s) w.r.t. t

  int ia, jb, kc, ld;
  for (int d = 0; d <= 1; ++d) {
    ld = l + d;
    for (int c = 0; c <= 1; ++c) {
      kc = k + c;
      for (int b = 0; b <= 1; ++b) {
        jb = j + b;
        for (int a = 0; a <= 1; ++a) {
          ia = i + a;
          dv = voxel_cast<RealType>(input->Get(ia, jb, kc, ld));
          dx += (wd[a] * wy[b] * wz[c] * wt[d]) * dv;
          dy += (wx[a] * wd[b] * wz[c] * wt[d]) * dv;
          dz += (wx[a] * wy[b] * wd[c] * wt[d]) * dv;
          dt += (wx[a] * wy[b] * wz[c] * wd[d]) * dv;
        }
      }
    }
  }

  for (int r = 0; r < input->N(); ++r) {
    jac(r, 0) = get(dx, r);
    jac(r, 1) = get(dy, r);
    jac(r, 2) = get(dz, r);
    jac(r, 3) = get(dt, r);
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
void GenericLinearInterpolateImageFunction<TImage>
::Jacobian(Matrix &jac, const TOtherImage *input,
           double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return Jacobian3D(jac, input, x, y, z, t);
    case 2:  return Jacobian2D(jac, input, x, y, z, t);
    default: return Jacobian4D(jac, input, x, y, z, t);
  }
}


} // namespace mirtk

#endif // MIRTK_LinearInterpolateImageFunction_HXX
