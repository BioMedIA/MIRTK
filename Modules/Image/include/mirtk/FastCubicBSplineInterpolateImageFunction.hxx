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

#ifndef MIRTK_FastCubicBSplineInterpolateImageFunction_HXX
#define MIRTK_FastCubicBSplineInterpolateImageFunction_HXX

#include "mirtk/FastCubicBSplineInterpolateImageFunction.h"

#include "mirtk/Math.h"
#include "mirtk/BSpline.h"
#include "mirtk/Matrix.h"
#include "mirtk/Vector3D.h" // get
#include "mirtk/Vector4D.h" // get

#include "mirtk/InterpolateImageFunction.hxx"
#include "mirtk/ImageToInterpolationCoefficients.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GenericFastCubicBSplineInterpolateImageFunction()
:
  _InfiniteCoefficient(NULL)
{
  // Default extrapolation mode is to apply the mirror boundary condition
  // which is also assumed when converting an input image to spline coefficients
  this->Extrapolator(ExtrapolatorType::New(Extrapolation_Mirror), true);
}

// -----------------------------------------------------------------------------
template <class TImage>
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::~GenericFastCubicBSplineInterpolateImageFunction()
{
  delete _InfiniteCoefficient;
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Initialize(bool coeff)
{
  // Initialize base class
  Superclass::Initialize(coeff);

  // Initialize kernel lookup table
  Kernel::Initialize();

  // Domain on which cubic B-spline interpolation is defined
  const double margin = 2.0;
  switch (this->NumberOfDimensions()) {
    case 4:
      this->_t1 = fdec(margin);
      this->_t2 = this->Input()->T() - margin - 1;
    case 3:
      this->_z1 = fdec(margin);
      this->_z2 = this->Input()->Z() - margin - 1;
    default:
      this->_y1 = fdec(margin);
      this->_y2 = this->Input()->Y() - margin - 1;
      this->_x1 = fdec(margin);
      this->_x2 = this->Input()->X() - margin - 1;
  }

  // Initialize coefficient image
  RealType *data = nullptr;
  _UseInputCoefficients = coeff;
  if (_UseInputCoefficients && this->Input()->GetDataType() == voxel_info<RealType>::type()) {
    data = reinterpret_cast<RealType *>(const_cast<void *>(this->Input()->GetDataPointer()));
  }
  _Coefficient.Initialize(this->Input()->Attributes(), data);
  this->Update();

  // Initialize infinite coefficient image (i.e., extrapolator)
  if (!_InfiniteCoefficient || _InfiniteCoefficient->ExtrapolationMode() != this->ExtrapolationMode()) {
    delete _InfiniteCoefficient;
    _InfiniteCoefficient = CoefficientExtrapolator::New(this->ExtrapolationMode(), &_Coefficient);
  }
  if (_InfiniteCoefficient) {
    _InfiniteCoefficient->Input(&_Coefficient);
    _InfiniteCoefficient->Initialize();
  }

  // Compute strides for fast iteration over coefficient image
  //_s1 = 1;
  _s2 =  this->Input()->X() - 4;
  _s3 = (this->Input()->Y() - 4) * this->Input()->X();
  _s4 = (this->Input()->Z() - 4) * this->Input()->Y() * this->Input()->X();
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>::Update()
{
  if (_Coefficient.GetDataPointer() != this->Input()->GetDataPointer()) {
    _Coefficient = *(this->Input());
    if (!_UseInputCoefficients) {
      Real poles[2];
      int  npoles;
      SplinePoles(3, poles, npoles);
      FillBackgroundBeforeConversionToSplineCoefficients(_Coefficient);
      switch (this->NumberOfDimensions()) {
        case 4:  ConvertToInterpolationCoefficientsT(_Coefficient, poles, npoles);
        case 3:  ConvertToInterpolationCoefficientsZ(_Coefficient, poles, npoles);
        default: ConvertToInterpolationCoefficientsY(_Coefficient, poles, npoles);
                 ConvertToInterpolationCoefficientsX(_Coefficient, poles, npoles);
      }
    }
  }
}

// =============================================================================
// Domain checks
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::BoundingInterval(double x, int &i, int &I) const
{
  i = ifloor(x), I = i + 3;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Get2D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  if (k < 0 || k >= _Coefficient.Z() ||
      l < 0 || l >= _Coefficient.T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));

  --i, --j;

  RealType val = voxel_cast<RealType>(0);
  Real     w;

  int ia, jb;
  for (int b = 0; b <= 3; ++b) {
    jb = j + b;
    DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
    for (int a = 0; a <= 3; ++a) {
      ia = i + a;
      DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
      w = Kernel::LookupTable[A][a] * Kernel::LookupTable[B][b];
      val += w * _Coefficient(ia, jb, k, l);
    }
  }

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetWithPadding2D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  if (k < 0 || k >= _Coefficient.Z() ||
      l < 0 || l >= _Coefficient.T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));

  --i, --j;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  int ia, jb;
  for (int b = 0; b <= 3; ++b) {
    jb = j + b;
    DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
    for (int a = 0; a <= 3; ++a) {
      ia = i + a;
      DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
      w = Kernel::LookupTable[A][a] * Kernel::LookupTable[B][b];
      val += w * _Coefficient(ia, jb, k, l);
      if (this->Input()->IsInsideForeground(i + a, j + b, k, l)) {
        fgw += w;
      } else {
        bgw += w;
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    val = voxel_cast<RealType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TCoefficient>
inline typename TCoefficient::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Get2D(const TCoefficient *coeff, double x, double y, double z, double t) const
{
  typedef typename TCoefficient::VoxelType VoxelType;
  typedef typename TCoefficient::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));

  --i, --j;

  RealType val = voxel_cast<RealType>(0);

  int jb;
  for (int b = 0; b <= 3; ++b) {
    jb   = j + b;
    val += Kernel::LookupTable[A][0] * Kernel::LookupTable[B][b] * coeff->Get(i,   jb, k, l);
    val += Kernel::LookupTable[A][1] * Kernel::LookupTable[B][b] * coeff->Get(i+1, jb, k, l);
    val += Kernel::LookupTable[A][2] * Kernel::LookupTable[B][b] * coeff->Get(i+2, jb, k, l);
    val += Kernel::LookupTable[A][3] * Kernel::LookupTable[B][b] * coeff->Get(i+3, jb, k, l);
  }

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage, class TCoefficient>
inline typename TCoefficient::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetWithPadding2D(const TOtherImage *input, const TCoefficient *coeff,
                   double x, double y, double z, double t) const
{
  typedef typename TCoefficient::VoxelType VoxelType;
  typedef typename TCoefficient::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));

  --i, --j;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), w;

  int ia, jb;
  for (int b = 0; b <= 3; ++b) {
    jb = j + b;
    for (int a = 0; a <= 3; ++a) {
      ia = i + a;
      w = Kernel::LookupTable[A][a] * Kernel::LookupTable[B][b];
      val += w * voxel_cast<RealType>(coeff->Get(ia, jb, k, l));
      if (input->IsForeground(ia, jb, k, l)) {
        fgw += w;
      } else {
        bgw += w;
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    val = voxel_cast<RealType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Get3D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  if (l < 0 || l >= _Coefficient.T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));

  --i, --j, --k;

  RealType val = voxel_cast<RealType>(0);
  Real     wyz, w;

  int ia, jb, kc;
  for (int c = 0; c <= 3; ++c) {
    kc = k + c;
    DefaultExtrapolator::Apply(kc, _Coefficient.Z() - 1);
    for (int b = 0; b <= 3; ++b) {
      jb = j + b;
      DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
      wyz = Kernel::LookupTable[B][b] * Kernel::LookupTable[C][c];
      for (int a = 0; a <= 3; ++a) {
        ia = i + a;
        DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
        w = Kernel::LookupTable[A][a] * wyz;
        val += w * _Coefficient(ia, jb, kc, l);
      }
    }
  }

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetWithPadding3D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  if (l < 0 || l >= _Coefficient.T()) {
    return voxel_cast<VoxelType>(this->DefaultValue());
  }

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));

  --i, --j, --k;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wyz, w;

  int ia, jb, kc;
  for (int c = 0; c <= 3; ++c) {
    kc = k + c;
    DefaultExtrapolator::Apply(kc, _Coefficient.Z() - 1);
    for (int b = 0; b <= 3; ++b) {
      jb = j + b;
      DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
      wyz = Kernel::LookupTable[B][b] * Kernel::LookupTable[C][c];
      for (int a = 0; a <= 3; ++a) {
        ia = i + a;
        DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
        w = Kernel::LookupTable[A][a] * wyz;
        val += w * _Coefficient(ia, jb, kc, l);
        if (this->Input()->IsInsideForeground(i + a, j + b, k + c, l)) {
          fgw += w;
        } else {
          bgw += w;
        }
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    val = voxel_cast<RealType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TCoefficient>
inline typename TCoefficient::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Get3D(const TCoefficient *coeff, double x, double y, double z, double t) const
{
  typedef typename TCoefficient::VoxelType VoxelType;
  typedef typename TCoefficient::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));

  --i, --j, --k;

  RealType val = voxel_cast<RealType>(0);

  int  jb, kc;
  Real wyz;
  for (int c = 0; c <= 3; ++c) {
    kc = k + c;
    for (int b = 0; b <= 3; ++b) {
      jb   = j + b;
      wyz  = Kernel::LookupTable[B][b] * Kernel::LookupTable[C][c];
      val += Kernel::LookupTable[A][0] * wyz * voxel_cast<RealType>(coeff->Get(i,   jb, kc, l));
      val += Kernel::LookupTable[A][1] * wyz * voxel_cast<RealType>(coeff->Get(i+1, jb, kc, l));
      val += Kernel::LookupTable[A][2] * wyz * voxel_cast<RealType>(coeff->Get(i+2, jb, kc, l));
      val += Kernel::LookupTable[A][3] * wyz * voxel_cast<RealType>(coeff->Get(i+3, jb, kc, l));
    }
  }

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage, class TCoefficient>
inline typename TCoefficient::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetWithPadding3D(const TOtherImage *input, const TCoefficient *coeff,
                   double x, double y, double z, double t) const
{
  typedef typename TCoefficient::VoxelType VoxelType;
  typedef typename TCoefficient::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));

  --i, --j, --k;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wyz, w;

  int ia, jb, kc;
  for (int c = 0; c <= 3; ++c) {
    kc = k + c;
    DefaultExtrapolator::Apply(kc, _Coefficient.Z() - 1);
    for (int b = 0; b <= 3; ++b) {
      jb = j + b;
      DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
      wyz = Kernel::LookupTable[B][b] * Kernel::LookupTable[C][c];
      for (int a = 0; a <= 3; ++a) {
        ia = i + a;
        DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
        w = Kernel::LookupTable[A][a] * wyz;
        val += w * voxel_cast<RealType>(coeff->Get(ia, jb, kc, l));
        if (input->IsForeground(i + a, j + b, k + c, l)) {
          fgw += w;
        } else {
          bgw += w;
        }
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    val = voxel_cast<RealType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Get4D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));
  const int D = Kernel::VariableToIndex(Real(t - l));

  --i, --j, --k, --l;

  RealType val = voxel_cast<RealType>(0);
  Real     wzt, wyzt, w;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 3; ++d) {
    ld = l + d;
    DefaultExtrapolator::Apply(ld, _Coefficient.T() - 1);
    for (int c = 0; c <= 3; ++c) {
      kc = k + c;
      DefaultExtrapolator::Apply(kc, _Coefficient.Z() - 1);
      wzt = Kernel::LookupTable[C][c] * Kernel::LookupTable[D][d];
      for (int b = 0; b <= 3; ++b) {
        jb = j + b;
        DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
        wyzt = Kernel::LookupTable[B][b] * wzt;
        for (int a = 0; a <= 3; ++a) {
          ia = i + a;
          DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
          w = Kernel::LookupTable[A][a] * wyzt;
          val += w * _Coefficient(ia, jb, kc, ld);
        }
      }
    }
  }

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetWithPadding4D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));
  const int D = Kernel::VariableToIndex(Real(t - l));

  --i, --j, --k, --l;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wzt, wyzt, w;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 3; ++d) {
    ld = l + d;
    DefaultExtrapolator::Apply(ld, _Coefficient.T() - 1);
    for (int c = 0; c <= 3; ++c) {
      kc = k + c;
      DefaultExtrapolator::Apply(kc, _Coefficient.Z() - 1);
      wzt = Kernel::LookupTable[C][c] * Kernel::LookupTable[D][d];
      for (int b = 0; b <= 3; ++b) {
        jb = j + b;
        DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
        wyzt = Kernel::LookupTable[B][b] * wzt;
        for (int a = 0; a <= 3; ++a) {
          ia = i + a;
          DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
          w = Kernel::LookupTable[A][a] * wyzt;
          val += w * _Coefficient(ia, jb, kc, ld);
          if (this->Input()->IsInsideForeground(i + a, j + b, k + c, l + d)) {
            fgw += w;
          } else {
            bgw += w;
          }
        }
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    val = voxel_cast<RealType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TCoefficient>
inline typename TCoefficient::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Get4D(const TCoefficient *coeff, double x, double y, double z, double t) const
{
  typedef typename TCoefficient::VoxelType VoxelType;
  typedef typename TCoefficient::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));
  const int D = Kernel::VariableToIndex(Real(t - l));

  --i, --j, --k, --l;

  RealType val = voxel_cast<RealType>(0);

  int  jb, kc, ld;
  Real wzt, wyzt;
  for (int d = 0; d <= 3; ++d) {
    ld = l + d;
    for (int c = 0; c <= 3; ++c) {
      kc  = k + c;
      wzt = Kernel::LookupTable[C][c] * Kernel::LookupTable[D][d];
      for (int b = 0; b <= 3; ++b) {
        jb   = j + b;
        wyzt = Kernel::LookupTable[B][b] * wzt;
        val += Kernel::LookupTable[A][0] * wyzt * voxel_cast<RealType>(coeff->Get(i,   jb, kc, ld));
        val += Kernel::LookupTable[A][1] * wyzt * voxel_cast<RealType>(coeff->Get(i+1, jb, kc, ld));
        val += Kernel::LookupTable[A][2] * wyzt * voxel_cast<RealType>(coeff->Get(i+2, jb, kc, ld));
        val += Kernel::LookupTable[A][3] * wyzt * voxel_cast<RealType>(coeff->Get(i+3, jb, kc, ld));
      }
    }
  }

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage, class TCoefficient>
inline typename TCoefficient::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetWithPadding4D(const TOtherImage *input, const TCoefficient *coeff,
                   double x, double y, double z, double t) const
{
  typedef typename TCoefficient::VoxelType VoxelType;
  typedef typename TCoefficient::RealType  RealType;

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));
  const int D = Kernel::VariableToIndex(Real(t - l));

  --i, --j, --k, --l;

  RealType val = voxel_cast<RealType>(0);
  Real     fgw(0), bgw(0), wzt, wyzt, w;

  int ia, jb, kc, ld;
  for (int d = 0; d <= 3; ++d) {
    ld = l + d;
    for (int c = 0; c <= 3; ++c) {
      kc = k + c;
      wzt = Kernel::LookupTable[C][c] * Kernel::LookupTable[D][d];
      for (int b = 0; b <= 3; ++b) {
        jb   = j + b;
        wyzt = Kernel::LookupTable[B][b] * wzt;
        for (int a = 0; a <= 3; ++a) {
          ia = i + a;
          w = Kernel::LookupTable[A][a] * wyzt;
          val += w * voxel_cast<RealType>(coeff->Get(ia, jb, kc, ld));
          if (input->IsForeground(ia, jb, kc, ld)) {
            fgw += w;
          } else {
            bgw += w;
          }
        }
      }
    }
  }

  if (bgw > fgw || AreEqual(bgw, fgw, 1e-3)) {
    val = voxel_cast<RealType>(this->DefaultValue());
  }
  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
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
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
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
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Get(const TOtherImage *coeff, double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return Get3D(coeff, x, y, z, t);
    case 2:  return Get2D(coeff, x, y, z, t);
    default: return Get4D(coeff, x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage, class TCoefficient>
inline typename TCoefficient::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetWithPadding(const TOtherImage *input, const TCoefficient *coeff,
                 double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return GetWithPadding3D(input, coeff, x, y, z, t);
    case 2:  return GetWithPadding2D(input, coeff, x, y, z, t);
    default: return GetWithPadding4D(input, coeff, x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetInside2D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));

  --i, --j;

  RealType        val   = voxel_cast<RealType>(0);
  const RealType *coeff = _Coefficient.Data(i, j, k, l);

  for (int b = 0; b <= 3; ++b, coeff += _s2) {
    val += Kernel::LookupTable[A][0] * Kernel::LookupTable[B][b] * (*coeff), ++coeff;
    val += Kernel::LookupTable[A][1] * Kernel::LookupTable[B][b] * (*coeff), ++coeff;
    val += Kernel::LookupTable[A][2] * Kernel::LookupTable[B][b] * (*coeff), ++coeff;
    val += Kernel::LookupTable[A][3] * Kernel::LookupTable[B][b] * (*coeff), ++coeff;
  }

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetInside3D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));

  --i, --j, --k;

  RealType        val   = voxel_cast<RealType>(0);
  const RealType *coeff = _Coefficient.Data(i, j, k, l);

  Real wyz;
  for (int c = 0; c <= 3; ++c, coeff += _s3) {
    for (int b = 0; b <= 3; ++b, coeff += _s2) {
      wyz  = Kernel::LookupTable[B][b] * Kernel::LookupTable[C][c];
      val += Kernel::LookupTable[A][0] * wyz * (*coeff), ++coeff;
      val += Kernel::LookupTable[A][1] * wyz * (*coeff), ++coeff;
      val += Kernel::LookupTable[A][2] * wyz * (*coeff), ++coeff;
      val += Kernel::LookupTable[A][3] * wyz * (*coeff), ++coeff;
    }
  }

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetInside4D(double x, double y, double z, double t) const
{
  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));
  const int D = Kernel::VariableToIndex(Real(t - l));

  --i, --j, --k, --l;

  RealType        val   = voxel_cast<RealType>(0);
  const RealType *coeff = _Coefficient.Data(i, j, k, l);

  Real wzt, wyzt;
  for (int d = 0; d <= 3; ++d, coeff += _s4) {
    for (int c = 0; c <= 3; ++c, coeff += _s3) {
      wzt = Kernel::LookupTable[C][c] * Kernel::LookupTable[D][d];
      for (int b = 0; b <= 3; ++b, coeff += _s2) {
        wyzt = Kernel::LookupTable[B][b] * wzt;
        val += Kernel::LookupTable[A][0] * wyzt * (*coeff), ++coeff;
        val += Kernel::LookupTable[A][1] * wyzt * (*coeff), ++coeff;
        val += Kernel::LookupTable[A][2] * wyzt * (*coeff), ++coeff;
        val += Kernel::LookupTable[A][3] * wyzt * (*coeff), ++coeff;
      }
    }
  }

  return voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetInside(double x, double y, double z, double t) const
{
  // Use faster coefficient iteration than Get(Coefficient(), x, y, z, t)
  switch (this->NumberOfDimensions()) {
    case 3:  return GetInside3D(x, y, z, t);
    case 2:  return GetInside2D(x, y, z, t);
    default: return GetInside4D(x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetOutside(double x, double y, double z, double t) const
{
  if (_InfiniteCoefficient) {
    return voxel_cast<VoxelType>(Get(_InfiniteCoefficient, x, y, z, t));
  } else {
    return Get(x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetWithPaddingInside(double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return voxel_cast<VoxelType>(GetWithPadding3D(this->Input(), &_Coefficient, x, y, z, t));
    case 2:  return voxel_cast<VoxelType>(GetWithPadding2D(this->Input(), &_Coefficient, x, y, z, t));
    default: return voxel_cast<VoxelType>(GetWithPadding4D(this->Input(), &_Coefficient, x, y, z, t));
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
inline typename GenericFastCubicBSplineInterpolateImageFunction<TImage>::VoxelType
GenericFastCubicBSplineInterpolateImageFunction<TImage>
::GetWithPaddingOutside(double x, double y, double z, double t) const
{
  if (this->Extrapolator() && _InfiniteCoefficient) {
    return voxel_cast<VoxelType>(GetWithPadding(this->Extrapolator(), _InfiniteCoefficient, x, y, z, t));
  } else {
    return GetWithPadding(x, y, z, t);
  }
}

// =============================================================================
// Evaluation of derivatives
// =============================================================================

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Jacobian2D(Matrix &jac, double x, double y, double z, double t) const
{
  jac.Initialize(_Coefficient.N(), 2);

  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  if (k < 0 || k >= _Coefficient.Z() ||
      l < 0 || l >= _Coefficient.T()) {
    return;
  }

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));

  --i, --j;

  RealType dx = voxel_cast<RealType>(0);
  RealType dy = voxel_cast<RealType>(0);
  Real     wx[2], wy[2];

  int ia, jb;
  for (int b = 0; b <= 3; ++b) {
    jb = j + b;
    DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
    wy[0] = Kernel::LookupTable  [B][b];
    wy[1] = Kernel::LookupTable_I[B][b];
    for (int a = 0; a <= 3; ++a) {
      ia = i + a;
      DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
      wx[0] = Kernel::LookupTable  [A][a];
      wx[1] = Kernel::LookupTable_I[A][a];
      dx += (wx[1] * wy[0]) * _Coefficient(ia, jb, k, l);
      dy += (wx[0] * wy[1]) * _Coefficient(ia, jb, k, l);
    }
  }

  for (int r = 0; r < _Coefficient.N(); ++r) {
    jac(r, 0) = get(dx, r);
    jac(r, 1) = get(dy, r);
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TCoefficient>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Jacobian2D(Matrix &jac, const TCoefficient *coeff,
             double x, double y, double z, double t) const
{
  jac.Initialize(coeff->N(), 2);

  int i = ifloor(x);
  int j = ifloor(y);
  int k = iround(z);
  int l = iround(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));

  --i, --j;

  RealType dx = voxel_cast<RealType>(0);
  RealType dy = voxel_cast<RealType>(0);
  Real     wx[2], wy[2];

  int ia, jb;
  for (int b = 0; b <= 3; ++b) {
    jb = j + b;
    wy[0] = Kernel::LookupTable  [B][b];
    wy[1] = Kernel::LookupTable_I[B][b];
    for (int a = 0; a <= 3; ++a) {
      ia = i + a;
      wx[0] = Kernel::LookupTable  [A][a];
      wx[1] = Kernel::LookupTable_I[A][a];
      dx += (wx[1] * wy[0]) * voxel_cast<RealType>(coeff->Get(ia, jb, k, l));
      dy += (wx[0] * wy[1]) * voxel_cast<RealType>(coeff->Get(ia, jb, k, l));
    }
  }

  for (int r = 0; r < coeff->N(); ++r) {
    jac(r, 0) = get(dx, r);
    jac(r, 1) = get(dy, r);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Jacobian3D(Matrix &jac, double x, double y, double z, double t) const
{
  jac.Initialize(_Coefficient.N(), 3);

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  if (l < 0 || l >= _Coefficient.T()) {
    return;
  }

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));

  --i, --j, --k;

  RealType dx = voxel_cast<RealType>(0);
  RealType dy = voxel_cast<RealType>(0);
  RealType dz = voxel_cast<RealType>(0);
  Real     wx[2], wy[2], wz[2];

  int ia, jb, kc;
  for (int c = 0; c <= 3; ++c) {
    kc = k + c;
    DefaultExtrapolator::Apply(kc, _Coefficient.Z() - 1);
    wz[0] = Kernel::LookupTable  [C][c];
    wz[1] = Kernel::LookupTable_I[C][c];
    for (int b = 0; b <= 3; ++b) {
      jb = j + b;
      DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
      wy[0] = Kernel::LookupTable  [B][b];
      wy[1] = Kernel::LookupTable_I[B][b];
      for (int a = 0; a <= 3; ++a) {
        ia = i + a;
        DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
        wx[0] = Kernel::LookupTable  [A][a];
        wx[1] = Kernel::LookupTable_I[A][a];
        dx += (wx[1] * wy[0] * wz[0]) * _Coefficient(ia, jb, kc, l);
        dy += (wx[0] * wy[1] * wz[0]) * _Coefficient(ia, jb, kc, l);
        dz += (wx[0] * wy[0] * wz[1]) * _Coefficient(ia, jb, kc, l);
      }
    }
  }

  for (int r = 0; r < _Coefficient.N(); ++r) {
    jac(r, 0) = get(dx, r);
    jac(r, 1) = get(dy, r);
    jac(r, 2) = get(dz, r);
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Jacobian3D(Matrix &jac, const TOtherImage *coeff,
             double x, double y, double z, double t) const
{
  jac.Initialize(coeff->N(), 3);

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = iround(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));

  --i, --j, --k;

  RealType dx = voxel_cast<RealType>(0);
  RealType dy = voxel_cast<RealType>(0);
  RealType dz = voxel_cast<RealType>(0);
  Real     wx[2], wy[2], wz[2];

  int ia, jb, kc;
  for (int c = 0; c <= 3; ++c) {
    kc = k + c;
    wz[0] = Kernel::LookupTable  [C][c];
    wz[1] = Kernel::LookupTable_I[C][c];
    for (int b = 0; b <= 3; ++b) {
      jb = j + b;
      wy[0] = Kernel::LookupTable  [B][b];
      wy[1] = Kernel::LookupTable_I[B][b];
      for (int a = 0; a <= 3; ++a) {
        ia = i + a;
        wx[0] = Kernel::LookupTable  [A][a];
        wx[1] = Kernel::LookupTable_I[A][a];
        dx += (wx[1] * wy[0] * wz[0]) * voxel_cast<RealType>(coeff->Get(ia, jb, kc, l));
        dy += (wx[0] * wy[1] * wz[0]) * voxel_cast<RealType>(coeff->Get(ia, jb, kc, l));
        dz += (wx[0] * wy[0] * wz[1]) * voxel_cast<RealType>(coeff->Get(ia, jb, kc, l));
      }
    }
  }

  for (int r = 0; r < coeff->N(); ++r) {
    jac(r, 0) = get(dx, r);
    jac(r, 1) = get(dy, r);
    jac(r, 2) = get(dz, r);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Jacobian4D(Matrix &jac, double x, double y, double z, double t) const
{
  jac.Initialize(_Coefficient.N(), 4);

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));
  const int D = Kernel::VariableToIndex(Real(t - l));

  --i, --j, --k, --l;

  RealType dx = voxel_cast<RealType>(0);
  RealType dy = voxel_cast<RealType>(0);
  RealType dz = voxel_cast<RealType>(0);
  RealType dt = voxel_cast<RealType>(0);
  Real     wx[2], wy[2], wz[2], wt[2];

  int ia, jb, kc, ld;
  for (int d = 0; d <= 3; ++d) {
    ld = l + d;
    DefaultExtrapolator::Apply(ld, _Coefficient.T() - 1);
    wt[0] = Kernel::LookupTable  [D][d];
    wt[1] = Kernel::LookupTable_I[D][d];
    for (int c = 0; c <= 3; ++c) {
      kc = k + c;
      DefaultExtrapolator::Apply(kc, _Coefficient.Z() - 1);
      wz[0] = Kernel::LookupTable  [C][c];
      wz[1] = Kernel::LookupTable_I[C][c];
      for (int b = 0; b <= 3; ++b) {
        jb = j + b;
        DefaultExtrapolator::Apply(jb, _Coefficient.Y() - 1);
        wy[0] = Kernel::LookupTable  [B][b];
        wy[1] = Kernel::LookupTable_I[B][b];
        for (int a = 0; a <= 3; ++a) {
          ia = i + a;
          DefaultExtrapolator::Apply(ia, _Coefficient.X() - 1);
          wx[0] = Kernel::LookupTable  [A][a];
          wx[1] = Kernel::LookupTable_I[A][a];
          dx += (wx[1] * wy[0] * wz[0] * wt[0]) * _Coefficient(ia, jb, kc, ld);
          dy += (wx[0] * wy[1] * wz[0] * wt[0]) * _Coefficient(ia, jb, kc, ld);
          dz += (wx[0] * wy[0] * wz[1] * wt[0]) * _Coefficient(ia, jb, kc, ld);
          dt += (wx[0] * wy[0] * wz[0] * wt[1]) * _Coefficient(ia, jb, kc, ld);
        }
      }
    }
  }

  for (int r = 0; r < _Coefficient.N(); ++r) {
    jac(r, 0) = get(dx, r);
    jac(r, 1) = get(dy, r);
    jac(r, 2) = get(dz, r);
    jac(r, 3) = get(dt, r);
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Jacobian4D(Matrix &jac, const TOtherImage *coeff,
             double x, double y, double z, double t) const
{
  jac.Initialize(coeff->N(), 4);

  int i = ifloor(x);
  int j = ifloor(y);
  int k = ifloor(z);
  int l = ifloor(t);

  const int A = Kernel::VariableToIndex(Real(x - i));
  const int B = Kernel::VariableToIndex(Real(y - j));
  const int C = Kernel::VariableToIndex(Real(z - k));
  const int D = Kernel::VariableToIndex(Real(t - l));

  --i, --j, --k, --l;

  RealType dx = voxel_cast<RealType>(0);
  RealType dy = voxel_cast<RealType>(0);
  RealType dz = voxel_cast<RealType>(0);
  RealType dt = voxel_cast<RealType>(0);
  Real     wx[2], wy[2], wz[2], wt[2];

  int ia, jb, kc, ld;
  for (int d = 0; d <= 3; ++d) {
    ld = l + d;
    wt[0] = Kernel::LookupTable  [D][d];
    wt[1] = Kernel::LookupTable_I[D][d];
    for (int c = 0; c <= 3; ++c) {
      kc = k + c;
      wz[0] = Kernel::LookupTable  [C][c];
      wz[1] = Kernel::LookupTable_I[C][c];
      for (int b = 0; b <= 3; ++b) {
        jb = j + b;
        wy[0] = Kernel::LookupTable  [B][b];
        wy[1] = Kernel::LookupTable_I[B][b];
        for (int a = 0; a <= 3; ++a) {
          ia = i + a;
          wx[0] = Kernel::LookupTable  [A][a];
          wx[1] = Kernel::LookupTable_I[A][a];
          dx += (wx[1] * wy[0] * wz[0] * wt[0]) * voxel_cast<RealType>(coeff->Get(ia, jb, kc, ld));
          dy += (wx[0] * wy[1] * wz[0] * wt[0]) * voxel_cast<RealType>(coeff->Get(ia, jb, kc, ld));
          dz += (wx[0] * wy[0] * wz[1] * wt[0]) * voxel_cast<RealType>(coeff->Get(ia, jb, kc, ld));
          dt += (wx[0] * wy[0] * wz[0] * wt[1]) * voxel_cast<RealType>(coeff->Get(ia, jb, kc, ld));
        }
      }
    }
  }

  for (int r = 0; r < _Coefficient.N(); ++r) {
    jac(r, 0) = get(dx, r);
    jac(r, 1) = get(dy, r);
    jac(r, 2) = get(dz, r);
    jac(r, 3) = get(dt, r);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Jacobian(Matrix &jac, double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return Jacobian3D(jac, x, y, z, t);
    case 2:  return Jacobian2D(jac, x, y, z, t);
    default: return Jacobian4D(jac, x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage> template <class TOtherImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::Jacobian(Matrix &jac, const TOtherImage *coeff,
           double x, double y, double z, double t) const
{
  switch (this->NumberOfDimensions()) {
    case 3:  return Jacobian3D(jac, coeff, x, y, z, t);
    case 2:  return Jacobian2D(jac, coeff, x, y, z, t);
    default: return Jacobian4D(jac, coeff, x, y, z, t);
  }
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::EvaluateJacobianInside(Matrix &jac, double x, double y, double z, double t) const
{
  Jacobian(jac, &_Coefficient, x, y, z, t);
}

// -----------------------------------------------------------------------------
template <class TImage>
void GenericFastCubicBSplineInterpolateImageFunction<TImage>
::EvaluateJacobianOutside(Matrix &jac, double x, double y, double z, double t) const
{
  Jacobian(jac, _InfiniteCoefficient, x, y, z, t);
}


} // namespace mirtk

#endif // MIRTK_FastCubicBSplineInterpolateImageFunction_HXX
