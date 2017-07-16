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

#ifndef MIRTK_ImageToInterpolationCoefficients_H
#define MIRTK_ImageToInterpolationCoefficients_H

#include "mirtk/Math.h"
#include "mirtk/Voxel.h"
#include "mirtk/GenericImage.h"
#include "mirtk/Parallel.h"
#include "mirtk/Stream.h"
#include "mirtk/Queue.h"


namespace mirtk {


// -----------------------------------------------------------------------------
template <class TData, class TReal>
TData InitialAntiCausalCoefficient(TData c[], int n, TReal z)
{
  // this initialization corresponds to mirror boundaries
  return ((z / (z * z - TReal(1))) * (z * c[n - 2] + c[n - 1]));
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
TData InitialCausalCoefficient(TData c[], int n, TReal z, TReal tol)
{
  TData sum;
  // this initialization corresponds to mirror boundaries
  int m;
  if (tol > 0) m = static_cast<int>(ceil(log(tol) / log(abs(z))));
  else         m = n;
  // accelerated loop
  if (m < n) {
    TReal zn = z;
    sum = c[0];
    for (int i = 1; i < m; ++i) {
      sum += zn * c[i];
      zn  *= z;
    }
  // full loop
  } else {
    const TReal iz  = TReal(1) / z;
    TReal       zn  = z;
    TReal       z2n = pow(z, TReal(n - 1));
    sum  = c[0] + z2n * c[n - 1];
    z2n *= z2n * iz;
    for (int i = 1; i <= n - 2; ++i) {
      sum += (zn + z2n) * c[i];
      zn  *= z;
      z2n *= iz;
    }
    sum /= (TReal(1) - zn * zn);
  }
  return sum;
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficients(TData *c, int n, const TReal *z, int npoles, TReal tol)
{
  // special case required by mirror boundaries
  if (n < 2) return;
  // compute the overall gain
  TReal lambda = TReal(1);
  for (int k = 0; k < npoles; ++k) {
    lambda *= (TReal(1) - z[k]) * (TReal(1) - TReal(1) / z[k]);
  }
  // apply the gain
  for (int i = 0; i < n; ++i) c[i] *= lambda;
  // loop over all poles
  for (int k = 0; k < npoles; ++k) {
    // causal initialization
    c[0] = InitialCausalCoefficient(c, n, z[k], tol);
    // causal recursion
    for (int i = 1; i < n; ++i) c[i] += z[k] * c[i - 1];
    // anticausal initialization
    c[n - 1] = InitialAntiCausalCoefficient(c, n, z[k]);
    // anticausal recursion
    for (int i = n - 2; i >= 0; --i) c[i] = z[k] * (c[i + 1] - c[i]);
  }
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
struct ConvertToInterpolationCoefficientsBase
{
  ConvertToInterpolationCoefficientsBase(GenericImage<TData> &image, const TReal *z, int npoles)
  :
    _image(image), _z(z), _npoles(npoles)
  {}

  void Process(TData *data, int n) const
  {
    ConvertToInterpolationCoefficients(data, n, _z, _npoles, static_cast<TReal>(DBL_EPSILON));
  }

protected:
  GenericImage<TData> &_image;
  const TReal         *_z;
  int                  _npoles;
};

// -----------------------------------------------------------------------------
template <class TData, class TReal>
struct ConvertToInterpolationCoefficientsXFunc : public ConvertToInterpolationCoefficientsBase<TData, TReal>
{
  ConvertToInterpolationCoefficientsXFunc(GenericImage<TData> &image, const TReal *z, int npoles)
  :
    ConvertToInterpolationCoefficientsBase<TData, TReal>(image, z, npoles)
  {}

  void operator()(const blocked_range3d<int> &re) const
  {
    TData *data = new TData[this->_image.X()];
    for (int l = re.pages().begin(); l != re.pages().end(); ++l)
    for (int k = re.rows ().begin(); k != re.rows ().end(); ++k)
    for (int j = re.cols ().begin(); j != re.cols ().end(); ++j) {
      for (int i = 0; i < this->_image.X(); ++i) {
        data[i] = this->_image(i, j, k, l);
      }
      this->Process(data, this->_image.X());
      for (int i = 0; i < this->_image.X(); ++i) {
        this->_image(i, j, k, l) = data[i];
      }
    }
    delete[] data;
  }

  void operator()(int l = -1)
  {
    int L = this->_image.T();
    if (0 <= l && l < this->_image.T()) L = l+1;
    else                                l = 0;
    blocked_range3d<int> re(l, L, 0, this->_image.Z(), 0, this->_image.Y());
    parallel_for(re, *this);
  }

  void operator()(int k, int l)
  {
    blocked_range3d<int> re(l, l+1, k, k+1, 0, this->_image.Y());
    parallel_for(re, *this);
  }
};

// -----------------------------------------------------------------------------
template <class TData, class TReal>
struct ConvertToInterpolationCoefficientsYFunc : public ConvertToInterpolationCoefficientsBase<TData, TReal>
{
  ConvertToInterpolationCoefficientsYFunc(GenericImage<TData> &image, const TReal *z, int npoles)
  :
    ConvertToInterpolationCoefficientsBase<TData, TReal>(image, z, npoles)
  {}

  void operator()(const blocked_range3d<int> &re) const
  {
    TData *data = new TData[this->_image.Y()];
    for (int l = re.pages().begin(); l != re.pages().end(); ++l)
    for (int k = re.rows ().begin(); k != re.rows ().end(); ++k)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      for (int j = 0; j < this->_image.Y(); ++j) {
        data[j] = this->_image(i, j, k, l);
      }
      this->Process(data, this->_image.Y());
      for (int j = 0; j < this->_image.Y(); ++j) {
        this->_image(i, j, k, l) = data[j];
      }
    }
    delete[] data;
  }

  void operator()(int l = -1)
  {
    int L = this->_image.T();
    if (0 <= l && l < this->_image.T()) L = l+1;
    else                                l = 0;
    blocked_range3d<int> re(l, L, 0, this->_image.Z(), 0, this->_image.X());
    parallel_for(re, *this);
  }

  void operator()(int k, int l)
  {
    blocked_range3d<int> re(l, l+1, k, k+1, 0, this->_image.X());
    parallel_for(re, *this);
  }
};

// -----------------------------------------------------------------------------
template <class TData, class TReal>
struct ConvertToInterpolationCoefficientsZFunc : public ConvertToInterpolationCoefficientsBase<TData, TReal>
{
  ConvertToInterpolationCoefficientsZFunc(GenericImage<TData> &image, const TReal *z, int npoles)
  :
    ConvertToInterpolationCoefficientsBase<TData, TReal>(image, z, npoles)
  {}

  void operator()(const blocked_range3d<int> &re) const
  {
    TData *data = new TData[this->_image.Z()];
    for (int l = re.pages().begin(); l != re.pages().end(); ++l)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      for (int k = 0; k < this->_image.Z(); ++k) {
        data[k] = this->_image(i, j, k, l);
      }
      this->Process(data, this->_image.Z());
      for (int k = 0; k < this->_image.Z(); ++k) {
        this->_image(i, j, k, l) = data[k];
      }
    }
    delete[] data;
  }

  void operator()(int l = -1)
  {
    int L = this->_image.T();
    if (0 <= l && l < this->_image.T()) L = l+1;
    else                                l = 0;
    blocked_range3d<int> re(l, L, 0, this->_image.Y(), 0, this->_image.X());
    parallel_for(re, *this);
  }
};

// -----------------------------------------------------------------------------
template <class TData, class TReal>
struct ConvertToInterpolationCoefficientsTFunc : public ConvertToInterpolationCoefficientsBase<TData, TReal>
{
  ConvertToInterpolationCoefficientsTFunc(GenericImage<TData> &image, const TReal *z, int npoles)
  :
    ConvertToInterpolationCoefficientsBase<TData, TReal>(image, z, npoles)
  {}

  void operator()(const blocked_range3d<int> &re) const
  {
    TData *data = new TData[this->_image.T()];
    for (int k = re.pages().begin(); k != re.pages().end(); ++k)
    for (int j = re.rows ().begin(); j != re.rows ().end(); ++j)
    for (int i = re.cols ().begin(); i != re.cols ().end(); ++i) {
      for (int l = 0; l < this->_image.T(); ++l) {
        data[l] = this->_image(i, j, k, l);
      }
      this->Process(data, this->_image.T());
      for (int l = 0; l < this->_image.T(); ++l) {
        this->_image(i, j, k, l) = data[l];
      }
    }
    delete[] data;
  }

  void operator()()
  {
    blocked_range3d<int> re(0, this->_image.Z(), 0, this->_image.Y(), 0, this->_image.X());
    parallel_for(re, *this);
  }
};

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficientsX(GenericImage<TData> &image, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsXFunc<TData, TReal> convert(image, z, npoles);
  convert();
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficientsY(GenericImage<TData> &image, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsYFunc<TData, TReal> convert(image, z, npoles);
  convert();
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficientsZ(GenericImage<TData> &image, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsZFunc<TData, TReal> convert(image, z, npoles);
  convert();
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficientsT(GenericImage<TData> &image, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsTFunc<TData, TReal> convert(image, z, npoles);
  convert();
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficients(GenericImage<TData> &image, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsX(image, z, npoles);
  ConvertToInterpolationCoefficientsY(image, z, npoles);
  ConvertToInterpolationCoefficientsZ(image, z, npoles);
  if (image.TSize() != .0 && image.T() > 1) {
    ConvertToInterpolationCoefficientsT(image, z, npoles);
  }
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficientsX(GenericImage<TData> &image, int l, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsXFunc<TData, TReal> convert(image, z, npoles);
  convert(l);
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficientsY(GenericImage<TData> &image, int l, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsYFunc<TData, TReal> convert(image, z, npoles);
  convert(l);
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficientsZ(GenericImage<TData> &image, int l, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsZFunc<TData, TReal> convert(image, z, npoles);
  convert(l);
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficients(GenericImage<TData> &image, int l, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsX(image, l, z, npoles);
  ConvertToInterpolationCoefficientsY(image, l, z, npoles);
  ConvertToInterpolationCoefficientsZ(image, l, z, npoles);
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficientsX(GenericImage<TData> &image, int k, int l, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsXFunc<TData, TReal> convert(image, z, npoles);
  convert(k, l);
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficientsY(GenericImage<TData> &image, int k, int l, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsYFunc<TData, TReal> convert(image, z, npoles);
  convert(k, l);
}

// -----------------------------------------------------------------------------
template <class TData, class TReal>
void ConvertToInterpolationCoefficients(GenericImage<TData> &image, int k, int l, const TReal *z, int npoles)
{
  ConvertToInterpolationCoefficientsX(image, k, l, z, npoles);
  ConvertToInterpolationCoefficientsY(image, k, l, z, npoles);
}

// -----------------------------------------------------------------------------
/// Recover spline poles from a lookup table
template <class TReal>
void SplinePoles(int degree, TReal pole[2], int &npoles)
{
  switch (degree) {
    case 2:
      npoles = 1;
      pole[0] = TReal(sqrt(8.0) - 3.0);
      break;
    case 3:
      npoles = 1;
      pole[0] = TReal(sqrt(3.0) - 2.0);
      break;
    case 4:
      npoles = 2;
      pole[0] = TReal(sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0);
      pole[1] = TReal(sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0);
      break;
    case 5:
      npoles = 2;
      pole[0] = TReal(sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0) - 13.0 / 2.0);
      pole[1] = TReal(sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0) - 13.0 / 2.0);
      break;
    default:
      cerr << "SplinePoles: Unsupported spline degree: " << degree << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
/// Convert 2D floating point scalar or vector image to spline coefficients
template <class TData>
void ConvertToSplineCoefficients(int degree, GenericImage<TData> &image, int k, int l)
{
  typename voxel_info<TData>::ScalarType poles[2];
  int                                    npoles;
  SplinePoles(degree, poles, npoles);
  ConvertToInterpolationCoefficients(image, k, l, poles, npoles);
}

// -----------------------------------------------------------------------------
/// Convert 3D floating point scalar or vector image to spline coefficients
template <class TData>
void ConvertToSplineCoefficients(int degree, GenericImage<TData> &image, int l)
{
  typename voxel_info<TData>::ScalarType poles[2];
  int                                    npoles;
  SplinePoles(degree, poles, npoles);
  ConvertToInterpolationCoefficients(image, l, poles, npoles);
}

// -----------------------------------------------------------------------------
/// Convert 4D floating point scalar or vector image to spline coefficients
template <class TData>
void ConvertToSplineCoefficients(int degree, GenericImage<TData> &image)
{
  typename voxel_info<TData>::ScalarType poles[2];
  int                                    npoles;
  SplinePoles(degree, poles, npoles);
  ConvertToInterpolationCoefficients(image, poles, npoles);
}

// -----------------------------------------------------------------------------
/// Convert 3D floating point scalar or vector image to cubic spline coefficients
template <class TData>
void ConvertToCubicBSplineCoefficients(GenericImage<TData> &image, int l)
{
  ConvertToSplineCoefficients(3, image, l);
}

// -----------------------------------------------------------------------------
/// Convert 4D floating point scalar or vector image to cubic spline coefficients
template <class TData>
void ConvertToCubicBSplineCoefficients(GenericImage<TData> &image)
{
  ConvertToSplineCoefficients(3, image);
}

// -----------------------------------------------------------------------------
/// Fill background by front propagation of foreground
template <class TData>
void FillBackgroundBeforeConversionToSplineCoefficients(GenericImage<TData> &image)
{
  int idx, nbr, i, j, k, l;
  Queue<int> active, next;
  BinaryImage bg(image.Attributes());
  for (idx = 0; idx < image.NumberOfVoxels(); ++idx) {
    if (image.IsBackground(idx)) {
      if (image.IsNextToForeground(idx)) {
        active.push(idx);
      }
      bg(idx) = 1;
    }
  }
  for (int gen = 0; gen < 3; ++gen) {
    while (!active.empty()) {
      idx = active.front();
      active.pop();
      if (bg(idx)) {
        int   count = 0;
        TData value = voxel_cast<TData>(0);
        image.IndexToVoxel(idx, i, j, k, l);
        for (int nl = l - 1; nl <= l + 1; ++nl)
        for (int nk = k - 1; nk <= k + 1; ++nk)
        for (int nj = j - 1; nj <= j + 1; ++nj)
        for (int ni = i - 1; ni <= i + 1; ++ni) {
          nbr = image.VoxelToIndex(ni, nj, nk, nl);
          if (nbr != idx && image.IsInside(nbr)) {
            if (bg(nbr)) {
              next.push(nbr);
            } else {
              value += image(nbr);
              count += 1;
            }
          }
        }
        value /= static_cast<float>(count);
        image.Put(idx, value);
        bg(idx) = 0;
      }
    }
    if (next.empty()) break;
    active = move(next);
  }
}


} // namespace mirtk

#endif // MIRTK_ImageToInterpolationCoefficients_H
