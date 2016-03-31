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

#include "mirtk/Assert.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Voxel.h"
#include "mirtk/BaseImage.h"
#include "mirtk/GenericImage.h"
#include "mirtk/VoxelFunction.h"
#include "mirtk/InterpolateImageFunction.h"
#include "mirtk/LieBracketImageFilter.h"
#include "mirtk/DifferenceOfCompositionLieBracketImageFilter3D.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/Memory.h"


namespace mirtk {

// =============================================================================
// Lie bracket voxel function (scalar-valued voxel type)
// =============================================================================

/**
 * Evaluates Lie bracket of 3D vector fields using difference of composition
 */
template <class ScalarType>
struct DifferenceOfCompositionLieBracket : public VoxelFunction
{
  typedef typename voxel_info<ScalarType>::RealType      RealType;
  typedef typename WorldCoordsImage::VoxelType           CoordType;
  typedef GenericImage<ScalarType>                       VectorImage;
  typedef GenericInterpolateImageFunction<VectorImage>   InputType;

  /// Constructor
  DifferenceOfCompositionLieBracket(InputType *X, InputType *Y)
  :
    _LeftField (X), _LeftScaling (1.0),
    _RightField(Y), _RightScaling(1.0),
    _Domain    (X->Input()),
    _y         (_Domain->NumberOfSpatialVoxels()),
    _z         (_y + _y)
  {
    if (!Y->Input()->HasSpatialAttributesOf(X->Input())) {
      cerr << "DifferenceOfCompositionLieBracket: ";
      cerr << "Vector fields must be defined on the same finite lattice" << endl;
      exit(1);
    }
    if (X->Input()->T() != 3 || Y->Input()->T() != 3) {
      cerr << "DifferenceOfCompositionLieBracket: ";
      cerr << "3D+t vector fields must have temporal dimension _t = 3" << endl;
      exit(1);
    }
  }

  /// Constructor
  DifferenceOfCompositionLieBracket(double sx, InputType *X,
                                    double sy, InputType *Y)
  :
    _LeftField (X), _LeftScaling (static_cast<RealType>(sx)),
    _RightField(Y), _RightScaling(static_cast<RealType>(sy)),
    _Domain    (X->Input()),
    _y         (_Domain->NumberOfSpatialVoxels()),
    _z         (_y + _y)
  {
    if (!Y->Input()->HasSpatialAttributesOf(X->Input())) {
      cerr << "DifferenceOfCompositionLieBracket: ";
      cerr << "Vector fields must be defined on the same finite lattice" << endl;
      exit(1);
    }
    if (X->Input()->T() != 3 || Y->Input()->T() != 3) {
      cerr << "DifferenceOfCompositionLieBracket: ";
      cerr << "3D+t vector fields must have temporal dimension _t = 3" << endl;
      exit(1);
    }
  }

  /// Calculate Lie bracket for single voxel
  void Evaluate(int i, int j, int k, CoordType wx, CoordType wy, CoordType wz,
                const ScalarType *lf, const ScalarType *rf,
                RealType &ox, RealType &oy, RealType &oz)
  {
    double d[3], x1, y1, z1, x2, y2, z2;
    // Transform point x by lf o rf
    _RightField->Evaluate(d, i, j, k);
    x1 = wx + _RightScaling * d[0];
    y1 = wy + _RightScaling * d[1];
    z1 = wz + _RightScaling * d[2];
    d[0] = x1, d[1] = y1, d[2] = z1;
    _Domain->WorldToImage(d[0], d[1], d[2]);
    _LeftField->Evaluate(d, d[0], d[1], d[2]);
    x1 += _LeftScaling * d[0];
    y1 += _LeftScaling * d[1];
    z1 += _LeftScaling * d[2];
    // Transform point x by rf o lf
    _LeftField->Evaluate(d, i, j, k);
    x2 = wx + _LeftScaling * d[0];
    y2 = wy + _LeftScaling * d[1];
    z2 = wz + _LeftScaling * d[2];
    d[0] = x2, d[1] = y2, d[2] = z2;
    _Domain->WorldToImage(d[0], d[1], d[2]);
    _RightField->Evaluate(d, d[0], d[1], d[2]);
    x2 += _RightScaling * d[0];
    y2 += _RightScaling * d[1];
    z2 += _RightScaling * d[2];
    // Compute the difference
    ox = static_cast<RealType>(x1 - x2);
    oy = static_cast<RealType>(y1 - y2);
    oz = static_cast<RealType>(z1 - z2);
  }

  /// Calculate Lie bracket for single voxel
  void Evaluate(int i, int j, int k, const ScalarType *lf, const ScalarType *rf,
                RealType &ox, RealType &oy, RealType &oz)
  {
    double wx = i, wy = j, wz = k;
    _Domain->ImageToWorld(wx, wy, wz);
    this->Evaluate(i, j, k, wx, wy, wz, lf, rf, ox, oy, oz);
  }

  /// Calculate Lie bracket for single voxel
  void Evaluate(int i, int j, int k, RealType &ox, RealType &oy, RealType &oz)
  {
    const ScalarType * const lf = _LeftField ->Input()->Data(i, j, k);
    const ScalarType * const rf = _RightField->Input()->Data(i, j, k);
    this->Evaluate(i, j, k, lf, rf, ox, oy, oz);
  }

  /// Calculate Lie bracket for single voxel without pre-computed world coordinates
  void operator ()(int i, int j, int k, int, const ScalarType *lf, const ScalarType *rf, RealType *out)
  {
    double wx = i, wy = j, wz = k;
    _Domain->ImageToWorld(wx, wy, wz);
    this->Evaluate(i, j, k, wx, wy, wz, lf, rf, out[_x], out[_y], out[_z]);
  }

  /// Calculate Lie bracket for single voxel using pre-computed world coordinates
  void operator ()(int i, int j, int k, int, const CoordType *w, const ScalarType *lf, const ScalarType *rf, RealType *out)
  {
    this->Evaluate(i, j, k, w[_x], w[_y], w[_z], lf, rf, out[_x], out[_y], out[_z]);
  }

  /// Calculate Lie bracket for single voxel using pre-computed world coordinates
  void operator ()(const VectorImage &img, int idx, const CoordType *w, const ScalarType *lf, const ScalarType *rf, RealType *out)
  {
    int i, j, k;
    img.IndexToLattice(idx, i, j, k);
    this->Evaluate(i, j, k, w[_x], w[_y], w[_z], lf, rf, out[_x], out[_y], out[_z]);
  }

protected:

  const InputType   *_LeftField;    ///< Continuous left  input vector field (X)
  RealType           _LeftScaling;  ///< Scaling of left  input vector field (X)
  const InputType   *_RightField;   ///< Continuous right input vector field (Y)
  RealType           _RightScaling; ///< Scaling of right input vector field (Y)
  const VectorImage *_Domain;       ///< Defines common attributes of vector fields

  static const int _x = 0;          ///< Temporal offset of first  vector component
  int              _y;              ///< Temporal offset of second vector component
  int              _z;              ///< Temporal offset of third  vector component

};

// =============================================================================
// Lie bracket voxel function (vector-valued voxel type)
// =============================================================================

/**
 * Evaluates Lie bracket of 3D vector fields using difference of composition
 */
template <class VectorType>
struct DifferenceOfCompositionLieBracket3D : public VoxelFunction
{
  // TODO: Change WorldCoordsImage to use real3 as voxel type
  typedef typename voxel_info<VectorType>::ScalarType    ScalarType;
  typedef typename voxel_info<ScalarType>::RealType      RealType;
  typedef typename WorldCoordsImage::VoxelType           CoordType;
  typedef Vector3D<CoordType>                            CoordVector;
  typedef GenericImage<VectorType>                       VectorImage;
  typedef GenericInterpolateImageFunction<VectorImage>   InputType;

  /// Constructor
  DifferenceOfCompositionLieBracket3D(const InputType *X, const InputType *Y)
  :
    _LeftField (X), _LeftScaling (RealType(1.0)),
    _RightField(Y), _RightScaling(RealType(1.0)),
    _Domain    (X->Input())
  {
    if (X->Input()->Attributes() != Y->Input()->Attributes()) {
      cerr << "DifferenceOfCompositionLieBracket3D: ";
      cerr << "Vector fields must be defined on the same finite lattice" << endl;
      exit(1);
    }
  }

  /// Constructor
  DifferenceOfCompositionLieBracket3D(double sx, const InputType *X,
                                      double sy, const InputType *Y)
  :
    _LeftField (X), _LeftScaling (static_cast<RealType>(sx)),
    _RightField(Y), _RightScaling(static_cast<RealType>(sy)),
    _Domain    (X->Input())
  {
    if (X->Input()->Attributes() != Y->Input()->Attributes()) {
      cerr << "DifferenceOfCompositionLieBracket3D: ";
      cerr << "Vector fields must be defined on the same finite lattice" << endl;
      exit(1);
    }
  }

  /// Evaluate Lie bracket at single voxel
  void Evaluate(int i, int j, int k, const CoordVector &wc, VectorType &lb)
  {
    CoordVector p1, p2, p;
    // Transform point x by lf o rf
    p1  = p = wc + _RightScaling * _RightField->Get(i, j, k);
    _Domain->WorldToImage(p._x, p._y, p._z);
    p1 += _LeftScaling * _LeftField->Get(p._x, p._y, p._z);
    // Transform point x by rf o lf
    p2  = p = wc + _LeftScaling * _LeftField->Get(i, j, k);
    _Domain ->WorldToImage(p._x, p._y, p._z);
    p2 += _RightScaling * _RightField->Get(p._x, p._y, p._z);
    // Compute the difference
    lb = p1 - p2;
  }

  /// Evaluate Lie bracket at single voxel
  void Evaluate(int i, int j, int k, VectorType &lb)
  {
    CoordVector wc(i, j, k);
    _Domain->ImageToWorld(wc._x, wc._y, wc._z);
    Evaluate(i, j, k, wc, lb);
  }

  /// Evaluate Lie bracket at single voxel without pre-computed world coordinates
  void operator ()(int i, int j, int k, int, VectorType *lb)
  {
    CoordVector wc(i, j, k);
    _Domain->ImageToWorld(wc._x, wc._y, wc._z);
    Evaluate(i, j, k, wc, *lb);
  }

  /// Evaluate Lie bracket at single voxel using pre-computed world coordinates
  void operator ()(const VectorImage &img, int idx, VectorType *lb)
  {
    int i, j, k;
    img.IndexToLattice(idx, i, j, k);
    CoordVector wc(i, j, k);
    _Domain->ImageToWorld(wc._x, wc._y, wc._z);
    Evaluate(i, j, k, wc, *lb);
  }

  /// Evaluate Lie bracket at single voxel using pre-computed world coordinates
  void operator ()(int i, int j, int k, int, const CoordVector *wc, VectorType *lb)
  {
    Evaluate(i, j, k, *wc, *lb);
  }

  /// Evaluate Lie bracket at single voxel using pre-computed world coordinates
  void operator ()(const VectorImage &img, int idx, const CoordVector *wc, VectorType *lb)
  {
    int i, j, k;
    img.IndexToLattice(idx, i, j, k);
    Evaluate(i, j, k, *wc, *lb);
  }

protected:

  const InputType   *_LeftField;    ///< Continuous left  input vector field (X)
  RealType           _LeftScaling;  ///< Scaling of left  input vector field (X)
  const InputType   *_RightField;   ///< Continuous right input vector field (Y)
  RealType           _RightScaling; ///< Scaling of right input vector field (Y)
  const VectorImage *_Domain;       ///< Defines common attributes of vector fields

};

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
DifferenceOfCompositionLieBracketImageFilter3D<VoxelType>
::DifferenceOfCompositionLieBracketImageFilter3D()
:
  _Interpolation(Interpolation_Linear),
  _Extrapolation(Extrapolation_NN),
  _ComputeInterpolationCoefficients(true)
{
  _Interpolator[0] = NULL, _Scaling[0] = 1.0;
  _Interpolator[1] = NULL, _Scaling[1] = 1.0;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
DifferenceOfCompositionLieBracketImageFilter3D<VoxelType>
::~DifferenceOfCompositionLieBracketImageFilter3D()
{
  Delete(_Interpolator[0]);
  Delete(_Interpolator[1]);
}

// =============================================================================
// Filter implementation
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void DifferenceOfCompositionLieBracketImageFilter3D<VoxelType>::Scaling(int i, double s)
{
  mirtkAssert(0 <= i && i < 2, "index is either 0 or 1");
  _Scaling[i] = s;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
double DifferenceOfCompositionLieBracketImageFilter3D<VoxelType>::Scaling(int i) const
{
  mirtkAssert(0 <= i && i < 2, "index is either 0 or 1");
  return _Scaling[i];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void DifferenceOfCompositionLieBracketImageFilter3D<VoxelType>::Initialize()
{
  // Initialize base class
  LieBracketImageFilter<VoxelType>::Initialize();

  // Initialize input field interpolators
  for (int i = 0; i < 2; ++i) {
    if (_Interpolator[i]) delete _Interpolator[i];
    _Interpolator[i] = InterpolatorType::New(_Interpolation, _Extrapolation, this->Input(i));
    _Interpolator[i]->Input(this->Input(i));
    _Interpolator[i]->Initialize(_ComputeInterpolationCoefficients == false);
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void DifferenceOfCompositionLieBracketImageFilter3D<VoxelType>
::Run(double out[3], int i, int j, int k)
{
  // FIXME: User cannot initialize interpolators
  typedef DifferenceOfCompositionLieBracket<VoxelType> LieBracketFunction;
  LieBracketFunction liebracket(_Scaling[0], _Interpolator[0],
                                _Scaling[1], _Interpolator[1]);
  VoxelType ox, oy, oz;
  liebracket.Evaluate(i, j, k, ox, oy, oz);
  out[0] = static_cast<double>(ox);
  out[1] = static_cast<double>(oy);
  out[2] = static_cast<double>(oz);
}

// -----------------------------------------------------------------------------
template <>
void DifferenceOfCompositionLieBracketImageFilter3D<Float3>
::Run(double out[3], int i, int j, int k)
{
  // FIXME: User cannot initialize interpolators
  typedef DifferenceOfCompositionLieBracket3D<Float3> LieBracketFunction;
  LieBracketFunction liebracket(_Scaling[0], _Interpolator[0],
                                _Scaling[1], _Interpolator[1]);
  Float3 tmp;
  liebracket.Evaluate(i, j, k, tmp);
  out[0] = static_cast<double>(tmp._x);
  out[1] = static_cast<double>(tmp._y);
  out[2] = static_cast<double>(tmp._z);
}

// -----------------------------------------------------------------------------
template <>
void DifferenceOfCompositionLieBracketImageFilter3D<Double3>
::Run(double out[3], int i, int j, int k)
{
  // FIXME: User cannot initialize interpolators
  typedef DifferenceOfCompositionLieBracket3D<Double3> LieBracketFunction;
  LieBracketFunction liebracket(_Scaling[0], _Interpolator[0],
                                _Scaling[1], _Interpolator[1]);
  Double3 tmp;
  liebracket.Evaluate(i, j, k, tmp);
  out[0] = tmp._x;
  out[1] = tmp._y;
  out[2] = tmp._z;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
double DifferenceOfCompositionLieBracketImageFilter3D<VoxelType>
::Run(int i, int j, int k, int t)
{
  // FIXME: User cannot initialize interpolators
  double out[3];
  this->Run(out, i, j, k);
  return out[t];
}

// ---------------------------------------------------------------------------
template <class VoxelType>
void DifferenceOfCompositionLieBracketImageFilter3D<VoxelType>::Run()
{
  MIRTK_START_TIMING();
  this->Initialize();
  typedef DifferenceOfCompositionLieBracket<VoxelType> LieBracketFunction;
  LieBracketFunction liebracket(_Scaling[0], _Interpolator[0],
                                _Scaling[1], _Interpolator[1]);
  blocked_range3d<int> roi(0, this->Output()->Z(),
                           0, this->Output()->Y(),
                           0, this->Output()->X());
  ParallelForEachVoxel(roi, this->GetInput(0), this->GetInput(1), this->Output(), liebracket);
  this->Finalize();
  MIRTK_DEBUG_TIMING(2, this->NameOfClass());
}

// ---------------------------------------------------------------------------
template <>
void DifferenceOfCompositionLieBracketImageFilter3D<Float3>::Run()
{
  MIRTK_START_TIMING();
  this->Initialize();
  typedef DifferenceOfCompositionLieBracket3D<Float3> LieBracketFunction;
  LieBracketFunction liebracket(_Scaling[0], _Interpolator[0],
                                _Scaling[1], _Interpolator[1]);
  blocked_range3d<int> roi(0, this->Output()->Z(),
                           0, this->Output()->Y(),
                           0, this->Output()->X());
  ParallelForEachVoxel(roi, this->Output(), liebracket);
  this->Finalize();
  MIRTK_DEBUG_TIMING(2, this->NameOfClass());
}

// ---------------------------------------------------------------------------
template <>
void DifferenceOfCompositionLieBracketImageFilter3D<Double3>::Run()
{
  MIRTK_START_TIMING();
  this->Initialize();
  typedef DifferenceOfCompositionLieBracket3D<Double3> LieBracketFunction;
  LieBracketFunction liebracket(_Scaling[0], _Interpolator[0],
                                _Scaling[1], _Interpolator[1]);
  blocked_range3d<int> roi(0, this->Output()->Z(),
                           0, this->Output()->Y(),
                           0, this->Output()->X());
  ParallelForEachVoxel(roi, this->Output(), liebracket);
  this->Finalize();
  MIRTK_DEBUG_TIMING(2, this->NameOfClass());
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class DifferenceOfCompositionLieBracketImageFilter3D<float>;
template class DifferenceOfCompositionLieBracketImageFilter3D<double>;
template class DifferenceOfCompositionLieBracketImageFilter3D<Float3>;
template class DifferenceOfCompositionLieBracketImageFilter3D<Double3>;


} // namespace mirtk
