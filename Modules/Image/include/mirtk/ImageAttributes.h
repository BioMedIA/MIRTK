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

#ifndef MIRTK_ImageAttributes_H
#define MIRTK_ImageAttributes_H

#include "mirtk/Indent.h"
#include "mirtk/Matrix.h"
#include "mirtk/Point.h"
#include "mirtk/Array.h"


namespace mirtk {


/**
 * Class which defines the attributes of the imaging geometry
 */
struct ImageAttributes
{
  // ---------------------------------------------------------------------------
  // Attributes
  int _x; ///< Image x-dimension (in voxels)
  int _y; ///< Image y-dimension (in voxels)
  int _z; ///< Image z-dimension (in voxels)
  int _t; ///< Image t-dimension (in voxels)

  double _dx; ///< Voxel x-dimensions (in mm)
  double _dy; ///< Voxel y-dimensions (in mm)
  double _dz; ///< Voxel z-dimensions (in mm)
  double _dt; ///< Voxel t-dimensions (in ms)

  double _xorigin; ///< Image x-origin (in mm)
  double _yorigin; ///< Image y-origin (in mm)
  double _zorigin; ///< Image z-origin (in mm)
  double _torigin; ///< Image t-origin (in ms)

  double _xaxis[3]; ///< Direction of x-axis
  double _yaxis[3]; ///< Direction of y-axis
  double _zaxis[3]; ///< Direction of z-axis

  Matrix _smat; ///< Affine transformation matrix

  const Matrix *_w2i;  ///< Pointer to pre-computed world to image matrix (cf. BaseImage::_matW2I)
  const Matrix *_i2w;  ///< Pointer to pre-computed image to world matrix (cf. BaseImage::_matI2W)

  /// Number of lattice points in image domain
  int NumberOfLatticePoints() const;

  /// Number of spatial lattice points
  int NumberOfSpatialPoints() const;

  /// Number of lattice points
  /// (i.e., NumberOfLatticePoints if _dt != 0, otherwise NumberOfSpatialPoints)
  int NumberOfPoints() const;

  /// Get number of lattice points in x dimension
  int X() const;

  /// Get number of lattice points in y dimension
  int Y() const;

  /// Get number of lattice points in z dimension
  int Z() const;

  /// Get number of lattice points in t dimension
  int T() const;

  /// Get number of lattice points in i-th dimension
  int N(int) const;

  /// Get spacing of lattice points in i-th dimension
  double Spacing(int) const;

  /// Get spacing of lattice points in x dimension
  double XSpacing() const;

  /// Get spacing of lattice points in y dimension
  double YSpacing() const;

  /// Get spacing of lattice points in z dimension
  double ZSpacing() const;

  /// Get spacing of lattice points in t dimension
  double TSpacing() const;

  /// Get voxel size/spacing of lattice points in x dimension
  double XSize() const;

  /// Get voxel size/spacing of lattice points in y dimension
  double YSize() const;

  /// Get voxel size/spacing of lattice points in z dimension
  double ZSize() const;

  /// Get voxel size/spacing of lattice points in t dimension
  double TSize() const;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  ImageAttributes();

  /// Constructor
  ImageAttributes(int, int, double = 1.0, double = 1.0);

  /// Constructor
  ImageAttributes(int, int, int, double = 1.0, double = 1.0, double = 1.0);

  /// Constructor
  ImageAttributes(int, int, int, int, double = 1.0, double = 1.0, double = 1.0, double = 1.0);

  /// Copy constructor
  ImageAttributes(const ImageAttributes &);

  /// Copy operator
  ImageAttributes& operator= (const ImageAttributes &);

  /// Check if attributes are valid
  operator bool() const;

  // ---------------------------------------------------------------------------
  // Comparison

  /// Whether given lattice is fully contained by this image lattice
  bool ContainsInSpace(const ImageAttributes &attr) const;

  /// Whether spatial attributes are equal
  bool EqualInSpace(const ImageAttributes &attr) const;

  /// Whether temporal attributes are equal
  bool EqualInTime(const ImageAttributes &attr) const;

  /// Equality operator
  bool operator==(const ImageAttributes &attr) const;

  /// Inequality operator
  bool operator!=(const ImageAttributes &attr) const;

  // ---------------------------------------------------------------------------
  // Coordinate conversion

  /// Get Index from Lattice
  int LatticeToIndex(int, int, int = 0, int = 0) const;

  /// Get Index from Lattice
  void IndexToLattice(int, int *, int *, int * = NULL, int * = NULL) const;

  /// Get Index from Lattice
  void IndexToLattice(int, int &, int &) const;

  /// Get Index from Lattice
  void IndexToLattice(int, int &, int &, int &) const;

  /// Get Index from Lattice
  void IndexToLattice(int, int &, int &, int &, int &) const;

  /// Get world coordinates (in mm) of lattice point
  void IndexToWorld(int, double &, double &) const;

  /// Get world coordinates (in mm) of lattice point
  void IndexToWorld(int, double &, double &, double &) const;

  /// Get world coordinates (in mm) of lattice point
  void IndexToWorld(int, Point &) const;

  /// Get world coordinates (in mm) of lattice point
  Point IndexToWorld(int) const;

  /// Convert lattice to world coordinate
  void LatticeToWorld(Point &) const;

  /// Convert lattice to world coordinate
  void LatticeToWorld(double &, double &, double &) const;

  /// Compute spatial world coordinates of lattice points (_x * _y * _z)
  void LatticeToWorld(double *, double *, double *) const;

  /// Compute world coordinates of all lattice points (_x * _y * _z * _t)
  void LatticeToWorld(double *, double *, double *, double *) const;

  /// Convert world to lattice coordinate
  void WorldToLattice(double &, double &, double &) const;

  /// Convert lattice to world coordinate
  void WorldToLattice(Point &) const;

  /// Convert lattice to time coordinate
  double LatticeToTime(double) const;

  /// Convert time to lattice coordinate
  double TimeToLattice(double) const;

  /// Put affine world coordinate transformation which is applied
  /// after the image to world coordinate transformation derived from the
  /// imaging geometry when mapping voxel indices to world coordinates.
  /// This transformation can be the inverse of the affine transformation
  /// obtained by an affine registration with this image as source.
  ///
  /// \param[in] m     Homogeneous transformation matrix.
  /// \param[in] apply Whether to apply the translation, rotation, and
  ///                  scaling directly to the image attributes and store
  ///                  only the shearing (if any) as additional transformation.
  void PutAffineMatrix(const Matrix &m, bool apply = false);

  /// Return transformation matrix for lattice to world coordinates
  Matrix GetLatticeToWorldMatrix() const;

  /// Return transformation matrix for world to lattice coordinates
  Matrix GetWorldToLatticeMatrix() const;

  /// Return orientation part of lattice to world coordinate transformation
  Matrix GetLatticeToWorldOrientation() const;

  /// Return orientation part of lattice to world coordinate transformation
  Matrix GetWorldToLatticeOrientation() const;

  /// Alias for GetLatticeToWorldMatrix
  inline Matrix GetImageToWorldMatrix()      const { return GetLatticeToWorldMatrix(); }
  /// Alias for GetWorldToLatticeMatrix
  inline Matrix GetWorldToImageMatrix()      const { return GetWorldToLatticeMatrix(); }
  /// Alias for GetLatticeToWorldOrientation
  inline Matrix GetImageToWorldOrientation() const { return GetLatticeToWorldOrientation(); }
  /// Alias for GetWorldToLatticeOrientation
  inline Matrix GetWorldToImageOrientation() const { return GetWorldToLatticeOrientation(); }

  // ---------------------------------------------------------------------------
  // Voxel index checks

  /// Whether voxel index is within finite image domain
  bool IsInside(int) const;

  /// Whether voxel indices are within finite 2D image domain
  bool IsInside(int, int) const;

  /// Whether voxel indices are within finite 3D image domain
  bool IsInside(int, int, int) const;

  /// Whether voxel indices are within finite 4D image domain
  bool IsInside(int, int, int, int) const;

  /// Whether voxel is index is outside finite image domain
  bool IsOutside(int) const;

  /// Whether voxel indices are outside finite 4D image domain
  bool IsOutside(int, int) const;

  /// Whether voxel indices are outside finite 4D image domain
  bool IsOutside(int, int, int) const;

  /// Whether voxel indices are outside finite 4D image domain
  bool IsOutside(int, int, int, int) const;

  /// Whether voxel index is at boundary of finite image domain
  bool IsBoundary(int) const;

  /// Whether voxel indices are at boundary of finite 2D image domain
  bool IsBoundary(int, int) const;

  /// Whether voxel indices are at boundary of finite 3D image domain
  bool IsBoundary(int, int, int) const;

  /// Whether voxel indices are at boundary of finite 4D image domain
  bool IsBoundary(int, int, int, int) const;

  // ---------------------------------------------------------------------------
  // Output

  /// Image slice area in world space
  double Area() const;

  /// Image volume in world space
  double Volume() const;

  /// Amount of space occupied by the image in n-D world space (excluding time dimension)
  double Space() const;

  /// Print attributes
  void Print(ostream &, Indent = 0) const;

  /// Print attributes
  void Print(Indent = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int ImageAttributes::X() const
{
  return _x;
}

// -----------------------------------------------------------------------------
inline int ImageAttributes::Y() const
{
  return _y;
}

// -----------------------------------------------------------------------------
inline int ImageAttributes::Z() const
{
  return _z;
}

// -----------------------------------------------------------------------------
inline int ImageAttributes::T() const
{
  return _t;
}

// -----------------------------------------------------------------------------
inline int ImageAttributes::N(int i) const
{
  if (i == 0) return _x;
  if (i == 1) return _y;
  if (i == 2) return _z;
  if (i == 3) return _t;
  return 0;
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::Spacing(int i) const
{
  if (i == 0) return _dx;
  if (i == 1) return _dy;
  if (i == 2) return _dz;
  if (i == 3) return _dt;
  return 0.;
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::XSpacing() const
{
  return _dx;
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::YSpacing() const
{
  return _dy;
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::ZSpacing() const
{
  return _dz;
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::TSpacing() const
{
  return _dt;
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::XSize() const
{
  return _dx;
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::YSize() const
{
  return _dy;
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::ZSize() const
{
  return _dz;
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::TSize() const
{
  return _dt;
}

// -----------------------------------------------------------------------------
inline ImageAttributes::operator bool() const
{
  // Note: _dz may be zero for 2D images, _dt may even be negative!
  return _x > 0 && _y > 0 && _z > 0 && _t > 0 && _dx > .0 && _dy > .0 && _dz >= .0;
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::operator==(const ImageAttributes &attr) const
{
  return EqualInSpace(attr) && EqualInTime(attr);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::operator!=(const ImageAttributes &attr) const
{
  return !(*this == attr);
}

// -----------------------------------------------------------------------------
inline int ImageAttributes::LatticeToIndex(int i, int j, int k, int l) const
{
  return i + _x * (j + _y * (k + _z * l));
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::IndexToLattice(int index, int *i, int *j, int *k, int *l) const
{
  int n = _x * _y * _z;
  if (l) *l = index / n;
  index = index % n;
  n = _x * _y;
  if (k) *k = index / n;
  *j = index % n / _x;
  *i = index % n % _x;
}

// -----------------------------------------------------------------------------
inline int ImageAttributes::NumberOfLatticePoints() const
{
  return _x * _y * _z * _t;
}

// -----------------------------------------------------------------------------
inline int ImageAttributes::NumberOfSpatialPoints() const
{
  return _x * _y * _z;
}

// -----------------------------------------------------------------------------
inline int ImageAttributes::NumberOfPoints() const
{
  return _dt ? NumberOfLatticePoints() : NumberOfSpatialPoints();
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::IndexToLattice(int index, int &i, int &j) const
{
  IndexToLattice(index, &i, &j);
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::IndexToLattice(int index, int &i, int &j, int &k) const
{
  IndexToLattice(index, &i, &j, &k);
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::IndexToLattice(int index, int &i, int &j, int &k, int &l) const
{
  IndexToLattice(index, &i, &j, &k, &l);
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::IndexToWorld(int idx, double &x, double &y) const
{
  int i, j, k;
  IndexToLattice(idx, i, j, k);
  x = i, y = j;
  double z = k;
  LatticeToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::IndexToWorld(int idx, double &x, double &y, double &z) const
{
  int i, j, k;
  IndexToLattice(idx, i, j, k);
  x = i, y = j, z = k;
  LatticeToWorld(x, y, z);
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::IndexToWorld(int idx, Point &p) const
{
  IndexToWorld(idx, p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline Point ImageAttributes::IndexToWorld(int idx) const
{
  Point p;
  IndexToWorld(idx, p);
  return p;
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::LatticeToWorld(double &x, double &y, double &z) const
{
  Transform(_i2w ? *_i2w : GetImageToWorldMatrix(), x, y, z);
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::LatticeToWorld(Point &p) const
{
  LatticeToWorld(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::LatticeToTime(double t) const
{
  return _torigin + t * _dt;
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::WorldToLattice(double &x, double &y, double &z) const
{
  Transform(_w2i ? *_w2i : GetWorldToImageMatrix(), x, y, z);
}

// -----------------------------------------------------------------------------
inline void ImageAttributes::WorldToLattice(Point &p) const
{
  WorldToLattice(p._x, p._y, p._z);
}

// -----------------------------------------------------------------------------
inline double ImageAttributes::TimeToLattice(double t) const
{
  return (_dt ? ((t - _torigin) / _dt) : .0);
}

// -----------------------------------------------------------------------------
inline Matrix ImageAttributes::GetLatticeToWorldOrientation() const
{
  Matrix R(3, 3);
  R(0, 0) = _xaxis[0];
  R(1, 0) = _xaxis[1];
  R(2, 0) = _xaxis[2];
  R(0, 1) = _yaxis[0];
  R(1, 1) = _yaxis[1];
  R(2, 1) = _yaxis[2];
  R(0, 2) = _zaxis[0];
  R(1, 2) = _zaxis[1];
  R(2, 2) = _zaxis[2];
  return R;
}

// -----------------------------------------------------------------------------
inline Matrix ImageAttributes::GetWorldToLatticeOrientation() const
{
  Matrix R(3, 3);
  R(0, 0) = _xaxis[0];
  R(0, 1) = _xaxis[1];
  R(0, 2) = _xaxis[2];
  R(1, 0) = _yaxis[0];
  R(1, 1) = _yaxis[1];
  R(1, 2) = _yaxis[2];
  R(2, 0) = _zaxis[0];
  R(2, 1) = _zaxis[1];
  R(2, 2) = _zaxis[2];
  return R;
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsInside(int idx) const
{
  return (0 <= idx && idx < _x * _y * _z * _t);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsInside(int i, int j) const
{
  return (0 <= i && i < _x) && (0 <= j && j < _y);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsInside(int i, int j, int k) const
{
  return IsInside(i, j) && (0 <= k && k < _z);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsInside(int i, int j, int k, int l) const
{
  return IsInside(i, j, k) && (0 <= l && l < _t);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsOutside(int idx) const
{
  return !IsInside(idx);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsOutside(int i, int j) const
{
  return !IsInside(i, j);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsOutside(int i, int j, int k) const
{
  return !IsInside(i, j, k);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsOutside(int i, int j, int k, int l) const
{
  return !IsInside(i, j, k, l);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsBoundary(int i, int j) const
{
  return (i == 0 || i == _x - 1) || (j == 0 || j == _y - 1);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsBoundary(int i, int j, int k) const
{
  return IsBoundary(i, j) || (k == 0 || k == _z - 1);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsBoundary(int i, int j, int k, int l) const
{
  return IsBoundary(i, j, k) || (l == 0 || l == _t - 1);
}

// -----------------------------------------------------------------------------
inline bool ImageAttributes::IsBoundary(int idx) const
{
  if (!IsInside(idx)) return false;
  int i, j, k, l;
  IndexToLattice(idx, i, j, k, l);
  if (_t == 1) {
    if (_z == 1) return IsBoundary(i, j);
    else         return IsBoundary(i, j, k);
  }
  return IsBoundary(i, j, k, l);
}

////////////////////////////////////////////////////////////////////////////////
// Image domain helpers
////////////////////////////////////////////////////////////////////////////////

/// Ensure that image domain has orthogonal basis vectors, i.e., no additional
/// affine transformation which contains any shearing
///
/// \returns Image grid which fully contains the input image domain but without
///          additional affine transformation (_smat is identity matrix).
///
/// \remarks The returned image domain need not be an axis-aligned bounding box!
///          When the Geometric Tools Engine (GTEngine) library is available,
///          the minimum-volume bounding box is returned.
ImageAttributes OrthogonalFieldOfView(const ImageAttributes &);

/// This method implements a method to find a common image grid which fully
/// contains all given image grids in world space
///
/// The voxel size of the resulting grid corresponds to the average voxel size
/// of the input images. The orientation and origin are chosen such that the
/// resulting image domain overlaps all input images in world space with a near
/// to minimal image volume. The current implementation only uses an heuristic
/// approach to find such minimal oriented bounding box of the input domains.
/// When the Geometric Tools Engine (GTEngine) library is available, the
/// minimum-volume bounding box is returned.
///
/// Additional information regarding algorithms to find the optimal
/// minimum-volume oriented bounding box (OBB) can be found at:
/// - http://stackoverflow.com/questions/7282805/where-can-i-find-a-c-c-implementation-of-the-minimum-bound-box-algorithm
/// - http://www.geometrictools.com/LibMathematics/Containment/Containment.html
/// - http://link.springer.com/article/10.1007%2FBF00991005?LI=true
///
/// \note The final overall image grid computed by this function must be
///       independent of the order in which the input attributes are given
///       in the input list. Otherwise an inverse consistent registration
///       might still depend on the order of the input images.
ImageAttributes OverallFieldOfView(const Array<ImageAttributes> &);


} // namespace mirtk

#endif // MIRTK_ImageAttributes_H
