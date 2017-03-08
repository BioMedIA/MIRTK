/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2017 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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

#include "mirtk/ImageConfig.h"

#include "mirtk/BaseImage.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Matrix3x3.h"

#include "mirtk/GenericImage.h"
#include "mirtk/ImageReader.h"

#if MIRTK_Image_WITH_VTK
  #include "vtkStructuredPoints.h"
#endif


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
BaseImage::BaseImage()
:
  _NumberOfVoxels(0),
  _mask(NULL),
  _maskOwner(false),
  _bg(-numeric_limits<double>::infinity()),
  _bgSet(false)
{
  _attr._i2w = &_matI2W;
  _attr._w2i = &_matW2I;
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
BaseImage::BaseImage(const ImageAttributes &attr, int n)
:
  _mask(NULL),
  _maskOwner(false),
  _bg(-numeric_limits<double>::infinity()),
  _bgSet(false)
{
  _attr = attr;
  if (n >= 1) {
    _attr._t  = n;
    _attr._dt = 0.;  // i.e., vector image with n components
  }
  PutAttributes(_attr);
}

// -----------------------------------------------------------------------------
BaseImage::BaseImage(const BaseImage &image)
:
  Object(image),
  _attr          (image._attr),
  _NumberOfVoxels(image._NumberOfVoxels), 
  _matI2W        (image._matI2W),
  _matW2I        (image._matW2I),
  _mask          (image._maskOwner ? new BinaryImage(*image._mask) : image._mask),
  _maskOwner     (image._maskOwner),
  _bg            (image._bg),
  _bgSet         (image._bgSet)
{
  _attr._i2w = &_matI2W;
  _attr._w2i = &_matW2I;
}

// -----------------------------------------------------------------------------
BaseImage::~BaseImage()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
BaseImage *BaseImage::New(const char *fname)
{
  UniquePtr<ImageReader> reader(ImageReader::New(fname));
  return reader->Run();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline BaseImage *NewImage(const BaseImage *base)
{
  const GenericImage<VoxelType> *image;
  if ((image = dynamic_cast<const GenericImage<VoxelType> *>(base))) {
    return new GenericImage<VoxelType>(*image);
  } else {
    cerr << "BaseImage::New: Input image is not of expected type" << endl;
    exit(1);
  }
}

// -----------------------------------------------------------------------------
BaseImage *BaseImage::New(const BaseImage *image)
{
  switch (image->GetDataType()) {
    case MIRTK_VOXEL_CHAR:           return NewImage<char>          (image);
    case MIRTK_VOXEL_UNSIGNED_CHAR:  return NewImage<unsigned char> (image);
    case MIRTK_VOXEL_SHORT:          return NewImage<short>         (image);
    case MIRTK_VOXEL_UNSIGNED_SHORT: return NewImage<unsigned short>(image);
    case MIRTK_VOXEL_INT:            return NewImage<int>           (image);
    case MIRTK_VOXEL_UNSIGNED_INT:   return NewImage<unsigned int>  (image);
    case MIRTK_VOXEL_FLOAT:          return NewImage<float>         (image);
    case MIRTK_VOXEL_DOUBLE:         return NewImage<double>        (image);
    default:
      cerr << "BaseImage::New: Cannot allocate image of unknown type: "
                << image->GetDataType() << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
BaseImage *BaseImage::New(int dtype)
{
  switch (dtype) {
    case MIRTK_VOXEL_CHAR:           return new GenericImage<char>;
    case MIRTK_VOXEL_UNSIGNED_CHAR:  return new GenericImage<unsigned char>;
    case MIRTK_VOXEL_SHORT:          return new GenericImage<short>;
    case MIRTK_VOXEL_UNSIGNED_SHORT: return new GenericImage<unsigned short>;
    case MIRTK_VOXEL_INT:            return new GenericImage<int>;
    case MIRTK_VOXEL_UNSIGNED_INT:   return new GenericImage<unsigned int>;
    case MIRTK_VOXEL_FLOAT:          return new GenericImage<float>;
    case MIRTK_VOXEL_DOUBLE:         return new GenericImage<double>;
    default:
      cerr << "BaseImage::New: Cannot allocate image of unknown type: " << dtype << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
BaseImage *BaseImage::Copy() const
{
  return BaseImage::New(this);
}

// -----------------------------------------------------------------------------
void BaseImage::UpdateMatrix()
{
  _matI2W = _attr.GetImageToWorldMatrix();
  _matW2I = _attr.GetWorldToImageMatrix();
}

// -----------------------------------------------------------------------------
void BaseImage::PutAttributes(const ImageAttributes &attr)
{
  _attr      = attr;
  _attr._i2w = &_matI2W;
  _attr._w2i = &_matW2I;
  _NumberOfVoxels = attr.NumberOfLatticePoints();
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
BaseImage& BaseImage::operator =(const BaseImage &image)
{
  // Copy image attributes
  this->Initialize(image.Attributes());
  // Copy/cast image data
  for (int idx = 0; idx < _NumberOfVoxels; idx++) {
    this->PutAsDouble(idx, image.GetAsDouble(idx));
  }
  // Copy foreground region
  if (image.OwnsMask()) _mask = new BinaryImage(*image.GetMask());
  else                  _mask = const_cast<BinaryImage *>(image.GetMask());
  _maskOwner = image.OwnsMask();
  if (image.HasBackgroundValue()) {
    this->PutBackgroundValueAsDouble(image.GetBackgroundValueAsDouble());
  }
  return *this;
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
// Based on nifti_mat44_to_orientation function from nifti1_io.c
void BaseImage::Orientation(OrientationCode &icod, OrientationCode &jcod, OrientationCode &kcod) const
{
  const Matrix &R = _matI2W; // FIXME: With or without affine transformation ?

  icod = jcod = kcod = Orientation_Unknown;

  // load column vectors for each (i,j,k) direction from matrix

  /*-- i axis --*/ /*-- j axis --*/ /*-- k axis --*/
  double xi = R(0, 0), xj = R(0, 1), xk = R(0, 2);
  double yi = R(1, 0), yj = R(1, 1), yk = R(1, 2);
  double zi = R(2, 0), zj = R(2, 1), zk = R(2, 2);

  // normalize column vectors to get unit vectors along each ijk-axis
  double val;

  // normalize i axis
  val = sqrt(xi*xi + yi*yi + zi*zi);
  if (val == .0) return; // bad orientation matrix
  xi /= val, yi /= val, zi /= val ;

  // normalize j axis
  val = sqrt(xj*xj + yj*yj + zj*zj);
  if (val == .0) return; // bad orientation matrix
  xj /= val, yj /= val, zj /= val ;

  // orthogonalize j axis to i axis, if needed
  val = xi*xj + yi*yj + zi*zj; // dot product between i and j
  if (abs(val) > 1.e-4) {
    xj -= val * xi;
    yj -= val * yi;
    zj -= val * zi;
    val = sqrt(xj*xj + yj*yj + zj*zj);  // must renormalize
    if (val == .0) return;              // j was parallel to i?
    xj /= val, yj /= val, zj /= val;
  }

  // normalize k axis; if it is zero, make it the cross product i x j
  val = sqrt(xk*xk + yk*yk + zk*zk);
  if (val == .0) {
    xk = yi * zj - zi * yj;
    yk = zi * xj - zj * xi;
    zk = xi * yj - yi * xj;
  } else {
    xk /= val, yk /= val, zk /= val;
  }

  // orthogonalize k to i
  val = xi*xk + yi*yk + zi*zk;    // dot product between i and k
  if (abs(val) > 1.e-4) {
    xk -= val * xi;
    yk -= val * yi;
    zk -= val * zi;
    val = sqrt(xk*xk + yk*yk + zk*zk);
    if (val == .0) return;      // bad
    xk /= val, yk /= val, zk /= val;
  }

  // orthogonalize k to j

  val = xj*xk + yj*yk + zj*zk;    // dot product between j and k
  if (abs(val) > 1.e-4) {
    xk -= val * xj;
    yk -= val * yj;
    zk -= val * zj;
    val = sqrt(xk*xk + yk*yk + zk*zk);
    if (val == .0) return;      // bad
    xk /= val, yk /= val, zk /= val;
  }

  Matrix3x3 Q(xi, xj, xk,  yi, yj, yk,  zi, zj, zk);

  // at this point, Q is the rotation matrix from the (i,j,k) to (x,y,z) axes
  double detQ = Q.Determinant();
  if (detQ == .0) return;  // shouldn't happen

  // Build and test all possible +1/-1 coordinate permutation matrices P;
  // then find the P such that the rotation matrix M=PQ is closest to the
  // identity, in the sense of M having the smallest total rotation angle.

  // Despite the formidable looking 6 nested loops, there are
  // only 3*3*3*2*2*2 = 216 passes, which will run very quickly.

  double vbest = -666.0;
  int    ibest = 1, pbest = 1, qbest = 1, rbest = 1, jbest = 2, kbest = 3;

  Matrix3x3 P, M;
  for (int i = 1; i <= 3; ++i)      // i = column number to use for row #1
  for (int j = 1; j <= 3; ++j) {    // j = column number to use for row #2
    if (i == j) continue;
    for (int k = 1; k <= 3; ++k) {  // k = column number to use for row #3
      if (i == k || j == k) continue;
      P = Matrix3x3::ZERO;
      for (int p = -1; p <= 1; p += 2)   // p,q,r are -1 or +1
      for (int q = -1; q <= 1; q += 2)   // and go into rows #1,2,3
      for (int r = -1; r <= 1; r += 2) {
        P[0][i-1] = p;
        P[1][j-1] = q;
        P[2][k-1] = r;
        // sign of permutation doesn't match sign of Q
        if (P.Determinant() * detQ <= .0) continue;
        M = P * Q;
        // angle of M rotation = 2.0*acos(0.5*sqrt(1.0+trace(M)))
        // we want largest trace(M) == smallest angle == M nearest to I
        val = M[0][0] + M[1][1] + M[2][2]; // trace
        if (val > vbest) {
          vbest = val;
          ibest = i, jbest = j, kbest = k;
          pbest = p, qbest = q, rbest = r;
        }
      }
    }
  }

  // At this point ibest is 1 or 2 or 3; pbest is -1 or +1; etc.

  // The matrix P that corresponds is the best permutation approximation
  // to Q-inverse; that is, P (approximately) takes (x,y,z) coordinates
  // to the (i,j,k) axes.
  //
  // For example, the first row of P (which contains pbest in column ibest)
  // determines the way the i axis points relative to the anatomical
  // (x,y,z) axes.  If ibest is 2, then the i axis is along the y axis,
  // which is direction P2A (if pbest > 0) or A2P (if pbest < 0).
  //
  // So, using ibest and pbest, we can assign the output code for
  // the i axis.  Mutatis mutandis for the j and k axes, of course.

  switch (ibest * pbest) {
    case  1: icod = L2R; break;
    case -1: icod = R2L; break;
    case  2: icod = P2A; break;
    case -2: icod = A2P; break;
    case  3: icod = I2S; break;
    case -3: icod = S2I; break;
  }

  switch (jbest * qbest) {
    case  1: jcod = L2R; break;
    case -1: jcod = R2L; break;
    case  2: jcod = P2A; break;
    case -2: jcod = A2P; break;
    case  3: jcod = I2S; break;
    case -3: jcod = S2I; break;
  }

  switch (kbest * rbest) {
    case  1: kcod = L2R; break;
    case -1: kcod = R2L; break;
    case  2: kcod = P2A; break;
    case -2: kcod = A2P; break;
    case  3: kcod = I2S; break;
    case -3: kcod = S2I; break;
  }
}

// -----------------------------------------------------------------------------
void BaseImage::ImageToWorld(WorldCoordsImage &i2w, bool _3D) const
{
  if (_attr._z == 1 && !_3D) {
    const int nvoxs = _attr._x * _attr._y;
    i2w.Initialize(_attr, 2);
    double *wx = i2w.GetPointerToVoxels();
    double *wy = wx + nvoxs;
    for (int j = 0; j < _attr._y; j++) {
      for (int i = 0; i < _attr._x; i++) {
        (*wx++) = _matI2W(0, 0) * i + _matI2W(0, 1) * j + _matI2W(0, 3);
        (*wy++) = _matI2W(1, 0) * i + _matI2W(1, 1) * j + _matI2W(1, 3);
      }
    }
  } else {
    const int nvoxs = _attr._x * _attr._y * _attr._z;
    i2w.Initialize(_attr, 3);
    double *wx = i2w.GetPointerToVoxels();
    double *wy = wx + nvoxs;
    double *wz = wy + nvoxs;
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          (*wx++) = _matI2W(0, 0) * i + _matI2W(0, 1) * j + _matI2W(0, 2) * k + _matI2W(0, 3);
          (*wy++) = _matI2W(1, 0) * i + _matI2W(1, 1) * j + _matI2W(1, 2) * k + _matI2W(1, 3);
          (*wz++) = _matI2W(2, 0) * i + _matI2W(2, 1) * j + _matI2W(2, 2) * k + _matI2W(2, 3);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
void BaseImage::ImageToWorld(double *i2w, bool _3D) const
{
  if (_attr._z == 1 && !_3D) {
    for (int j = 0; j < _attr._y; j++) {
      for (int i = 0; i < _attr._x; i++) {
        (*i2w++) = _matI2W(0, 0) * i + _matI2W(0, 1) * j + _matI2W(0, 3);
        (*i2w++) = _matI2W(1, 0) * i + _matI2W(1, 1) * j + _matI2W(1, 3);
      }
    }
  } else {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          (*i2w++) = _matI2W(0, 0) * i + _matI2W(0, 1) * j + _matI2W(0, 2) * k + _matI2W(0, 3);
          (*i2w++) = _matI2W(1, 0) * i + _matI2W(1, 1) * j + _matI2W(1, 2) * k + _matI2W(1, 3);
          (*i2w++) = _matI2W(2, 0) * i + _matI2W(2, 1) * j + _matI2W(2, 2) * k + _matI2W(2, 3);
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
void BaseImage::ImageToWorld(PointSet &points) const
{
  Point p;
  points.Reserve(points.Size() + NumberOfSpatialVoxels());
  for (int k = 0; k < Z(); ++k)
  for (int j = 0; j < Y(); ++j)
  for (int i = 0; i < X(); ++i) {
    p._x = i, p._y = j, p._z = k;
    ImageToWorld(p);
    points.Add(p);
  }
}

// =============================================================================
// Common image statistics
// =============================================================================

// -----------------------------------------------------------------------------
void BaseImage::GetMinMaxAsDouble(double &min, double &max) const
{
  min = max = .0;

  bool   first = true;
  double value;

  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      value = this->GetAsDouble(idx);
      if (first) {
        min   = max = value;
        first = false;
      } else {
        if (value < min) min = value;
        if (value > max) max = value;
      }
    }
  }
}

// -----------------------------------------------------------------------------
void BaseImage::PutMinMaxAsDouble(double min, double max)
{
  // Check if full range can be represented by data type
  if (min < this->GetDataTypeMin()) {
    cerr << this->NameOfClass() << "::PutMinMaxAsDouble: "
              << "Requested minimum value is out of range for the voxel type: "
              << this->GetDataType() << endl;
    exit(1);
  }
  if (max > this->GetDataTypeMax()) {
    cerr << this->NameOfClass() << "::PutMinMaxAsDouble: "
              << "Requested maximum value is out of range for the voxel type: "
              << this->GetDataType() << endl;
    exit(1);
  }

  // Get current min and max values
  double min_val, max_val;
  this->GetMinMaxAsDouble(min_val, max_val);

  // Rescale foreground to desired [min, max] range
  const double slope = (max - min) / (max_val - min_val);
  const double inter = min - slope * min_val;
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      this->PutAsDouble(idx, inter + slope * this->GetAsDouble(idx));
    }
  }
}

// =============================================================================
// Foreground region
// =============================================================================

// -----------------------------------------------------------------------------
void BaseImage::PutMask(BinaryImage *mask, bool owner)
{
  if (_maskOwner) delete _mask;
  _mask      = mask;
  _maskOwner = owner;
}

// -----------------------------------------------------------------------------
void BaseImage::InitializeMask(int t, bool force)
{
  if (!_mask) {
    _mask      = new BinaryImage(_attr, t);
    _maskOwner = 2;
    force      = true;
  }
  if (force) {
    if (!_bgSet) {
      cerr << "BaseImage::InitializeMask: Background value not set!" << endl;
      exit(1);
    }
    if (_mask->GetT() == _attr._t) {
      BinaryPixel *ptr2msk = _mask->GetPointerToVoxels();
      for (int l = 0; l < _attr._t; l++) {
        for (int k = 0; k < _attr._z; k++) {
          for (int j = 0; j < _attr._y; j++) {
            for (int i = 0; i < _attr._x; i++) {
              (*ptr2msk++) = (this->GetAsDouble(i, j, k, l) == _bg ? false : true);
            }
          }
        }
      }
    } else {
      BinaryPixel *ptr2msk = _mask->GetPointerToVoxels();
      for (int k = 0; k < _attr._z; k++) {
        for (int j = 0; j < _attr._y; j++) {
          for (int i = 0; i < _attr._x; i++) {
            (*ptr2msk++) = (this->GetAsDouble(i, j, k, 0) == _bg ? false : true);
          }
        }
      }
      for (int l = 1; l < _attr._t; l++) {
        BinaryPixel *ptr2msk = _mask->GetPointerToVoxels();
        for (int k = 0; k < _attr._z; k++) {
          for (int j = 0; j < _attr._y; j++) {
            for (int i = 0; i < _attr._x; i++) {
              if (*ptr2msk != BinaryPixel(0)) {
                (*ptr2msk) = (this->GetAsDouble(i, j, k, l) == _bg ? false : true);
              }
              ++ptr2msk;
            }
          }
        }
      }
    }
  }
  if (_maskOwner > 1) {
    if (_maskOwner == numeric_limits<int>::max()) {
      cerr << "BaseImage::InitializeMask: Too many nested calls! Do also not forget to call ClearMask at some point." << endl;
      exit(1);
    }
    _maskOwner++;
  }
}

// -----------------------------------------------------------------------------
void BaseImage::ClearMask(bool force)
{
  if (_maskOwner > 2) {
    _maskOwner = (force ? 2 : _maskOwner-1);
    if (_maskOwner == 2) delete _mask, _mask = NULL;
  }
}

// -----------------------------------------------------------------------------
void BaseImage::ResetBackgroundValueAsDouble(double bg)
{
  if (this->HasBackgroundValue()) {
    const int nvox = this->NumberOfVoxels();
    for (int vox = 0; vox < nvox; ++vox) {
      if (AreEqualOrNaN(this->GetAsDouble(vox), _bg, 1e-6)) {
        this->PutAsDouble(vox, bg);
      }
    }
  }
  this->PutBackgroundValueAsDouble(bg);
}

// -----------------------------------------------------------------------------
void BaseImage::BoundingBox(int &i1, int &j1,
                            int &i2, int &j2) const
{
  // Determine lower bound along x axis: i1
  i1 = _attr._x;
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    if (this->IsForeground(i, j, k, l)) {
      if (i < i1) i1 = i;
      break;
    }
  }
  // Determine upper bound along x axis: i2
  i2 = -1;
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = _attr._x - 1; i >= i1; --i) {
    if (this->IsForeground(i, j, k, l)) {
      if (i > i2) i2 = i;
      break;
    }
  }
  // Determine lower bound along y axis: j1
  j1 = _attr._y;
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int i = i1; i <= i2; ++i)
  for (int j = 0; j < _attr._y; ++j) {
    if (this->IsForeground(i, j, k, l)) {
      if (j < j1) j1 = j;
      break;
    }
  }
  // Determine upper bound along y axis: j2
  j2 = -1;
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int i = i1; i <= i2; ++i)
  for (int j = _attr._y - 1; j >= j1; --j) {
    if (this->IsForeground(i, j, k, l)) {
      if (j > j2) j2 = j;
      break;
    }
  }
}

// -----------------------------------------------------------------------------
void BaseImage::BoundingBox(int &i1, int &j1, int &k1,
                            int &i2, int &j2, int &k2) const
{
  // Get 2D bounding box
  this->BoundingBox(i1, j1, i2, j2);
  // Determine lower bound along z axis: k1
  k1 = _attr._z;
  for (int l = 0; l < _attr._t; ++l)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i)
  for (int k = 0; k < _attr._z; ++k) {
    if (this->IsForeground(i, j, k, l)) {
      if (k < k1) k1 = k;
      break;
    }
  }
  // Determine upper bound along z axis: k2
  k2 = -1;
  for (int l = 0; l < _attr._t; ++l)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i)
  for (int k = _attr._z - 1; k >= k1; --k) {
    if (this->IsForeground(i, j, k, l)) {
      if (k > k2) k2 = k;
      break;
    }
  }
}

// -----------------------------------------------------------------------------
void BaseImage::BoundingBox(int &i1, int &j1, int &k1, int &l1,
                            int &i2, int &j2, int &k2, int &l2) const
{
  // Get spatial bounding box
  this->BoundingBox(i1, j1, k1, i2, j2, k2);
  // Determine lower bound along t axis: l1
  l1 = _attr._t;
  for (int k = k1; k <= k2; ++k)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i)
  for (int l = 0; l < _attr._t; ++l) {
    if (this->IsForeground(i, j, k, l)) {
      if (l < l1) l1 = l;
      break;
    }
  }
  // Determine upper bound along t axis: l2
  l2 = -1;
  for (int k = k1; k <= k2; ++k)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i)
  for (int l = _attr._t - 1; l >= l1; --l) {
    if (this->IsForeground(i, j, k, l)) {
      if (l > l2) l2 = l;
      break;
    }
  }
}

// -----------------------------------------------------------------------------
int BaseImage::CenterOfForeground(Point &center) const
{
  int n = 0;
  center._x = .0, center._y = .0, center._z = .0;
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    if (IsForeground(i, j, k)) {
      center._x += i;
      center._y += j;
      center._z += k;
      ++n;
    }
  }
  if (n > 0) center /= n;
  ImageToWorld(center);
  return n;
}

// -----------------------------------------------------------------------------
int BaseImage::CenterOfForeground(Point &center, double padding) const
{
  int n = 0;
  center._x = .0, center._y = .0, center._z = .0;
  if (IsNaN(padding)) {
    for (int k = 0; k < _attr._z; ++k)
    for (int j = 0; j < _attr._y; ++j)
    for (int i = 0; i < _attr._x; ++i) {
      if (!IsNaN(this->GetAsDouble(i, j, k))) {
        center._x += i;
        center._y += j;
        center._z += k;
        ++n;
      }
    }
  } else {
    for (int k = 0; k < _attr._z; ++k)
    for (int j = 0; j < _attr._y; ++j)
    for (int i = 0; i < _attr._x; ++i) {
      if (this->GetAsDouble(i, j, k) > padding) {
        center._x += i;
        center._y += j;
        center._z += k;
        ++n;
      }
    }
  }
  if (n > 0) center /= n;
  ImageToWorld(center);
  return n;
}

// -----------------------------------------------------------------------------
// Auxiliary function for BaseImage::ForegroundDomain overloads
ImageAttributes BaseImage::ForegroundDomain(int i1, int j1, int k1,
                                            int i2, int j2, int k2,
                                            bool orthogonal) const
{
  // Copy image attributes
  ImageAttributes attr = Attributes();
  // Adjust image attributes
  if (attr._z == 1) attr._dz = .0;
  if (attr._t == 1) attr._dt = .0;
  if (i1 <= i2 && j1 <= j2 && k1 <= k2) {
    attr._x = i2 - i1 + 1;
    attr._y = j2 - j1 + 1;
    attr._z = k2 - k1 + 1;
    attr._xorigin = i1 + .5 * (i2 - i1);
    attr._yorigin = j1 + .5 * (j2 - j1);
    attr._zorigin = k1 + .5 * (k2 - k1);
    ImageToWorld(attr._xorigin, attr._yorigin, attr._zorigin);
  }
  // Orthogonalize coordinate system; required in case of input images where
  // a previous affine (12 DoFs) alignment has been applied to the attributes
  return (orthogonal ? OrthogonalFieldOfView(attr) : attr);
}

// -----------------------------------------------------------------------------
ImageAttributes BaseImage::ForegroundDomain(bool orthogonal) const
{
  ImageAttributes attr = Attributes();
  // Determine lower bound along x axis: i1
  int i1 = attr._x;
  for (int k = 0; k < attr._z; ++k)
  for (int j = 0; j < attr._y; ++j)
  for (int i = 0; i < attr._x; ++i) {
    if (IsForeground(i, j, k)) {
      if (i < i1) i1 = i;
      break;
    }
  }
  // Determine upper bound along x axis: i2
  int i2 = -1;
  for (int k = 0; k < attr._z; ++k)
  for (int j = 0; j < attr._y; ++j)
  for (int i = attr._x - 1; i >= i1; --i) {
    if (IsForeground(i, j, k)) {
      if (i > i2) i2 = i;
      break;
    }
  }
  // Determine lower bound along y axis: j1
  int j1 = attr._y;
  for (int k = 0; k < attr._z; ++k)
  for (int i = i1; i <= i2; ++i)
  for (int j = 0; j < attr._y; ++j) {
    if (IsForeground(i, j, k)) {
      if (j < j1) j1 = j;
      break;
    }
  }
  // Determine upper bound along y axis: j2
  int j2 = -1;
  for (int k = 0; k < attr._z; ++k)
  for (int i = i1; i <= i2; ++i)
  for (int j = attr._y - 1; j >= j1; --j) {
    if (IsForeground(i, j, k)) {
      if (j > j2) j2 = j;
      break;
    }
  }
  // Determine lower bound along z axis: k1
  int k1 = attr._z;
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i)
  for (int k = 0; k < attr._z; ++k) {
    if (IsForeground(i, j, k)) {
      if (k < k1) k1 = k;
      break;
    }
  }
  // Determine upper bound along z axis: k2
  int k2 = -1;
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i)
  for (int k = attr._z - 1; k >= k1; --k) {
    if (IsForeground(i, j, k)) {
      if (k > k2) k2 = k;
      break;
    }
  }
  return ForegroundDomain(i1, j1, k1, i2, j2, k2, orthogonal);
}

// -----------------------------------------------------------------------------
ImageAttributes BaseImage::ForegroundDomain(double padding, bool orthogonal) const
{
  ImageAttributes attr = Attributes();

  int i1 = attr._x;
  int j1 = attr._y;
  int k1 = attr._z;

  int i2 = -1;
  int j2 = -1;
  int k2 = -1;

  if (IsNaN(padding)) {

    // Determine lower bound along x axis: i1
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      if (!IsNaN(this->GetAsDouble(i, j, k))) {
        if (i < i1) i1 = i;
        break;
      }
    }
    // Determine upper bound along x axis: i2
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = attr._x - 1; i >= i1; --i) {
      if (!IsNaN(this->GetAsDouble(i, j, k))) {
        if (i > i2) i2 = i;
        break;
      }
    }
    // Determine lower bound along y axis: j1
    for (int k = 0; k < attr._z; ++k)
    for (int i = i1; i <= i2; ++i)
    for (int j = 0; j < attr._y; ++j) {
      if (!IsNaN(this->GetAsDouble(i, j, k))) {
        if (j < j1) j1 = j;
        break;
      }
    }
    // Determine upper bound along y axis: j2
    for (int k = 0; k < attr._z; ++k)
    for (int i = i1; i <= i2; ++i)
    for (int j = attr._y - 1; j >= j1; --j) {
      if (!IsNaN(this->GetAsDouble(i, j, k))) {
        if (j > j2) j2 = j;
        break;
      }
    }
    // Determine lower bound along z axis: k1
    for (int j = j1; j <= j2; ++j)
    for (int i = i1; i <= i2; ++i)
    for (int k = 0; k < attr._z; ++k) {
      if (!IsNaN(this->GetAsDouble(i, j, k))) {
        if (k < k1) k1 = k;
        break;
      }
    }
    // Determine upper bound along z axis: k2
    for (int j = j1; j <= j2; ++j)
    for (int i = i1; i <= i2; ++i)
    for (int k = attr._z - 1; k >= k1; --k) {
      if (!IsNaN(this->GetAsDouble(i, j, k))) {
        if (k > k2) k2 = k;
        break;
      }
    }

  } else {

    // Determine lower bound along x axis: i1
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      if (this->GetAsDouble(i, j, k) > padding) {
        if (i < i1) i1 = i;
        break;
      }
    }
    // Determine upper bound along x axis: i2
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = attr._x - 1; i >= i1; --i) {
      if (this->GetAsDouble(i, j, k) > padding) {
        if (i > i2) i2 = i;
        break;
      }
    }
    // Determine lower bound along y axis: j1
    for (int k = 0; k < attr._z; ++k)
    for (int i = i1; i <= i2; ++i)
    for (int j = 0; j < attr._y; ++j) {
      if (this->GetAsDouble(i, j, k) > padding) {
        if (j < j1) j1 = j;
        break;
      }
    }
    // Determine upper bound along y axis: j2
    for (int k = 0; k < attr._z; ++k)
    for (int i = i1; i <= i2; ++i)
    for (int j = attr._y - 1; j >= j1; --j) {
      if (this->GetAsDouble(i, j, k) > padding) {
        if (j > j2) j2 = j;
        break;
      }
    }
    // Determine lower bound along z axis: k1
    for (int j = j1; j <= j2; ++j)
    for (int i = i1; i <= i2; ++i)
    for (int k = 0; k < attr._z; ++k) {
      if (this->GetAsDouble(i, j, k) > padding) {
        if (k < k1) k1 = k;
        break;
      }
    }
    // Determine upper bound along z axis: k2
    for (int j = j1; j <= j2; ++j)
    for (int i = i1; i <= i2; ++i)
    for (int k = attr._z - 1; k >= k1; --k) {
      if (this->GetAsDouble(i, j, k) > padding) {
        if (k > k2) k2 = k;
        break;
      }
    }

  }
  return ForegroundDomain(i1, j1, k1, i2, j2, k2, orthogonal);
}

// =============================================================================
// VTK interface
// =============================================================================
#if MIRTK_Image_WITH_VTK

// -----------------------------------------------------------------------------
int BaseImage::ImageToVTKScalarType() const
{
  return ToVTKDataType(this->GetDataType());
}

// -----------------------------------------------------------------------------
void BaseImage::ImageToVTK(vtkStructuredPoints *vtk) const
{
  if (this->ImageToVTKScalarType() == VTK_VOID) {
    cerr << this->NameOfType() << "::ImageToVTK: Cannot convert image to VTK structured points" << endl;
    exit(1);
  }
  double x = .0, y = .0, z = .0;
  this->ImageToWorld(x, y, z);
  vtk->SetOrigin    (x, y, z);
  vtk->SetDimensions(_attr._x,  _attr._y,  _attr._z);
  vtk->SetSpacing   (_attr._dx, _attr._dy, _attr._dz);
  vtk->AllocateScalars(this->ImageToVTKScalarType(), _attr._t);
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    vtk->SetScalarComponentFromDouble(i, j, k, l, this->GetAsDouble(i, j, k, l));
  }
}

// -----------------------------------------------------------------------------
void BaseImage::VTKToImage(vtkStructuredPoints *)
{
  cerr << this->NameOfClass() << "::VTKToImage: Not implemented" << endl;
  exit(1);
}

#endif // MIRTK_Image_WITH_VTK
// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void BaseImage::Print(Indent indent) const
{
  cout << indent << "Image lattice:\n";
  _attr.Print(indent + 1);
  cout << indent << "Foreground mask:  " << ToString(_mask != NULL) << "\n";
  cout << indent << "Background value: ";
  if (_bgSet) cout << _bg;
  else        cout << "n/a";
  cout << "\n";
}


} // namespace mirtk
