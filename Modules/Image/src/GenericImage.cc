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

#include "mirtk/ImageConfig.h"
#include "mirtk/GenericImage.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Path.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/VoxelCast.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Point.h"

#include "mirtk/ImageReader.h"
#include "mirtk/ImageWriter.h"

#if MIRTK_Image_WITH_VTK
#  include "vtkStructuredPoints.h"
#endif

// Default output image file name extension used by GenericImage::Write
// if none was provided (e.g., when called by debugging library functions
// such as overridden EnergyTerm::WriteDataSets implementations).
#ifndef MIRTK_Image_DEFAULT_EXT
#  define MIRTK_Image_DEFAULT_EXT ".gipl"
#endif


namespace mirtk {

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
// Note: Base class BaseImage must be initialized before calling this function!
template <class VoxelType>
void GenericImage<VoxelType>::AllocateImage(VoxelType *data)
{
  // Delete existing mask (if any)
  if (_maskOwner) Delete(_mask);
  // Free previously allocated memory
  Deallocate(_matrix, _data);
  if (_dataOwner) Deallocate(_data);
  _dataOwner = false;
  // Initialize memory
  const int nvox = _attr.NumberOfLatticePoints();
  if (nvox > 0) {
    if (data) {
      _data      = data;
      _dataOwner = false;
    } else {
      _data      = CAllocate<VoxelType>(nvox);
      _dataOwner = true;
    }
    Allocate(_matrix, _attr._x, _attr._y, _attr._z, _attr._t, _data);
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::GenericImage()
:
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::GenericImage(const char *fname)
:
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  Read(fname);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::GenericImage(int x, int y, int z, int t, VoxelType *data)
:
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  ImageAttributes attr;
  attr._x = x;
  attr._y = y;
  attr._z = z;
  attr._t = t;
  PutAttributes(attr);
  AllocateImage(data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::GenericImage(int x, int y, int z, int t, int n, VoxelType *data)
:
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  if (t > 1 && n > 1) {
    cerr << "GenericImage::GenericImage: 5D images not supported! Use 4D image with vector voxel type instead." << endl;
    exit(1);
  }
  ImageAttributes attr;
  if (n > 1) t = n, attr._dt = .0; // i.e., vector image with n components
  attr._x = x;
  attr._y = y;
  attr._z = z;
  attr._t = t;
  PutAttributes(attr);
  AllocateImage(data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::GenericImage(const ImageAttributes &attr, VoxelType *data)
:
  BaseImage(attr),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  AllocateImage(data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::GenericImage(const ImageAttributes &attr, int n, VoxelType *data)
:
  BaseImage(attr, n),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  AllocateImage(data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::GenericImage(const BaseImage &image)
:
  BaseImage(image),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  // Initialize image
  AllocateImage();
  // Copy/cast data
  VoxelType *ptr = _data;
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
    (*ptr) = voxel_cast<VoxelType>(image.GetAsVector(idx));
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::GenericImage(const GenericImage &image)
:
  BaseImage(image),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  if (image._dataOwner) {
    AllocateImage();
    memcpy(_data, image._data, _NumberOfVoxels * sizeof(VoxelType));
  } else {
    AllocateImage(const_cast<VoxelType *>(image.Data()));
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
GenericImage<VoxelType>::GenericImage(const GenericImage<VoxelType2> &image)
:
  BaseImage(image),
  _matrix   (NULL),
  _data     (NULL),
  _dataOwner(false)
{
  AllocateImage();
  VoxelType        *ptr1 = this->Data();
  const VoxelType2 *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    ptr1[idx] = voxel_cast<VoxelType>(ptr2[idx]);
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::~GenericImage()
{
  Deallocate(_matrix, _data);
  if (_dataOwner) Deallocate(_data);
  if (_maskOwner) Delete(_mask);
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
BaseImage *GenericImage<VoxelType>::Copy() const
{
  return new GenericImage<VoxelType>(*this);
}

// -----------------------------------------------------------------------------
template <class VoxelType> void GenericImage<VoxelType>::Initialize()
{
  if (_matrix) *this = VoxelType();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::Initialize(const ImageAttributes &a, int n, VoxelType *data)
{
  // Initialize attributes
  ImageAttributes attr(a);
  if (n >= 1) {
    attr._t  = n;
    attr._dt = 0;  // i.e., vector image with n components
  }
  // Initialize memory
  if (_attr._x != attr._x || _attr._y != attr._y || _attr._z != attr._z || _attr._t != attr._t) {
    PutAttributes(attr);
    AllocateImage(data);
  } else {
    PutAttributes(attr);
    if (_dataOwner) *this = VoxelType();
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::Initialize(const ImageAttributes &attr, int n)
{
  this->Initialize(attr, n, nullptr);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::Initialize(const ImageAttributes &attr, VoxelType *data)
{
  this->Initialize(attr, -1, data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::Initialize(int x, int y, int z, int t, int n, VoxelType *data)
{
  ImageAttributes attr(_attr);
  attr._x = x;
  attr._y = y;
  attr._z = z;
  attr._t = t;
  this->Initialize(attr, n, data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::Initialize(int x, int y, int z, int t, VoxelType *data)
{
  this->Initialize(x, y, z, t, 1, data);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::CopyFrom(const VoxelType *data)
{
  if (_data != data) {
    memcpy(_data, data, _NumberOfVoxels * sizeof(VoxelType));
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::CopyFrom(const BaseImage &image)
{
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    _matrix[l][k][j][i] = voxel_cast<VoxelType>(image.GetAsVector(i, j, k, l));
  }
  if (_maskOwner) delete _mask;
  if (image.OwnsMask()) {
    _mask      = new BinaryImage(*image.GetMask());
    _maskOwner = true;
  } else {
    _mask      = const_cast<BinaryImage *>(image.GetMask());
    _maskOwner = false;
  }
  if (image.HasBackgroundValue()) {
    this->PutBackgroundValueAsDouble(image.GetBackgroundValueAsDouble());
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::CopyFrom(const GenericImage &image)
{
  CopyFrom(image.Data());
  if (_maskOwner) delete _mask;
  if (image.OwnsMask()) {
    _mask      = new BinaryImage(*image.GetMask());
    _maskOwner = true;
  } else {
    _mask      = const_cast<BinaryImage *>(image.GetMask());
    _maskOwner = false;
  }
  if (image.HasBackgroundValue()) {
    this->PutBackgroundValueAsDouble(image.GetBackgroundValueAsDouble());
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator=(VoxelType scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    ptr[idx] = scalar;
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator=(const BaseImage &image)
{
  if (this != &image) {
    this->Initialize(image.Attributes());
    this->CopyFrom(image);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator=(const GenericImage &image)
{
  if (this != &image) {
    this->Initialize(image.Attributes());
    this->CopyFrom(image);
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType> void GenericImage<VoxelType>::Clear()
{
  Deallocate(_matrix, _data);
  if (_dataOwner) Deallocate(_data);
  if (_maskOwner) Delete(_mask);
  _attr = ImageAttributes();
}

// =============================================================================
// Region-of-interest extraction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>
::GetRegion(GenericImage<VoxelType> &image, int k, int m) const
{
  int i, j;
  double x1, y1, z1, t1, x2, y2, z2, t2;

  if ((k < 0) || (k >= _attr._z) || (m < 0) || (m >= _attr._t)) {
    cerr << "GenericImage<VoxelType>::GetRegion: Parameter out of range" << endl;
    exit(1);
  }

  // Initialize
  ImageAttributes attr = this->Attributes();
  attr._z = 1;
  attr._t = 1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  image.Initialize(attr);

  // Calculate position of first voxel in roi in original image
  x1 = 0;
  y1 = 0;
  z1 = k;
  this->ImageToWorld(x1, y1, z1);
  t1 = this->ImageToTime(m);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  t2 = 0;
  image.ImageToWorld(x2, y2, z2);
  t2 = image.ImageToTime(0);

  // Shift origin of new image accordingly
  image.PutOrigin(x1 - x2, y1 - y2, z1 - z2, t1 - t2);

  // Copy region
  for (j = 0; j < _attr._y; j++) {
    for (i = 0; i < _attr._x; i++) {
      image._matrix[0][0][j][i] = _matrix[m][k][j][i];
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>
::GetRegion(int k, int m) const
{
  GenericImage<VoxelType> image;
  this->GetRegion(image, k, m);
  return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>
::GetRegion(BaseImage *&base, int k, int m) const
{
  GenericImage<VoxelType> *image = dynamic_cast<GenericImage<VoxelType> *>(base);
  if (image == NULL) {
    delete base;
    image = new GenericImage<VoxelType>();
    base  = image;
  }
  this->GetRegion(*image, k, m);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>
::GetRegion(GenericImage<VoxelType> &image, int i1, int j1, int k1,
                                            int i2, int j2, int k2) const
{
  int i, j, k, l;
  double x1, y1, z1, x2, y2, z2;

  if ((i1 < 0) || (i1 >= i2) ||
      (j1 < 0) || (j1 >= j2) ||
      (k1 < 0) || (k1 >= k2) ||
      (i2 > _attr._x) || (j2 > _attr._y) || (k2 > _attr._z)) {
    cerr << "GenericImage<VoxelType>::GetRegion: Parameter out of range\n";
    exit(1);
  }

  // Initialize
  ImageAttributes attr = this->Attributes();
  attr._x = i2 - i1;
  attr._y = j2 - j1;
  attr._z = k2 - k1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  image.Initialize(attr);

  // Calculate position of first voxel in roi in original image
  x1 = i1;
  y1 = j1;
  z1 = k1;
  this->ImageToWorld(x1, y1, z1);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  image.ImageToWorld(x2, y2, z2);

  // Shift origin of new image accordingly
  image.PutOrigin(x1 - x2, y1 - y2, z1 - z2);

  // Copy region
  for (l = 0; l < _attr._t; l++) {
    for (k = k1; k < k2; k++) {
      for (j = j1; j < j2; j++) {
        for (i = i1; i < i2; i++) {
          image._matrix[l][k-k1][j-j1][i-i1] = _matrix[l][k][j][i];
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>
::GetRegion(int i1, int j1, int k1, int i2, int j2, int k2) const
{
  GenericImage<VoxelType> image;
  this->GetRegion(image, i1, j1, k1, i2, j2, k2);
  return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>
::GetRegion(BaseImage *&base, int i1, int j1, int k1, int i2, int j2, int k2) const
{
  GenericImage<VoxelType> *image = dynamic_cast<GenericImage<VoxelType> *>(base);
  if (image == NULL) {
    delete base;
    image = new GenericImage<VoxelType>();
    base  = image;
  }
  this->GetRegion(*image, i1, j1, k1, i2, j2, k2);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>
::GetRegion(GenericImage<VoxelType> &image, int i1, int j1, int k1, int l1,
                                            int i2, int j2, int k2, int l2) const
{
  int i, j, k, l;
  double x1, y1, z1, x2, y2, z2;

  if ((i1 < 0) || (i1 >= i2) ||
      (j1 < 0) || (j1 >= j2) ||
      (k1 < 0) || (k1 >= k2) ||
      (l1 < 0) || (l1 >= l2) ||
      (i2 > _attr._x) || (j2 > _attr._y) || (k2 > _attr._z) || (l2 > _attr._t)) {
    cerr << "GenericImage<VoxelType>::GetRegion: Parameter out of range\n";
    exit(1);
  }

  // Initialize
  ImageAttributes attr = this->Attributes();
  attr._x = i2 - i1;
  attr._y = j2 - j1;
  attr._z = k2 - k1;
  attr._t = l2 - l1;
  attr._xorigin = 0;
  attr._yorigin = 0;
  attr._zorigin = 0;
  image.Initialize(attr);

  // Calculate position of first voxel in roi in original image
  x1 = i1;
  y1 = j1;
  z1 = k1;
  this->ImageToWorld(x1, y1, z1);

  // Calculate position of first voxel in roi in new image
  x2 = 0;
  y2 = 0;
  z2 = 0;
  image.ImageToWorld(x2, y2, z2);

  // Shift origin of new image accordingly
  image.PutOrigin(x1 - x2, y1 - y2, z1 - z2);

  // Copy region
  for (l = l1; l < l2; l++) {
    for (k = k1; k < k2; k++) {
      for (j = j1; j < j2; j++) {
        for (i = i1; i < i2; i++) {
          image._matrix[l-l1][k-k1][j-j1][i-i1] = _matrix[l][k][j][i];
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>
::GetRegion(int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2) const
{
  GenericImage<VoxelType> image;
  this->GetRegion(image, i1, j1, k1, l1, i2, j2, k2, l2);
  return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>
::GetRegion(BaseImage *&base, int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2) const
{
  GenericImage<VoxelType> *image = dynamic_cast<GenericImage<VoxelType> *>(base);
  if (image == NULL) {
    delete base;
    image = new GenericImage<VoxelType>();
    base  = image;
  }
  this->GetRegion(*image, i1, j1, k1, l1, i2, j2, k2, l2);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::GetFrame(GenericImage<VoxelType> &image, int l1, int l2) const
{
  if (l2 < 0) l2 = l1;

  if ((l2 < 0) || (l1 >= _attr._t)) {
    cerr << "GenericImage<VoxelType>::GetFrame: Parameter out of range\n";
    exit(1);
  }

  if (l1 < 0) l1 = 0;
  if (l2 >= _attr._t) l2 = _attr._t - 1;

  // Initialize
  ImageAttributes attr = this->Attributes();
  attr._t       = l2 - l1 + 1;
  attr._torigin = this->ImageToTime(l1);
  image.Initialize(attr);

  // Copy region
  for (int l = l1; l <= l2; l++) {
    for (int k = 0; k < _attr._z; k++) {
      for (int j = 0; j < _attr._y; j++) {
        for (int i = 0; i < _attr._x; i++) {
          image._matrix[l-l1][k][j][i] = _matrix[l][k][j][i];
        }
      }
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::GetFrame(int l1, int l2) const
{
  GenericImage<VoxelType> image;
  this->GetFrame(image, l1, l2);
  return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::GetFrame(BaseImage *&base, int l1, int l2) const
{
  GenericImage<VoxelType> *image = dynamic_cast<GenericImage<VoxelType> *>(base);
  if (image == NULL) {
    delete base;
    image = new GenericImage<VoxelType>();
    base  = image;
  }
  this->GetFrame(*image, l1, l2);
}

// =============================================================================
// Image arithmetic
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator+=(const GenericImage &image)
{
  if (image.Attributes() != this->Attributes()) {
    cerr << "GenericImage<VoxelType>::operator+=: Size mismatch in images" << endl;
    this->Attributes().Print();
    image.Attributes().Print();
    exit(1);
  }
  VoxelType       *ptr1 = this->Data();
  const VoxelType *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx)) ptr1[idx] += ptr2[idx];
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator-=(const GenericImage &image)
{
  if (image.Attributes() != this->Attributes()) {
    cerr << "GenericImage<VoxelType>::operator-=: Size mismatch in images" << endl;
    this->Attributes().Print();
    image.Attributes().Print();
    exit(1);
  }
  VoxelType       *ptr1 = this->Data();
  const VoxelType *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx)) ptr1[idx] -= ptr2[idx];
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator*=(const GenericImage &image)
{
  if (image.Attributes() != this->Attributes()) {
    cerr << "GenericImage<VoxelType>::operator*=: Size mismatch in images" << endl;
    this->Attributes().Print();
    image.Attributes().Print();
    exit(1);
  }
  VoxelType       *ptr1 = this->Data();
  const VoxelType *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx)) ptr1[idx] *= ptr2[idx];
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator/=(const GenericImage &image)
{
  if (image.Attributes() != this->Attributes()) {
    cerr << "GenericImage<VoxelType>::operator/=: Size mismatch in images" << endl;
    this->Attributes().Print();
    image.Attributes().Print();
    exit(1);
  }
  VoxelType       *ptr1 = this->Data();
  const VoxelType *ptr2 = image.Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx)) {
      if (ptr2[idx] == VoxelType()) {
        if (HasBackgroundValue()) {
          ptr1[idx] = voxel_cast<VoxelType>(GetBackgroundValueAsDouble());
        } else {
          ptr1[idx] = VoxelType();
        }
      } else {
        ptr1[idx] /= ptr2[idx];
      }
    }
  }
  return *this;
}

template <> GenericImage<float3x3 > &GenericImage<float3x3 >::operator/=(const GenericImage &)
{
  cerr << "GenericImage<float3x3>::operator /=: Not implemented" << endl;
  exit(1);
}

template <> GenericImage<double3x3> &GenericImage<double3x3>::operator/=(const GenericImage &)
{
  cerr << "GenericImage<double3x3>::operator /=: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator+=(double scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      ptr[idx] += voxel_cast<VoxelType>(scalar);
    }
  }
  return *this;
}

template <> GenericImage<float3x3 > &GenericImage<float3x3 >::operator+=(double scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      ptr[idx] += static_cast<float>(scalar);
    }
  }
  return *this;
}

template <> GenericImage<double3x3> &GenericImage<double3x3>::operator+=(double scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      ptr[idx] += scalar;
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator-=(double scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      ptr[idx] -= voxel_cast<VoxelType>(scalar);
    }
  }
  return *this;
}

template <> GenericImage<float3x3 > &GenericImage<float3x3 >::operator-=(double scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      ptr[idx] -= static_cast<float>(scalar);
    }
  }
  return *this;
}

template <> GenericImage<double3x3> &GenericImage<double3x3>::operator-=(double scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      ptr[idx] -= scalar;
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator*=(double scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      ptr[idx] *= voxel_cast<VoxelType>(scalar);
    }
  }
  return *this;
}

template <> GenericImage<float3x3 > &GenericImage<float3x3 >::operator*=(double scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      ptr[idx] *= static_cast<float>(scalar);
    }
  }
  return *this;
}

template <> GenericImage<double3x3> &GenericImage<double3x3>::operator*=(double scalar)
{
  VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (IsForeground(idx)) {
      ptr[idx] *= scalar;
    }
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator/=(double scalar)
{
  if (scalar != .0) {
    VoxelType *ptr = this->Data();
    for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
      if (IsForeground(idx)) {
        ptr[idx] /= voxel_cast<VoxelType>(scalar);
      }
    }
  } else {
    cerr << "GenericImage<VoxelType>::operator/=: Division by zero" << endl;
  }
  return *this;
}

template <> GenericImage<float3x3 > &GenericImage<float3x3 >::operator/=(double scalar)
{
  if (scalar != 0.) {
    VoxelType *ptr = this->Data();
    for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
      if (IsForeground(idx)) {
        ptr[idx] /= static_cast<float>(scalar);
      }
    }
  } else {
    cerr << "GenericImage<VoxelType>::operator/=: Division by zero" << endl;
  }
  return *this;
}

template <> GenericImage<double3x3> &GenericImage<double3x3>::operator/=(double scalar)
{
  if (scalar != 0.) {
    VoxelType *ptr = this->Data();
    for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
      if (IsForeground(idx)) {
        ptr[idx] /= scalar;
      }
    }
  } else {
    cerr << "GenericImage<VoxelType>::operator/=: Division by zero" << endl;
  }
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator+(const GenericImage &image) const
{
  GenericImage<VoxelType> tmp(*this);
  tmp += image;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator-(const GenericImage &image) const
{
  GenericImage<VoxelType> tmp(*this);
  tmp -= image;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator*(const GenericImage &image) const
{
  GenericImage<VoxelType> tmp(*this);
  tmp *= image;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator/(const GenericImage &image) const
{
  GenericImage<VoxelType> tmp(*this);
  tmp /= image;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator+(double scalar) const
{
  GenericImage<VoxelType> tmp(*this);
  tmp += scalar;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator-(double scalar) const
{
  GenericImage<VoxelType> tmp(*this);
  tmp -= scalar;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator*(double scalar) const
{
  GenericImage<VoxelType> tmp(*this);
  tmp *= scalar;
  return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator/(double scalar) const
{
  GenericImage<VoxelType> tmp(*this);
  tmp /= scalar;
  return tmp;
}

// =============================================================================
// Thresholding
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::PutBackgroundValueAsDouble(double value, bool threshold)
{
  BaseImage::PutBackgroundValueAsDouble(value);
  if (threshold && !IsNaN(value)) {
    const VoxelType bg = voxel_cast<VoxelType>(this->_bg);
    VoxelType *ptr = this->GetPointerToVoxels();
    for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
      if (*ptr < bg) *ptr = bg;
    }
  }
}

template <> void GenericImage<float3x3>::PutBackgroundValueAsDouble(double value, bool threshold)
{
  BaseImage::PutBackgroundValueAsDouble(value, threshold);
}

template <> void GenericImage<double3x3>::PutBackgroundValueAsDouble(double value, bool threshold)
{
  BaseImage::PutBackgroundValueAsDouble(value, threshold);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> &GenericImage<VoxelType>::operator>=(VoxelType pixel)
{
  VoxelType *ptr = this->GetPointerToVoxels();
  for (int idx = 0; idx < _NumberOfVoxels; idx++) {
    if (IsForeground(idx) && ptr[idx] > pixel) ptr[idx] = pixel;
  }
  return *this;
}

template <> GenericImage<float3x3> &GenericImage<float3x3>::operator>=(float3x3)
{
  cerr << "GenericImage<float3x3 >::operator >=: Not implemented" << endl;
  exit(1);
}

template <> GenericImage<double3x3> &GenericImage<double3x3>::operator>=(double3x3)
{
  cerr << "GenericImage<double3x3 >::operator >=: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator<=(VoxelType pixel)
{
  VoxelType *ptr = this->GetPointerToVoxels();
  for (int idx = 0; idx < _NumberOfVoxels; idx++) {
    if (IsForeground(idx) && ptr[idx] < pixel) ptr[idx] = pixel;
  }
  return *this;
}

template <> GenericImage<float3x3> &GenericImage<float3x3>::operator<=(float3x3)
{
  cerr << "GenericImage<float3x3>::operator <=: Not implemented" << endl;
  exit(1);
}

template <> GenericImage<double3x3> &GenericImage<double3x3>::operator<=(double3x3)
{
  cerr << "GenericImage<double3x3>::operator <=: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator>(VoxelType pixel) const
{
  GenericImage<VoxelType> image(*this);
  image >= pixel;
  return image;
}

template <> GenericImage<float3x3> GenericImage<float3x3>::operator>(float3x3) const
{
  cerr << "GenericImage<float3x3>::operator >: Not implemented" << endl;
  exit(1);
}

template <> GenericImage<double3x3> GenericImage<double3x3>::operator>(double3x3) const
{
  cerr << "GenericImage<double3x3>::operator >: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType> GenericImage<VoxelType>::operator<(VoxelType pixel) const
{
  GenericImage<VoxelType> image(*this);
  image <= pixel;
  return image;
}

template <> GenericImage<float3x3> GenericImage<float3x3>::operator<(float3x3 ) const
{
  cerr << "GenericImage<float3x3>::operator <: Not implemented" << endl;
  exit(1);
}

template <> GenericImage<double3x3> GenericImage<double3x3>::operator<(double3x3) const
{
  cerr << "GenericImage<double3x3>::operator <: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
BinaryImage GenericImage<VoxelType>::operator!=(VoxelType pixel) const
{
  BinaryImage mask(_attr);
  const VoxelType *ptr1 = this->GetPointerToVoxels();
  BinaryPixel *ptr2 = mask .GetPointerToVoxels();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr1, ++ptr2) {
    *ptr2 = (*ptr1 != pixel);
  }
  return mask;
}

// =============================================================================
// Common image statistics
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::GetMinMax(VoxelType &min, VoxelType &max) const
{
  min = max = VoxelType();
  
  const VoxelType *ptr   = this->Data();
  bool             first = true;
  
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
    if (IsForeground(idx)) {
      if (first) {
        min   = max = *ptr;
        first = false;
      } else {
        if (*ptr < min) min = *ptr;
        if (*ptr > max) max = *ptr;
      }
    }
  }
}

template <> void GenericImage<float3x3>::GetMinMax(VoxelType &, VoxelType &) const
{
  cerr << "GenericImage<float3x3>::GetMinMax: Not implemented" << endl;
  exit(1);
}

template <> void GenericImage<double3x3>::GetMinMax(VoxelType &, VoxelType &) const
{
  cerr << "GenericImage<double3x3>::GetMinMax: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::GetMinMax(VoxelType &min, VoxelType &max, VoxelType pad) const
{
  min = max = VoxelType();

  const VoxelType *ptr   = this->Data();
  bool             first = true;
  
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
    if (*ptr != pad) {
      if (first) {
        min   = max = *ptr;
        first = false;
      } else {
        if (*ptr < min) min = *ptr;
        if (*ptr > max) max = *ptr;
      }
    }
  }
}

template <> void GenericImage<float3x3>::GetMinMax(VoxelType &, VoxelType &, VoxelType) const
{
  cerr << "GenericImage<float3x3>::GetMinMax: Not implemented" << endl;
  exit(1);
}

template <> void GenericImage<double3x3>::GetMinMax(VoxelType &, VoxelType &, VoxelType) const
{
  cerr << "GenericImage<double3x3>::GetMinMax: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::PutMinMax(VoxelType min, VoxelType max)
{
  VoxelType min_val, max_val;
  this->GetMinMax(min_val, max_val);
  VoxelType *ptr = this->Data();
  RealType slope = voxel_cast<RealType>(max  - min)  / voxel_cast<RealType>(max_val - min_val);
  RealType inter = voxel_cast<RealType>(min) - slope * voxel_cast<RealType>(min_val);
  for (int idx = 0; idx < _NumberOfVoxels; ++idx, ++ptr) {
    if (IsForeground(idx)) *ptr = static_cast<VoxelType>(inter + slope * static_cast<RealType>(*ptr));
  }
}

template <> void GenericImage<float3x3>::PutMinMax(VoxelType, VoxelType)
{
  cerr << "GenericImage<float3x3>::PutMinMax: Not implemented" << endl;
  exit(1);
}

template <> void GenericImage<double3x3>::PutMinMax(VoxelType, VoxelType)
{
  cerr << "GenericImage<double3x3>::PutMinMax: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
typename GenericImage<VoxelType>::RealType
GenericImage<VoxelType>::Mean(bool fg) const
{
  int num = 0;
  RealType sum = voxel_cast<RealType>(0);
  const VoxelType *ptr = this->Data();
  for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
    if (!fg || IsForeground(idx)) {
      sum += voxel_cast<RealType>(*ptr);
      num += 1;
    }
    ++ptr;
  }
  return (num > 0 ? sum / static_cast<RealScalarType>(num) : sum);
}

template <> typename GenericImage<float3x3>::RealType GenericImage<float3x3 >::Mean(bool) const
{
  cerr << "GenericImage<float3x3>::Mean: Not implemented" << endl;
  exit(1);
}

template <> typename GenericImage<double3x3>::RealType GenericImage<double3x3>::Mean(bool) const
{
  cerr << "GenericImage<double3x3>::Mean: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
typename GenericImage<VoxelType>::RealType
GenericImage<VoxelType>::GetAverage(int toggle) const
{
  return Mean(toggle == 0 ? false : true);
}

template <> typename GenericImage<float3x3 >::RealType GenericImage<float3x3 >::GetAverage(int) const
{
  cerr << "GenericImage<float3x3>::GetAverage: Not implemented" << endl;
  exit(1);
}

template <> typename GenericImage<double3x3>::RealType GenericImage<double3x3>::GetAverage(int) const
{
  cerr << "GenericImage<double3x3>::GetAverage: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
typename GenericImage<VoxelType>::RealType
GenericImage<VoxelType>::GetSD(int toggle) const
{
  const VoxelType  zero = voxel_cast<VoxelType>(0);
  RealType         std  = voxel_cast<RealType >(0);
  const RealType   avg  = this->GetAverage(toggle);
  const VoxelType *ptr;

  if (toggle) {
    int n = 0;
    ptr = this->Data();
    for (int i = 0; i < _NumberOfVoxels; i++) {
      if (IsForeground(i) && (*ptr) > zero) n++;
      ++ptr;
    }
    const RealType norm = voxel_cast<RealType>(1.0 / n);
    ptr = this->Data();
    for (int i = 0; i < _NumberOfVoxels; i++) {
      if (IsForeground(i) && (*ptr) > zero) {
        std += norm * pow(voxel_cast<RealType>(*ptr) - avg, 2);
      }
      ++ptr;
    }
  } else {
    const RealType norm = voxel_cast<RealType>(1.0 / _NumberOfVoxels);
    ptr = this->Data();
    for (int i = 0; i < _NumberOfVoxels; i++) {
      if (IsForeground(i)) {
        std += norm * pow(voxel_cast<RealType>(*ptr) - avg, 2);
      }
      ++ptr;
    }
  }

  return sqrt(std);
}

template <> typename GenericImage<float3x3>::RealType GenericImage<float3x3 >::GetSD(int) const
{
  cerr << "GenericImage<float3x3>::GetSD: Not implemented" << endl;
  exit(1);
}

template <> typename GenericImage<double3x3>::RealType GenericImage<double3x3>::GetSD(int) const
{
  cerr << "GenericImage<double3x3>::GetSD: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::GetMaxPosition(Point& p, int ds, int t) const
{
  if (ds < 0) ds = -ds;
  
  this->WorldToImage(p);
  const int x = iround(p._x);
  const int y = iround(p._y);
  const int z = iround(p._z);
  
  const VoxelType *ptr = this->Data(0, 0, 0, t);
  VoxelType        max = this->Get(x, y, z, t);
  
  p._z = static_cast<double>(z);
  for (int j = y - ds; j <= y + ds; j++) {
    for (int i = x - ds; i <= x + ds; i++) {
      if (IsForeground(i, j, z, t) && *ptr > max) {
        p._x = static_cast<double>(i);
        p._y = static_cast<double>(j);
        max  = *ptr;
      }
      ++ptr;
    }
  }
  
  this->ImageToWorld(p);
}

template <> void GenericImage<float3x3 >::GetMaxPosition(Point &, int, int) const
{
  cerr << "GenericImage<float3x3>::GetMaxPosition: Not implemented" << endl;
  exit(1);
}

template <> void GenericImage<double3x3>::GetMaxPosition(Point &, int, int) const
{
  cerr << "GenericImage<double3x3>::GetMaxPosition: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::GravityCenter(Point& p, int ds, int t) const
{
  if (ds < 0) ds = -ds;

  this->WorldToImage(p);
  const int x = iround(p._x);
  const int y = iround(p._y);
  const int z = iround(p._z);

  double si = .0, sj = .0, sk = .0, sw = .0, ssw = .0;
  for (int k = z - ds; k <= z + ds; k++) {
    for (int j = y - ds; j <= y + ds; j++) {
      for (int i = x - ds; i <= x + ds; i++) {
        if (IsForeground(i, j, k, t)) {
          sw   = this->GetAsDouble(i, j, k, t);
          si  += i * sw;
          sj  += j * sw;
          sk  += k * sw;
          ssw +=     sw;
        }
      }
    }
  }

  if (ssw) p._x = si/ssw, p._y = sj/ssw, p._z = sk/ssw;
  this->ImageToWorld(p);
}

template <> void GenericImage<float3x3 >::GravityCenter(Point &, int, int) const
{
  cerr << "GenericImage<float3x3>::GravityCenter: Not implemented" << endl;
  exit(1);
}

template <> void GenericImage<double3x3>::GravityCenter(Point &, int, int) const
{
  cerr << "GenericImage<double3x3>::GravityCenter: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Common image manipulations
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::ReflectX(bool modify_axes)
{
  for (int t = 0; t < _attr._t; ++t)
  for (int z = 0; z < _attr._z; ++z)
  for (int y = 0; y < _attr._y; ++y)
  for (int x = 0; x < _attr._x / 2; ++x) {
    swap(_matrix[t][z][y][x], _matrix[t][z][y][_attr._x-(x+1)]);
  }
  if (modify_axes) {
    _attr._xaxis[0] = -_attr._xaxis[0];
    _attr._xaxis[1] = -_attr._xaxis[1];
    _attr._xaxis[2] = -_attr._xaxis[2];
    UpdateMatrix();
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::ReflectY(bool modify_axes)
{
  for (int t = 0; t < _attr._t; ++t)
  for (int z = 0; z < _attr._z; ++z)
  for (int y = 0; y < _attr._y / 2; ++y)
  for (int x = 0; x < _attr._x; ++x) {
    swap(_matrix[t][z][y][x], _matrix[t][z][_attr._y-(y+1)][x]);
  }
  if (modify_axes) {
    _attr._yaxis[0] = -_attr._yaxis[0];
    _attr._yaxis[1] = -_attr._yaxis[1];
    _attr._yaxis[2] = -_attr._yaxis[2];
    UpdateMatrix();
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::ReflectZ(bool modify_axes)
{
  for (int t = 0; t < _attr._t; ++t)
  for (int z = 0; z < _attr._z / 2; ++z)
  for (int y = 0; y < _attr._y; ++y)
  for (int x = 0; x < _attr._x; ++x) {
    swap(_matrix[t][z][y][x], _matrix[t][_attr._z-(z+1)][y][x]);
  }
  if (modify_axes) {
    _attr._zaxis[0] = -_attr._zaxis[0];
    _attr._zaxis[1] = -_attr._zaxis[1];
    _attr._zaxis[2] = -_attr._zaxis[2];
    UpdateMatrix();
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::ReflectT(bool modify_axes)
{
  for (int t = 0; t < _attr._t / 2; ++t)
  for (int z = 0; z < _attr._z; ++z)
  for (int y = 0; y < _attr._y; ++y)
  for (int x = 0; x < _attr._x; ++x) {
    swap(_matrix[t][z][y][x], _matrix[_attr._t-(t+1)][z][y][x]);
  }
  if (modify_axes) {
    _attr._dt = -_attr._dt;
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::FlipXY(bool modify_origin)
{
  // TODO: Implement BaseImage::FlipXY which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._y, _attr._x, _attr._z, _attr._t);

  // Flip image in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[l][k][i][j] = _matrix[l][k][j][i];
  }

  // Swap image dimensions
  swap(_attr._x, _attr._y);

  // Swap voxel dimensions
  swap(_attr._dx, _attr._dy);

  // Swap origin coordinates
  if (modify_origin) swap(_attr._xorigin, _attr._yorigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);

  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::FlipXZ(bool modify_origin)
{
  // TODO: Implement BaseImage::FlipXZ which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._z, _attr._y, _attr._x, _attr._t);

  // Flip image in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[l][i][j][k] = _matrix[l][k][j][i];
  }

  // Swap image dimensions
  swap(_attr._x, _attr._z);

  // Swap voxel dimensions
  swap(_attr._dx, _attr._dz);

  // Swap origin coordinates
  if (modify_origin) swap(_attr._xorigin, _attr._zorigin);

  // Reshape image data
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);

  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::FlipYZ(bool modify_origin)
{
  // TODO: Implement BaseImage::FlipYZ which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory for flipped image
  VoxelType ****matrix = Allocate<VoxelType>(_attr._x, _attr._z, _attr._y, _attr._t);

  // Flip image in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[l][j][k][i] = _matrix[l][k][j][i];
  }

  // Swap image dimensions
  swap(_attr._y, _attr._z);

  // Swap voxel dimensions
  swap(_attr._dy, _attr._dz);

  // Swap origin coordinates
  if (modify_origin) swap(_attr._yorigin, _attr._zorigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);

  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::FlipXT(bool modify_origin)
{
  // TODO: Implement BaseImage::FlipXT which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._t, _attr._y, _attr._z, _attr._x);

  // Flip image in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[i][k][j][l] = _matrix[l][k][j][i];
  }

  // Swap image dimensions
  swap(_attr._x, _attr._t);

  // Swap voxel dimensions
  swap(_attr._dx, _attr._dt);

  // Swap origin coordinates
  if (modify_origin) swap(_attr._xorigin, _attr._torigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);

  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::FlipYT(bool modify_origin)
{
  // TODO: Implement BaseImage::FlipYT which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._x, _attr._t, _attr._z, _attr._y);

  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[j][k][l][i] = _matrix[l][k][j][i];
  }

  // Swap image dimensions
  swap(_attr._y, _attr._t);

  // Swap voxel dimensions
  swap(_attr._dy, _attr._dt);

  // Swap origin coordinates
  if (modify_origin) swap(_attr._yorigin, _attr._torigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);

  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::FlipZT(bool modify_origin)
{
  // TODO: Implement BaseImage::FlipZT which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._x, _attr._y, _attr._t, _attr._z);

  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[k][l][j][i] = _matrix[l][k][j][i];
  }

  // Swap image dimensions
  swap(_attr._z, _attr._t);

  // Swap voxel dimensions
  swap(_attr._dz, _attr._dt);

  // Swap origin coordinates
  if (modify_origin) swap(_attr._zorigin, _attr._torigin);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);

  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::SwapXY(bool modify_axes)
{
  // TODO: Implement BaseImage::FlipXY which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._y, _attr._x, _attr._z, _attr._t);

  // Swap image dimensions in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[l][k][i][j] = _matrix[l][k][j][i];
  }
  swap(_attr._x, _attr._y);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);

  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  if (modify_axes) {
    swap(_attr._dx, _attr._dy);
    for (int d = 0; d < 3; ++d) {
      swap(_attr._xaxis[d], _attr._yaxis[d]);
    }
  }
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::SwapXZ(bool modify_axes)
{
  // TODO: Implement BaseImage::FlipXZ which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._z, _attr._y, _attr._x, _attr._t);

  // Swap image dimensions in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[l][i][j][k] = _matrix[l][k][j][i];
  }
  swap(_attr._x, _attr._z);

  // Reshape image data
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);

  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  if (modify_axes) {
    swap(_attr._dx, _attr._dz);
    for (int d = 0; d < 3; ++d) {
      swap(_attr._xaxis[d], _attr._zaxis[d]);
    }
  }
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::SwapYZ(bool modify_axes)
{
  // TODO: Implement BaseImage::FlipYZ which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory for flipped image
  VoxelType ****matrix = Allocate<VoxelType>(_attr._x, _attr._z, _attr._y, _attr._t);

  // Swap image dimensions in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[l][j][k][i] = _matrix[l][k][j][i];
  }
  swap(_attr._y, _attr._z);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);
  
  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  if (modify_axes) {
    swap(_attr._dy, _attr._dz);
    for (int d = 0; d < 3; ++d) {
      swap(_attr._yaxis[d], _attr._zaxis[d]);
    }
  }
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::SwapXT(bool modify_axes)
{
  // TODO: Implement BaseImage::FlipXT which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._t, _attr._y, _attr._z, _attr._x);

  // Swap image dimensions in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[i][k][j][l] = _matrix[l][k][j][i];
  }
  swap(_attr._x, _attr._t);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);
  
  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  if (modify_axes) {
    swap(_attr._dx, _attr._dt);
    _attr._xaxis[0] = 1.;
    _attr._xaxis[1] = 0.;
    _attr._xaxis[2] = 0.;
  }
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::SwapYT(bool modify_axes)
{
  // TODO: Implement BaseImage::FlipYT which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._x, _attr._t, _attr._z, _attr._y);

  // Swap image dimensions in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[j][k][l][i] = _matrix[l][k][j][i];
  }
  swap(_attr._y, _attr._t);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);
  
  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  if (modify_axes) {
    swap(_attr._dy, _attr._dt);
    _attr._yaxis[0] = 1.;
    _attr._yaxis[1] = 0.;
    _attr._yaxis[2] = 0.;
  }
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::SwapZT(bool modify_axes)
{
  // TODO: Implement BaseImage::FlipZT which flips the foreground mask (if any),
  //       adjusts the attributes, and updates the coordinate transformation matrices.
  //       The subclass then only needs to reshape the image _matrix data itself.

  // Allocate memory
  VoxelType ****matrix = Allocate<VoxelType>(_attr._x, _attr._y, _attr._t, _attr._z);

  // Swap image dimensions in memory
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i) {
    matrix[k][l][j][i] = _matrix[l][k][j][i];
  }
  swap(_attr._z, _attr._t);

  // Reshape image matrix
  //
  // Attention: DO NOT just swap the pointers to the data elements as this
  //            changes the memory location of the image data. This is not
  //            predictable by users of the class which may still hold a pointer
  //            to the old memory and in particular complicates the synchronization
  //            of host and device memory in CUGenericImage used by CUDA code.
  _matrix = Reshape(_matrix, _attr._x, _attr._y, _attr._z, _attr._t);
  
  // Copy flipped image
  CopyFrom(matrix[0][0][0]);

  // Deallocate memory
  Deallocate(matrix);

  // Update coordinate transformation
  if (modify_axes) {
    swap(_attr._dz, _attr._dt);
    _attr._zaxis[0] = 1.;
    _attr._zaxis[1] = 0.;
    _attr._zaxis[2] = 0.;
  }
  UpdateMatrix();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
bool GenericImage<VoxelType>::CropPad(int margin)
{
  // Get bounding box of image foreground
  int i1, j1, k1, l1, i2, j2, k2, l2;
  this->BoundingBox(i1, j1, k1, l1, i2, j2, k2, l2);
  // Do nothing if all voxels are background, but report it
  if (i1 > i2 || j1 > j2 || k1 > k2 || l1 > l2) return false;
  // Convert upper index bounds to margin widths
  i2 = (_attr._x - 1) - i2;
  j2 = (_attr._y - 1) - j2;
  k2 = (_attr._z - 1) - k2;
  l2 = (_attr._t - 1) - l2;
  // Leave a margin of background voxels with specified width
  // Note: Negative value gives the number of voxels to add.
  if (_attr._x > 1) i1 -= margin, i2 -= margin;
  if (_attr._y > 1) j1 -= margin, j2 -= margin;
  if (_attr._z > 1) k1 -= margin, k2 -= margin;
  if (_attr._t > 1) l1 -= margin, l2 -= margin;
  // Do nothing, if nothing to be done
  if (i1 == 0 && i2 == 0 && j1 == 0 && j2 == 0 &&
      k1 == 0 && k2 == 0 && l1 == 0 && l2 == 0) return true;
  // Adjust image lattice
  ImageAttributes attr = _attr;
  attr._x -= i1 + i2;
  attr._y -= j1 + j2;
  attr._z -= k1 + k2;
  attr._t -= l1 + l2;
  attr._xorigin = 0.5 * ((_attr._x - 1) + (i1 - i2));
  attr._yorigin = 0.5 * ((_attr._y - 1) + (j1 - j2));
  attr._zorigin = 0.5 * ((_attr._z - 1) + (k1 - k2));
  _attr.LatticeToWorld(attr._xorigin, attr._yorigin, attr._zorigin);
  attr._torigin = _attr.LatticeToTime(l1);
  // Convert upper margin widths to index bounds
  i2 = (_attr._x - 1) - i2;
  j2 = (_attr._y - 1) - j2;
  k2 = (_attr._z - 1) - k2;
  l2 = (_attr._t - 1) - l2;
  // Copy remaining voxels and pad lattice where needed
  const int nvoxels = attr.NumberOfLatticePoints();
  VoxelType *data       = Allocate<VoxelType>(nvoxels);
  VoxelType *data_iter  = data;
  for (int l = l1; l <= l2; ++l)
  for (int k = k1; k <= k2; ++k)
  for (int j = j1; j <= j2; ++j)
  for (int i = i1; i <= i2; ++i, ++data_iter) {
    if (0 <= i && i < _attr._x &&
        0 <= j && j < _attr._y &&
        0 <= k && k < _attr._z &&
        0 <= l && l < _attr._t) {
      (*data_iter) = Get(i, j, k, l);
    } else {
      // Padded voxel to extend margin
      (*data_iter) = voxel_cast<VoxelType>(_bg);
    }
  }
  // Initialize new image lattice
  this->Initialize(attr);
  data_iter = data;
  for (int l = 0; l < _attr._t; ++l)
  for (int k = 0; k < _attr._z; ++k)
  for (int j = 0; j < _attr._y; ++j)
  for (int i = 0; i < _attr._x; ++i, ++data_iter) {
    Put(i, j, k, l, (*data_iter));
  }
  // Free temporary allocated memory
  Deallocate(data);
  return true;
}

// =============================================================================
// VTK interface
// =============================================================================
#if MIRTK_Image_WITH_VTK

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::ImageToVTK(vtkStructuredPoints *vtk) const
{
  if (this->ImageToVTKScalarType() == VTK_VOID) {
    cerr << "GenericImage::ImageToVTK: Cannot convert image to VTK structured points" << endl;
    exit(1);
  }
  double x = 0, y = 0, z = 0;
  this->ImageToWorld(x, y, z);
  vtk->SetOrigin    (x, y, z);
  vtk->SetDimensions(_attr._x,  _attr._y,  _attr._z);
  vtk->SetSpacing   (_attr._dx, _attr._dy, _attr._dz);
#if VTK_MAJOR_VERSION >= 6
  vtk->AllocateScalars(this->ImageToVTKScalarType(), 1);
#else
  vtk->SetScalarType(this->ImageToVTKScalarType());
  vtk->AllocateScalars();
#endif
  const int        nvox = _attr._x * _attr._y * _attr._z;
  const VoxelType *ptr1 = this->Data();
  VoxelType       *ptr2 = reinterpret_cast<VoxelType *>(vtk->GetScalarPointer());
  for (int i = 0; i < nvox; ++i, ++ptr1) {
    for (int l = 0; l < _attr._t; ++l, ++ptr2) *ptr2 = ptr1[l * nvox];
  }
}

// -----------------------------------------------------------------------------
template <class Type>
void GenericImage<Type>::VTKToImage(vtkStructuredPoints *)
{
  cerr << this->NameOfClass() << "::VTKToImage: Not implemented" << endl;
  exit(1);
}

#endif // MIRTK_Image_WITH_VTK
// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::Read(const char *fname)
{
  // Read image
  UniquePtr<ImageReader> reader(ImageReader::New(fname));
  UniquePtr<BaseImage>   image(reader->Run());
  // Convert image
  switch (image->GetDataType()) {
    case MIRTK_VOXEL_CHAR:           { *this = *(dynamic_cast<GenericImage<char>           *>(image.get())); } break;
    case MIRTK_VOXEL_UNSIGNED_CHAR:  { *this = *(dynamic_cast<GenericImage<unsigned char>  *>(image.get())); } break;
    case MIRTK_VOXEL_SHORT:          { *this = *(dynamic_cast<GenericImage<short>          *>(image.get())); } break;
    case MIRTK_VOXEL_UNSIGNED_SHORT: { *this = *(dynamic_cast<GenericImage<unsigned short> *>(image.get())); } break;
    case MIRTK_VOXEL_INT:            { *this = *(dynamic_cast<GenericImage<int>            *>(image.get())); } break;
    case MIRTK_VOXEL_FLOAT:          { *this = *(dynamic_cast<GenericImage<float>          *>(image.get())); } break;
    case MIRTK_VOXEL_DOUBLE:         { *this = *(dynamic_cast<GenericImage<double>         *>(image.get())); } break;
    default:
      cerr << this->NameOfClass() << "::Read: Unknown data type: " << image->GetDataType() << endl;
      exit(1);
  }
  // Apply rescaling function
  if (reader->Slope() != .0 && reader->Slope() != 1.0) {
    *this *= reader->Slope();
  }
  if (reader->Intercept() != .0) {
    *this += reader->Intercept();
  }
}

template <> void GenericImage<float3x3>::Read(const char *)
{
  cerr << "GenericImage<float3x3>::Read: Not implemented" << endl;
  exit(1);
}

template <> void GenericImage<double3x3>::Read(const char *)
{
  cerr << "GenericImage<double3x3>::Read: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void GenericImage<VoxelType>::Write(const char *fname) const
{
  string name(fname);
  if (Extension(fname).empty()) name += MIRTK_Image_DEFAULT_EXT;
  UniquePtr<ImageWriter> writer(ImageWriter::New(name.c_str()));
  writer->Input(this);
  writer->Run();
}

template <> void GenericImage<float3x3>::Write(const char *) const
{
  cerr << "GenericImage<float3x3>::Write: Not implemented" << endl;
  exit(1);
}

template <> void GenericImage<double3x3>::Write(const char *) const
{
  cerr << "GenericImage<double3x3>::Read: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// Explicit template instantiations
// =============================================================================

template class GenericImage<char>;
template class GenericImage<unsigned char>;
template class GenericImage<short>;
template class GenericImage<unsigned short>;
template class GenericImage<int>;
template class GenericImage<unsigned int>;
template class GenericImage<float>;
template class GenericImage<float2>;
template class GenericImage<float3>;
template class GenericImage<Float3>;
template class GenericImage<float4>;
template class GenericImage<Float9>;
template class GenericImage<double>;
template class GenericImage<double2>;
template class GenericImage<double3>;
template class GenericImage<Double3>;
template class GenericImage<double4>;
template class GenericImage<Double9>;
template class GenericImage<float3x3>;
template class GenericImage<double3x3>;

template GenericImage<char>::GenericImage(const GenericImage<unsigned char> &);
template GenericImage<char>::GenericImage(const GenericImage<short> &);
template GenericImage<char>::GenericImage(const GenericImage<unsigned short> &);
template GenericImage<char>::GenericImage(const GenericImage<int> &);
template GenericImage<char>::GenericImage(const GenericImage<float> &);
template GenericImage<char>::GenericImage(const GenericImage<double> &);

template GenericImage<unsigned char>::GenericImage(const GenericImage<char> &);
template GenericImage<unsigned char>::GenericImage(const GenericImage<short> &);
template GenericImage<unsigned char>::GenericImage(const GenericImage<unsigned short> &);
template GenericImage<unsigned char>::GenericImage(const GenericImage<int> &);
template GenericImage<unsigned char>::GenericImage(const GenericImage<float> &);
template GenericImage<unsigned char>::GenericImage(const GenericImage<double> &);

template GenericImage<short>::GenericImage(const GenericImage<char> &);
template GenericImage<short>::GenericImage(const GenericImage<unsigned char> &);
template GenericImage<short>::GenericImage(const GenericImage<unsigned short> &);
template GenericImage<short>::GenericImage(const GenericImage<int> &);
template GenericImage<short>::GenericImage(const GenericImage<float> &);
template GenericImage<short>::GenericImage(const GenericImage<double> &);

template GenericImage<unsigned short>::GenericImage(const GenericImage<char> &);
template GenericImage<unsigned short>::GenericImage(const GenericImage<unsigned char> &);
template GenericImage<unsigned short>::GenericImage(const GenericImage<short> &);
template GenericImage<unsigned short>::GenericImage(const GenericImage<int> &);
template GenericImage<unsigned short>::GenericImage(const GenericImage<float> &);
template GenericImage<unsigned short>::GenericImage(const GenericImage<double> &);

template GenericImage<int>::GenericImage(const GenericImage<char> &);
template GenericImage<int>::GenericImage(const GenericImage<unsigned char> &);
template GenericImage<int>::GenericImage(const GenericImage<short> &);
template GenericImage<int>::GenericImage(const GenericImage<unsigned short> &);
template GenericImage<int>::GenericImage(const GenericImage<float> &);
template GenericImage<int>::GenericImage(const GenericImage<double> &);

template GenericImage<float>::GenericImage(const GenericImage<char> &);
template GenericImage<float>::GenericImage(const GenericImage<unsigned char> &);
template GenericImage<float>::GenericImage(const GenericImage<short> &);
template GenericImage<float>::GenericImage(const GenericImage<unsigned short> &);
template GenericImage<float>::GenericImage(const GenericImage<int> &);
template GenericImage<float>::GenericImage(const GenericImage<double> &);

template GenericImage<double>::GenericImage(const GenericImage<char> &);
template GenericImage<double>::GenericImage(const GenericImage<unsigned char> &);
template GenericImage<double>::GenericImage(const GenericImage<short> &);
template GenericImage<double>::GenericImage(const GenericImage<unsigned short> &);
template GenericImage<double>::GenericImage(const GenericImage<int> &);
template GenericImage<double>::GenericImage(const GenericImage<float> &);

template GenericImage<float2  >::GenericImage(const GenericImage<double2  > &);
template GenericImage<float3  >::GenericImage(const GenericImage<double3  > &);
template GenericImage<Float3  >::GenericImage(const GenericImage<Double3  > &);
template GenericImage<float4  >::GenericImage(const GenericImage<double4  > &);
template GenericImage<Float9  >::GenericImage(const GenericImage<Double9  > &);
template GenericImage<float3x3>::GenericImage(const GenericImage<double3x3> &);

template GenericImage<double2  >::GenericImage(const GenericImage<float2  > &);
template GenericImage<double3  >::GenericImage(const GenericImage<float3  > &);
template GenericImage<Double3  >::GenericImage(const GenericImage<Float3  > &);
template GenericImage<double4  >::GenericImage(const GenericImage<float4  > &);
template GenericImage<Double9  >::GenericImage(const GenericImage<Float9  > &);
template GenericImage<double3x3>::GenericImage(const GenericImage<float3x3> &);


} // namespace mirtk
