/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Antonios Makropoulos
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
#include "mirtk/HashImage.h"

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

// Default output image file name extension used by HashImage::Write
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
void HashImage<VoxelType>::AllocateImage()
{
    // Delete existing mask (if any)
    if (_maskOwner) Delete(_mask);
    // Free previously allocated memory
    _Data.clear();
    _DefaultValue = VoxelType();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>::HashImage(){}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>::HashImage(const char *fname)
{
    Read(fname);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>::HashImage(int x, int y, int z, int t)
{
    ImageAttributes attr;
    attr._x = x;
    attr._y = y;
    attr._z = z;
    attr._t = t;
    PutAttributes(attr);
    AllocateImage();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>::HashImage(int x, int y, int z, int t, int n)
{
    if (t > 1 && n > 1) {
        cerr << "HashImage::HashImage: 5D images not supported! Use 4D image with vector voxel type instead." << endl;
        exit(1);
    }
    ImageAttributes attr;
    if (n > 1) t = n, attr._dt = .0; // i.e., vector image with n components
    attr._x = x;
    attr._y = y;
    attr._z = z;
    attr._t = t;
    PutAttributes(attr);
    AllocateImage();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>::HashImage(const ImageAttributes &attr)
    :
      BaseImage(attr)
{
    AllocateImage();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>::HashImage(const ImageAttributes &attr, int n)
    :
      BaseImage(attr, n)
{
    AllocateImage();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>::HashImage(const BaseImage &image)
    :
      BaseImage(image)
{
    // Initialize image
    AllocateImage();
    // Copy/cast data
    CopyFrom(image);
}


// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>::HashImage(const HashImage<VoxelType> &image)
    :
      BaseImage(image)
{
    AllocateImage();
    CopyFrom(image);
}

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
HashImage<VoxelType>::HashImage(const HashImage<VoxelType2> &image)
    : BaseImage(image)
{
    AllocateImage();
    CopyFrom(image);
}

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
HashImage<VoxelType>::HashImage(const GenericImage<VoxelType2> &image)
    : BaseImage(image)
{
    AllocateImage();
    CopyFrom(image);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>::~HashImage()
{
    Clear();
}


// -----------------------------------------------------------------------------
template <class VoxelType> void HashImage<VoxelType>::Clear()
{
    _Data.clear();
    if (_maskOwner) Delete(_mask);
    _attr = ImageAttributes();
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
BaseImage *HashImage<VoxelType>::Copy() const
{
    return new HashImage<VoxelType>(*this);
}

// -----------------------------------------------------------------------------
template <class VoxelType> void HashImage<VoxelType>::Initialize()
{
    *this = VoxelType();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::Initialize(const ImageAttributes &a, int n)
{
    // Initialize attributes
    ImageAttributes attr(a);
    if (n >= 1) attr._t = n, attr._dt = .0; // i.e., vector image with n components
    // Initialize memory
    if (_attr._x != attr._x || _attr._y != attr._y || _attr._z != attr._z || _attr._t != attr._t) {
        PutAttributes(attr);
        AllocateImage();
    } else {
        PutAttributes(attr);
        _DefaultValue = VoxelType();
        *this = VoxelType();
    }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::Initialize(int x, int y, int z, int t, int n)
{
    ImageAttributes attr(_attr);
    attr._x = x;
    attr._y = y;
    attr._z = z;
    attr._t = t;
    this->Initialize(attr, n);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::Initialize(int x, int y, int z, int t)
{
    this->Initialize(x, y, z, t, 1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::CopyFrom(const BaseImage &image)
{
    if (this != &image) {
        VoxelType value;
        _Data.clear();
        _DefaultValue  = VoxelType();

        for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
            Put(idx, voxel_cast<VoxelType>(image.GetAsVector(idx) ));
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
}

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
void HashImage<VoxelType>::CopyFrom(const GenericImage<VoxelType2> &image)
{
    _Data.clear();
    _DefaultValue  = VoxelType();

    for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
        Put(idx, voxel_cast<VoxelType>(image.Get(idx)));
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
template <class VoxelType> template <class VoxelType2>
void HashImage<VoxelType>::CopyFrom(const HashImage<VoxelType2> &image)
{
    _Data.clear();
    _DefaultValue = voxel_cast<VoxelType>(image.DefaultValue());

    VoxelType value;
    for ( auto it = image.Begin(); it != image.End(); ++it ){
        Put(it->first, voxel_cast<VoxelType>(it->second));
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
template <class VoxelType> template <class VoxelType2>
void HashImage<VoxelType>::CopyTo(GenericImage<VoxelType2> &image) const{
    image.Initialize(_attr);
    image = voxel_cast<VoxelType2>(_DefaultValue);
    for ( auto it = Begin(); it != End(); ++it ){
        image.Put(it->first, voxel_cast<VoxelType2>(it->second));
    }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator=(VoxelType scalar)
{
    _Data.clear();
    _DefaultValue=scalar;
    return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator=(const BaseImage &image)
{
    if (this != &image) {
        this->Initialize(image.Attributes());
        this->CopyFrom(image);
    }
    return *this;
}
// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
HashImage<VoxelType>& HashImage<VoxelType>::operator=(const GenericImage<VoxelType2> &image)
{
    this->Initialize(image.Attributes());
    this->CopyFrom(image);

    return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator=(const HashImage &image)
{
    if (this != &image) {
        this->Initialize(image.Attributes());
        this->CopyFrom(image);
    }
    return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
HashImage<VoxelType>& HashImage<VoxelType>::operator=(const HashImage<VoxelType2> &image)
{
    this->Initialize(image.GetImageAttributes());
    CopyFrom(image);
}



// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
bool HashImage<VoxelType>::operator==(const HashImage<VoxelType2> &image) const
{
    if (this->GetImageAttributes() != image.GetImageAttributes()) return false;
    if (_DefaultValue != image.DefaultValue()) return false;
    if (_Data.size() != image.NumberOfNonDefaultVoxels()) return false;

    int idx;
    for (int idx = 0; idx < _NumberOfVoxels; ++idx) {
        if (IsForeground(idx) && image.IsForeground(idx)){
            if(Get(idx)!= voxel_cast<VoxelType>(image.Get(idx)))
                return false;
        }
    }

    return true;
}



// =============================================================================
// Region-of-interest extraction
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>
::GetRegion(HashImage<VoxelType> &image, int k, int m) const
{
    int i, j;
    double x1, y1, z1, t1, x2, y2, z2, t2;

    if ((k < 0) || (k >= _attr._z) || (m < 0) || (m >= _attr._t)) {
        cerr << "HashImage<VoxelType>::GetRegion: Parameter out of range" << endl;
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
    VoxelType value, empty=image.DefaultValue();
    for (j = 0; j < _attr._y; j++) {
        for (i = 0; i < _attr._x; i++) {
            value=Get(i,j,k,m);
            if(value!=empty) image.Put(i,j,0,0, value );
        }
    }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>
::GetRegion(int k, int m) const
{
    HashImage<VoxelType> image;
    this->GetRegion(image, k, m);
    return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>
::GetRegion(BaseImage *&base, int k, int m) const
{
    HashImage<VoxelType> *image = dynamic_cast<HashImage<VoxelType> *>(base);
    if (image == NULL) {
        delete base;
        image = new HashImage<VoxelType>();
        base  = image;
    }
    this->GetRegion(*image, k, m);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>
::GetRegion(HashImage<VoxelType> &image, int i1, int j1, int k1,
            int i2, int j2, int k2) const
{
    int i, j, k, l;
    double x1, y1, z1, x2, y2, z2;

    if ((i1 < 0) || (i1 >= i2) ||
            (j1 < 0) || (j1 >= j2) ||
            (k1 < 0) || (k1 >= k2) ||
            (i2 > _attr._x) || (j2 > _attr._y) || (k2 > _attr._z)) {
        cerr << "HashImage<VoxelType>::GetRegion: Parameter out of range\n";
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
    VoxelType value, empty=image.DefaultValue();
    for (l = 0; l < _attr._t; l++) {
        for (k = k1; k < k2; k++) {
            for (j = j1; j < j2; j++) {
                for (i = i1; i < i2; i++) {
                    value= Get(i,j,k,l);
                    if(value!=empty) image.Put(i-i1,j-j1,k-k1,l,  value );
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>
::GetRegion(int i1, int j1, int k1, int i2, int j2, int k2) const
{
    HashImage<VoxelType> image;
    this->GetRegion(image, i1, j1, k1, i2, j2, k2);
    return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>
::GetRegion(BaseImage *&base, int i1, int j1, int k1, int i2, int j2, int k2) const
{
    HashImage<VoxelType> *image = dynamic_cast<HashImage<VoxelType> *>(base);
    if (image == NULL) {
        delete base;
        image = new HashImage<VoxelType>();
        base  = image;
    }
    this->GetRegion(*image, i1, j1, k1, i2, j2, k2);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>
::GetRegion(HashImage<VoxelType> &image, int i1, int j1, int k1, int l1,
            int i2, int j2, int k2, int l2) const
{
    int i, j, k, l;
    double x1, y1, z1, x2, y2, z2;

    if ((i1 < 0) || (i1 >= i2) ||
            (j1 < 0) || (j1 >= j2) ||
            (k1 < 0) || (k1 >= k2) ||
            (l1 < 0) || (l1 >= l2) ||
            (i2 > _attr._x) || (j2 > _attr._y) || (k2 > _attr._z) || (l2 > _attr._t)) {
        cerr << "HashImage<VoxelType>::GetRegion: Parameter out of range\n";
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
    VoxelType value, empty=image.DefaultValue();
    for (l = l1; l < l2; l++) {
        for (k = k1; k < k2; k++) {
            for (j = j1; j < j2; j++) {
                for (i = i1; i < i2; i++) {
                    value= Get(i,j,k,l);
                    if(value!=empty) image.Put(i-i1,j-j1,k-k1,l-l1,  value );
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>
::GetRegion(int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2) const
{
    HashImage<VoxelType> image;
    this->GetRegion(image, i1, j1, k1, l1, i2, j2, k2, l2);
    return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>
::GetRegion(BaseImage *&base, int i1, int j1, int k1, int l1, int i2, int j2, int k2, int l2) const
{
    HashImage<VoxelType> *image = dynamic_cast<HashImage<VoxelType> *>(base);
    if (image == NULL) {
        delete base;
        image = new HashImage<VoxelType>();
        base  = image;
    }
    this->GetRegion(*image, i1, j1, k1, l1, i2, j2, k2, l2);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::GetFrame(HashImage<VoxelType> &image, int l1, int l2) const
{
    if (l2 < 0) l2 = l1;

    if ((l2 < 0) || (l1 >= _attr._t)) {
        cerr << "HashImage<VoxelType>::GetFrame: Parameter out of range\n";
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
    VoxelType value, empty=image.DefaultValue();
    for (int l = l1; l <= l2; l++) {
        for (int k = 0; k < _attr._z; k++) {
            for (int j = 0; j < _attr._y; j++) {
                for (int i = 0; i < _attr._x; i++) {
                    value=  Get(i,j,k,l);
                    if(value!=empty) image.Put(i,j,k,l-l1, value );
                }
            }
        }
    }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>::GetFrame(int l1, int l2) const
{
    HashImage<VoxelType> image;
    this->GetFrame(image, l1, l2);
    return image;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::GetFrame(BaseImage *&base, int l1, int l2) const
{
    HashImage<VoxelType> *image = dynamic_cast<HashImage<VoxelType> *>(base);
    if (image == NULL) {
        delete base;
        image = new HashImage<VoxelType>();
        base  = image;
    }
    this->GetFrame(*image, l1, l2);
}

// =============================================================================
// Image arithmetic
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator+=(const HashImage &image)
{
    if (image.Attributes() != this->Attributes()) {
        cerr << "HashImage<VoxelType>::operator+=: Size mismatch in images" << endl;
        this->Attributes().Print();
        image.Attributes().Print();
        exit(1);
    }

    for (int i = 0; i < _NumberOfVoxels; i++) {
        if (IsForeground(i) && image.IsForeground(i))
            Put(i, Get(i)+image.Get(i));
    }
    return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator-=(const HashImage &image)
{
    if (image.Attributes() != this->Attributes()) {
        cerr << "HashImage<VoxelType>::operator-=: Size mismatch in images" << endl;
        this->Attributes().Print();
        image.Attributes().Print();
        exit(1);
    }

    for (int i = 0; i < _NumberOfVoxels; i++) {
        if (IsForeground(i) && image.IsForeground(i))
            Put(i, Get(i)-image.Get(i));
    }
    return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator*=(const HashImage &image)
{
    if (image.Attributes() != this->Attributes()) {
        cerr << "HashImage<VoxelType>::operator*=: Size mismatch in images" << endl;
        this->Attributes().Print();
        image.Attributes().Print();
        exit(1);
    }

    for (int i = 0; i < _NumberOfVoxels; i++) {
        if (IsForeground(i) && image.IsForeground(i))
            Put(i, Get(i)*image.Get(i));
    }
    return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator/=(const HashImage &image)
{
    if (image.Attributes() != this->Attributes()) {
        cerr << "HashImage<VoxelType>::operator/=: Size mismatch in images" << endl;
        this->Attributes().Print();
        image.Attributes().Print();
        exit(1);
    }

    VoxelType value;
    for (int i = 0; i < _NumberOfVoxels; i++) {
        if (IsForeground(i) && image.IsForeground(i)){
            value=image.Get(i);
            if(value != VoxelType()){
                value=Get(i)/value;
            }
            Put(i, value);
        }
    }
    return *this;
}

template <> HashImage<float3x3 > &HashImage<float3x3 >::operator/=(const HashImage &)
{
    cerr << "HashImage<float3x3>::operator /=: Not implemented" << endl;
    exit(1);
}

template <> HashImage<double3x3> &HashImage<double3x3>::operator/=(const HashImage &)
{
    cerr << "HashImage<double3x3>::operator /=: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator+=(ScalarType scalar)
{
    _DefaultValue += scalar;
    for ( auto it = Begin(); it != End(); ++it ){
        if (IsForeground(it->first)){
            _Data[it->first] +=scalar;
        }
    }
    return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator-=(ScalarType scalar)
{
    _DefaultValue -= scalar;
    for ( auto it = Begin(); it != End(); ++it ){
        if (IsForeground(it->first)){
            _Data[it->first] -=scalar;
        }
    }
    return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator*=(ScalarType scalar)
{
    _DefaultValue *= scalar;
    if(scalar==0){
        _Data.clear();
    }else{
        for ( auto it = Begin(); it != End(); ++it ){
            if (IsForeground(it->first)){
                _Data[it->first] *=scalar;
            }
        }
    }
    return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType>& HashImage<VoxelType>::operator/=(ScalarType scalar)
{
    if (scalar) {
        _DefaultValue /= scalar;
        for ( auto it = Begin(); it != End(); ++it ){
            if (IsForeground(it->first)){
                _Data[it->first] /=scalar;
            }
        }
    } else {
        cerr << "HashImage<VoxelType>::operator/=: Division by zero" << endl;
    }
    return *this;
}

// -----------------------------------------------------------------------------

template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>::operator+(const HashImage &image) const
{
    HashImage<VoxelType> tmp(*this);
    tmp += image;
    return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>::operator-(const HashImage &image) const
{
    HashImage<VoxelType> tmp(*this);
    tmp -= image;
    return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>::operator*(const HashImage &image) const
{
    HashImage<VoxelType> tmp(*this);
    tmp *= image;
    return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>::operator/(const HashImage &image) const
{
    HashImage<VoxelType> tmp(*this);
    tmp /= image;
    return tmp;
}
// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>::operator+(ScalarType scalar) const
{
    HashImage<VoxelType> tmp(*this);
    tmp += scalar;
    return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>::operator-(ScalarType scalar) const
{
    HashImage<VoxelType> tmp(*this);
    tmp -= scalar;
    return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>::operator*(ScalarType scalar) const
{
    HashImage<VoxelType> tmp(*this);
    tmp *= scalar;
    return tmp;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
HashImage<VoxelType> HashImage<VoxelType>::operator/(ScalarType scalar) const
{
    HashImage<VoxelType> tmp(*this);
    tmp /= scalar;
    return tmp;
}

// =============================================================================
// Thresholding
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::PutBackgroundValueAsDouble(double value, bool threshold)
{
    _bg    = value;
    _bgSet = true;
    if (threshold) {
        const VoxelType bg = voxel_cast<VoxelType>(this->_bg);

        bool changedValue=false;
        if( _DefaultValue < bg){
            changedValue=true;
            _DefaultValue=bg;
        }

        bool modify=false;
        for ( auto it = Begin(); it != End();){
            modify=(it->second < bg);
            if(modify && changedValue){
                // this erase also increases the iterator to the next element...
                it = _Data.erase( it );
            }else {
                if(modify){
                    _Data[it->first]=bg;
                }
                it++;
            }
        }

    }
}

template <> void HashImage<float3x3>::PutBackgroundValueAsDouble(double value, bool threshold)
{
    BaseImage::PutBackgroundValueAsDouble(value, threshold);
}

template <> void HashImage<double3x3>::PutBackgroundValueAsDouble(double value, bool threshold)
{
    BaseImage::PutBackgroundValueAsDouble(value, threshold);
}
// =============================================================================
// Common image statistics
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::GetMinMax(VoxelType &min, VoxelType &max) const
{
    min = max = _DefaultValue;

    VoxelType value;
    for ( auto it = Begin(); it != End(); ++it ){
        if (IsForeground(it->first)) {
            value=it->second;
            if (value < min) min = value;
            if (value > max) max = value;
        }
    }
}

template <> void HashImage<float3x3 >::GetMinMax(VoxelType &, VoxelType &) const
{
    cerr << "HashImage<float3x3>::GetMinMax: Not implemented" << endl;
    exit(1);
}

template <> void HashImage<double3x3>::GetMinMax(VoxelType &, VoxelType &) const
{
    cerr << "HashImage<double3x3>::GetMinMax: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::GetMinMax(VoxelType &min, VoxelType &max, VoxelType pad) const
{
    min = max = _DefaultValue;
    VoxelType value;
    for ( auto it = Begin(); it != End(); ++it ){
        value=it->second;
        if (value != pad) {
            if (value < min) min = value;
            if (value > max) max = value;
        }
    }
}

template <> void HashImage<float3x3 >::GetMinMax(VoxelType &, VoxelType &, VoxelType) const
{
    cerr << "HashImage<float3x3>::GetMinMax: Not implemented" << endl;
    exit(1);
}

template <> void HashImage<double3x3>::GetMinMax(VoxelType &, VoxelType &, VoxelType) const
{
    cerr << "HashImage<double3x3>::GetMinMax: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::PutMinMax(VoxelType min, VoxelType max)
{
    VoxelType min_val, max_val;
    this->GetMinMax(min_val, max_val);
    RealType slope = voxel_cast<RealType>(max  - min)  / voxel_cast<RealType>(max_val - min_val);
    RealType inter = voxel_cast<RealType>(min) - slope * voxel_cast<RealType>(min_val);

    _DefaultValue = inter + slope * _DefaultValue;
    for ( auto it = Begin(); it != End(); ++it ){
        if (IsForeground(it->first)) {
            _Data[it->first]=static_cast<VoxelType>(inter + slope *it->second);
        }
    }
}

template <> void HashImage<float3x3 >::PutMinMax(VoxelType, VoxelType)
{
    cerr << "HashImage<float3x3>::PutMinMax: Not implemented" << endl;
    exit(1);
}

template <> void HashImage<double3x3>::PutMinMax(VoxelType, VoxelType)
{
    cerr << "HashImage<double3x3>::PutMinMax: Not implemented" << endl;
    exit(1);
}


// =============================================================================
// Common image manipulations
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::ReflectX()
{
    VoxelType value;
    for (int t = 0; t < _attr._t; ++t)
        for (int z = 0; z < _attr._z; ++z)
            for (int y = 0; y < _attr._y; ++y)
                for (int x = 0; x < _attr._x / 2; ++x) {
                    value=Get(x,y,z,t);
                    Put(x,y,z,t,  Get(_attr._x-(x+1),y,z,t));
                    Put(_attr._x-(x+1),y,z,t,  value);
                }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::ReflectY()
{
    VoxelType value;
    for (int t = 0; t < _attr._t; ++t)
        for (int z = 0; z < _attr._z; ++z)
            for (int y = 0; y < _attr._y / 2; ++y)
                for (int x = 0; x < _attr._x; ++x) {
                    value=Get(x,y,z,t);
                    Put(x,y,z,t,  Get(x,_attr._y-(y+1),z,t));
                    Put(x,_attr._y-(y+1),z,t,  value);
                }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::ReflectZ()
{
    VoxelType value;
    for (int t = 0; t < _attr._t; ++t)
        for (int z = 0; z < _attr._z / 2; ++z)
            for (int y = 0; y < _attr._y; ++y)
                for (int x = 0; x < _attr._x; ++x) {
                    value=Get(x,y,z,t);
                    Put(x,y,z,t,  Get(x,y,_attr._z-(z+1),t));
                    Put(x,y,_attr._z-(z+1),t,  value);
                }
}


// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::FlipXY(bool modifyOrigin)
{
    // TODO: Implement BaseImage::FlipXY which flips the foreground mask (if any),
    //       adjusts the attributes, and updates the coordinate transformation matrices.
    //       The subclass then only needs to reshape the image _matrix data itself.

    cerr << "HashImage<VoxelType>::FlipXY: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::FlipXZ(bool modifyOrigin)
{
    // TODO: Implement BaseImage::FlipXZ which flips the foreground mask (if any),
    //       adjusts the attributes, and updates the coordinate transformation matrices.
    //       The subclass then only needs to reshape the image _matrix data itself.

    cerr << "HashImage<VoxelType>::FlipXZ: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::FlipYZ(bool modifyOrigin)
{
    // TODO: Implement BaseImage::FlipYZ which flips the foreground mask (if any),
    //       adjusts the attributes, and updates the coordinate transformation matrices.
    //       The subclass then only needs to reshape the image _matrix data itself.

    cerr << "HashImage<VoxelType>::FlipYZ: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::FlipXT(bool modifyOrigin)
{
    // TODO: Implement BaseImage::FlipXT which flips the foreground mask (if any),
    //       adjusts the attributes, and updates the coordinate transformation matrices.
    //       The subclass then only needs to reshape the image _matrix data itself.

    cerr << "HashImage<VoxelType>::FlipXT: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::FlipYT(bool modifyOrigin)
{
    // TODO: Implement BaseImage::FlipYT which flips the foreground mask (if any),
    //       adjusts the attributes, and updates the coordinate transformation matrices.
    //       The subclass then only needs to reshape the image _matrix data itself.

    cerr << "HashImage<VoxelType>::FlipYT: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::FlipZT(bool modifyOrigin)
{
    // TODO: Implement BaseImage::FlipZT which flips the foreground mask (if any),
    //       adjusts the attributes, and updates the coordinate transformation matrices.
    //       The subclass then only needs to reshape the image _matrix data itself.

    cerr << "HashImage<VoxelType>::FlipZT: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
/*template <class VoxelType>
bool HashImage<VoxelType>::CropPad(int margin)
{
    cerr << "HashImage<VoxelType>::CropPad: Not implemented" << endl;
    exit(1);
}*/

// =============================================================================
// VTK interface
// =============================================================================
#if MIRTK_Image_WITH_VTK

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::ImageToVTK(vtkStructuredPoints *vtk) const
{
    if (this->ImageToVTKScalarType() == VTK_VOID) {
        cerr << "HashImage::ImageToVTK: Cannot convert image to VTK structured points" << endl;
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
    VoxelType       *ptr2 = reinterpret_cast<VoxelType *>(vtk->GetScalarPointer());
    for (int i = 0; i < nvox; ++i) {
        for (int l = 0; l < _attr._t; ++l, ++ptr2) *ptr2 = Get(l * nvox);
    }
}
template <>
void HashImage<Matrix3x3>::ImageToVTK(vtkStructuredPoints *) const
{
    cerr << "HashImage<Matrix3x3>::VTKToImage: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class Type>
void HashImage<Type>::VTKToImage(vtkStructuredPoints *)
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
void HashImage<VoxelType>::Read(const char *fname)
{
    // Read image
    GenericImage<VoxelType> img(fname);
    Initialize(img.Attributes());
    CopyFrom(img);
}

template <> void HashImage<float3x3>::Read(const char *)
{
    cerr << "HashImage<float3x3>::Read: Not implemented" << endl;
    exit(1);
}

template <> void HashImage<double3x3>::Read(const char *)
{
    cerr << "HashImage<double3x3>::Read: Not implemented" << endl;
    exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
void HashImage<VoxelType>::Write(const char *fname) const
{
    string name(fname);
    if (Extension(fname).empty()) name += MIRTK_Image_DEFAULT_EXT;
    unique_ptr<ImageWriter> writer(ImageWriter::New(name.c_str()));
    GenericImage<VoxelType> img=this->ToGenericImage();
    writer->Input(&img);
    writer->Run();
}

template <> void HashImage<float3x3>::Write(const char *) const
{
    cerr << "HashImage<float3x3>::Write: Not implemented" << endl;
    exit(1);
}

template <> void HashImage<double3x3>::Write(const char *) const
{
    cerr << "HashImage<double3x3>::Read: Not implemented" << endl;
    exit(1);
}


template <class VoxelType>
GenericImage<VoxelType> HashImage<VoxelType>::ToGenericImage() const
{
    GenericImage<VoxelType> dense_image;
    this->CopyTo(dense_image);
    return dense_image;
}


template <class VoxelType> template <class VoxelType2>
void HashImage<VoxelType>::ToGenericImage(GenericImage<VoxelType2>& image) const
{
    this->CopyTo(image);
}



// =============================================================================
// Explicit template instantiations
// =============================================================================

template class HashImage<char>;
template class HashImage<unsigned char>;
template class HashImage<short>;
template class HashImage<unsigned short>;
template class HashImage<int>;
template class HashImage<unsigned int>;
template class HashImage<float>;
template class HashImage<float2>;
template class HashImage<float3>;
template class HashImage<float4>;
template class HashImage<double>;
template class HashImage<double2>;
template class HashImage<double3>;
template class HashImage<double4>;
template class HashImage<float3x3>;
template class HashImage<double3x3>;

template HashImage<char>::HashImage(const HashImage<unsigned char> &);
template HashImage<char>::HashImage(const HashImage<short> &);
template HashImage<char>::HashImage(const HashImage<unsigned short> &);
template HashImage<char>::HashImage(const HashImage<int> &);
template HashImage<char>::HashImage(const HashImage<float> &);
template HashImage<char>::HashImage(const HashImage<double> &);

template HashImage<unsigned char>::HashImage(const HashImage<char> &);
template HashImage<unsigned char>::HashImage(const HashImage<short> &);
template HashImage<unsigned char>::HashImage(const HashImage<unsigned short> &);
template HashImage<unsigned char>::HashImage(const HashImage<int> &);
template HashImage<unsigned char>::HashImage(const HashImage<float> &);
template HashImage<unsigned char>::HashImage(const HashImage<double> &);

template HashImage<short>::HashImage(const HashImage<char> &);
template HashImage<short>::HashImage(const HashImage<unsigned char> &);
template HashImage<short>::HashImage(const HashImage<unsigned short> &);
template HashImage<short>::HashImage(const HashImage<int> &);
template HashImage<short>::HashImage(const HashImage<float> &);
template HashImage<short>::HashImage(const HashImage<double> &);

template HashImage<unsigned short>::HashImage(const HashImage<char> &);
template HashImage<unsigned short>::HashImage(const HashImage<unsigned char> &);
template HashImage<unsigned short>::HashImage(const HashImage<short> &);
template HashImage<unsigned short>::HashImage(const HashImage<int> &);
template HashImage<unsigned short>::HashImage(const HashImage<float> &);
template HashImage<unsigned short>::HashImage(const HashImage<double> &);

template HashImage<int>::HashImage(const HashImage<char> &);
template HashImage<int>::HashImage(const HashImage<unsigned char> &);
template HashImage<int>::HashImage(const HashImage<short> &);
template HashImage<int>::HashImage(const HashImage<unsigned short> &);
template HashImage<int>::HashImage(const HashImage<float> &);
template HashImage<int>::HashImage(const HashImage<double> &);

template HashImage<float>::HashImage(const HashImage<char> &);
template HashImage<float>::HashImage(const HashImage<unsigned char> &);
template HashImage<float>::HashImage(const HashImage<short> &);
template HashImage<float>::HashImage(const HashImage<unsigned short> &);
template HashImage<float>::HashImage(const HashImage<int> &);
template HashImage<float>::HashImage(const HashImage<double> &);

template HashImage<double>::HashImage(const HashImage<char> &);
template HashImage<double>::HashImage(const HashImage<unsigned char> &);
template HashImage<double>::HashImage(const HashImage<short> &);
template HashImage<double>::HashImage(const HashImage<unsigned short> &);
template HashImage<double>::HashImage(const HashImage<int> &);
template HashImage<double>::HashImage(const HashImage<float> &);

template HashImage<char>::HashImage(const GenericImage<char> &);
template HashImage<char>::HashImage(const GenericImage<unsigned char> &);
template HashImage<char>::HashImage(const GenericImage<short> &);
template HashImage<char>::HashImage(const GenericImage<unsigned short> &);
template HashImage<char>::HashImage(const GenericImage<int> &);
template HashImage<char>::HashImage(const GenericImage<float> &);
template HashImage<char>::HashImage(const GenericImage<double> &);

template HashImage<unsigned char>::HashImage(const GenericImage<char> &);
template HashImage<unsigned char>::HashImage(const GenericImage<unsigned char> &);
template HashImage<unsigned char>::HashImage(const GenericImage<short> &);
template HashImage<unsigned char>::HashImage(const GenericImage<unsigned short> &);
template HashImage<unsigned char>::HashImage(const GenericImage<int> &);
template HashImage<unsigned char>::HashImage(const GenericImage<float> &);
template HashImage<unsigned char>::HashImage(const GenericImage<double> &);

template HashImage<short>::HashImage(const GenericImage<char> &);
template HashImage<short>::HashImage(const GenericImage<unsigned char> &);
template HashImage<short>::HashImage(const GenericImage<short> &);
template HashImage<short>::HashImage(const GenericImage<unsigned short> &);
template HashImage<short>::HashImage(const GenericImage<int> &);
template HashImage<short>::HashImage(const GenericImage<float> &);
template HashImage<short>::HashImage(const GenericImage<double> &);

template HashImage<unsigned short>::HashImage(const GenericImage<char> &);
template HashImage<unsigned short>::HashImage(const GenericImage<unsigned char> &);
template HashImage<unsigned short>::HashImage(const GenericImage<short> &);
template HashImage<unsigned short>::HashImage(const GenericImage<unsigned short> &);
template HashImage<unsigned short>::HashImage(const GenericImage<int> &);
template HashImage<unsigned short>::HashImage(const GenericImage<float> &);
template HashImage<unsigned short>::HashImage(const GenericImage<double> &);

template HashImage<int>::HashImage(const GenericImage<char> &);
template HashImage<int>::HashImage(const GenericImage<unsigned char> &);
template HashImage<int>::HashImage(const GenericImage<short> &);
template HashImage<int>::HashImage(const GenericImage<unsigned short> &);
template HashImage<int>::HashImage(const GenericImage<int> &);
template HashImage<int>::HashImage(const GenericImage<float> &);
template HashImage<int>::HashImage(const GenericImage<double> &);

template HashImage<float>::HashImage(const GenericImage<char> &);
template HashImage<float>::HashImage(const GenericImage<unsigned char> &);
template HashImage<float>::HashImage(const GenericImage<short> &);
template HashImage<float>::HashImage(const GenericImage<unsigned short> &);
template HashImage<float>::HashImage(const GenericImage<int> &);
template HashImage<float>::HashImage(const GenericImage<float> &);
template HashImage<float>::HashImage(const GenericImage<double> &);

template HashImage<double>::HashImage(const GenericImage<char> &);
template HashImage<double>::HashImage(const GenericImage<unsigned char> &);
template HashImage<double>::HashImage(const GenericImage<short> &);
template HashImage<double>::HashImage(const GenericImage<unsigned short> &);
template HashImage<double>::HashImage(const GenericImage<int> &);
template HashImage<double>::HashImage(const GenericImage<float> &);
template HashImage<double>::HashImage(const GenericImage<double> &);

template HashImage<float2  >::HashImage(const HashImage<double2  > &);
template HashImage<float3  >::HashImage(const HashImage<double3  > &);
template HashImage<float4  >::HashImage(const HashImage<double4  > &);
template HashImage<float3x3  >::HashImage(const HashImage<double3x3  > &);

template HashImage<double2  >::HashImage(const HashImage<float2  > &);
template HashImage<double3  >::HashImage(const HashImage<float3  > &);
template HashImage<double4  >::HashImage(const HashImage<float4  > &);
template HashImage<double3x3  >::HashImage(const HashImage<float3x3  > &);

template void HashImage<char>::CopyFrom(const HashImage<char> &image);
template void HashImage<char>::CopyFrom(const HashImage<unsigned char> &image);
template void HashImage<char>::CopyFrom(const HashImage<short> &image);
template void HashImage<char>::CopyFrom(const HashImage<unsigned short> &image);
template void HashImage<char>::CopyFrom(const HashImage<int> &image);
template void HashImage<char>::CopyFrom(const HashImage<float> &image);
template void HashImage<char>::CopyFrom(const HashImage<double> &image);

template void HashImage<unsigned char>::CopyFrom(const HashImage<char> &image);
template void HashImage<unsigned char>::CopyFrom(const HashImage<unsigned char> &image);
template void HashImage<unsigned char>::CopyFrom(const HashImage<short> &image);
template void HashImage<unsigned char>::CopyFrom(const HashImage<unsigned short> &image);
template void HashImage<unsigned char>::CopyFrom(const HashImage<int> &image);
template void HashImage<unsigned char>::CopyFrom(const HashImage<float> &image);
template void HashImage<unsigned char>::CopyFrom(const HashImage<double> &image);

template void HashImage<short>::CopyFrom(const HashImage<char> &image);
template void HashImage<short>::CopyFrom(const HashImage<unsigned char> &image);
template void HashImage<short>::CopyFrom(const HashImage<short> &image);
template void HashImage<short>::CopyFrom(const HashImage<unsigned short> &image);
template void HashImage<short>::CopyFrom(const HashImage<int> &image);
template void HashImage<short>::CopyFrom(const HashImage<float> &image);
template void HashImage<short>::CopyFrom(const HashImage<double> &image);

template void HashImage<unsigned short>::CopyFrom(const HashImage<char> &image);
template void HashImage<unsigned short>::CopyFrom(const HashImage<unsigned char> &image);
template void HashImage<unsigned short>::CopyFrom(const HashImage<short> &image);
template void HashImage<unsigned short>::CopyFrom(const HashImage<unsigned short> &image);
template void HashImage<unsigned short>::CopyFrom(const HashImage<int> &image);
template void HashImage<unsigned short>::CopyFrom(const HashImage<float> &image);
template void HashImage<unsigned short>::CopyFrom(const HashImage<double> &image);

template void HashImage<int>::CopyFrom(const HashImage<char> &image);
template void HashImage<int>::CopyFrom(const HashImage<unsigned char> &image);
template void HashImage<int>::CopyFrom(const HashImage<short> &image);
template void HashImage<int>::CopyFrom(const HashImage<unsigned short> &image);
template void HashImage<int>::CopyFrom(const HashImage<int> &image);
template void HashImage<int>::CopyFrom(const HashImage<float> &image);
template void HashImage<int>::CopyFrom(const HashImage<double> &image);

template void HashImage<float>::CopyFrom(const HashImage<char> &image);
template void HashImage<float>::CopyFrom(const HashImage<unsigned char> &image);
template void HashImage<float>::CopyFrom(const HashImage<short> &image);
template void HashImage<float>::CopyFrom(const HashImage<unsigned short> &image);
template void HashImage<float>::CopyFrom(const HashImage<int> &image);
template void HashImage<float>::CopyFrom(const HashImage<float> &image);
template void HashImage<float>::CopyFrom(const HashImage<double> &image);

template void HashImage<double>::CopyFrom(const HashImage<char> &image);
template void HashImage<double>::CopyFrom(const HashImage<unsigned char> &image);
template void HashImage<double>::CopyFrom(const HashImage<short> &image);
template void HashImage<double>::CopyFrom(const HashImage<unsigned short> &image);
template void HashImage<double>::CopyFrom(const HashImage<int> &image);
template void HashImage<double>::CopyFrom(const HashImage<float> &image);
template void HashImage<double>::CopyFrom(const HashImage<double> &image);

template void HashImage<float2>::CopyFrom(const HashImage<double2  > &);
template void HashImage<float3>::CopyFrom(const HashImage<double3  > &);
template void HashImage<float4>::CopyFrom(const HashImage<double4  > &);
template void HashImage<float3x3>::CopyFrom(const HashImage<double3x3  > &);

template void HashImage<double2>::CopyFrom(const HashImage<float2  > &);
template void HashImage<double3>::CopyFrom(const HashImage<float3  > &);
template void HashImage<double4>::CopyFrom(const HashImage<float4  > &);
template void HashImage<double3x3>::CopyFrom(const HashImage<float3x3  > &);

template void HashImage<char>::CopyFrom(const GenericImage<char> &image);
template void HashImage<char>::CopyFrom(const GenericImage<unsigned char> &image);
template void HashImage<char>::CopyFrom(const GenericImage<short> &image);
template void HashImage<char>::CopyFrom(const GenericImage<unsigned short> &image);
template void HashImage<char>::CopyFrom(const GenericImage<int> &image);
template void HashImage<char>::CopyFrom(const GenericImage<float> &image);
template void HashImage<char>::CopyFrom(const GenericImage<double> &image);

template void HashImage<unsigned char>::CopyFrom(const GenericImage<char> &image);
template void HashImage<unsigned char>::CopyFrom(const GenericImage<unsigned char> &image);
template void HashImage<unsigned char>::CopyFrom(const GenericImage<short> &image);
template void HashImage<unsigned char>::CopyFrom(const GenericImage<unsigned short> &image);
template void HashImage<unsigned char>::CopyFrom(const GenericImage<int> &image);
template void HashImage<unsigned char>::CopyFrom(const GenericImage<float> &image);
template void HashImage<unsigned char>::CopyFrom(const GenericImage<double> &image);

template void HashImage<short>::CopyFrom(const GenericImage<char> &image);
template void HashImage<short>::CopyFrom(const GenericImage<unsigned char> &image);
template void HashImage<short>::CopyFrom(const GenericImage<short> &image);
template void HashImage<short>::CopyFrom(const GenericImage<unsigned short> &image);
template void HashImage<short>::CopyFrom(const GenericImage<int> &image);
template void HashImage<short>::CopyFrom(const GenericImage<float> &image);
template void HashImage<short>::CopyFrom(const GenericImage<double> &image);

template void HashImage<unsigned short>::CopyFrom(const GenericImage<char> &image);
template void HashImage<unsigned short>::CopyFrom(const GenericImage<unsigned char> &image);
template void HashImage<unsigned short>::CopyFrom(const GenericImage<short> &image);
template void HashImage<unsigned short>::CopyFrom(const GenericImage<unsigned short> &image);
template void HashImage<unsigned short>::CopyFrom(const GenericImage<int> &image);
template void HashImage<unsigned short>::CopyFrom(const GenericImage<float> &image);
template void HashImage<unsigned short>::CopyFrom(const GenericImage<double> &image);

template void HashImage<int>::CopyFrom(const GenericImage<char> &image);
template void HashImage<int>::CopyFrom(const GenericImage<unsigned char> &image);
template void HashImage<int>::CopyFrom(const GenericImage<short> &image);
template void HashImage<int>::CopyFrom(const GenericImage<unsigned short> &image);
template void HashImage<int>::CopyFrom(const GenericImage<int> &image);
template void HashImage<int>::CopyFrom(const GenericImage<float> &image);
template void HashImage<int>::CopyFrom(const GenericImage<double> &image);

template void HashImage<float>::CopyFrom(const GenericImage<char> &image);
template void HashImage<float>::CopyFrom(const GenericImage<unsigned char> &image);
template void HashImage<float>::CopyFrom(const GenericImage<short> &image);
template void HashImage<float>::CopyFrom(const GenericImage<unsigned short> &image);
template void HashImage<float>::CopyFrom(const GenericImage<int> &image);
template void HashImage<float>::CopyFrom(const GenericImage<float> &image);
template void HashImage<float>::CopyFrom(const GenericImage<double> &image);

template void HashImage<double>::CopyFrom(const GenericImage<char> &image);
template void HashImage<double>::CopyFrom(const GenericImage<unsigned char> &image);
template void HashImage<double>::CopyFrom(const GenericImage<short> &image);
template void HashImage<double>::CopyFrom(const GenericImage<unsigned short> &image);
template void HashImage<double>::CopyFrom(const GenericImage<int> &image);
template void HashImage<double>::CopyFrom(const GenericImage<float> &image);
template void HashImage<double>::CopyFrom(const GenericImage<double> &image);

template void HashImage<char>::CopyTo(GenericImage<char> &image) const;
template void HashImage<char>::CopyTo(GenericImage<unsigned char> &image) const;
template void HashImage<char>::CopyTo(GenericImage<short> &image) const;
template void HashImage<char>::CopyTo(GenericImage<unsigned short> &image) const;
template void HashImage<char>::CopyTo(GenericImage<int> &image) const;
template void HashImage<char>::CopyTo(GenericImage<float> &image) const;
template void HashImage<char>::CopyTo(GenericImage<double> &image) const;

template void HashImage<unsigned char>::CopyTo(GenericImage<char> &image) const;
template void HashImage<unsigned char>::CopyTo(GenericImage<unsigned char> &image) const;
template void HashImage<unsigned char>::CopyTo(GenericImage<short> &image) const;
template void HashImage<unsigned char>::CopyTo(GenericImage<unsigned short> &image) const;
template void HashImage<unsigned char>::CopyTo(GenericImage<int> &image) const;
template void HashImage<unsigned char>::CopyTo(GenericImage<float> &image) const;
template void HashImage<unsigned char>::CopyTo(GenericImage<double> &image) const;

template void HashImage<short>::CopyTo(GenericImage<char> &image) const;
template void HashImage<short>::CopyTo(GenericImage<unsigned char> &image) const;
template void HashImage<short>::CopyTo(GenericImage<short> &image) const;
template void HashImage<short>::CopyTo(GenericImage<unsigned short> &image) const;
template void HashImage<short>::CopyTo(GenericImage<int> &image) const;
template void HashImage<short>::CopyTo(GenericImage<float> &image) const;
template void HashImage<short>::CopyTo(GenericImage<double> &image) const;

template void HashImage<unsigned short>::CopyTo(GenericImage<char> &image) const;
template void HashImage<unsigned short>::CopyTo(GenericImage<unsigned char> &image) const;
template void HashImage<unsigned short>::CopyTo(GenericImage<short> &image) const;
template void HashImage<unsigned short>::CopyTo(GenericImage<unsigned short> &image) const;
template void HashImage<unsigned short>::CopyTo(GenericImage<int> &image) const;
template void HashImage<unsigned short>::CopyTo(GenericImage<float> &image) const;
template void HashImage<unsigned short>::CopyTo(GenericImage<double> &image) const;

template void HashImage<int>::CopyTo(GenericImage<char> &image) const;
template void HashImage<int>::CopyTo(GenericImage<unsigned char> &image) const;
template void HashImage<int>::CopyTo(GenericImage<short> &image) const;
template void HashImage<int>::CopyTo(GenericImage<unsigned short> &image) const;
template void HashImage<int>::CopyTo(GenericImage<int> &image) const;
template void HashImage<int>::CopyTo(GenericImage<float> &image) const;
template void HashImage<int>::CopyTo(GenericImage<double> &image) const;

template void HashImage<float>::CopyTo(GenericImage<char> &image) const;
template void HashImage<float>::CopyTo(GenericImage<unsigned char> &image) const;
template void HashImage<float>::CopyTo(GenericImage<short> &image) const;
template void HashImage<float>::CopyTo(GenericImage<unsigned short> &image) const;
template void HashImage<float>::CopyTo(GenericImage<int> &image) const;
template void HashImage<float>::CopyTo(GenericImage<float> &image) const;
template void HashImage<float>::CopyTo(GenericImage<double> &image) const;

template void HashImage<double>::CopyTo(GenericImage<char> &image) const;
template void HashImage<double>::CopyTo(GenericImage<unsigned char> &image) const;
template void HashImage<double>::CopyTo(GenericImage<short> &image) const;
template void HashImage<double>::CopyTo(GenericImage<unsigned short> &image) const;
template void HashImage<double>::CopyTo(GenericImage<int> &image) const;
template void HashImage<double>::CopyTo(GenericImage<float> &image) const;
template void HashImage<double>::CopyTo(GenericImage<double> &image) const;

template void HashImage<char>::ToGenericImage(GenericImage<char>& image) const;
template void HashImage<char>::ToGenericImage(GenericImage<unsigned char>& image) const;
template void HashImage<char>::ToGenericImage(GenericImage<short>& image) const;
template void HashImage<char>::ToGenericImage(GenericImage<unsigned short>& image) const;
template void HashImage<char>::ToGenericImage(GenericImage<int>& image) const;
template void HashImage<char>::ToGenericImage(GenericImage<float>& image) const;
template void HashImage<char>::ToGenericImage(GenericImage<double>& image) const;

template void HashImage<unsigned char>::ToGenericImage(GenericImage<char>& image) const;
template void HashImage<unsigned char>::ToGenericImage(GenericImage<unsigned char>& image) const;
template void HashImage<unsigned char>::ToGenericImage(GenericImage<short>& image) const;
template void HashImage<unsigned char>::ToGenericImage(GenericImage<unsigned short>& image) const;
template void HashImage<unsigned char>::ToGenericImage(GenericImage<int>& image) const;
template void HashImage<unsigned char>::ToGenericImage(GenericImage<float>& image) const;
template void HashImage<unsigned char>::ToGenericImage(GenericImage<double>& image) const;

template void HashImage<short>::ToGenericImage(GenericImage<char>& image) const;
template void HashImage<short>::ToGenericImage(GenericImage<unsigned char>& image) const;
template void HashImage<short>::ToGenericImage(GenericImage<short>& image) const;
template void HashImage<short>::ToGenericImage(GenericImage<unsigned short>& image) const;
template void HashImage<short>::ToGenericImage(GenericImage<int>& image) const;
template void HashImage<short>::ToGenericImage(GenericImage<float>& image) const;
template void HashImage<short>::ToGenericImage(GenericImage<double>& image) const;

template void HashImage<unsigned short>::ToGenericImage(GenericImage<char>& image) const;
template void HashImage<unsigned short>::ToGenericImage(GenericImage<unsigned char>& image) const;
template void HashImage<unsigned short>::ToGenericImage(GenericImage<short>& image) const;
template void HashImage<unsigned short>::ToGenericImage(GenericImage<unsigned short>& image) const;
template void HashImage<unsigned short>::ToGenericImage(GenericImage<int>& image) const;
template void HashImage<unsigned short>::ToGenericImage(GenericImage<float>& image) const;
template void HashImage<unsigned short>::ToGenericImage(GenericImage<double>& image) const;

template void HashImage<int>::ToGenericImage(GenericImage<char>& image) const;
template void HashImage<int>::ToGenericImage(GenericImage<unsigned char>& image) const;
template void HashImage<int>::ToGenericImage(GenericImage<short>& image) const;
template void HashImage<int>::ToGenericImage(GenericImage<unsigned short>& image) const;
template void HashImage<int>::ToGenericImage(GenericImage<int>& image) const;
template void HashImage<int>::ToGenericImage(GenericImage<float>& image) const;
template void HashImage<int>::ToGenericImage(GenericImage<double>& image) const;

template void HashImage<float>::ToGenericImage(GenericImage<char>& image) const;
template void HashImage<float>::ToGenericImage(GenericImage<unsigned char>& image) const;
template void HashImage<float>::ToGenericImage(GenericImage<short>& image) const;
template void HashImage<float>::ToGenericImage(GenericImage<unsigned short>& image) const;
template void HashImage<float>::ToGenericImage(GenericImage<int>& image) const;
template void HashImage<float>::ToGenericImage(GenericImage<float>& image) const;
template void HashImage<float>::ToGenericImage(GenericImage<double>& image) const;

template void HashImage<double>::ToGenericImage(GenericImage<char>& image) const;
template void HashImage<double>::ToGenericImage(GenericImage<unsigned char>& image) const;
template void HashImage<double>::ToGenericImage(GenericImage<short>& image) const;
template void HashImage<double>::ToGenericImage(GenericImage<unsigned short>& image) const;
template void HashImage<double>::ToGenericImage(GenericImage<int>& image) const;
template void HashImage<double>::ToGenericImage(GenericImage<float>& image) const;
template void HashImage<double>::ToGenericImage(GenericImage<double>& image) const;



template HashImage<char>& HashImage<char>::operator= (const HashImage<unsigned char> &);
template HashImage<char>& HashImage<char>::operator= (const HashImage<short> &);
template HashImage<char>& HashImage<char>::operator= (const HashImage<unsigned short> &);
template HashImage<char>& HashImage<char>::operator= (const HashImage<int> &);
template HashImage<char>& HashImage<char>::operator= (const HashImage<float> &);
template HashImage<char>& HashImage<char>::operator= (const HashImage<double> &);

template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const HashImage<char> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const HashImage<short> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const HashImage<unsigned short> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const HashImage<int> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const HashImage<float> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const HashImage<double> &);

template HashImage<short>& HashImage<short>::operator= (const HashImage<char> &);
template HashImage<short>& HashImage<short>::operator= (const HashImage<unsigned char> &);
template HashImage<short>& HashImage<short>::operator= (const HashImage<unsigned short> &);
template HashImage<short>& HashImage<short>::operator= (const HashImage<int> &);
template HashImage<short>& HashImage<short>::operator= (const HashImage<float> &);
template HashImage<short>& HashImage<short>::operator= (const HashImage<double> &);

template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const HashImage<char> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const HashImage<unsigned char> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const HashImage<short> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const HashImage<int> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const HashImage<float> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const HashImage<double> &);

template HashImage<int>& HashImage<int>::operator= (const HashImage<char> &);
template HashImage<int>& HashImage<int>::operator= (const HashImage<unsigned char> &);
template HashImage<int>& HashImage<int>::operator= (const HashImage<short> &);
template HashImage<int>& HashImage<int>::operator= (const HashImage<unsigned short> &);
template HashImage<int>& HashImage<int>::operator= (const HashImage<float> &);
template HashImage<int>& HashImage<int>::operator= (const HashImage<double> &);

template HashImage<float>& HashImage<float>::operator= (const HashImage<char> &);
template HashImage<float>& HashImage<float>::operator= (const HashImage<unsigned char> &);
template HashImage<float>& HashImage<float>::operator= (const HashImage<short> &);
template HashImage<float>& HashImage<float>::operator= (const HashImage<unsigned short> &);
template HashImage<float>& HashImage<float>::operator= (const HashImage<int> &);
template HashImage<float>& HashImage<float>::operator= (const HashImage<double> &);

template HashImage<double>& HashImage<double>::operator= (const HashImage<char> &);
template HashImage<double>& HashImage<double>::operator= (const HashImage<unsigned char> &);
template HashImage<double>& HashImage<double>::operator= (const HashImage<short> &);
template HashImage<double>& HashImage<double>::operator= (const HashImage<unsigned short> &);
template HashImage<double>& HashImage<double>::operator= (const HashImage<int> &);
template HashImage<double>& HashImage<double>::operator= (const HashImage<float> &);

template HashImage<char>& HashImage<char>::operator= (const GenericImage<char> &);
template HashImage<char>& HashImage<char>::operator= (const GenericImage<unsigned char> &);
template HashImage<char>& HashImage<char>::operator= (const GenericImage<short> &);
template HashImage<char>& HashImage<char>::operator= (const GenericImage<unsigned short> &);
template HashImage<char>& HashImage<char>::operator= (const GenericImage<int> &);
template HashImage<char>& HashImage<char>::operator= (const GenericImage<float> &);
template HashImage<char>& HashImage<char>::operator= (const GenericImage<double> &);

template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const GenericImage<char> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const GenericImage<unsigned char> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const GenericImage<short> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const GenericImage<unsigned short> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const GenericImage<int> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const GenericImage<float> &);
template HashImage<unsigned char>& HashImage<unsigned char>::operator= (const GenericImage<double> &);

template HashImage<short>& HashImage<short>::operator= (const GenericImage<char> &);
template HashImage<short>& HashImage<short>::operator= (const GenericImage<unsigned char> &);
template HashImage<short>& HashImage<short>::operator= (const GenericImage<short> &);
template HashImage<short>& HashImage<short>::operator= (const GenericImage<unsigned short> &);
template HashImage<short>& HashImage<short>::operator= (const GenericImage<int> &);
template HashImage<short>& HashImage<short>::operator= (const GenericImage<float> &);
template HashImage<short>& HashImage<short>::operator= (const GenericImage<double> &);

template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const GenericImage<char> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const GenericImage<unsigned char> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const GenericImage<short> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const GenericImage<unsigned short> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const GenericImage<int> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const GenericImage<float> &);
template HashImage<unsigned short>& HashImage<unsigned short>::operator= (const GenericImage<double> &);

template HashImage<int>& HashImage<int>::operator= (const GenericImage<char> &);
template HashImage<int>& HashImage<int>::operator= (const GenericImage<unsigned char> &);
template HashImage<int>& HashImage<int>::operator= (const GenericImage<short> &);
template HashImage<int>& HashImage<int>::operator= (const GenericImage<unsigned short> &);
template HashImage<int>& HashImage<int>::operator= (const GenericImage<int> &);
template HashImage<int>& HashImage<int>::operator= (const GenericImage<float> &);
template HashImage<int>& HashImage<int>::operator= (const GenericImage<double> &);

template HashImage<float>& HashImage<float>::operator= (const GenericImage<char> &);
template HashImage<float>& HashImage<float>::operator= (const GenericImage<unsigned char> &);
template HashImage<float>& HashImage<float>::operator= (const GenericImage<short> &);
template HashImage<float>& HashImage<float>::operator= (const GenericImage<unsigned short> &);
template HashImage<float>& HashImage<float>::operator= (const GenericImage<int> &);
template HashImage<float>& HashImage<float>::operator= (const GenericImage<float> &);
template HashImage<float>& HashImage<float>::operator= (const GenericImage<double> &);

template HashImage<double>& HashImage<double>::operator= (const GenericImage<char> &);
template HashImage<double>& HashImage<double>::operator= (const GenericImage<unsigned char> &);
template HashImage<double>& HashImage<double>::operator= (const GenericImage<short> &);
template HashImage<double>& HashImage<double>::operator= (const GenericImage<unsigned short> &);
template HashImage<double>& HashImage<double>::operator= (const GenericImage<int> &);
template HashImage<double>& HashImage<double>::operator= (const GenericImage<float> &);
template HashImage<double>& HashImage<double>::operator= (const GenericImage<double> &);

template bool HashImage<char>::operator==(const HashImage<char> &) const;
template bool HashImage<char>::operator==(const HashImage<unsigned char> &) const;
template bool HashImage<char>::operator==(const HashImage<short> &) const;
template bool HashImage<char>::operator==(const HashImage<unsigned short> &) const;
template bool HashImage<char>::operator==(const HashImage<int> &) const;
template bool HashImage<char>::operator==(const HashImage<float> &) const;
template bool HashImage<char>::operator==(const HashImage<double> &) const;

template bool HashImage<unsigned char>::operator==(const HashImage<char> &) const;
template bool HashImage<unsigned char>::operator==(const HashImage<unsigned char> &) const;
template bool HashImage<unsigned char>::operator==(const HashImage<short> &) const;
template bool HashImage<unsigned char>::operator==(const HashImage<unsigned short> &) const;
template bool HashImage<unsigned char>::operator==(const HashImage<int> &) const;
template bool HashImage<unsigned char>::operator==(const HashImage<float> &) const;
template bool HashImage<unsigned char>::operator==(const HashImage<double> &) const;

template bool HashImage<short>::operator==(const HashImage<char> &) const;
template bool HashImage<short>::operator==(const HashImage<unsigned char> &) const;
template bool HashImage<short>::operator==(const HashImage<short> &) const;
template bool HashImage<short>::operator==(const HashImage<unsigned short> &) const;
template bool HashImage<short>::operator==(const HashImage<int> &) const;
template bool HashImage<short>::operator==(const HashImage<float> &) const;
template bool HashImage<short>::operator==(const HashImage<double> &) const;

template bool HashImage<unsigned short>::operator==(const HashImage<char> &) const;
template bool HashImage<unsigned short>::operator==(const HashImage<unsigned char> &) const;
template bool HashImage<unsigned short>::operator==(const HashImage<short> &) const;
template bool HashImage<unsigned short>::operator==(const HashImage<unsigned short> &) const;
template bool HashImage<unsigned short>::operator==(const HashImage<int> &) const;
template bool HashImage<unsigned short>::operator==(const HashImage<float> &) const;
template bool HashImage<unsigned short>::operator==(const HashImage<double> &) const;

template bool HashImage<int>::operator==(const HashImage<char> &) const;
template bool HashImage<int>::operator==(const HashImage<unsigned char> &) const;
template bool HashImage<int>::operator==(const HashImage<short> &) const;
template bool HashImage<int>::operator==(const HashImage<unsigned short> &) const;
template bool HashImage<int>::operator==(const HashImage<int> &) const;
template bool HashImage<int>::operator==(const HashImage<float> &) const;
template bool HashImage<int>::operator==(const HashImage<double> &) const;

template bool HashImage<float>::operator==(const HashImage<char> &) const;
template bool HashImage<float>::operator==(const HashImage<unsigned char> &) const;
template bool HashImage<float>::operator==(const HashImage<short> &) const;
template bool HashImage<float>::operator==(const HashImage<unsigned short> &) const;
template bool HashImage<float>::operator==(const HashImage<int> &) const;
template bool HashImage<float>::operator==(const HashImage<float> &) const;
template bool HashImage<float>::operator==(const HashImage<double> &) const;

template bool HashImage<double>::operator==(const HashImage<char> &) const;
template bool HashImage<double>::operator==(const HashImage<unsigned char> &) const;
template bool HashImage<double>::operator==(const HashImage<short> &) const;
template bool HashImage<double>::operator==(const HashImage<unsigned short> &) const;
template bool HashImage<double>::operator==(const HashImage<int> &) const;
template bool HashImage<double>::operator==(const HashImage<float> &) const;
template bool HashImage<double>::operator==(const HashImage<double> &) const;


// TODO: Remove deprecated template instantiations below
template class HashImage<Vector3D<float>  >;
template class HashImage<Vector3D<double> >;
template HashImage<Vector3D<double> >::HashImage(const HashImage<Vector3D<float> > &);


} // namespace mirtk
