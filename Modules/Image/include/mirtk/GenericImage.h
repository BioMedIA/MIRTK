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

#include "mirtk/BaseImage.h" // MUST be before include guard because of
                            // cyclic dependency between BaseImage and
                            // GenericImage

#ifndef MIRTK_GenericImage_H
#define MIRTK_GenericImage_H

#include "mirtk/VoxelCast.h"


namespace mirtk {


/**
 * Generic class for 2D or 3D images
 *
 * This class implements generic 2D and 3D images. It provides functions
 * for accessing, reading, writing and manipulating images. This class can
 * be used for images with arbitrary voxel types using templates.
 */
template <class TVoxel>
class GenericImage : public BaseImage
{
  mirtkObjectMacro(GenericImage);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Voxel type
  typedef TVoxel VoxelType;

  /// Floating point type corresponding to voxel type
  /// \note The VoxelType as well as the RealType may be a matrix/vector type!
  typedef typename voxel_info<VoxelType>::RealType RealType;

  /// Scalar type corresponding to voxel type
  typedef typename voxel_info<VoxelType>::ScalarType ScalarType;

  /// Floating point type corresponding to scalar type of voxel type
  typedef typename voxel_info<ScalarType>::RealType RealScalarType;

  // ---------------------------------------------------------------------------
  // Data members

protected:

  /// Pointer array for access to image data
  ///
  /// \note The image data is stored in a contiguous memory block which can
  ///       be alternatively referred to as 1D data array using \c _data.
  VoxelType ****_matrix;

  /// Pointer to image data
  VoxelType *_data;

  /// Whether image data memory itself is owned by this instance
  bool _dataOwner;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Allocate image memory
  void AllocateImage(VoxelType * = NULL);

public:

  /// Default constructor
  GenericImage();

  /// Constructor from image file
  GenericImage(const char *);

  /// Constructor for given image size
  explicit GenericImage(int, int, int = 1, int = 1, VoxelType *data = NULL);

  /// Constructor for given image size
  explicit GenericImage(int, int, int, int, int, VoxelType *data = NULL);

  /// Constructor for given image attributes
  explicit GenericImage(const ImageAttributes &, VoxelType *data = NULL);

  /// Constructor for given image attributes
  explicit GenericImage(const ImageAttributes &, int, VoxelType *data = NULL);

  /// Copy constructor for image
  explicit GenericImage(const BaseImage &);

  /// Copy constructor for image
  GenericImage(const GenericImage &);

  /// Copy constructor for image of different type
  template <class TVoxel2>
  GenericImage(const GenericImage<TVoxel2> &);

  /// Destructor
  ~GenericImage();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Create copy of this image
  virtual BaseImage *Copy() const;

  /// Initialize a previously allocated image
  virtual void Initialize();

  /// Initialize an image
  virtual void Initialize(const ImageAttributes &, int, VoxelType *data);

  /// Initialize an image
  void Initialize(const ImageAttributes &, int);

  /// Initialize an image
  void Initialize(const ImageAttributes &, VoxelType *data = NULL);

  /// Initialize an image
  void Initialize(int, int, int, int, int, VoxelType *data = NULL);

  /// Initialize an image
  void Initialize(int, int, int = 1, int = 1, VoxelType *data = NULL);

  /// Copy image data from 1D array
  void CopyFrom(const VoxelType *);

  /// Copy image data from other image of same size
  void CopyFrom(const BaseImage &);

  /// Copy image data from other image of same size
  void CopyFrom(const GenericImage &);

  /// Assign constant value to each voxel
  GenericImage& operator= (VoxelType);

  /// Assignment operator with implicit cast to double and then VoxelType
  GenericImage<VoxelType>& operator= (const BaseImage &);

  /// Assignment operator
  GenericImage<VoxelType>& operator= (const GenericImage &);

  /// Assignment operator with implicit cast
  template <class TVoxel2>
  GenericImage<VoxelType>& operator= (const GenericImage<TVoxel2> &);

  /// Cast to bool, checks if this image has been initialized and ready for use
  operator bool() const;

  /// Clear an image
  void Clear();

  // ---------------------------------------------------------------------------
  // Lattice

  /// Number of vector components per voxel
  int N() const;

  // ---------------------------------------------------------------------------
  // Image data access

  /// Function for pixel get access
  VoxelType Get(int) const;

  /// Function for pixel get access
  VoxelType Get(int, int, int = 0, int = 0) const;

  /// Function for pixel put access
  void Put(int, VoxelType);

  /// Function for pixel put access
  void Put(int, int, VoxelType);

  /// Function for pixel put access
  void Put(int, int, int, VoxelType);

  /// Function for pixel put access
  void Put(int, int, int, int, VoxelType);

  /// Function for pixel access from via operators
  VoxelType& operator()(int);

  /// Function for pixel access from via operators
  const VoxelType& operator()(int) const;

  /// Function for pixel access from via operators
  VoxelType& operator()(int, int, int = 0, int = 0);

  /// Function for pixel access from via operators
  const VoxelType& operator()(int, int, int = 0, int = 0) const;

  // ---------------------------------------------------------------------------
  // Type independent access to scalar image data

  /// Function for pixel get access as double
  virtual double GetAsDouble(int) const;

  /// Function for pixel get access as double
  virtual double GetAsDouble(int, int, int = 0, int = 0) const;

  /// Function for pixel put access
  virtual void PutAsDouble(int, double);
  
  /// Function for pixel put access
  virtual void PutAsDouble(int, int, double);
  
  /// Function for pixel put access
  virtual void PutAsDouble(int, int, int, double);
  
  /// Function for pixel put access
  virtual void PutAsDouble(int, int, int, int, double);

  /// Function for pixel get access as double
  virtual void GetAsVector(Vector &, int) const;

  /// Function for pixel get access as double
  virtual void GetAsVector(Vector &, int, int, int = 0, int = 0) const;

  /// Function for pixel get access as double
  virtual Vector GetAsVector(int) const;

  /// Function for pixel get access as double
  virtual Vector GetAsVector(int, int, int = 0, int = 0) const;

  /// Function for pixel put access
  virtual void PutAsVector(int, const Vector &);

  /// Function for pixel put access
  virtual void PutAsVector(int, int, const Vector &);

  /// Function for pixel put access
  virtual void PutAsVector(int, int, int, const Vector &);

  /// Function for pixel put access
  virtual void PutAsVector(int, int, int, int, const Vector &);

  // ---------------------------------------------------------------------------
  // Access to raw image data
  
  /// Get raw pointer to contiguous image data
  VoxelType *Data(int = 0);

  /// Get raw pointer to contiguous image data
  VoxelType *Data(int, int, int = 0, int = 0);

  /// Get raw pointer to contiguous image data
  const VoxelType *Data(int = 0) const;

  /// Get raw pointer to contiguous image data
  const VoxelType *Data(int, int, int = 0, int = 0) const;

  /// Get raw pointer to contiguous image data
  virtual void *GetDataPointer(int = 0);

  /// Get raw pointer to contiguous image data
  virtual const void *GetDataPointer(int = 0) const;

  /// Get raw pointer to contiguous image data
  virtual void *GetDataPointer(int, int, int = 0, int = 0);

  /// Get raw pointer to contiguous image data
  virtual const void *GetDataPointer(int, int, int = 0, int = 0) const;

  /// Get enumeration value corresponding to voxel type
  virtual int GetDataType() const;

  /// Get size of each voxel in bytes
  virtual int GetDataTypeSize() const;

  /// Minimum value a pixel can hold without overflowing
  virtual double GetDataTypeMin() const;
  
  /// Maximum value a pixel can hold without overflowing
  virtual double GetDataTypeMax() const;

  // ---------------------------------------------------------------------------
  // Region-of-interest extraction

  /// Get image consisting of specified 2D slice
  GenericImage GetRegion(int, int) const;

  /// Get image consisting of specified 2D slice
  void GetRegion(GenericImage &, int, int) const;

  /// Get image consisting of specified 2D slice
  virtual void GetRegion(BaseImage *&, int, int) const;

  /// Get image consisting of specified 3D subregion
  GenericImage GetRegion(int, int, int,
                         int, int, int) const;

  /// Get image consisting of specified 3D subregion
  void GetRegion(GenericImage &, int, int, int,
                                 int, int, int) const;

  /// Get image consisting of specified 3D subregion
  virtual void GetRegion(BaseImage *&, int, int, int,
                                       int, int, int) const;

  /// Get image consisting of specified 4D subregion
  GenericImage GetRegion(int, int, int, int,
                         int, int, int, int) const;

  /// Get image consisting of specified 4D subregion
  void GetRegion(GenericImage &, int, int, int, int,
                                 int, int, int, int) const;

  /// Get image consisting of specified 4D subregion
  virtual void GetRegion(BaseImage *&, int, int, int, int,
                                       int, int, int, int) const;

  /// Get time instance (i.e., frame) or channel of image
  GenericImage GetFrame(int, int = -1) const;

  /// Get time instance (i.e., frame) or channel of image
  void GetFrame(GenericImage &, int, int = -1) const;

  /// Get time instance (i.e., frame) or channel of image
  virtual void GetFrame(BaseImage *&, int, int = -1) const;

  // ---------------------------------------------------------------------------
  // Image arithmetic

  /// Equality operator
  /// \note Use explicit negation for inequality comparison.
  ///       The overloaded != operator is used for binarization of the image.
  template <class TVoxel2>
  bool operator== (const GenericImage<TVoxel2> &) const;

  GenericImage& operator+=(const GenericImage &); ///< Add image
  GenericImage& operator-=(const GenericImage &); ///< Subtract image
  GenericImage& operator*=(const GenericImage &); ///< Multipy voxels
  GenericImage& operator/=(const GenericImage &); ///< Divide voxels

  GenericImage& operator+=(double); ///< Add scalar
  GenericImage& operator-=(double); ///< Subtract scalar
  GenericImage& operator*=(double); ///< Multiply by scalar
  GenericImage& operator/=(double); ///< Divide by scalar

  GenericImage  operator+ (const GenericImage &) const; ///< Add images
  GenericImage  operator- (const GenericImage &) const; ///< Subtract images
  GenericImage  operator* (const GenericImage &) const; ///< Multiply images voxel-wise
  GenericImage  operator/ (const GenericImage &) const; ///< Divide images voxel-wise

  GenericImage  operator+ (double) const; ///< Add scalar to image
  GenericImage  operator- (double) const; ///< Subtract scalar from image
  GenericImage  operator* (double) const; ///< Multiply image by scalar
  GenericImage  operator/ (double) const; ///< Divide image by scalar

  // ---------------------------------------------------------------------------
  // Thresholding

  // Import other overload
  using BaseImage::PutBackgroundValueAsDouble;

  /// Put background value
  virtual void PutBackgroundValueAsDouble(double, bool);

  GenericImage& operator>=(VoxelType);       ///< Clamp image given upper threshold
  GenericImage& operator<=(VoxelType);       ///< Clamp image given lower threshold

  GenericImage  operator> (VoxelType) const; ///< Clamp image given upper threshold
  GenericImage  operator< (VoxelType) const; ///< Clamp image given lower threshold

  /// Get binary mask for voxels which are not equal the scalar
  BinaryImage operator!=(VoxelType) const;

  // ---------------------------------------------------------------------------
  // Common image statistics

  /// Minimum and maximum pixel values get accessor
  void GetMinMax(VoxelType &, VoxelType &) const;

  /// Minimum and maximum pixel values get accessor with padding
  void GetMinMax(VoxelType &, VoxelType &, VoxelType) const;

  /// Linearly rescale intensities
  void PutMinMax(VoxelType, VoxelType);

  /// Mean pixel value
  ///
  /// \param[in] fg Calculate mean of foreground only.
  RealType Mean(bool fg = true) const;

  /// Average pixel values get accessor
  /// \deprecated Use Mean instead.
  RealType GetAverage(int = 1) const;

  /// Standard Deviation of the pixels
  RealType GetSD(int = 1) const;

  /// Get Max Intensity position around the point
  void GetMaxPosition(Point &, int = 1, int = 0) const;
  
  /// Get Gravity center position of a given window
  void GravityCenter(Point &, int = 1, int = 0) const;

  // ---------------------------------------------------------------------------
  // Common image manipulations

  virtual void ReflectX();  ///< Reflect image along x
  virtual void ReflectY();  ///< Reflect image along y
  virtual void ReflectZ();  ///< Reflect image along z

  virtual void FlipXY(bool); ///< Flip x and y axis
  virtual void FlipXZ(bool); ///< Flip x and z axis
  virtual void FlipYZ(bool); ///< Flip y and z axis
  virtual void FlipXT(bool); ///< Flip x and t axis
  virtual void FlipYT(bool); ///< Flip y and t axis
  virtual void FlipZT(bool); ///< Flip z and t axis

  bool CropPad(int margin = 0); ///< Crop/pad image background

  // ---------------------------------------------------------------------------
  // VTK interface
  #if MIRTK_Image_WITH_VTK

  /// Convert image to VTK structured points
  ///
  /// \note Use only when MIRTK_Image_WITH_VTK is 1.
  virtual void ImageToVTK(vtkStructuredPoints *) const;

  /// Convert VTK structured points to image
  ///
  /// \note Use only when MIRTK_Image_WITH_VTK is 1.
  virtual void VTKToImage(vtkStructuredPoints *);

  #endif // MIRTK_Image_WITH_VTK
  // ---------------------------------------------------------------------------
  // I/O

  /// Read image from file
  virtual void Read(const char *);
  
  /// Write image to file
  virtual void Write(const char *) const;

  // ---------------------------------------------------------------------------
  // Deprecated

  /// Minimum and maximum pixel values get accessor
  /// \deprecated Use respective overloaded method of GetMinMax instead.
  void GetMinMax(VoxelType *, VoxelType *) const;
  
  /// Minimum and maximum pixel values get accessor with padding
  /// \deprecated Use respective overloaded method of GetMinMax instead.
  void GetMinMax(VoxelType *, VoxelType *, VoxelType) const;

  /// Minimum and maximum pixel values get accessor with padding
  /// \deprecated Use respective overloaded method of GetMinMax instead.
  void GetMinMaxPad(VoxelType *, VoxelType *, VoxelType) const;

  /// \returns Raw pointer to contiguous image data.
  /// \deprecated Use Data instead.
  VoxelType *GetPointerToVoxels(int = 0, int = 0, int = 0, int = 0);
  
  /// \returns Raw pointer to contiguous image data.
  /// \deprecated Use Data instead.
  const VoxelType *GetPointerToVoxels(int = 0, int = 0, int = 0, int = 0) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
GenericImage<VoxelType>& GenericImage<VoxelType>::operator=(const GenericImage<VoxelType2> &image)
{
  this->Initialize(image.GetImageAttributes());
  VoxelType        *ptr1 = this->GetPointerToVoxels();
  const VoxelType2 *ptr2 = image.GetPointerToVoxels();
  for (int idx = 0; idx < _NumberOfVoxels; idx++) {
    ptr1[idx] = voxel_cast<VoxelType>(ptr2[idx]);
  }
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
  return *this;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
GenericImage<VoxelType>::operator bool() const
{
  return bool(_attr) && _matrix != nullptr;
}

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int GenericImage<VoxelType>::N() const
{
  return voxel_info<VoxelType>::vector_size();
}

// =============================================================================
// Image data access
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::Put(int index, VoxelType val)
{
  _data[index] = val;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::Put(int x, int y, VoxelType val)
{
  _matrix[0][0][y][x] = val;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::Put(int x, int y, int z, VoxelType val)
{
  _matrix[0][z][y][x] = val;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::Put(int x, int y, int z, int t, VoxelType val)
{
  _matrix[t][z][y][x] = val;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType GenericImage<VoxelType>::Get(int index) const
{
  return _data[index];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType GenericImage<VoxelType>::Get(int x, int y, int z, int t) const
{
  return _matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &GenericImage<VoxelType>::operator ()(int index)
{
  return _data[index];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType &GenericImage<VoxelType>::operator ()(int index) const
{
  return _data[index];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType& GenericImage<VoxelType>::operator()(int x, int y, int z, int t)
{
  return _matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType& GenericImage<VoxelType>::operator()(int x, int y, int z, int t) const
{
  return _matrix[t][z][y][x];
}

// =============================================================================
// Image arithmetics
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType> template <class VoxelType2>
bool GenericImage<VoxelType>::operator==(const GenericImage<VoxelType2> &image) const
{
  if (this->GetImageAttributes() != image.GetImageAttributes()) return false;
  const VoxelType  *ptr1 = this->GetPointerToVoxels();
  const VoxelType2 *ptr2 = image.GetPointerToVoxels();
  for (int idx = 0; idx < image; ++idx) {
    if (IsForeground(idx) && image.IsForeground(idx) && ptr1[idx] != ptr2[idx]) {
      return false;
    }
  }
  return true;
}

// =============================================================================
// Type independent access to scalar image data
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::PutAsDouble(int index, double val)
{
  _data[index] = voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::PutAsDouble(int x, int y, double val)
{
  _matrix[0][0][y][x] = voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::PutAsDouble(int x, int y, int z, double val)
{
  _matrix[0][z][y][x] = voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::PutAsDouble(int x, int y, int z, int t, double val)
{
  _matrix[t][z][y][x] = voxel_cast<VoxelType>(val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double GenericImage<VoxelType>::GetAsDouble(int index) const
{
  return voxel_cast<double>(_data[index]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double GenericImage<VoxelType>::GetAsDouble(int x, int y, int z, int t) const
{
  return voxel_cast<double>(_matrix[t][z][y][x]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::PutAsVector(int index, const Vector &value)
{
  _data[index] = voxel_cast<VoxelType>(value);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::PutAsVector(int x, int y, const Vector &value)
{
  _matrix[0][0][y][x] = voxel_cast<VoxelType>(value);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::PutAsVector(int x, int y, int z, const Vector &value)
{
  _matrix[0][z][y][x] = voxel_cast<VoxelType>(value);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::PutAsVector(int x, int y, int z, int t, const Vector &value)
{
  _matrix[t][z][y][x] = voxel_cast<VoxelType>(value);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::GetAsVector(Vector &value, int index) const
{
  value = voxel_cast<Vector>(_data[index]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::GetAsVector(Vector &value, int x, int y, int z, int t) const
{
  value = voxel_cast<Vector>(_matrix[t][z][y][x]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline Vector GenericImage<VoxelType>::GetAsVector(int index) const
{
  return voxel_cast<Vector>(_data[index]);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline Vector GenericImage<VoxelType>::GetAsVector(int x, int y, int z, int t) const
{
  return voxel_cast<Vector>(_matrix[t][z][y][x]);
}

// =============================================================================
// Access to raw image data
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *GenericImage<VoxelType>::Data(int i)
{
  return _data + i;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *GenericImage<VoxelType>::Data(int i) const
{
  return _data + i;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *GenericImage<VoxelType>::Data(int x, int y, int z, int t)
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *GenericImage<VoxelType>::Data(int x, int y, int z, int t) const
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void *GenericImage<VoxelType>::GetDataPointer(int i)
{
  return _data + i;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const void *GenericImage<VoxelType>::GetDataPointer(int i) const
{
  return _data + i;
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void *GenericImage<VoxelType>::GetDataPointer(int x, int y, int z, int t)
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const void *GenericImage<VoxelType>::GetDataPointer(int x, int y, int z, int t) const
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int GenericImage<VoxelType>::GetDataType() const
{
  return voxel_info<VoxelType>::type();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int GenericImage<VoxelType>::GetDataTypeSize() const
{
  return static_cast<int>(sizeof(VoxelType));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double GenericImage<VoxelType>::GetDataTypeMin() const
{
  return voxel_limits<VoxelType>::min();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double GenericImage<VoxelType>::GetDataTypeMax() const
{
  return voxel_limits<VoxelType>::max();
}

// =============================================================================
// Deprecated
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::GetMinMax(VoxelType *min, VoxelType *max) const
{
  this->GetMinMax(*min, *max);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::GetMinMax(VoxelType *min, VoxelType *max, VoxelType pad) const
{
  this->GetMinMax(*min, *max, pad);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void GenericImage<VoxelType>::GetMinMaxPad(VoxelType *min, VoxelType *max, VoxelType pad) const
{
  this->GetMinMax(*min, *max, pad);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *GenericImage<VoxelType>::GetPointerToVoxels(int x, int y, int z, int t)
{
  return &_matrix[t][z][y][x];
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *GenericImage<VoxelType>::GetPointerToVoxels(int x, int y, int z, int t) const
{
  return &_matrix[t][z][y][x];
}

////////////////////////////////////////////////////////////////////////////////
// Common specializations
////////////////////////////////////////////////////////////////////////////////

typedef GenericImage<BytePixel> ByteImage;
typedef GenericImage<GreyPixel> GreyImage;
typedef GenericImage<RealPixel> RealImage;


} // namespace mirtk

#endif // MIRTK_GenericImage_H
