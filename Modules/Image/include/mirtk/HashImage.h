/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2015 Antonios Makropoulos
 * Copyright 2017 Andreas Schuh
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *   http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef MIRTK_HashImage_H
#define MIRTK_HashImage_H

#include "mirtk/BaseImage.h" 
#include "mirtk/VoxelCast.h"
#include "mirtk/UnorderedMap.h"

namespace mirtk {


/**
 * Hash class for N-D images
 *
 * This class implements N-D images based on hashmaps. It provides functions
 * for accessing, reading, writing and manipulating images. This class can
 * be used for images with arbitrary voxel types using templates.
 */

template <class TVoxel>
class HashImage : public BaseImage
{
  mirtkObjectMacro(HashImage);

  // ---------------------------------------------------------------------------
  // Types

public:

  /// Voxel type
  typedef TVoxel VoxelType;

  /// Floating point type corresponding to voxel type
  /// \note The VoxelType as well as the RealType may be a matrix/vector type!
  typedef typename voxel_info<VoxelType>::RealType RealType;

  /// Data map (hashmap)
  typedef UnorderedMap<int, VoxelType> DataMap;

  /// Data map iterator
  typedef typename DataMap::const_iterator DataIterator;

  /// Scalar type corresponding to voxel type
  typedef typename voxel_info<VoxelType>::ScalarType ScalarType;

  // ---------------------------------------------------------------------------
  // Data members

protected:

  /// Pointer array for access to image data
  ///
  /// \note The image data is stored in a hash map
  mirtkPublicAttributeMacro(DataMap, Data);

  /// \note Voxels that their value==DefaultValue are not stored
  mirtkPublicAttributeMacro(VoxelType, DefaultValue);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Allocate image memory
  void AllocateImage();

  /// Function for not const pixel get access
  VoxelType Access(int);

public:

  /// Default constructor
  HashImage();

  /// Constructor from image file
  HashImage(const char *);

  /// Constructor for given image size
  explicit HashImage(int, int, int = 1, int = 1);

  /// Constructor for given image size
  explicit HashImage(int, int, int, int, int);

  /// Constructor for given image attributes
  explicit HashImage(const ImageAttributes &, int = -1);

  /// Copy constructor for image
  explicit HashImage(const BaseImage &);

  /// Copy constructor for image
  HashImage(const HashImage<VoxelType> &);

  /// Copy constructor for image of different type
  template <class TVoxel2>
  HashImage(const GenericImage<TVoxel2> &);

  /// Copy constructor for image of different type
  template <class TVoxel2>
  HashImage(const HashImage<TVoxel2> &);

  /// Destructor
  ~HashImage();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Create copy of this image
  virtual BaseImage *Copy() const;

  /// Initialize a previously allocated image
  virtual void Initialize();

  /// Initialize an image
  virtual void Initialize(const ImageAttributes &, int = -1);

  /// Initialize an image
  void Initialize(int, int, int, int, int);

  /// Initialize an image
  void Initialize(int, int, int = 1, int = 1);

  /// Copy image data from other image of same size
  void CopyFrom(const BaseImage &);

  /// Copy image data from other image of same size
  template <class TVoxel2>
  void CopyFrom(const GenericImage<TVoxel2> &);

  /// Copy image data from other image of same size
  template <class TVoxel2>
  void CopyFrom(const HashImage<TVoxel2> &);

  /// Copy image data to GenericImage
  template <class TVoxel2>
  void CopyTo(GenericImage<TVoxel2> &) const;

  /// Assign constant value to each voxel
  HashImage& operator= (VoxelType);

  /// Assignment operator with implicit cast to double and then VoxelType
  HashImage<VoxelType>& operator= (const BaseImage &);

  /// Assignment operator with implicit cast to double and then VoxelType
  template <class TVoxel2>
  HashImage<VoxelType>& operator= (const GenericImage<TVoxel2> &);

  /// Assignment operator
  HashImage<VoxelType>& operator= (const HashImage &);

  /// Assignment operator with implicit cast
  template <class TVoxel2>
  HashImage<VoxelType>& operator= (const HashImage<TVoxel2> &);

  /// Clear an image
  void Clear();

  // ---------------------------------------------------------------------------
  // Lattice

  /// Number of vector components per voxel
  int N() const;

  /// Function to convert pixel to index
  /// more efficient than overwritten base class implementation
  int VoxelToIndex(int, int, int = 0, int = 0) const;

  // ---------------------------------------------------------------------------
  // Image data access

  DataIterator Begin() const;
  DataIterator End() const;

  int NumberOfNonDefaultVoxels() const;

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
  HashImage GetRegion(int, int) const;

  /// Get image consisting of specified 2D slice
  void GetRegion(HashImage &, int, int) const;

  /// Get image consisting of specified 2D slice
  virtual void GetRegion(BaseImage *&, int, int) const;

  /// Get image consisting of specified 3D subregion
  HashImage GetRegion(int, int, int,
             int, int, int) const;

  /// Get image consisting of specified 3D subregion
  void GetRegion(HashImage &, int, int, int,
                 int, int, int) const;

  /// Get image consisting of specified 3D subregion
  virtual void GetRegion(BaseImage *&, int, int, int,
                     int, int, int) const;

  /// Get image consisting of specified 4D subregion
  HashImage GetRegion(int, int, int, int,
             int, int, int, int) const;

  /// Get image consisting of specified 4D subregion
  void GetRegion(HashImage &, int, int, int, int,
                 int, int, int, int) const;

  /// Get image consisting of specified 4D subregion
  virtual void GetRegion(BaseImage *&, int, int, int, int,
                     int, int, int, int) const;

  /// Get time instance (i.e., frame) or channel of image
  HashImage GetFrame(int, int = -1) const;

  /// Get time instance (i.e., frame) or channel of image
  void GetFrame(HashImage &, int, int = -1) const;

  /// Get time instance (i.e., frame) or channel of image
  virtual void GetFrame(BaseImage *&, int, int = -1) const;

  // ---------------------------------------------------------------------------
  // Image arithmetic
  HashImage& operator+=(ScalarType); ///< Add scalar
  HashImage& operator-=(ScalarType); ///< Subtract scalar
  HashImage& operator*=(ScalarType); ///< Multiply by scalar
  HashImage& operator/=(ScalarType); ///< Divide by scalar

  /// Equality operator
  /// \note Use explicit negation for inequality comparison.
  ///     The overloaded != operator is used for binarization of the image.
  template <class TVoxel2>
  bool operator== (const HashImage<TVoxel2> &) const;

  HashImage& operator+=(const HashImage &); ///< Add image
  HashImage& operator-=(const HashImage &); ///< Subtract image
  HashImage& operator*=(const HashImage &); ///< Multipy voxels
  HashImage& operator/=(const HashImage &); ///< Divide voxels


  HashImage  operator+ (const HashImage &) const; ///< Add images
  HashImage  operator- (const HashImage &) const; ///< Subtract images
  HashImage  operator* (const HashImage &) const; ///< Multiply images voxel-wise
  HashImage  operator/ (const HashImage &) const; ///< Divide images voxel-wise

  HashImage  operator+ (ScalarType) const; ///< Add scalar to image
  HashImage  operator- (ScalarType) const; ///< Subtract scalar from image
  HashImage  operator* (ScalarType) const; ///< Multiply image by scalar
  HashImage  operator/ (ScalarType) const; ///< Divide image by scalar
  
  // ---------------------------------------------------------------------------
  // Thresholding

  // Import other overload
  using BaseImage::PutBackgroundValueAsDouble;

  /// Put background value
  virtual void PutBackgroundValueAsDouble(double, bool);

  // ---------------------------------------------------------------------------
  // Common image statistics

  /// Minimum and maximum pixel values get accessor
  void GetMinMax(VoxelType &, VoxelType &) const;

  /// Minimum and maximum pixel values get accessor with padding
  void GetMinMax(VoxelType &, VoxelType &, VoxelType) const;

  /// Linearly rescale intensities
  void PutMinMax(VoxelType, VoxelType);

  // ---------------------------------------------------------------------------
  // Common image manipulations

  virtual void ReflectX(bool modify_axes = false);  ///< Reflect image along x
  virtual void ReflectY(bool modify_axes = false);  ///< Reflect image along y
  virtual void ReflectZ(bool modify_axes = false);  ///< Reflect image along z
  virtual void ReflectT(bool modify_axes = false);  ///< Reflect image along t

  virtual void FlipXY(bool modify_origin = false); ///< Flip x and y axis, always also swaps voxel size
  virtual void FlipXZ(bool modify_origin = false); ///< Flip x and z axis, always also swaps voxel size
  virtual void FlipYZ(bool modify_origin = false); ///< Flip y and z axis, always also swaps voxel size
  virtual void FlipXT(bool modify_origin = false); ///< Flip x and t axis, always also swaps voxel size
  virtual void FlipYT(bool modify_origin = false); ///< Flip y and t axis, always also swaps voxel size
  virtual void FlipZT(bool modify_origin = false); ///< Flip z and t axis, always also swaps voxel size

  virtual void SwapXY(bool modify_axes = true); ///< Swap x and y axis
  virtual void SwapXZ(bool modify_axes = true); ///< Swap x and z axis
  virtual void SwapYZ(bool modify_axes = true); ///< Swap y and z axis
  virtual void SwapXT(bool modify_axes = true); ///< Swap x and t axis
  virtual void SwapYT(bool modify_axes = true); ///< Swap y and t axis
  virtual void SwapZT(bool modify_axes = true); ///< Swap z and t axis

  //bool CropPad(int margin = 0); ///< Crop/pad image background

  // ---------------------------------------------------------------------------
  // VTK interface
  #if defined(HAVE_VTK) && MIRTK_Image_WITH_VTK

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
  // Convertors
  GenericImage<VoxelType> ToGenericImage() const;

  template <class TVoxel2>
  void ToGenericImage(GenericImage<TVoxel2>& image) const;

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

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------

// =============================================================================
// Lattice
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int HashImage<VoxelType>::N() const
{
  return voxel_info<VoxelType>::vector_size();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int HashImage<VoxelType>::VoxelToIndex(int x, int y, int z, int t) const
{
  return x + y*_attr._x + z*_attr._x*_attr._y + t*_attr._x*_attr._y*_attr._z;
}

// =============================================================================
// Image data access
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::Put(int index, VoxelType val)
{
  if(index>=0){
    if(val==_DefaultValue){
      _Data.erase(index);
     }else{
      _Data[index]=val;
    }
  }
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::Put(int x, int y, VoxelType val)
{
  Put(x,y,0,0,val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::Put(int x, int y, int z, VoxelType val)
{
  Put(x,y,z,0,val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::Put(int x, int y, int z, int t, VoxelType val)
{
  Put(VoxelToIndex(x, y, z, t) ,val);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType HashImage<VoxelType>::Access(int index)
{
  cerr << "HashImage::Access: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------

template <class VoxelType>
inline typename HashImage<VoxelType>::DataIterator HashImage<VoxelType>::Begin() const
{
  return _Data.begin();
}

// -----------------------------------------------------------------------------

template <class VoxelType>
inline typename HashImage<VoxelType>::DataIterator HashImage<VoxelType>::End() const
{
  return _Data.end();
}

// -----------------------------------------------------------------------------

template <class VoxelType>
int HashImage<VoxelType>::NumberOfNonDefaultVoxels() const
{
  return _Data.size();
}

// -----------------------------------------------------------------------------

template <class VoxelType>
inline VoxelType HashImage<VoxelType>::Get(int index) const
{
  auto pos=_Data.find (index);
  if ( pos == End() ) return _DefaultValue;
  return voxel_cast<VoxelType>(pos->second);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType HashImage<VoxelType>::Get(int x, int y, int z, int t) const
{
  return Get(VoxelToIndex(x, y, z, t));
}

// =============================================================================
// Type independent access to scalar image data
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::PutAsDouble(int index, double val)
{
  Put(index, voxel_cast<VoxelType>(val));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::PutAsDouble(int x, int y, double val)
{
  Put(VoxelToIndex(x, y, 0, 0), voxel_cast<VoxelType>(val));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::PutAsDouble(int x, int y, int z, double val)
{
  Put(VoxelToIndex(x, y, z, 0), voxel_cast<VoxelType>(val));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::PutAsDouble(int x, int y, int z, int t, double val)
{
  Put(VoxelToIndex(x, y, z, t), voxel_cast<VoxelType>(val));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double HashImage<VoxelType>::GetAsDouble(int index) const
{
  return voxel_cast<double>(Get(index));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double HashImage<VoxelType>::GetAsDouble(int x, int y, int z, int t) const
{
  return voxel_cast<double>(Get(VoxelToIndex(x, y, z, t)));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::PutAsVector(int index, const Vector &value)
{
  Put(index, voxel_cast<VoxelType>(value));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::PutAsVector(int x, int y, const Vector &value)
{
  Put(VoxelToIndex(x, y, 0, 0), voxel_cast<VoxelType>(value));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::PutAsVector(int x, int y, int z, const Vector &value)
{
  Put(VoxelToIndex(x, y, z, 0), voxel_cast<VoxelType>(value));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::PutAsVector(int x, int y, int z, int t, const Vector &value)
{
  Put(VoxelToIndex(x, y, z, t), voxel_cast<VoxelType>(value));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::GetAsVector(Vector &value, int index) const
{
  value = voxel_cast<Vector>(Get(index));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::GetAsVector(Vector &value, int x, int y, int z, int t) const
{
  value = voxel_cast<Vector>(Get(VoxelToIndex(x, y, z, t)));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline Vector HashImage<VoxelType>::GetAsVector(int index) const
{
  return voxel_cast<Vector>(Get(index));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline Vector HashImage<VoxelType>::GetAsVector(int x, int y, int z, int t) const
{
  return voxel_cast<Vector>(Get(VoxelToIndex(x, y, z, t)));
}

// =============================================================================
// Access to raw image data
// =============================================================================

// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void *HashImage<VoxelType>::GetDataPointer(int i)
{
  cerr << "HashImage<VoxelType>::GetDataPointer: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const void *HashImage<VoxelType>::GetDataPointer(int i) const
{
  cerr << "HashImage<VoxelType>::GetDataPointer: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void *HashImage<VoxelType>::GetDataPointer(int x, int y, int z, int t)
{
  cerr << "HashImage<VoxelType>::GetDataPointer: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const void *HashImage<VoxelType>::GetDataPointer(int x, int y, int z, int t) const
{
  cerr << "HashImage<VoxelType>::GetDataPointer: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int HashImage<VoxelType>::GetDataType() const
{
  return voxel_info<VoxelType>::type();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline int HashImage<VoxelType>::GetDataTypeSize() const
{
  return static_cast<int>(sizeof(VoxelType));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double HashImage<VoxelType>::GetDataTypeMin() const
{
  return voxel_limits<VoxelType>::min();
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline double HashImage<VoxelType>::GetDataTypeMax() const
{
  return voxel_limits<VoxelType>::max();
}

// =============================================================================
// Deprecated
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::GetMinMax(VoxelType *min, VoxelType *max) const
{
  this->GetMinMax(*min, *max);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::GetMinMax(VoxelType *min, VoxelType *max, VoxelType pad) const
{
  this->GetMinMax(*min, *max, pad);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline void HashImage<VoxelType>::GetMinMaxPad(VoxelType *min, VoxelType *max, VoxelType pad) const
{
  this->GetMinMax(*min, *max, pad);
}
////////////////////////////////////////////////////////////////////////////////
// Common specializations
////////////////////////////////////////////////////////////////////////////////

typedef HashImage<BytePixel> HashByteImage;
typedef HashImage<GreyPixel> HashGreyImage;
typedef HashImage<RealPixel> HashRealImage;

} // namespace mirtk

#include "mirtk/HashImage.hxx"

#endif // MIRTK_HashImage_H
