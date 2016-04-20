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

#ifndef MIRTK_ConstImageIterator_H
#define MIRTK_ConstImageIterator_H

#include "mirtk/Voxel.h"
#include "mirtk/Vector4D.h"
#include "mirtk/ImageAttributes.h"
#include "mirtk/BaseImage.h"


namespace mirtk {


/**
 * Base class of const image iterator
 */
class ConstImageIterator
{
  // ---------------------------------------------------------------------------
  // Initialization
protected:

  /// Calculate pointer increments
  void CalculateStride();

  /// Initialize iterator, called by first GoTo command
  void Initialize();

public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  ConstImageIterator(const ImageAttributes &, int);

  /// Constructor
  ConstImageIterator(const ImageAttributes &, const void * = NULL, int = MIRTK_VOXEL_UNKNOWN);

  /// Constructor
  ConstImageIterator(const BaseImage &);

  /// Constructor
  ConstImageIterator(const BaseImage *);

  /// Copy constructor
  ConstImageIterator(const ConstImageIterator &);

  /// Assignment operator
  ConstImageIterator &operator =(const ConstImageIterator &);

  /// Destructor
  virtual ~ConstImageIterator();

  /// Get stride between columns in number of voxels
  int ColumnStride() const;

  /// Get stride between rows/lines in number of voxels
  int LineStride() const;

  /// Get stride between slices in number of voxels
  int SliceStride() const;

  /// Get stride between frames in number of voxels
  int FrameStride() const;

  // ---------------------------------------------------------------------------
  // Image attributes

  /// Set data type - defines number of bytes per voxel
  void SetDataType(int);

  /// Set raw data pointer to start of entire image
  void SetData(const void *, int = MIRTK_VOXEL_UNKNOWN);

  /// Whether the image is a 3D+t image sequence
  /// (i.e., number of voxels in t dimension is greater 1 and dt is not 0)
  bool IsImageSequence() const;

  /// Whether the image is a scalar image
  /// (i.e., number of voxels in t dimension is 1 or dt equal 0)
  bool IsScalarImage() const;

  /// Number of voxels in the entire image
  int NumberOfImageVoxels() const;

  /// Get number of image channels
  /// (i.e., number of voxels in t dimension if dt <= 0 or 1 otherwise)
  int NumberOfImageChannels() const;

  /// Get number of vector components
  /// (i.e., number of voxels in t dimension if dt <= 0 or 1 otherwise)
  int NumberOfVectorComponents() const;

  /// Get number of iterated frames
  /// (i.e., number of voxels in t dimension if dt > 0 or 1 otherwise)
  int NumberOfSequenceFrames() const;

  // ---------------------------------------------------------------------------
  // Region attributes

  /// Whether the image region is a 3D+t image sequence
  /// (i.e., number of voxels in t dimension > 1 and dt > 0)
  bool IsSequence() const;

  /// Whether the image region is scalar
  /// (i.e., number of voxels in t dimension is 1)
  bool IsScalar() const;

  /// Number of voxels in image region
  int NumberOfVoxels() const;

  /// Get number of channels in image region
  int NumberOfChannels() const;

  /// Get number of vector components considered
  int NumberOfComponents() const;

  /// Get number of frames in image region
  int NumberOfFrames() const;

  // ---------------------------------------------------------------------------
  // Region/Neighborhood

  /// Set 2D image region (start and size)
  void SetRegion(int, int, int, int);

  /// Set 2D image region (start and end)
  void SetRegion(const blocked_range2d<int> &);

  /// Set 3D image region (start and size)
  void SetRegion(int, int, int, int, int, int);

  /// Set 3D image region (start and end)
  void SetRegion(const blocked_range3d<int> &);

  /// Set 2D neighborhood (center and radius)
  void SetNeighborhood(int, int, int, int);

  /// Set 3D neighborhood (center and radius)
  void SetNeighborhood(int, int, int, int, int, int);

  /// Set temporal region (start and size)
  void SetFrame(int, int = 1);

  /// Set temporal region (start and end)
  void SetFrame(const blocked_range<int> &);

  /// Set temporal region (start and size)
  void SetChannel(int, int = 1);

  /// Set temporal region (start and end)
  void SetChannel(const blocked_range<int> &);

  /// Set temporal region (start and size)
  void SetComponent(int, int = 1);

  /// Set temporal region (start and end)
  void SetComponent(const blocked_range<int> &);

  /// Set 4D image region (start and size)
  void SetRegion(int, int, int, int, int, int, int, int);

  /// Set 4D neighborhood (center and radius)
  void SetNeighborhood(int, int, int, int, int, int, int, int);

  // ---------------------------------------------------------------------------
  // Position

  /// Convert iterator position to voxel coordinates
  void IndexToVoxel(int, int &, int &) const;

  /// Convert iterator position to voxel coordinates
  void IndexToVoxel(int, int &, int &, int &) const;

  /// Convert iterator position to voxel coordinates
  void IndexToVoxel(int, int &, int &, int &, int &) const;

  /// Convert iterator position to voxel coordinates
  void IndexToVoxel(int, Vector4D<int> &) const;

  /// Convert voxel coordinates to index
  int VoxelToIndex(int, int, int = 0, int = 0) const;

  /// Convert voxel coordinates to index
  int VoxelToIndex(const Vector4D<int> &) const;

  /// Convert iterator position to voxel coordinates
  void PosToVoxel(int, int &, int &) const;

  /// Convert iterator position to voxel coordinates
  void PosToVoxel(int, int &, int &, int &) const;

  /// Convert iterator position to voxel coordinates
  void PosToVoxel(int, int &, int &, int &, int &) const;

  /// Convert iterator position to voxel coordinates
  void PosToVoxel(int, Vector4D<int> &) const;

  /// Convert iterator position to voxel index
  int PosToIndex(int) const;

  /// Convert voxel index to iterator position
  int IndexToPos(int) const;

  /// Convert voxel cooridnates to iterator position
  int VoxelToPos(int, int, int = 0, int = 0) const;

  /// Convert voxel cooridnates to iterator position
  int VoxelToPos(const Vector4D<int> &) const;

  /// Current iterator position within image region
  int Pos() const;

  /// Index of voxel at current iterator position
  int Index() const;

  /// Coordinates of voxel at current iterator position
  void Voxel(int &, int &) const;

  /// Coordinates of voxel at current iterator position
  void Voxel(int &, int &, int &) const;

  /// Coordinates of voxel at current iterator position
  void Voxel(int &, int &, int &, int &) const;

  // ---------------------------------------------------------------------------
  // Iteration

  /// Go to voxel with specified coordinates relative to entire image
  void GoToVoxel(int, int, int = -1, int = -1);

  /// Go to voxel with specified coordinates relative to entire image
  void GoToVoxel(const Vector4D<int> &);

  /// Go to voxel with specified index relative to entire image
  void GoToIndex(int);

  /// Go to begin of region
  void GoToBegin();

  /// Go to center of region
  void GoToCenter();

  /// Go to end of region
  void GoToEnd();

  /// Go to specified position within image region
  void GoToPos(int);

  /// Whether iterator is valid and not yet at end of region
  operator bool() const;

  /// Whether iterator reached the start of the region and is invalid now
  bool IsAtBegin() const;

  /// Whether iterator reached the end of the region and is invalid now
  bool IsAtEnd() const;

  /// Pre-decrement operator
  void operator --();

  /// Pre-increment operator
  void operator ++();

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  template <class VoxelType>
  const VoxelType *Current() const;

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame stride relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  template <class VoxelType>
  const VoxelType *Current(int) const;

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  template <class VoxelType>
  const VoxelType *Next();

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame stride relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  template <class VoxelType>
  const VoxelType *Next(int);

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  template <class VoxelType>
  const VoxelType &Value() const;

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame stride relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  template <class VoxelType>
  const VoxelType &Value(int t) const;

  /// Get current voxel value casted to double
  virtual double ValueAsDouble() const;

  /// Get current voxel value casted to double
  ///
  /// \param[in] t Channel/Component/Frame stride relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  virtual double ValueAsDouble(int) const;

  /// Move given raw image data pointer to the current iterator position assuming
  /// it is at the voxel position at which the iterator has been before
  ///
  /// This method is useful when iterating over multiple images with common
  /// image attributes simultaneously. An iterator instance associated with one
  /// of the images or initialized with the image attributes corresponding to all
  /// images is then sufficient to reduce the associated overhead. Raw pointers
  /// to the voxel data of the other images can then be moved to the next position
  /// using this method after the iterator itself has been incremented or decremented.
  ///
  /// Note that after a GoTo method of the iterator has been called, this method
  /// moves the raw image data pointer returned by BaseImage::GetScalarPointer
  /// to the same position as the iterator was moved to by the GoTo method.
  ///
  /// \code
  /// ConstImageIterator it(gradient.GetImageAttributes(), gradient.GetScalarType());
  /// double *gx = gradient.GetPointerToVoxels(0, 0, 0, 0);
  /// double *gy = gradient.GetPointerToVoxels(0, 0, 0, 1);
  /// double *gz = gradient.GetPointerToVoxels(0, 0, 0, 2);
  ///
  /// it.GoToBegin();
  /// while (!it.IsAtEnd()) {
  ///   // do something with the data
  ///   ++it; // must be done before moving the other data pointers
  ///   it.Move(gx), it.Move(gy), it.Move(gz);
  /// }
  ///
  /// it.GoToEnd();
  /// while (!it.IsAtBegin()) {
  ///   // do something with the data
  ///   --it; // must be done before moving the other data pointers
  ///   it.Move(gx), it.Move(gy), it.Move(gz);
  /// }
  /// \endcode
  template <class VoxelType>
  void Move(const VoxelType *&) const;

  /// Move given raw image data pointer to the current iterator position assuming
  /// it is at the voxel position at which the iterator has been before
  template <class VoxelType>
  void Move(VoxelType *&) const;

  // ---------------------------------------------------------------------------
  // Members
protected:

  const char      *_Data;              ///< Pointer to begin of entire image
  Vector4D<int>    _DataSize;          ///< Size of the entire image
  bool             _IsImageSequence;   ///< Whether the image is a sequence

  Vector4D<int>    _Index;             ///< Start index of the image region
  Vector4D<int>    _Begin;             ///< First index of the image region
  Vector4D<int>    _End;               ///< Last index of the image region
  Vector4D<int>    _Size;              ///< Size of image region
  int              _XYZ;               ///< Number of voxels per frame of the image
  int              _ColumnStride;      ///< Increment in number of bytes (size of voxel type)
  int              _LineStride;        ///< Increment at end of line  in number of bytes
  int              _SliceStride;       ///< Increment at end of slice in number of bytes
  int              _FrameStride;       ///< Increment at end of frame in number of bytes
  int              _Inc;               ///< Previous pointer increment in number of bytes
  const char      *_Next;              ///< Pointer to next image voxel

};

//////////////////////////////////////////////////////////////////////////////
// Inline definitions
//////////////////////////////////////////////////////////////////////////////

// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
inline ConstImageIterator::ConstImageIterator(const ImageAttributes &attr, int type)
:
  _Data             (NULL),
  _DataSize         (attr._x, attr._y, attr._z, attr._t),
  _IsImageSequence  (attr._t > 1 && attr._dt > .0),
  _Index            (0, 0, 0, 0),
  _Begin            (0, 0, 0, 0),
  _End              (attr._x-1, attr._y-1, attr._z-1, attr._t-1),
  _Size             (attr._x,   attr._y,   attr._z,   attr._t),
  _XYZ              (attr._x * attr._y * attr._z),
  _ColumnStride     (DataTypeSize(type)),
  _LineStride       (0),
  _SliceStride      (0),
  _FrameStride      (0),
  _Inc              (0),
  _Next             (NULL)
{
}

// ---------------------------------------------------------------------------
inline ConstImageIterator::ConstImageIterator(const ImageAttributes &attr, const void *data, int type)
:
  _Data             (reinterpret_cast<const char *>(data)),
  _DataSize         (attr._x, attr._y, attr._z, attr._t),
  _IsImageSequence  (attr._t > 1 && attr._dt > .0),
  _Index            (0, 0, 0, 0),
  _Begin            (0, 0, 0, 0),
  _End              (attr._x-1, attr._y-1, attr._z-1, attr._t-1),
  _Size             (attr._x,   attr._y,   attr._z,   attr._t),
  _XYZ              (attr._x * attr._y * attr._z),
  _ColumnStride     (DataTypeSize(type)),
  _LineStride       (0),
  _SliceStride      (0),
  _FrameStride      (0),
  _Inc              (0),
  _Next             (NULL)
{
}

// ---------------------------------------------------------------------------
inline ConstImageIterator::ConstImageIterator(const BaseImage &image)
:
  _Data             (reinterpret_cast<const char *>(image.GetScalarPointer())),
  _DataSize         (image.GetX(), image.GetY(), image.GetZ(), image.GetT()),
  _IsImageSequence  (image.GetT() > 1 && image.GetTSize() > .0),
  _Index            (0, 0, 0, 0),
  _Begin            (0, 0, 0, 0),
  _End              (image.GetX()-1, image.GetY()-1, image.GetZ()-1, image.GetT()-1),
  _Size             (image.GetX(),   image.GetY(),   image.GetZ(),   image.GetT()),
  _XYZ              (image.GetX() * image.GetY() * image.GetZ()),
  _ColumnStride     (image.GetScalarTypeSize()),
  _LineStride       (0),
  _SliceStride      (0),
  _FrameStride      (0),
  _Inc              (0),
  _Next             (NULL)
{
}

// ---------------------------------------------------------------------------
inline ConstImageIterator::ConstImageIterator(const BaseImage *image)
:
  _Data             (reinterpret_cast<const char *>(image->GetScalarPointer())),
  _DataSize         (image->GetX(), image->GetY(), image->GetZ(), image->GetT()),
  _IsImageSequence  (image->GetT() > 1 && image->GetTSize() > .0),
  _Index            (0, 0, 0, 0),
  _Begin            (0, 0, 0, 0),
  _End              (image->GetX()-1, image->GetY()-1, image->GetZ()-1, image->GetT()-1),
  _Size             (image->GetX(),   image->GetY(),   image->GetZ(),   image->GetT()),
  _XYZ              (image->GetX() * image->GetY() * image->GetZ()),
  _ColumnStride     (image->GetScalarTypeSize()),
  _LineStride       (0),
  _SliceStride      (0),
  _FrameStride      (0),
  _Inc              (0),
  _Next             (NULL)
{
}

// ---------------------------------------------------------------------------
inline ConstImageIterator::ConstImageIterator(const ConstImageIterator &other)
:
  _Data             (other._Data),
  _DataSize         (other._DataSize),
  _IsImageSequence  (other._IsImageSequence),
  _Index            (other._Index),
  _Begin            (other._Begin),
  _End              (other._End),
  _Size             (other._Size),
  _XYZ              (other._XYZ),
  _ColumnStride     (other._ColumnStride),
  _LineStride       (other._LineStride),
  _SliceStride      (other._SliceStride),
  _FrameStride      (other._FrameStride),
  _Inc              (other._Inc),
  _Next             (other._Next)
{
}

// ---------------------------------------------------------------------------
inline ConstImageIterator &ConstImageIterator::operator =(const ConstImageIterator &rhs)
{
  _Data              = rhs._Data;
  _DataSize          = rhs._DataSize;
  _IsImageSequence   = rhs._IsImageSequence;
  _Index             = rhs._Index;
  _Begin             = rhs._Begin;
  _End               = rhs._End;
  _Size              = rhs._Size;
  _XYZ               = rhs._XYZ;
  _ColumnStride      = rhs._ColumnStride;
  _LineStride        = rhs._LineStride;
  _SliceStride       = rhs._SliceStride;
  _FrameStride       = rhs._FrameStride;
  _Inc               = rhs._Inc;
  _Next              = rhs._Next;
  return *this;
}

// ---------------------------------------------------------------------------
inline ConstImageIterator::~ConstImageIterator()
{
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::CalculateStride()
{
  // Pointer increments at end of line/slice/frame in number of voxels
  _LineStride  =   _DataSize._x - _Size._x + 1;
  _SliceStride =  (_DataSize._y - _Size._y + 1) * _DataSize._x - _Size._x + 1;
  _FrameStride = ((_DataSize._z - _Size._z + 1) * _DataSize._y - _Size._y + 1) * _DataSize._x - _Size._x + 1;
  // Multiply by number of bytes per voxel
  if (_ColumnStride == 0) _ColumnStride = 1;
  _LineStride  *= _ColumnStride;
  _SliceStride *= _ColumnStride;
  _FrameStride *= _ColumnStride;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::ColumnStride() const
{
  return 1;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::LineStride() const
{
  return _LineStride / _ColumnStride;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::SliceStride() const
{
  return _SliceStride / _ColumnStride;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::FrameStride() const
{
  return _FrameStride / _ColumnStride;
}

// ===========================================================================
// Image attributes
// ===========================================================================

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetDataType(int type)
{
  const int size = DataTypeSize(type);
  if (size == 0) {
    cerr << "ConstImageIterator::SetDataType: Unknown data type: " << type << endl;
    exit(1);
  }
  if (_ColumnStride != size) {
    _ColumnStride = size;
    _LineStride   = 0;
  }
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetData(const void *data, int type)
{
  if (type != MIRTK_VOXEL_UNKNOWN) {
    SetDataType(type);
  }
  _Data = reinterpret_cast<const char *>(data);
  _Next = NULL;
}

// ---------------------------------------------------------------------------
inline bool ConstImageIterator::IsImageSequence() const
{
  return _IsImageSequence;
}

// ---------------------------------------------------------------------------
inline bool ConstImageIterator::IsScalarImage() const
{
  return _DataSize._t == 1 || _IsImageSequence;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::NumberOfImageVoxels() const
{
  return _DataSize._x * _DataSize._y * _DataSize._z * _DataSize._t;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::NumberOfImageChannels() const
{
  return _IsImageSequence ? 1 : _DataSize._t;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::NumberOfVectorComponents() const
{
  return NumberOfImageChannels();
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::NumberOfSequenceFrames() const
{
  return _IsImageSequence ? _DataSize._t : 1;
}

// ===========================================================================
// Region attributes
// ===========================================================================

// ---------------------------------------------------------------------------
inline bool ConstImageIterator::IsSequence() const
{
  return _Size._t > 1 && _IsImageSequence;
}

// ---------------------------------------------------------------------------
inline bool ConstImageIterator::IsScalar() const
{
  return _Size._t == 1;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::NumberOfVoxels() const
{
  return _Size._x * _Size._y * _Size._z * _Size._t;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::NumberOfChannels() const
{
  return IsSequence() ? 1 : _Size._t;
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::NumberOfComponents() const
{
  return NumberOfChannels();
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::NumberOfFrames() const
{
  return IsSequence() ? _Size._t : 1;
}

// ===========================================================================
// Region
// ===========================================================================

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetRegion(int i, int j, int k, int ni, int nj, int nk)
{
  // Allow negative size value which means all voxels up to the image boundary
  if (ni < 0) ni = _DataSize._x - i;
  if (nj < 0) nj = _DataSize._y - j;
  if (nk < 0) nk = _DataSize._z - k;
  // Note: DO NOT change frame/channel selection here!
  //       Done by SetFrame/SetChannel/SetComponent.
  _Begin._x = i;
  _Begin._y = j;
  _Begin._z = k;
  _End  ._x = i + ni - 1;
  _End  ._y = j + nj - 1;
  _End  ._z = k + nk - 1;
  _Size ._x = ni;
  _Size ._y = nj;
  _Size ._z = nk;
  // Check image region
  if (_Begin._x < 0 || _Begin._x >= _DataSize._x || _End._x >= _DataSize._x ||
      _Begin._y < 0 || _Begin._y >= _DataSize._y || _End._y >= _DataSize._y ||
      _Begin._z < 0 || _Begin._z >= _DataSize._z || _End._z >= _DataSize._z ||
      _Begin._t < 0 || _Begin._t >= _DataSize._t || _End._t >= _DataSize._t) {
    cerr << "ConstImageIterator::SetRegion: Region at least partially outside of image domain" << endl;
    exit(1);
  }
  // Invalidate current strides
  _LineStride = 0;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetRegion(int i, int j, int ni, int nj)
{
  SetRegion(i, j, 0, ni, nj, 1);
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetRegion(const blocked_range2d<int> &r)
{
  SetRegion(r.cols().begin(),
            r.rows().begin(),
            r.cols().end() - r.cols().begin(),
            r.rows().end() - r.rows().begin());
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetRegion(const blocked_range3d<int> &r)
{
  SetRegion(r.cols ().begin(),
            r.rows ().begin(),
            r.pages().begin(),
            r.cols ().end() - r.cols ().begin(),
            r.rows ().end() - r.rows ().begin(),
            r.pages().end() - r.pages().begin());
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetChannel(int c, int nc)
{
  if (!_IsImageSequence) {
    if (c < 0 || c >= _DataSize._t) {
      cerr << "ConstImageIterator::SetChannel/Component: Channel index out of bounds (c=" << c << ")" << endl;
      exit(1);
    }
    // Allow negative value to select all channels starting at c
    if (nc < 1) nc = _DataSize._t - c;
    if ((c + nc) > _DataSize._t) {
      cerr << "ConstImageIterator::SetChannel/Component: Channel range out of bounds (c=" << c << ", nc=" << nc << ")" << endl;
      exit(1);
    }
    // Set range of channels
    _Begin._t = c;
    _End  ._t = c + nc - 1;
    _Size ._t = nc;
    // Invalidate current strides
    _LineStride = 0;
  } else if (c != 0 || nc != 1) {
    cerr << "ConstImageIterator::SetChannel/Component: Channel index/range out of bounds (c=" << c << ", nc=" << nc << ")" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetChannel(const blocked_range<int> &r)
{
  SetChannel(r.begin(), r.end() - r.begin());
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetComponent(int c, int nc)
{
  SetChannel(c, nc);
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetComponent(const blocked_range<int> &r)
{
  SetFrame(r.begin(), r.end() - r.begin());
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetFrame(int l, int nl)
{
  if (_IsImageSequence) {
    if (l < 0 || l >= _DataSize._t) {
      cerr << "ConstImageIterator::SetFrame: Frame index out of bounds (l=" << l << ")" << endl;
      exit(1);
    }
    // Allow negative value to select all frames starting at c
    if (nl < 1) nl = _DataSize._t - l;
    if ((l + nl) > _DataSize._t) {
      cerr << "ConstImageIterator::SetFrame: Frame range out of bounds (l=" << l << ", nl=" << nl << ")" << endl;
      exit(1);
    }
    // Set range of frames
    _Begin._t = l;
    _End  ._t = l + nl - 1;
    _Size ._t = nl;
    // Invalidate current strides
    _LineStride = 0;
  } else if (l != 0 || nl != 1) {
    cerr << "ConstImageIterator::SetFrame: Frame index/range out of bounds (l=" << l << ", nl=" << nl << ")" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetFrame(const blocked_range<int> &r)
{
  SetFrame(r.begin(), r.end() - r.begin());
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetRegion(int i, int j, int k, int l, int ni, int nj, int nk, int nl)
{
  SetRegion(i, j, k, ni, nj, nk);
  if (IsSequence()) SetFrame  (l, nl);
  else              SetChannel(l, nl);
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetNeighborhood(int i, int j, int k, int ri, int rj, int rk)
{
  if (ri < 0) ri = 0;
  if (rj < 0) rj = 0;
  if (rk < 0) rk = 0;
  if (0 <= i && i < _DataSize._x && 0 <= j && j < _DataSize._y && 0 <= k && k < _DataSize._z) {
    int i1 = i - ri;
    int j1 = j - rj;
    int k1 = k - rk;
    int i2 = i + ri;
    int j2 = j + rj;
    int k2 = k + rk;

    if (i1 <               0) i1 = 0;
    if (j1 <               0) j1 = 0;
    if (k1 <               0) k1 = 0;
    if (i2 >= _DataSize._x) i2 = _DataSize._x - 1;
    if (j2 >= _DataSize._y) j2 = _DataSize._y - 1;
    if (k2 >= _DataSize._z) k2 = _DataSize._z - 1;

    SetRegion(i1, j1, k1, i2-i1+1, j2-j1+1, k2-k1+1);
  } else {
    cerr << "ConstImageIterator::SetNeighborhood: Neighborhood out of bounds" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetNeighborhood(int i, int j, int ri, int rj)
{
  SetNeighborhood(i, j, 0, ri, rj, 0);
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::SetNeighborhood(int i, int j, int k, int l, int ri, int rj, int rk, int rl)
{
  if (ri < 0) ri = 0;
  if (rj < 0) rj = 0;
  if (rk < 0) rk = 0;
  if (rl < 0) rl = 0;
  SetNeighborhood(i, j, k, ri, rj, rk);
  if (_IsImageSequence) {
    if (0 <= l && l < _DataSize._t) {
      int l1 = l - rl;
      int l2 = l + rl;
      if (l1 <               0) l1 = 0;
      if (l2 >= _DataSize._t) l2 = _DataSize._t - 1;
      SetFrame(l1, l2-l1+1);
    } else {
      cerr << "ConstImageIterator::SetNeighborhood: Neighborhood out of bounds" << endl;
      exit(1);
    }
  } else if (l != 0 || rl != 0) {
    cerr << "ConstImageIterator::SetNeighborhood: 4D neighborhood only valid for image sequences (i.e., dt > 0)" << endl;
    exit(1);
  }
}

// ===========================================================================
// Position
// ===========================================================================

// ---------------------------------------------------------------------------
inline void ConstImageIterator::IndexToVoxel(int idx, int &i, int &j) const
{
  const int n = _DataSize._x * _DataSize._y;
  idx = idx % (n * _DataSize._z) % n;
  j   = idx / _DataSize._x;
  i   = idx % _DataSize._x;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::IndexToVoxel(int idx, int &i, int &j, int &k) const
{
  int n;
  n   = _DataSize._x * _DataSize._y * _DataSize._z;
  idx = idx % n;
  n   = _DataSize._x * _DataSize._y;
  k   = idx / n;
  j   = idx % n / _DataSize._x;
  i   = idx % n % _DataSize._x;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::IndexToVoxel(int idx, int &i, int &j, int &k, int &l) const
{
  int n;
  n   = _DataSize._x * _DataSize._y * _DataSize._z;
  l   = idx / n;
  idx = idx % n;
  n   = _DataSize._x * _DataSize._y;
  k   = idx / n;
  j   = idx % n / _DataSize._x;
  i   = idx % n % _DataSize._x;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::IndexToVoxel(int idx, Vector4D<int> &v) const
{
  IndexToVoxel(v._x, v._y, v._z, v._t);
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::VoxelToIndex(int i, int j, int k, int l) const
{
  return i + _DataSize._x * (j + _DataSize._y * (k + _DataSize._z * l));
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::VoxelToIndex(const Vector4D<int> &v) const
{
  return VoxelToIndex(v._x, v._y, v._z, v._t);
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::PosToVoxel(int pos, int &i, int &j) const
{
  const int n = _Size._x * _Size._y;
  pos = pos % (n * _Size._z) % n;
  j   = pos / _Size._x;
  i   = pos % _Size._x;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::PosToVoxel(int pos, int &i, int &j, int &k) const
{
  int n;
  n   = _Size._x * _Size._y * _Size._z;
  pos = pos % n;
  n   = _Size._x * _Size._y;
  k   = pos / n;
  j   = pos % n / _Size._x;
  i   = pos % n % _Size._x;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::PosToVoxel(int pos, int &i, int &j, int &k, int &l) const
{
  int n;
  n   = _Size._x * _Size._y * _Size._z;
  l   = pos / n;
  pos = pos % n;
  n   = _Size._x * _Size._y;
  k   = pos / n;
  j   = pos % n / _Size._x;
  i   = pos % n % _Size._x;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::PosToVoxel(int pos, Vector4D<int> &v) const
{
  return PosToVoxel(v._x, v._y, v._z, v._t);
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::PosToIndex(int pos) const
{
  int i, j, k, l;
  PosToVoxel(pos, i, j, k, l);
  return VoxelToIndex(i, j, k, l);
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::VoxelToPos(int i, int j, int k, int l) const
{
  return (i - _Begin._x) + _Size._x * ((j - _Begin._y) + _Size._y * ((k - _Begin._z) + _Size._z * (l - _Begin._t)));
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::VoxelToPos(const Vector4D<int> &v) const
{
  return VoxelToPos(v._x, v._y, v._z, v._t);
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::IndexToPos(int idx) const
{
  int i, j, k, l;
  IndexToVoxel(idx, i, j, k, l);
  return VoxelToPos(i, j, k, l);
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::Pos() const
{
  return VoxelToPos(_Index);
}

// ---------------------------------------------------------------------------
inline int ConstImageIterator::Index() const
{
  return static_cast<int>(_Next - reinterpret_cast<const char *>(_Data)) / _ColumnStride;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::Voxel(int &i, int &j) const
{
  i = _Index._x;
  j = _Index._y;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::Voxel(int &i, int &j, int &k) const
{
  i = _Index._x;
  j = _Index._y;
  k = _Index._z;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::Voxel(int &i, int &j, int &k, int &l) const
{
  i = _Index._x;
  j = _Index._y;
  k = _Index._z;
  l = _Index._t;
}

// ===========================================================================
// Iteration
// ===========================================================================

// ---------------------------------------------------------------------------
inline void ConstImageIterator::Initialize()
{
  if (_LineStride == 0) CalculateStride();
  _Inc = VoxelToIndex(_Index) * _ColumnStride;
  if (_Data) _Next = reinterpret_cast<const char *>(_Data) + _Inc;
  else       _Next = nullptr;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::GoToVoxel(const Vector4D<int> &v)
{
  _Index = v;
  Initialize();
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::GoToVoxel(int i, int j, int k, int l)
{
  _Index = Vector4D<int>(i, j, k, l);
  Initialize();
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::GoToIndex(int idx)
{
  IndexToVoxel(idx, _Index);
  Initialize();
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::GoToPos(int pos)
{
  PosToVoxel(pos, _Index);
  Initialize();
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::GoToBegin()
{
  _Index = _Begin;
  Initialize();
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::GoToCenter()
{
  _Index = _Begin + _Size/2;
  Initialize();
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::GoToEnd()
{
  _Index = _End;
  Initialize();
}

// ---------------------------------------------------------------------------
inline bool ConstImageIterator::IsAtBegin() const
{
  return _Index._x < _Begin._x;
}

// ---------------------------------------------------------------------------
inline bool ConstImageIterator::IsAtEnd() const
{
  return _Index._t > _End._t;
}

// ---------------------------------------------------------------------------
inline ConstImageIterator::operator bool() const
{
  return !(IsAtBegin() || IsAtEnd());
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::operator --()
{
  // Note: Store number of bytes skipped in member variable so that the
  //       Move method can make use of it to increment external pointers.
  if (_Index._t > _Begin._t) {
    _Index._x = _End._x;
    _Index._y = _End._y;
    _Index._z = _End._z;
    --_Index._t;
    _Inc = -_FrameStride;
  } else if (_Index._z > _Begin._z) {
    _Index._x = _End._x;
    _Index._y = _End._y;
    --_Index._z;
    _Inc = -_SliceStride;
  } else if (_Index._y > _Begin._y) {
    _Index._x = _End._x;
    --_Index._y;
    _Inc = -_LineStride;
  } else {
    --_Index._x;
    _Inc = -_ColumnStride;
  }
  if (_Next) _Next += _Inc;
}

// ---------------------------------------------------------------------------
inline void ConstImageIterator::operator ++()
{
  // Note: Store number of bytes skipped in member variable so that the
  //       Move method can make use of it to increment external pointers.
  if (_Index._x < _End._x) {
    ++_Index._x;
    _Inc = _ColumnStride;
  } else if (_Index._y < _End._y) {
    _Index._x = _Begin._x;
    ++_Index._y;
    _Inc = _LineStride;
  } else if (_Index._z < _End._z) {
    _Index._x = _Begin._x;
    _Index._y = _Begin._y;
    ++_Index._z;
    _Inc = _SliceStride;
  } else {
    _Index._x = _Begin._x;
    _Index._y = _Begin._y;
    _Index._z = _Begin._z;
    ++_Index._t;
    _Inc = _FrameStride;
  }
  if (_Next) _Next += _Inc;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ConstImageIterator::Current() const
{
  return reinterpret_cast<const VoxelType *>(_Next);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ConstImageIterator::Current(int t) const
{
  return reinterpret_cast<const VoxelType *>(_Next + t * _XYZ);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ConstImageIterator::Next()
{
  const VoxelType *current = reinterpret_cast<const VoxelType *>(_Next);
  this->operator ++();
  return current;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ConstImageIterator::Next(int t)
{
  const VoxelType *current = reinterpret_cast<const VoxelType *>(_Next + t * _XYZ);
  this->operator ++();
  return current;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType &ConstImageIterator::Value() const
{
  return *reinterpret_cast<const VoxelType *>(_Next);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType &ConstImageIterator::Value(int t) const
{
  return *reinterpret_cast<const VoxelType *>(_Next + t * _XYZ);
}

// ---------------------------------------------------------------------------
inline double ConstImageIterator::ValueAsDouble() const
{
  return Value<double>();
}

// ---------------------------------------------------------------------------
inline double ConstImageIterator::ValueAsDouble(int t) const
{
  return Value<double>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline void ConstImageIterator::Move(const VoxelType *&p) const
{
  p += _Inc / _ColumnStride;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline void ConstImageIterator::Move(VoxelType *&p) const
{
  p += _Inc / _ColumnStride;
}


} // namespace mirtk

#endif // MIRTK_ConstImageIterator_H 
