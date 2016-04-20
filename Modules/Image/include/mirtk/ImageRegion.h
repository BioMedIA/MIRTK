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

#ifndef MIRTK_ImageRegion_H
#define MIRTK_ImageRegion_H

#include "mirtk/Voxel.h"
#include "mirtk/Parallel.h"
#include "mirtk/Vector4D.h"
#include "mirtk/GenericImage.h"
#include "mirtk/Stream.h"


namespace mirtk {


////////////////////////////////////////////////////////////////////////////////
// Local macros (undefined again at end of this file)
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
// All get Next/Prev pointer method overloads
#define _NextPrevMethod(name, stride)                                          \
  template <class T>                                                           \
  inline const T *Next##name(const T *&p) const { return p + stride; }         \
  template <class T>                                                           \
  inline T *Next##name(T *&p) const { return p + stride; }                     \
  template <class T>                                                           \
  inline const T *Prev##name(const T *&p) const { return p - stride; }         \
  template <class T>                                                           \
  inline T *Prev##name(T *&p) const { return p - stride; }

// -----------------------------------------------------------------------------
// ToNext/ToPrev method with one image pointer argument
#define _ToMethod1(name, op, stride)                                           \
  template <class T>                                                           \
  inline void To##name(const T *&p) const { p = p op stride; }                 \
  template <class T>                                                           \
  inline void To##name(T *&p) const { p = p op stride; }

// -----------------------------------------------------------------------------
// ToNext/ToPrev method with two image pointer arguments
#define _ToMethod2(name, op, stride)                                           \
  template <class T1, class T2>                                                \
  inline void To##name(const T1 *&p1, const T2 *&p2) const                     \
    { p1 = p1 op stride; p2 = p2 op stride; }                                  \
  template <class T1, class T2>                                                \
  inline void To##name(const T1 *&p1, T2 *&p2) const                           \
    { p1 = p1 op stride; p2 = p2 op stride; }                                  \
  template <class T1, class T2>                                                \
  inline void To##name(T1 *&p1, const T2 *&p2) const                           \
    { p1 = p1 op stride; p2 = p2 op stride; }                                  \
  template <class T1, class T2>                                                \
  inline void To##name(T1 *&p1, T2 *&p2) const                                 \
    { p1 = p1 op stride; p2 = p2 op stride; }

// -----------------------------------------------------------------------------
// ToNext/ToPrev method with three image pointer arguments
#define _ToMethod3(name, op, stride)                                           \
  template <class T1, class T2, class T3>                                      \
  inline void To##name(const T1 *&p1, const T2 *&p2, const T3 *&p3) const      \
    { p1 = p1 op stride; p2 = p2 op stride; p3 = p3 op stride; }               \
  template <class T1, class T2, class T3>                                      \
  inline void To##name(const T1 *&p1, const T2 *&p2, T3 *&p3) const            \
    { p1 = p1 op stride; p2 = p2 op stride; p3 = p3 op stride; }               \
  template <class T1, class T2, class T3>                                      \
  inline void To##name(const T1 *&p1, T2 *&p2, const T3 *&p3) const            \
    { p1 = p1 op stride; p2 = p2 op stride; p3 = p3 op stride; }               \
  template <class T1, class T2, class T3>                                      \
  inline void To##name(const T1 *&p1, T2 *&p2, T3 *&p3) const                  \
    { p1 = p1 op stride; p2 = p2 op stride; p3 = p3 op stride; }               \
  template <class T1, class T2, class T3>                                      \
  inline void To##name(T1 *&p1, const T2 *&p2, const T3 *&p3) const            \
    { p1 = p1 op stride; p2 = p2 op stride; p3 = p3 op stride; }               \
  template <class T1, class T2, class T3>                                      \
  inline void To##name(T1 *&p1, const T2 *&p2, T3 *&p3) const                  \
    { p1 = p1 op stride; p2 = p2 op stride; p3 = p3 op stride; }               \
  template <class T1, class T2, class T3>                                      \
  inline void To##name(T1 *&p1, T2 *&p2, const T3 *&p3) const                  \
    { p1 = p1 op stride; p2 = p2 op stride; p3 = p3 op stride; }               \
  template <class T1, class T2, class T3>                                      \
  inline void To##name(T1 *&p1, T2 *&p2, T3 *&p3) const                        \
    { p1 = p1 op stride; p2 = p2 op stride; p3 = p3 op stride; }

// -----------------------------------------------------------------------------
// All ToNext/ToPrev method overloads
#define _ToMethod(name, stride)     \
  _ToMethod1(Next##name, +, stride) \
  _ToMethod2(Next##name, +, stride) \
  _ToMethod3(Next##name, +, stride) \
  _ToMethod1(Prev##name, -, stride) \
  _ToMethod2(Prev##name, -, stride) \
  _ToMethod3(Prev##name, -, stride) \

////////////////////////////////////////////////////////////////////////////////
// Fast image iterator declaration
////////////////////////////////////////////////////////////////////////////////

/**
 * Helper for iterating over an image region using raw image pointers
 *
 * This class makes it convenient to set a specific image region using one or
 * more of the many available setter methods and calculates the start and end
 * index of the first and last voxel in the region. Moreover, it computes how
 * many voxels have to be skipped when moving from one line of the image to
 * another, from one slice to another, or from one frame of an image sequence
 * to another. It therefore helps to iterate over the set image region using
 * a pointer to the image data. This makes the low-level iteration of an image
 * using such fast image pointers more convenient to implement.
 *
 * An instance of this class does not keep track of the current position of the
 * image pointer. It can therefore not tell when it reached the end of a line,
 * slice, or image frame/channel. This lies yet in the responsibility of the user.
 *
 * Note that a single instance can (and for the sake of speed should be) used
 * to move pointers to more than one image if needed. Herefore, all images must
 * have the same size, however. Otherwise, use different instances.
 *
 * The following example demonstrates how to iterate over a 5x5x3 neighborhood
 * of a 3D image centered at voxel (128, 128, 64).
 * \code
 * GreyImage   im(256, 256, 128);
 * ImageRegion it(image);
 *
 * // Set image region to iterate over
 * it.SetCenter(128, 128, 64);
 * it.SetRadius(  5,   5,  3);
 *
 * // Get pointer to start of image region
 * GreyPixel *p = it.GetPointerToBegin(image);
 *
 * // Use image iterator to advance pointer from one voxel to the other
 * // while iterating over the voxels of the set image region
 * GreyPixel min = MIRTK_MAX_GREY;
 * GreyPixel max = MIRTK_MIN_GREY;
 * for (int k = it.BeginZ(); k != it.EndZ(); ++k) {
 *   for (int j = it.BeginY(); j != it.EndY(); ++j) {
 *     for (int i = it.BeginX(); i != it.EndX(); ++i) {
 *       if (*p < min) min = *p;
 *       if (*p > max) max = *p;
 *       it.ToNextColumn(p);
 *     }
 *     it.ToNextLine(p);
 *   }
 *   it.ToNextSlice(p);
 * }
 * \endcode
 *
 * \sa See ConstImageIterator, ImageIterator, ConstGenericImageIterator,
 *     and GenericImageIterator for more advanced image iterators which do not
 *     require the use of up to four for-loops, but only one while-loop as these
 *     keep track of the position of the iterator themselves and decide with each
 *     increment whether to move the iterator to the next column, line, slice, or
 *     frame/channel. This, however, comes with the expense of more decisions and
 *     operations to be made per voxel and results in a noticeably slower iteration.
 *
 * \sa ForEachVoxel template functions for another convenient and fast way of
 *     performing a single operation per voxel. An example of such operation is
 *     the GetMin class which can be used as Operation template argument
 *     of the ForEachVoxel functions. These template functions and basic
 *     voxel-wise operators are defined in the mirtkVoxelFunction.h header file.
 */

class ImageRegion
{
public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  ImageRegion(const ImageAttributes &);

  /// Constructor
  ImageRegion(const ImageAttributes &, const blocked_range2d<int> &);

  /// Constructor
  ImageRegion(const ImageAttributes &, const blocked_range3d<int> &);

  /// Constructor
  ImageRegion(const BaseImage &);

  /// Constructor
  ImageRegion(const BaseImage &, const blocked_range2d<int> &);

  /// Constructor
  ImageRegion(const BaseImage &, const blocked_range3d<int> &);

  /// Constructor
  ImageRegion(const BaseImage *);

  /// Constructor
  ImageRegion(const BaseImage *, const blocked_range2d<int> &);

  /// Constructor
  ImageRegion(const BaseImage *, const blocked_range3d<int> &);

  /// Copy constructor
  ImageRegion(const ImageRegion &);

  /// Assignment operator
  ImageRegion &operator =(const ImageRegion &);

  /// Destructor
  ~ImageRegion();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Initialize iteration range and calculate strides
  ///
  /// This method determines the begin and end indices of the image region and
  /// reduces the size of the region in each dimension individually if it is
  /// outside the image domain. Use IsOutside to check if the resulting image
  /// region is completely outside the image domain and NumberOfVoxels to get
  /// the actual number of voxels that are within the overlap of the set image
  /// region and the image domain.
  ///
  /// The Initialize method is called by GetPointerToBegin or GetPointerToEnd
  /// when necessary, i.e., when the image region has been modified. If the
  /// image pointer is initialized otherwise, the Initialize method has to be
  /// called explicitly after the region has been specified. Note that the
  /// constructor initializes the iterator already. Thus, if the region is
  /// not adjusted after construction of the iterator, there is no need to
  /// re-initialize it.
  void Initialize();

  // ---------------------------------------------------------------------------
  // Image attributes

  /// Whether the image is a 3D+t image sequence
  /// (i.e., number of voxels in t dimension > 1 and dt > 0)
  bool IsImageSequence() const;

  /// Whether the image is a scalar image
  /// (i.e., number of voxels in t dimension is 1)
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

  /// Number of voxels in set image region
  int MaxNumberOfVoxels() const;

  /// Get number of channels in set image region
  int MaxNumberOfChannels() const;

  /// Get number of vector components in set image region
  int MaxNumberOfComponents() const;

  /// Get number of frames in set image region
  int MaxNumberOfFrames() const;

  /// Number of voxels in overlap of set image region and image domain
  int NumberOfVoxels() const;

  /// Get number of channels in overlap of set image region and image domain
  int NumberOfChannels() const;

  /// Get number of vector components in overlap of set image region and image domain
  int NumberOfComponents() const;

  /// Get number of frames in overlap of set image region and image domain
  int NumberOfFrames() const;

  // ---------------------------------------------------------------------------
  // Region/Neighborhood

  /// Set start of image region
  void SetStart(int, int, int = 0, int = -1);

  /// Set size of image region
  void SetSize(int);

  /// Set size of image region
  void SetSize(int, int, int = 1, int = 1);

  /// Set center of image region
  void SetCenter(int, int, int = 0, int = -1);

  /// Set radius of image region
  void SetRadius(int);

  /// Set radius of image region
  void SetRadius(int, int, int = 0, int = 0);

  /// Set 2D image region
  void SetRegion(int, int, int, int);

  /// Set 2D image region
  void SetRegion(const blocked_range2d<int> &);

  /// Set 3D image region
  void SetRegion(int, int, int, int, int, int);

  /// Set 3D image region
  void SetRegion(const blocked_range3d<int> &);

  /// Set 4D image region
  void SetRegion(int, int, int, int, int, int, int, int);

  /// Set 2D neighborhood
  void SetNeighborhood(int, int, int, int);

  /// Set 3D neighborhood
  void SetNeighborhood(int, int, int, int, int, int);

  /// Set 4D neighborhood
  void SetNeighborhood(int, int, int, int, int, int, int, int);

  /// Set temporal region
  void SetFrame(int);

  /// Set temporal region
  void SetFrame(int, int);

  /// Set temporal region
  void SetFrame(const blocked_range<int> &);

  /// Set temporal region
  void SetChannel(int);

  /// Set temporal region
  void SetChannel(int, int);

  /// Set temporal region
  void SetChannel(const blocked_range<int> &);

  /// Set temporal region
  void SetComponent(int);

  /// Set temporal region
  void SetComponent(int, int);

  /// Set temporal region
  void SetComponent(const blocked_range<int> &);

  // ---------------------------------------------------------------------------
  // Iteration range

  /// Get start of image region to iterate over in x dimension
  inline int BeginX() const { return _Begin._x; }

  /// Get start of image region to iterate over in y dimension
  inline int BeginY() const { return _Begin._y; }

  /// Get start of image region to iterate over in z dimension
  inline int BeginZ() const { return _Begin._z; }

  /// Get start of image region to iterate over in t dimension
  inline int BeginT() const { return _Begin._t; }

  /// Get end of image region to iterate over in x dimension
  inline int EndX()   const { return _End  ._x; }

  /// Get end of image region to iterate over in y dimension
  inline int EndY()   const { return _End  ._y; }

  /// Get end of image region to iterate over in z dimension
  inline int EndZ()   const { return _End  ._z; }

  /// Get end of image region to iterate over in t dimension
  inline int EndT()   const { return _End  ._t; }

  // ---------------------------------------------------------------------------
  // Strides

  /// Get stride between columns in number of voxels
  inline int ColumnStride() const { return 1; }

  /// Get stride between rows/lines in number of voxels
  inline int LineStride() const { return _LineStride; }

  /// Get stride between slices in number of voxels
  inline int SliceStride() const { return _SliceStride; }

  /// Get stride between frames in number of voxels
  inline int FrameStride() const { return _FrameStride; }

  // ---------------------------------------------------------------------------
  // Get initial image pointer

  /// Get raw image pointer to first voxel of of region
  template <class VoxelType>
  const VoxelType *GetPointerToBegin(const BaseImage &) const;

  /// Get raw image pointer to first voxel of of region
  template <class VoxelType>
  const VoxelType *GetPointerToBegin(const BaseImage *) const;

  /// Get raw image pointer to first voxel of of region
  template <class VoxelType>
  const VoxelType *GetPointerToBegin(const GenericImage<VoxelType> &) const;

  /// Get raw image pointer to first voxel of of region
  template <class VoxelType>
  VoxelType *GetPointerToBegin(GenericImage<VoxelType> &) const;

  /// Get raw image pointer to first voxel of of region
  template <class VoxelType>
  const VoxelType *GetPointerToBegin(const GenericImage<VoxelType> *) const;

  /// Get raw image pointer to first voxel of of region
  template <class VoxelType>
  VoxelType *GetPointerToBegin(GenericImage<VoxelType> *) const;

  /// Get raw image pointer to last voxel image of region
  template <class VoxelType>
  const VoxelType *GetPointerToEnd(const BaseImage &) const;

  /// Get raw image pointer to last voxel image of region
  template <class VoxelType>
  const VoxelType *GetPointerToEnd(const BaseImage *) const;

  /// Get raw image pointer to last voxel image of region
  template <class VoxelType>
  const VoxelType *GetPointerToEnd(const GenericImage<VoxelType> &) const;

  /// Get raw image pointer to last voxel image of region
  template <class VoxelType>
  VoxelType *GetPointerToEnd(GenericImage<VoxelType> &) const;

  /// Get raw image pointer to last voxel image of region
  template <class VoxelType>
  const VoxelType *GetPointerToEnd(const GenericImage<VoxelType> *) const;

  /// Get raw image pointer to last voxel image of region
  template <class VoxelType>
  VoxelType *GetPointerToEnd(GenericImage<VoxelType> *) const;

  // ---------------------------------------------------------------------------
  // Get next/previous image pointer

  // Declaration and inline definition of Next/Prev methods
  _NextPrevMethod(Column,               1)
  _NextPrevMethod(Row,        _LineStride)
  _NextPrevMethod(Line,       _LineStride)
  _NextPrevMethod(Slice,     _SliceStride)
  _NextPrevMethod(Page,      _SliceStride)
  _NextPrevMethod(Channel,   _FrameStride)
  _NextPrevMethod(Component, _FrameStride)
  _NextPrevMethod(Frame,     _FrameStride)

  // ---------------------------------------------------------------------------
  // Increment/decrement image pointer

  // Declaration and inline definition of ToNext/ToPrev methods
  _ToMethod(Column,               1)
  _ToMethod(Row,        _LineStride)
  _ToMethod(Line,       _LineStride)
  _ToMethod(Slice,     _SliceStride)
  _ToMethod(Page,      _SliceStride)
  _ToMethod(Channel,   _FrameStride)
  _ToMethod(Component, _FrameStride)
  _ToMethod(Frame,     _FrameStride)

  // ---------------------------------------------------------------------------
  // Members
protected:

  // Image attributes
  Vector4D<int> _DataSize;        ///< Size of the entire image
  bool          _IsImageSequence; ///< Whether the image is a sequence

  // Region attributes
  Vector4D<int> _Index; ///< Start of set image region
  Vector4D<int> _Size;  ///< Size of  set image region
  Vector4D<int> _Begin; ///< First index of actual image region
  Vector4D<int> _End;   ///< Last  index of actual image region

  // Strides in number of voxels
  int _LineStride;   ///< Increment at end of line  in number of voxels
  int _SliceStride;  ///< Increment at end of slice in number of voxels
  int _FrameStride;  ///< Increment at end of frame in number of voxels

};

//////////////////////////////////////////////////////////////////////////////
// Inline definitions
//////////////////////////////////////////////////////////////////////////////

// ===========================================================================
// Initialization
// ===========================================================================

// ---------------------------------------------------------------------------
inline void ImageRegion::Initialize()
{
  // Set iteration range
  _Begin = _Index;
  _End   = _Index + _Size;
  // Handle boundary conditions
  if (_Begin._x <            0) _Begin._x = 0;
  if (_Begin._y <            0) _Begin._y = 0;
  if (_Begin._z <            0) _Begin._z = 0;
  if (_Begin._t <            0) _Begin._t = 0;
  if (_End  ._x > _DataSize._x) _End  ._x = _DataSize._x;
  if (_End  ._y > _DataSize._y) _End  ._y = _DataSize._y;
  if (_End  ._z > _DataSize._z) _End  ._z = _DataSize._z;
  if (_End  ._t > _DataSize._t) _End  ._t = _DataSize._t;
  // Calculate strides in number of voxels
  const int nx = _End._x - _Begin._x;
  const int ny = _End._y - _Begin._y;
  const int nz = _End._z - _Begin._z;
  _LineStride  =  (_DataSize._x - nx);
  _SliceStride =  (_DataSize._y - ny) * _DataSize._x - nx;
  _FrameStride = ((_DataSize._z - nz) * _DataSize._y - ny) * _DataSize._x - nx;
}

// ===========================================================================
// Image attributes
// ===========================================================================

// ---------------------------------------------------------------------------
inline bool ImageRegion::IsImageSequence() const
{
  return _IsImageSequence;
}

// ---------------------------------------------------------------------------
inline bool ImageRegion::IsScalarImage() const
{
  return _DataSize._t == 1;
}

// ---------------------------------------------------------------------------
inline int ImageRegion::NumberOfImageVoxels() const
{
  return _DataSize._x * _DataSize._y * _DataSize._z * _DataSize._t;
}

// ---------------------------------------------------------------------------
inline int ImageRegion::NumberOfImageChannels() const
{
  return _IsImageSequence ? 1 : _DataSize._t;
}

// ---------------------------------------------------------------------------
inline int ImageRegion::NumberOfVectorComponents() const
{
  return NumberOfImageChannels();
}

// ---------------------------------------------------------------------------
inline int ImageRegion::NumberOfSequenceFrames() const
{
  return _IsImageSequence ? _DataSize._t : 1;
}

// ===========================================================================
// Region attributes
// ===========================================================================

// ---------------------------------------------------------------------------
inline bool ImageRegion::IsSequence() const
{
  return _Size._t > 1 && _IsImageSequence;
}

// ---------------------------------------------------------------------------
inline bool ImageRegion::IsScalar() const
{
  return _Size._t == 1;
}

// ---------------------------------------------------------------------------
inline int ImageRegion::MaxNumberOfVoxels() const
{
  return _Size._x * _Size._y * _Size._z * _Size._t;
}

// ---------------------------------------------------------------------------
inline int ImageRegion::MaxNumberOfChannels() const
{
  return IsSequence() ? 1 : _Size._t;
}

// ---------------------------------------------------------------------------
inline int ImageRegion::MaxNumberOfComponents() const
{
  return MaxNumberOfChannels();
}

// ---------------------------------------------------------------------------
inline int ImageRegion::MaxNumberOfFrames() const
{
  return IsSequence() ? _Size._t : 1;
}

// ---------------------------------------------------------------------------
inline int ImageRegion::NumberOfVoxels() const
{
  return (_End._x - _Begin._x) * (_End._y - _Begin._y) * (_End._z - _Begin._z) * (_End._t - _Begin._t);
}

// ---------------------------------------------------------------------------
inline int ImageRegion::NumberOfChannels() const
{
  return IsSequence() ? 1 : (_End._t - _Begin._t);
}

// ---------------------------------------------------------------------------
inline int ImageRegion::NumberOfComponents() const
{
  return NumberOfChannels();
}

// ---------------------------------------------------------------------------
inline int ImageRegion::NumberOfFrames() const
{
  return IsSequence() ? (_End._t - _Begin._t) : 1;
}

// ===========================================================================
// Region
// ===========================================================================

// ---------------------------------------------------------------------------
inline void ImageRegion::SetStart(int i, int j, int k, int l)
{
  _Index._x = i;
  _Index._y = j;
  _Index._z = k;
  if (l >= 0) _Index._t = l;
  _Begin._x = -1; // mark as invalid
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetCenter(int i, int j, int k, int l)
{
  _Index._x = i - _Size._x/2;
  _Index._y = j - _Size._y/2;
  _Index._z = k - _Size._z/2;
  if (l >= 0) _Index._t = l;
  _Begin._x = -1; // mark as invalid
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetSize(int s)
{
  _Size = s;
  _Begin._x = -1; // mark as invalid
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetSize(int ni, int nj, int nk, int nl)
{
  _Size._x = ni;
  _Size._y = nj;
  _Size._z = nk;
  _Size._t = nl;
  _Begin._x = -1; // mark as invalid
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetRadius(int r)
{
  _Size = 2 * r + 1;
  _Begin._x = -1; // mark as invalid
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetRadius(int ri, int rj, int rk, int rl)
{
  _Size._x = 2 * ri + 1;
  _Size._y = 2 * rj + 1;
  _Size._z = 2 * rk + 1;
  _Size._t = 2 * rl + 1;
  _Begin._x = -1; // mark as invalid
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetRegion(int bi, int bj, int ei, int ej)
{
  _Index._x = bi;
  _Index._y = bj;
  _Index._z = 0;
  _Size ._x = ei - bi;
  _Size ._y = ej - bj;
  _Size ._z = 1;
  _Begin._x = -1; // mark as invalid
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetRegion(const blocked_range2d<int> &r)
{
  SetRegion(r.cols().begin(), r.rows().begin(), r.cols().end(), r.rows().end());
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetRegion(int bi, int bj, int bk, int ei, int ej, int ek)
{
  _Index._x = bi;
  _Index._y = bj;
  _Index._z = bk;
  _Size ._x = ei - bi;
  _Size ._y = ej - bj;
  _Size ._z = ek - bk;
  _Begin._x = -1; // mark as invalid
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetRegion(int bi, int bj, int bk, int bl, int ei, int ej, int ek, int el)
{
  _Index._x = bi;
  _Index._y = bj;
  _Index._z = bk;
  _Index._t = bl;
  _Size ._x = ei - bi;
  _Size ._y = ej - bj;
  _Size ._z = ek - bk;
  _Size ._t = el - bl;
  _Begin._x = -1; // mark as invalid
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetRegion(const blocked_range3d<int> &r)
{
  SetRegion(r.cols().begin(), r.rows().begin(), r.pages().begin(),
            r.cols().end(),   r.rows().end(),   r.pages().end());
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetChannel(int l)
{
  if (!_IsImageSequence) {
    if (l < 0 || l >= _DataSize._t) {
      cerr << "ImageRegion::SetChannel: Index out of bounds" << endl;
      exit(1);
    }
    _Index._t = l;
    _Begin._x = -1; // mark as invalid
  } else if (l != 0) {
    cerr << "ImageRegion::SetChannel: Index out of bounds (an image sequence can only have one channel)" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetChannel(int bl, int el)
{
  if (!_IsImageSequence) {
    if (bl < 0 || bl >= _DataSize._t || el < 0 || el >= _DataSize._t) {
      cerr << "ImageRegion::SetChannel: Index out of bounds" << endl;
      exit(1);
    }
    _Index._t = bl;
    _Size ._t = el - bl;
    _Begin._x = -1; // mark as invalid
  } else if (bl != 0 || el != 0) {
    cerr << "ImageRegion::SetChannel: Index out of bounds (an image sequence can only have one channel)" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetChannel(const blocked_range<int> &r)
{
  SetChannel(r.begin(), r.end());
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetComponent(int l)
{
  if (!_IsImageSequence) {
    if (l < 0 || l >= _DataSize._t) {
      cerr << "ImageRegion::SetComponent: Index out of bounds" << endl;
      exit(1);
    }
    _Index._t = l;
    _Begin._x = -1; // mark as invalid
  } else if (l != 0) {
    cerr << "ImageRegion::SetComponent: Index out of bounds (an image sequence can only have one component)" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetComponent(int bl, int el)
{
  if (!_IsImageSequence) {
    if (bl < 0 || bl >= _DataSize._t || el < 0 || el >= _DataSize._t) {
      cerr << "ImageRegion::SetComponent: Index out of bounds" << endl;
      exit(1);
    }
    _Index._t = bl;
    _Size ._t = el - bl;
    _Begin._x = -1; // mark as invalid
  } else if (bl != 0 || el != 0) {
    cerr << "ImageRegion::SetComponent: Index out of bounds (an image sequence can only have one component)" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetComponent(const blocked_range<int> &r)
{
  SetComponent(r.begin(), r.end());
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetFrame(int l)
{
  if (_IsImageSequence) {
    if (l < 0 || l >= _DataSize._t) {
      cerr << "ImageRegion::SetFrame: Index out of bounds" << endl;
      exit(1);
    }
    _Index._t = l;
    _Begin._x = -1; // mark as invalid
  } else if (l != 0) {
    cerr << "ImageRegion::SetFrame: Index out of bounds (image seems to have multiple channels/vector components instead)" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetFrame(int bl, int el)
{
  if (_IsImageSequence) {
    if (bl < 0 || bl >= _DataSize._t || el < 0 || el >= _DataSize._t) {
      cerr << "ImageRegion::SetFrame: Index out of bounds" << endl;
      exit(1);
    }
    _Index._t = bl;
    _Size ._t = el - bl;
    _Begin._x = -1; // mark as invalid
  } else if (bl != 0 || el != 0) {
    cerr << "ImageRegion::SetFrame: Index out of bounds (image seems to have multiple channels/vector components instead)" << endl;
    exit(1);
  }
}

// ---------------------------------------------------------------------------
inline void ImageRegion::SetFrame(const blocked_range<int> &r)
{
  SetFrame(r.begin(), r.end());
}

// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const ImageAttributes &attr)
:
  _DataSize       (attr._x, attr._y, attr._z, attr._t),
  _IsImageSequence(attr._t > 1 && attr._dt > .0),
  _Index          (0, 0, 0, 0),
  _Size           (attr._x, attr._y, attr._z, attr._t)
{
  Initialize();
}

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const ImageAttributes &attr, const blocked_range2d<int> &region)
:
  _DataSize       (attr._x, attr._y, attr._z, attr._t),
  _IsImageSequence(attr._t > 1 && attr._dt > .0),
  _Index          (0, 0, 0, 0),
  _Size           (region.cols ().end() - region.cols ().begin(),
                   region.rows ().end() - region.rows ().begin(), 1, 1)
{
  Initialize();
}

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const ImageAttributes &attr, const blocked_range3d<int> &region)
:
  _DataSize       (attr._x, attr._y, attr._z, attr._t),
  _IsImageSequence(attr._t > 1 && attr._dt > .0),
  _Index          (0, 0, 0, 0),
  _Size           (region.cols ().end() - region.cols ().begin(),
                   region.rows ().end() - region.rows ().begin(),
                   region.pages().end() - region.pages().begin(), 1)
{
  Initialize();
}

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const BaseImage &image)
:
  _DataSize       (image.X(), image.Y(), image.Z(), image.T()),
  _IsImageSequence(image.T() > 1 && image.GetTSize() > .0),
  _Index          (0, 0, 0, 0),
  _Size           (image.X(), image.Y(), image.GetZ(), image.T())
{
  Initialize();
}

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const BaseImage &image, const blocked_range2d<int> &region)
:
  _DataSize       (image.X(), image.Y(), image.Z(), image.T()),
  _IsImageSequence(image.T() > 1 && image.GetTSize() > .0),
  _Index          (0, 0, 0, 0),
  _Size           (region.cols ().end() - region.cols ().begin(),
                   region.rows ().end() - region.rows ().begin(), 1, 1)
{
  Initialize();
}

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const BaseImage &image, const blocked_range3d<int> &region)
:
  _DataSize       (image.X(), image.Y(), image.Z(), image.T()),
  _IsImageSequence(image.T() > 1 && image.GetTSize() > .0),
  _Index          (0, 0, 0, 0),
  _Size           (region.cols ().end() - region.cols ().begin(),
                   region.rows ().end() - region.rows ().begin(),
                   region.pages().end() - region.pages().begin(), 1)
{
  Initialize();
}

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const BaseImage *image)
:
  _DataSize       (image->X(), image->Y(), image->Z(), image->T()),
  _IsImageSequence(image->T() > 1 && image->GetTSize() > .0),
  _Index          (0, 0, 0, 0),
  _Size           (image->X(), image->Y(), image->Z(), image->T())
{
  Initialize();
}

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const BaseImage *image, const blocked_range2d<int> &region)
:
  _DataSize       (image->X(), image->Y(), image->Z(), image->T()),
  _IsImageSequence(image->T() > 1 && image->GetTSize() > .0),
  _Index          (0, 0, 0, 0),
  _Size           (region.cols ().end() - region.cols ().begin(),
                   region.rows ().end() - region.rows ().begin(), 1, 1)
{
  Initialize();
}

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const BaseImage *image, const blocked_range3d<int> &region)
:
  _DataSize       (image->X(), image->Y(), image->Z(), image->T()),
  _IsImageSequence(image->T() > 1 && image->GetTSize() > .0),
  _Index          (0, 0, 0, 0),
  _Size           (region.cols ().end() - region.cols ().begin(),
                   region.rows ().end() - region.rows ().begin(),
                   region.pages().end() - region.pages().begin(), 1)
{
  Initialize();
}

// ---------------------------------------------------------------------------
inline ImageRegion::ImageRegion(const ImageRegion &other)
:
  _DataSize       (other._DataSize),
  _IsImageSequence(other._IsImageSequence),
  _Index          (other._Index),
  _Size           (other._Size),
  _Begin          (other._Begin),
  _End            (other._End),
  _LineStride     (other._LineStride),
  _SliceStride    (other._SliceStride),
  _FrameStride    (other._FrameStride)
{
}

// ---------------------------------------------------------------------------
inline ImageRegion &ImageRegion::operator =(const ImageRegion &rhs)
{
  _DataSize        = rhs._DataSize;
  _IsImageSequence = rhs._IsImageSequence;
  _Index           = rhs._Index;
  _Size            = rhs._Size;
  _Begin           = rhs._Begin;
  _End             = rhs._End;
  _LineStride      = rhs._LineStride;
  _SliceStride     = rhs._SliceStride;
  _FrameStride     = rhs._FrameStride;
  return *this;
}

// ---------------------------------------------------------------------------
inline ImageRegion::~ImageRegion()
{
}

// =============================================================================
// Initialize image pointer
// =============================================================================

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ImageRegion::GetPointerToBegin(const BaseImage &image) const
{
  if (image.GetScalarType() != voxel_info<VoxelType>::type()) {
    cerr << "ImageRegion::GetPointerToBegin: Type of image differs from VoxelType template argument" << endl;
    exit(1);
  }
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return reinterpret_cast<const VoxelType *>(image.GetScalarPointer(_Begin._x, _Begin._y, _Begin._z, _Begin._t));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ImageRegion::GetPointerToBegin(const BaseImage *image) const
{
  if (image->GetScalarType() != voxel_info<VoxelType>::type()) {
    cerr << "ImageRegion::GetPointerToBegin: Type of image differs from VoxelType template argument" << endl;
    exit(1);
  }
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return reinterpret_cast<const VoxelType *>(image->GetScalarPointer(_Begin._x, _Begin._y, _Begin._z, _Begin._t));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ImageRegion::GetPointerToBegin(const GenericImage<VoxelType> &image) const
{
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return image.GetPointerToVoxels(_Begin._x, _Begin._y, _Begin._z, _Begin._t);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *ImageRegion::GetPointerToBegin(GenericImage<VoxelType> &image) const
{
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return image.GetPointerToVoxels(_Begin._x, _Begin._y, _Begin._z, _Begin._t);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ImageRegion::GetPointerToBegin(const GenericImage<VoxelType> *image) const
{
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return image->GetPointerToVoxels(_Begin._x, _Begin._y, _Begin._z, _Begin._t);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *ImageRegion::GetPointerToBegin(GenericImage<VoxelType> *image) const
{
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return image->GetPointerToVoxels(_Begin._x, _Begin._y, _Begin._z, _Begin._t);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ImageRegion::GetPointerToEnd(const BaseImage &image) const
{
  if (image.GetScalarType() != voxel_info<VoxelType>::type()) {
    cerr << "ImageRegion::GetPointerToEnd: Type of image differs from VoxelType template argument" << endl;
    exit(1);
  }
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return reinterpret_cast<const VoxelType *>(image.GetScalarPointer(_End._x - 1, _End._y - 1, _End._z - 1, _End._t - 1));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ImageRegion::GetPointerToEnd(const BaseImage *image) const
{
  if (image->GetScalarType() != voxel_info<VoxelType>::type()) {
    cerr << "ImageRegion::GetPointerToEnd: Type of image differs from VoxelType template argument" << endl;
    exit(1);
  }
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return reinterpret_cast<const VoxelType *>(image->GetScalarPointer(_End._x - 1, _End._y - 1, _End._z - 1, _End._t - 1));
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ImageRegion::GetPointerToEnd(const GenericImage<VoxelType> &image) const
{
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return image.GetPointerToVoxels(_End._x - 1, _End._y - 1, _End._z - 1, _End._t - 1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *ImageRegion::GetPointerToEnd(GenericImage<VoxelType> &image) const
{
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return image.GetPointerToVoxels(_End._x - 1, _End._y - 1, _End._z - 1, _End._t - 1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ImageRegion::GetPointerToEnd(const GenericImage<VoxelType> *image) const
{
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return image->GetPointerToVoxels(_End._x - 1, _End._y - 1, _End._z - 1, _End._t - 1);
}

// -----------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *ImageRegion::GetPointerToEnd(GenericImage<VoxelType> *image) const
{
  if (_Begin._x < 0) const_cast<ImageRegion *>(this)->Initialize();
  return image->GetPointerToVoxels(_End._x - 1, _End._y - 1, _End._z - 1, _End._t - 1);
}

//////////////////////////////////////////////////////////////////////////////
// Undefine local macros
//////////////////////////////////////////////////////////////////////////////

#undef _NextPrevMethod
#undef _ToMethod1
#undef _ToMethod2
#undef _ToMethod3
#undef _ToMethod


} // namespace mirtk

#endif // MIRTK_ImageRegion_H
