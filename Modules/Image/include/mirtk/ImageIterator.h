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

#ifndef MIRTK_ImageIterator_H
#define MIRTK_ImageIterator_H

#include "mirtk/Voxel.h"
#include "mirtk/BaseImage.h"
#include "mirtk/ConstImageIterator.h"


namespace mirtk {


/**
 * Base class of non-const image iterator
 */
class ImageIterator : public ConstImageIterator
{
public:

  // ---------------------------------------------------------------------------
  // Construction/Destruction

  /// Constructor
  ImageIterator(const ImageAttributes &, int);

  /// Constructor
  ImageIterator(const ImageAttributes &, void * = NULL, int = MIRTK_VOXEL_UNKNOWN);

  /// Constructor
  ImageIterator(BaseImage &);

  /// Constructor
  ImageIterator(BaseImage *);

  /// Copy constructor
  ImageIterator(const ConstImageIterator &);

  /// Assignment operator
  ImageIterator &operator =(const ImageIterator &);

  /// Destructor
  virtual ~ImageIterator();

  // ---------------------------------------------------------------------------
  // Iteration

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  template <class VoxelType>
  VoxelType *Current() const;

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  template <class VoxelType>
  VoxelType *Current(int) const;

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  template <class VoxelType>
  VoxelType *Next();

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  template <class VoxelType>
  VoxelType *Next(int);

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  template <class VoxelType>
  VoxelType &Value() const;

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  template <class VoxelType>
  VoxelType &Value(int t) const;

};

//////////////////////////////////////////////////////////////////////////////
// Inline definitions
//////////////////////////////////////////////////////////////////////////////

// ===========================================================================
// Construction/Destruction
// ===========================================================================

// ---------------------------------------------------------------------------
inline ImageIterator::ImageIterator(const ImageAttributes &attr, int type)
:
  ConstImageIterator(attr, type)
{
}

// ---------------------------------------------------------------------------
inline ImageIterator::ImageIterator(const ImageAttributes &attr, void *data, int type)
:
  ConstImageIterator(attr, data, type)
{
}

// ---------------------------------------------------------------------------
inline ImageIterator::ImageIterator(BaseImage &image)
:
  ConstImageIterator(image)
{
}

// ---------------------------------------------------------------------------
inline ImageIterator::ImageIterator(BaseImage *image)
:
  ConstImageIterator(image)
{
}

// ---------------------------------------------------------------------------
inline ImageIterator::ImageIterator(const ConstImageIterator &other)
:
  ConstImageIterator(other)
{
}

// ---------------------------------------------------------------------------
inline ImageIterator &ImageIterator::operator =(const ImageIterator &rhs)
{
  ConstImageIterator::operator =(rhs);
  return *this;
}

// ---------------------------------------------------------------------------
inline ImageIterator::~ImageIterator()
{
}

// ===========================================================================
// Iteration
// ===========================================================================

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *ImageIterator::Current() const
{
  return reinterpret_cast<VoxelType *>(const_cast<char *>(_Next));
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *ImageIterator::Current(int t) const
{
  return reinterpret_cast<VoxelType *>(const_cast<char *>(_Next) + t * _XYZ);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *ImageIterator::Next()
{
  VoxelType *current = reinterpret_cast<VoxelType *>(const_cast<char *>(_Next));
  this->operator ++();
  return current;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *ImageIterator::Next(int t)
{
  VoxelType *current = reinterpret_cast<VoxelType *>(const_cast<char *>(_Next) + t * _XYZ);
  this->operator ++();
  return current;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &ImageIterator::Value() const
{
  return *reinterpret_cast<VoxelType *>(const_cast<char *>(_Next));
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &ImageIterator::Value(int t) const
{
  return *reinterpret_cast<VoxelType *>(const_cast<char *>(_Next) + t * _XYZ);
}


} // namespace mirtk

#endif // MIRTK_ImageIterator_H
