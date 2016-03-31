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

#ifndef MIRTK_GenericImageIterator_H
#define MIRTK_GenericImageIterator_H

#include "mirtk/ImageIterator.h"
#include "mirtk/ConstImageIterator.h"
#include "mirtk/ImageAttributes.h"


namespace mirtk {


/**
 * Non-const image iterator
 */
template <class VoxelType>
class GenericImageIterator : public ImageIterator
{
public:

  /// Constructor
  GenericImageIterator(const ImageAttributes &, VoxelType * = NULL);

  /// Constructor
  GenericImageIterator(GenericImage<VoxelType> &);

  /// Constructor
  GenericImageIterator(GenericImage<VoxelType> *);

  /// Copy constructor
  GenericImageIterator(const ConstImageIterator &);

  /// Assignment operator
  GenericImageIterator &operator =(const GenericImageIterator &);

  /// Destructor
  ~GenericImageIterator();

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  VoxelType *Current() const;

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  VoxelType *Current(int) const;

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  VoxelType *Next();

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  VoxelType *Next(int);

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  VoxelType &Value() const;

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  VoxelType &Value(int t) const;

  /// Get current voxel value casted to double
  double ValueAsDouble() const;

  /// Get current voxel value casted to double
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  double ValueAsDouble(int) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// ---------------------------------------------------------------------------
template <class VoxelType>
inline GenericImageIterator<VoxelType>
::GenericImageIterator(const ImageAttributes &attr, VoxelType *data)
:
  ImageIterator(attr, data, voxel_info<VoxelType>::type())
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline GenericImageIterator<VoxelType>::GenericImageIterator(GenericImage<VoxelType> &image)
:
  ImageIterator(image)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline GenericImageIterator<VoxelType>::GenericImageIterator(GenericImage<VoxelType> *image)
:
  ImageIterator(image)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline GenericImageIterator<VoxelType>::GenericImageIterator(const ConstImageIterator &other)
:
  ImageIterator(other)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline GenericImageIterator<VoxelType> &GenericImageIterator<VoxelType>
::operator =(const GenericImageIterator &rhs)
{
  ImageIterator::operator =(rhs);
  return *this;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline GenericImageIterator<VoxelType>::~GenericImageIterator()
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *GenericImageIterator<VoxelType>::Current() const
{
  return ImageIterator::Current<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *GenericImageIterator<VoxelType>::Current(int t) const
{
  return ImageIterator::Current<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *GenericImageIterator<VoxelType>::Next()
{
  return ImageIterator::Next<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType *GenericImageIterator<VoxelType>::Next(int t)
{
  return ImageIterator::Next<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &GenericImageIterator<VoxelType>::Value() const
{
  return ImageIterator::Value<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline VoxelType &GenericImageIterator<VoxelType>::Value(int t) const
{
  return ImageIterator::Value<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline double GenericImageIterator<VoxelType>::ValueAsDouble() const
{
  return static_cast<double>(Value());
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline double GenericImageIterator<VoxelType>::ValueAsDouble(int t) const
{
  return static_cast<double>(Value(t));
}


} // namespace mirtk

#endif // MIRTK_GenericImageIterator_H
