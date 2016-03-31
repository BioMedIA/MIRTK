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

#ifndef MIRTK_ConstGenericImageIterator_H
#define MIRTK_ConstGenericImageIterator_H

#include "mirtk/Voxel.h"
#include "mirtk/ImageAttributes.h"
#include "mirtk/GenericImage.h"
#include "mirtk/ConstImageIterator.h"


namespace mirtk {


/**
 * Const image iterator
 */
template <class VoxelType>
class ConstGenericImageIterator : public ConstImageIterator
{
public:

  /// Constructor
  ConstGenericImageIterator(const ImageAttributes &, const VoxelType * = NULL);

  /// Constructor
  ConstGenericImageIterator(GenericImage<VoxelType> &);

  /// Constructor
  ConstGenericImageIterator(GenericImage<VoxelType> *);

  /// Copy constructor
  ConstGenericImageIterator(const ConstImageIterator &);

  /// Assignment operator
  ConstGenericImageIterator &operator =(const ConstGenericImageIterator &);

  /// Destructor
  ~ConstGenericImageIterator();

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  const VoxelType *Current() const;

  /// Get pointer to current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  const VoxelType *Current(int) const;

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  const VoxelType *Next();

  /// Get pointer to current iterator position and post-increment iterator
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  const VoxelType *Next(int);

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  const VoxelType &Value() const;

  /// Get reference to voxel value at current iterator position
  ///
  /// The VoxelType template argument must match the actual scalar type of the image.
  ///
  /// \param[in] t Channel/Component/Frame offset relative to current iterator position.
  ///              For example, set iterator region to only the first channel/frame and
  ///              then access other channels/vector components/frames using this method.
  const VoxelType &Value(int t) const;

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
inline ConstGenericImageIterator<VoxelType>
::ConstGenericImageIterator(const ImageAttributes &attr, const VoxelType *data)
:
  ConstImageIterator(attr, data, voxel_info<VoxelType>::type())
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline ConstGenericImageIterator<VoxelType>::ConstGenericImageIterator(GenericImage<VoxelType> &image)
:
  ConstImageIterator(image)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline ConstGenericImageIterator<VoxelType>::ConstGenericImageIterator(GenericImage<VoxelType> *image)
:
  ConstImageIterator(image)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline ConstGenericImageIterator<VoxelType>::ConstGenericImageIterator(const ConstImageIterator &other)
:
  ConstImageIterator(other)
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline ConstGenericImageIterator<VoxelType> &ConstGenericImageIterator<VoxelType>
::operator =(const ConstGenericImageIterator &rhs)
{
  ConstImageIterator::operator =(rhs);
  return *this;
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline ConstGenericImageIterator<VoxelType>::~ConstGenericImageIterator()
{
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ConstGenericImageIterator<VoxelType>::Current() const
{
  return ConstImageIterator::Current<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ConstGenericImageIterator<VoxelType>::Current(int t) const
{
  return ConstImageIterator::Current<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ConstGenericImageIterator<VoxelType>::Next()
{
  return ConstImageIterator::Next<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType *ConstGenericImageIterator<VoxelType>::Next(int t)
{
  return ConstImageIterator::Next<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType &ConstGenericImageIterator<VoxelType>::Value() const
{
  return ConstImageIterator::Value<VoxelType>();
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline const VoxelType &ConstGenericImageIterator<VoxelType>::Value(int t) const
{
  return ConstImageIterator::Value<VoxelType>(t);
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline double ConstGenericImageIterator<VoxelType>::ValueAsDouble() const
{
  return static_cast<double>(Value());
}

// ---------------------------------------------------------------------------
template <class VoxelType>
inline double ConstGenericImageIterator<VoxelType>::ValueAsDouble(int t) const
{
  return static_cast<double>(Value(t));
}


} // namespace mirtk

#endif // MIRTK_ConstGenericImageIterator_H
