/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
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

#ifndef MIRTK_NeighborhoodOffsets_H
#define MIRTK_NeighborhoodOffsets_H

#include "mirtk/BaseImage.h"


namespace mirtk {


/// Type of image connectivity, i.e., number of neighbors for each voxel
enum ConnectivityType
{
  CONNECTIVITY_4  = 4,
  CONNECTIVITY_6  = 6,
  CONNECTIVITY_18 = 18,
  CONNECTIVITY_26 = 26
};


/**
 * Image neighborhood offsets.
 */
class NeighborhoodOffsets : public Object
{
  mirtkObjectMacro(NeighborhoodOffsets);

  // ---------------------------------------------------------------------------
  // Types

public:

  // ---------------------------------------------------------------------------
  // Attributes

private:

  /// Image connectivity
  mirtkPublicAttributeMacro(ConnectivityType, Connectivity);

  /// Size of neighborhood
  mirtkReadOnlyAttributeMacro(int, Size);

  /// Neighborhood offsets, set to maximum possible size
  int _Offsets[26];

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  NeighborhoodOffsets();

  /// Constructor with image and connectivity specified
  NeighborhoodOffsets(const BaseImage *, ConnectivityType = CONNECTIVITY_26);

  /// Initializer with image and connectivity specified
  void Initialize(const BaseImage *, ConnectivityType = CONNECTIVITY_26);

  /// Initializer with slice dimensions and connectivity specified
  void Initialize(int, int, ConnectivityType);

  /// Destructor
  virtual ~NeighborhoodOffsets();

  /// Get i-th neighborhood offset
  ///
  /// \param[in] i Index of neighboring voxel.
  ///
  /// \returns Offset of neighboring voxel.
  int operator()(int i) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int NeighborhoodOffsets::operator()(int i) const
{
  return _Offsets[i];
}

// =============================================================================
// Connectivity type / string conversion
// =============================================================================

// -----------------------------------------------------------------------------
/// Convert image connectivity type to string
template <>
inline string ToString(const ConnectivityType &value, int w, char c, bool left)
{
  return ToString(static_cast<int>(value), w, c, left);
}

// -----------------------------------------------------------------------------
/// Convert string to image connectivity type
template <>
inline bool FromString(const char *str, ConnectivityType &value)
{
  int n;
  if (!FromString(str, n)) return false;
  switch (n) {
    case 4:  value = CONNECTIVITY_4;  break;
    case 6:  value = CONNECTIVITY_6;  break;
    case 18: value = CONNECTIVITY_18; break;
    case 26: value = CONNECTIVITY_26; break;
    default: return false;
  }
  return true;
}


} // namespace mirtk

#endif // MIRTK_NeighborhoodOffsets_H
