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

#ifndef MIRTK_ConnectedComponents_H
#define MIRTK_ConnectedComponents_H

#include "mirtk/ImageToImage.h"

#include "mirtk/Array.h"
#include "mirtk/OrderedSet.h"
#include "mirtk/NeighborhoodOffsets.h"


namespace mirtk {


/// Enumeration of possible orderings of connected components
enum ConnectedComponentsOrdering
{
  CC_NoOrdering,   ///< No sorting of output components
  CC_LargestFirst, ///< Sort by decreasing size
  CC_SmallestFirst ///< Sort by increasing size
};


/**
 * Label the connected components of a segmentation.
 *
 * The components are sorted by decreasing size, i.e., the first
 * component is the largest connected component.
 */
template <class TVoxel = GreyPixel>
class ConnectedComponents : public ImageToImage<TVoxel>
{
  mirtkImageFilterMacro(ConnectedComponents, TVoxel);

  // ---------------------------------------------------------------------------
  // Attributes

protected:

  /// Ordering of output components
  mirtkPublicAttributeMacro(ConnectedComponentsOrdering, Ordering);

  /// What connectivity to assume when running the filter.
  mirtkPublicAttributeMacro(ConnectivityType, Connectivity);

  /// List of voxel offsets of the neighborhood
  mirtkAttributeMacro(NeighborhoodOffsets, Offsets);

  /// Number of connected components
  mirtkReadOnlyAttributeMacro(int, NumberOfComponents);

  /// Sizes of connected components
  mirtkReadOnlyAttributeMacro(Array<int>, ComponentSize);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  ConnectedComponents(ConnectedComponentsOrdering = CC_LargestFirst, ConnectivityType = CONNECTIVITY_26);

  /// Destructor
  virtual ~ConnectedComponents();

  // ---------------------------------------------------------------------------
  // Execution

  /// Run erosion
  virtual void Run();

  /// Remove specified component from the output image
  ///
  /// This method must be called after Run(). Another execution of
  /// Run() will add the component again to the output.
  virtual void DeleteComponent(VoxelType);

  /// Size of the specified component
  ///
  /// \param[in] label Component label (1-based).
  int ComponentSize(VoxelType label) const;

protected:

  /// Initialize the filter execution
  virtual void Initialize();

  /// Finalize the filter execution
  virtual void Finalize();

};

// =============================================================================
// Ordering / string conversion
// =============================================================================

// -----------------------------------------------------------------------------
/// Convert connected components ordering to string
template <>
inline string ToString(const ConnectedComponentsOrdering &value, int w, char c, bool left)
{
  const char *str;
  switch (value) {
    case CC_NoOrdering:    str = "none";           break;
    case CC_LargestFirst:  str = "largest";  break;
    case CC_SmallestFirst: str = "smallest"; break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
/// Convert string to connected components ordering
template <>
inline bool FromString(const char *str, ConnectedComponentsOrdering &value)
{
  const string lstr = ToLower(str);
  if (lstr == "none" || lstr == "no ordering") {
    value = CC_NoOrdering;
  } else if (lstr == "largest" || lstr == "largestfirst" || lstr == "largest first") {
    value = CC_LargestFirst;
  } else if (lstr == "smallest" || lstr == "smallestfirst" || lstr == "smallest first") {
    value = CC_SmallestFirst;
  } else return false;
  return true;
}

// =============================================================================
// Inline definitions
// =============================================================================

// -----------------------------------------------------------------------------
template <class TVoxel>
inline int ConnectedComponents<TVoxel>::ComponentSize(VoxelType label) const
{
  return _ComponentSize[label-1];
}


} // namespace mirtk

#endif // MIRTK_ConnectedComponents_H
