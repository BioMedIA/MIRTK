/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#ifndef MIRTK_SurfaceBoundary_H
#define MIRTK_SurfaceBoundary_H

#include "mirtk/Object.h"

#include "mirtk/Memory.h"
#include "mirtk/Array.h"
#include "mirtk/UnorderedMap.h"
#include "mirtk/Point.h"

#include "mirtk/EdgeTable.h"
#include "mirtk/BoundarySegment.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"


namespace mirtk {


/**
 * Boundary points and boundary segments of surface mesh
 */
class SurfaceBoundary : public Object
{
  mirtkObjectMacro(SurfaceBoundary);

  // ---------------------------------------------------------------------------
  // Types

  /// Type of shared edge table pointer
  typedef SharedPtr<class EdgeTable> EdgeTablePointer;

  /// Type of surface point ID to boundary point index map
  typedef UnorderedMap<int, int> PointIdToIndexMap;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Surface mesh
  mirtkPublicAttributeMacro(vtkSmartPointer<vtkPolyData>, Surface);

  /// Pre-computed edge table (optional)
  mirtkPublicAttributeMacro(EdgeTablePointer, EdgeTable);

  /// IDs of boundary points
  mirtkReadOnlyAttributeMacro(Array<int>, PointIds);

  /// Map from surface point ID to boundary point index
  mirtkAttributeMacro(PointIdToIndexMap, Index);

  /// Closed boundary segments
  mirtkReadOnlyAttributeMacro(Array<BoundarySegment>, Segments);

  /// Copy attributes of this class from another instance
  void CopyAttributes(const SurfaceBoundary &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Constructor
  SurfaceBoundary(vtkPolyData *, EdgeTablePointer = nullptr);

  /// Copy constructor
  SurfaceBoundary(const SurfaceBoundary &);

  /// Assignment operator
  SurfaceBoundary &operator =(const SurfaceBoundary &);

  /// Destructor
  virtual ~SurfaceBoundary();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Pre-compute maps of surface point ID to boundary (segment) point index
  void InitializeIndex();

  /// Pre-compute boundary edge lengths and length of each segment
  void InitializeLengths();

  // ---------------------------------------------------------------------------
  // Boundary points

  /// Number of surface boundary points
  int NumberOfPoints() const;

  /// Surface point index of boundary point
  ///
  /// \param[in] i SurfaceBoundary point index.
  ///
  /// \returns Surface point index of i-th boundary point.
  int PointId(int i) const;

  /// Get coordinates of boundary point
  ///
  /// \param[in]  i Index of boundary point.
  /// \param[out] p Boundary point coordinates.
  void GetPoint(int i, double p[3]) const;

  /// Get coordinates of boundary point
  ///
  /// \param[in] i Index of boundary point.
  ///
  /// \return Boundary point coordinates.
  class Point Point(int i) const;

  /// Find surface point ID in list of boundary points
  ///
  /// \param[in] ptId Surface point ID.
  ///
  /// \returns Index of surface boundary point or -1 if point is not on the boundary.
  int Find(int ptId) const;
 
  /// Check if surface point is on the boundary
  ///
  /// \param[in] ptId Surface point ID.
  ///
  /// \returns Whether a given surface point is on the boundary.
  bool Contains(int ptId) const;

  // ---------------------------------------------------------------------------
  // Boundary segments

  /// Number of boundary segments
  int NumberOfSegments() const;

  /// Get index of boundary segment that a surface point belongs to
  ///
  /// \param[in]  ptId Index of surface point.
  /// \param[out] i    Corresponding boundary segment point index.
  ///
  /// \returns Index of first encountered boundary segment.
  /// \retval -1 if point is not a boundary point.
  int FindSegment(int ptId, int *i = nullptr) const;

  /// Get index of longest boundary segment
  int FindLongestSegment() const;

  /// Get index of boundary segment which contains the most boundary points
  int FindLargestSegment() const;

  /// Get n-th boundary segment
  ///
  /// \param[in] n Index of boundary segment.
  ///
  /// \returns Reference to n-th boundary segment.
  BoundarySegment &Segment(int n);

  /// Get n-th boundary segment
  ///
  /// \param[in] n Index of boundary segment.
  ///
  /// \returns Reference to n-th boundary segment.
  const BoundarySegment &Segment(int n) const;

  /// Get longest boundary segment
  const BoundarySegment &LongestSegment() const;

  /// Get boundary segment which contains the most boundary points
  const BoundarySegment &LargestSegment() const;

  /// Number of boundary segment points
  int NumberOfPoints(int n) const;

  /// Surface point IDs of points making up a boundary point segment
  ///
  /// \param[in] n Index of boundary segment.
  ///
  /// \returns Indices of surface points making up the n-th boundary segment.
  const Array<int> &PointIds(int n) const;

  /// Surface point ID of boundary segment point
  ///
  /// \param[in] n Index of boundary segment.
  /// \param[in] i Index of boundary segment point.
  ///
  /// \returns Surface point ID of i-th point of n-th boundary segment.
  int PointId(int n, int i) const;

  /// Get indices of boundary points making up a boundary segment
  ///
  /// \param[in]  n Index of boundary segment.
  /// \param[out] i Indices of boundary points making up the n-th boundary segment.
  void PointIndices(int n, Array<int> &i) const;

  /// Get indices of boundary points making up a boundary segment
  ///
  /// \param[in] n Index of boundary segment.
  ///
  /// \returns Indices of boundary points making up the n-th boundary segment.
  Array<int> PointIndices(int n) const;

  /// Boundary point index of boundary segment point
  ///
  /// \param[in] n Index of boundary segment.
  /// \param[in] i Index of boundary segment point.
  ///
  /// \returns Boundary point index of i-th point of n-th boundary segment.
  int PointIndex(int n, int i) const;

  // ---------------------------------------------------------------------------
  // Selected points

  /// Select i-th boundary point
  ///
  /// \param[in] i Boundary point index.
  void SelectPoint(int i);

  /// Deselect i-th boundary point
  ///
  /// \param[in] i Boundary point index.
  void DeselectPoint(int i);

  // ---------------------------------------------------------------------------
  // Debug

  /// Write boundary lines to polygonal data set file
  bool Write(const char *) const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Boundary points
// =============================================================================

// -----------------------------------------------------------------------------
inline int SurfaceBoundary::NumberOfPoints() const
{
  return static_cast<int>(_PointIds.size());
}

// -----------------------------------------------------------------------------
inline int SurfaceBoundary::PointId(int i) const
{
  return _PointIds[i];
}

// -----------------------------------------------------------------------------
inline void SurfaceBoundary::GetPoint(int i, double p[3]) const
{
  const int ptId = PointId(i);
  _Surface->GetPoint(static_cast<vtkIdType>(ptId), p);
}

// -----------------------------------------------------------------------------
inline Point SurfaceBoundary::Point(int i) const
{
  double p[3];
  GetPoint(i, p);
  return p;
}

// -----------------------------------------------------------------------------
inline bool SurfaceBoundary::Contains(int ptId) const
{
  return Find(ptId) != -1;
}

// =============================================================================
// Boundary segments
// =============================================================================

// -----------------------------------------------------------------------------
inline int SurfaceBoundary::NumberOfSegments() const
{
  return static_cast<int>(_Segments.size());
}

// -----------------------------------------------------------------------------
inline BoundarySegment &SurfaceBoundary::Segment(int n)
{
  return _Segments[n];
}

// -----------------------------------------------------------------------------
inline const BoundarySegment &SurfaceBoundary::Segment(int n) const
{
  return _Segments[n];
}

// -----------------------------------------------------------------------------
inline const BoundarySegment &SurfaceBoundary::LongestSegment() const
{
  return Segment(FindLongestSegment());
}

// -----------------------------------------------------------------------------
inline const BoundarySegment &SurfaceBoundary::LargestSegment() const
{
  return Segment(FindLargestSegment());
}

// -----------------------------------------------------------------------------
inline int SurfaceBoundary::NumberOfPoints(int n) const
{
  return Segment(n).NumberOfPoints();
}

// -----------------------------------------------------------------------------
inline const Array<int> &SurfaceBoundary::PointIds(int n) const
{
  return Segment(n).PointIds();
}

// -----------------------------------------------------------------------------
inline int SurfaceBoundary::PointId(int n, int i) const
{
  return Segment(n).PointId(i);
}

// -----------------------------------------------------------------------------
inline void SurfaceBoundary::PointIndices(int n, Array<int> &i) const
{
  i = PointIds(n);
  for (auto it = i.begin(); it != i.end(); ++it) {
    *it = Find(*it);
  }
}

// -----------------------------------------------------------------------------
inline Array<int> SurfaceBoundary::PointIndices(int n) const
{
  Array<int> i;
  PointIndices(n, i);
  return i;
}

// -----------------------------------------------------------------------------
inline int SurfaceBoundary::PointIndex(int n, int i) const
{
  return Find(Segment(n).PointId(i));
}


} // namespace mirtk

#endif // MIRTK_SurfaceBoundary_H
