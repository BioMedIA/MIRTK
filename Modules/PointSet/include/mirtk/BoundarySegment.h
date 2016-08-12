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

#ifndef MIRTK_BoundarySegment_H
#define MIRTK_BoundarySegment_H

#include "mirtk/Object.h"
#include "mirtk/Array.h"
#include "mirtk/Vector.h"
#include "mirtk/Point.h"
#include "mirtk/UnorderedMap.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"


namespace mirtk {


/**
 * Closed boundary segment of surface mesh
 *
 * A boundary segment is a closed line strip defined by an ordered list of point
 * indices with an edge between each consecutive pair of points. The boundary
 * loop is closed by the edge connecting the last point with the first point.
 * A surface mesh may have any number of boundary segments, including no boundary
 * at all in case of a closed genus-0 surface mesh which is topologically
 * equivalent to a sphere. This class represents a single boundary segment.
 */
class BoundarySegment : public Object
{
  mirtkObjectMacro(BoundarySegment);

  // ---------------------------------------------------------------------------
  // Types

  /// Type of surface point ID to boundary segment point index map
  typedef UnorderedMap<int, int> PointIdToIndexMap;

  // ---------------------------------------------------------------------------
  // Attributes

  /// Surface mesh which this boundary segment belongs to
  mirtkReadOnlyAttributeMacro(vtkSmartPointer<vtkPolyData>, Surface);

  /// IDs of boundary segment points
  mirtkReadOnlyAttributeMacro(Array<int>, PointIds);

  /// Pre-computed map from surface point ID to boundary segment point index
  mirtkAttributeMacro(PointIdToIndexMap, Index);

  /// Indices of selected points
  mirtkAttributeMacro(Array<int>, Selection);

  /// Pre-computed boundary segment edge lengths
  Vector _EdgeLengths;

  /// Length of boundary segment
  double _Length;

  /// Copy attributes of this class from another instance
  void CopyAttributes(const BoundarySegment &);

  // ---------------------------------------------------------------------------
  // Construction/Destruction

public:

  /// Default constructor
  BoundarySegment();

  /// Constructor
  ///
  /// \param[in] surface Surface mesh.
  /// \param[in] ptIds   Ordered list of IDs of surface points making up this
  ///                    boundary segment.
  BoundarySegment(vtkPolyData *surface, const Array<int> &ptIds);

  /// Copy constructor
  BoundarySegment(const BoundarySegment &);

  /// Assignment operator
  BoundarySegment &operator =(const BoundarySegment &);

  /// Destructor
  virtual ~BoundarySegment();

  // ---------------------------------------------------------------------------
  // Initialization

  /// Pre-compute map of surface point ID to boundary segment point index
  void InitializeIndex();

  // Pre-compute edge lengths and total length
  void InitializeLengths();

  /// Compute edge lengths
  Vector ComputeEdgeLengths() const;

  /// Compute length of boundary segment
  double ComputeLength() const;

  // ---------------------------------------------------------------------------
  // Boundary points

  /// Number of boundary segment points
  int NumberOfPoints() const;
 
  /// Get boundary segment point index in range [0, N)
  int IndexModuloNumberOfPoints(int i) const;

  /// Get surface point ID of i-th boundary segment point
  ///
  /// \param[in] i Index of boundary segment point.
  ///
  /// \returns Surface point ID of i-th boundary segment point.
  int PointId(int i) const;

  /// Get boundary segment point coordinates
  ///
  /// \param[in]  i Index of boundary segment point.
  /// \param[out] p Coordinates of i-th boundary segment point.
  void GetPoint(int i, double p[3]) const;

  /// Set boundary segment point coordinates
  ///
  /// \param[in] i Index of boundary segment point.
  /// \param[in] p New coordinates of i-th boundary segment point.
  void SetPoint(int i, const double p[3]);

  /// Get boundary segment point coordinates
  ///
  /// \param[in] i Index of boundary segment point.
  ///
  /// \returns Coordinates of i-th boundary segment point.
  class Point Point(int i) const;

  /// Set boundary segment point coordinates
  ///
  /// \param[in] i Index of boundary segment point.
  /// \param[in] p New coordinates of i-th boundary segment point.
  void Point(int i, const class Point &p);

  /// Find surface point in boundary segment
  ///
  /// \param[in] ptId Surface point ID.
  ///
  /// \returns Boundary segment point index or -1 if segment does not contain this point.
  int Find(int ptId) const;

  /// Find closest boundary point
  ///
  /// \param[in] x     Point coordinates.
  /// \param[in] dist2 Squared distance of closest boundary point.
  ///
  /// \returns Index of closest boundary segment point.
  int FindClosestPoint(const class Point &x, double *dist2 = nullptr) const;

  /// Whether this boundary segment contains a given surface point
  ///
  /// \param[in] ptId Surface point ID.
  ///
  /// \returns Whether surface point ID is part of this boundary segment.
  bool Contains(int ptId) const;

  // ---------------------------------------------------------------------------
  // Boundary edges

  /// Get total length of boundary segment
  double Length() const;

  /// Get lengths of boundary segment edges
  Vector EdgeLengths() const;

  /// Length of boundary segment edge
  ///
  /// \param[in] i  Index of first boundary segment point.
  /// \param[in] di Index increment +1 or -1, i.e., traversal direction.
  ///
  /// \returns Boundary segment edge lengths.
  double EdgeLength(int i, int di = +1) const;

  // ---------------------------------------------------------------------------
  // Selected points

  /// Reserve space for n selected points
  ///
  /// \param[in] n Number of points to be selected.
  void ReserveSelection(int n);

  /// Remove the i-th selected point
  ///
  /// \param[in] i Selected point index.
  void RemoveSelection(int i);

  /// Deselect all points
  void ClearSelection();

  /// Add i-th boundary segment point to selection
  ///
  /// \param[in] i Boundary segment point index.
  void SelectPoint(int i);

  /// Remove i-th boundary segment point from selection
  ///
  /// \param[in] i Boundary segment point index.
  void DeselectPoint(int i);

  /// Number of selected curve points
  int NumberOfSelectedPoints() const;

  /// Get boundary segment point index of i-th selected boundary point
  ///
  /// \param[in] i Index of selected point.
  ///
  /// \return Index of corresponding boundary segment point.
  int SelectedPointIndex(int i) const;

  /// Get surface point ID of i-th selected boundary point
  ///
  /// \param[in] i Index of selected point.
  ///
  /// \return Index of corresponding surface point.
  int SelectedPointId(int i) const;

  /// Get selected boundary segment point coordinates
  ///
  /// \param[in]  i Index of selected boundary segment point.
  /// \param[out] p Coordinates of i-th selected boundary segment point.
  void GetSelectedPoint(int i, double p[3]) const;

  /// Get selected boundary segment point coordinates
  ///
  /// \param[in] i Index of selected boundary segment point.
  ///
  /// \returns Coordinates of i-th selected boundary segment point.
  class Point SelectedPoint(int i) const;

  /// Check if boundary segment point is selected
  ///
  /// \param[in] i Index of boundary segment point.
  ///              The index is taken modulo the number of boundary
  ///              points and negative values count from the end of
  ///              the list of boundary points.
  ///
  /// \return Whether boundary segment point is selected.
  bool IsSelected(int i) const;

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
inline double BoundarySegment::ComputeLength() const
{
  if (_EdgeLengths) return _EdgeLengths.Sum();
  return ComputeEdgeLengths().Sum();
}

// =============================================================================
// Boundary points
// =============================================================================

// -----------------------------------------------------------------------------
inline int BoundarySegment::NumberOfPoints() const
{
  return static_cast<int>(_PointIds.size());
}

// -----------------------------------------------------------------------------
inline int BoundarySegment::IndexModuloNumberOfPoints(int i) const
{
  const int n = NumberOfPoints();
  if      (i <  0) i = (i + 1) % n + n - 1;
  else if (i >= n) i =  i      % n;
  return i;
}

// -----------------------------------------------------------------------------
inline int BoundarySegment::PointId(int i) const
{
  i = IndexModuloNumberOfPoints(i);
  return _PointIds[i];
}

// -----------------------------------------------------------------------------
inline void BoundarySegment::GetPoint(int i, double p[3]) const
{
  i = IndexModuloNumberOfPoints(i);
  _Surface->GetPoint(static_cast<vtkIdType>(PointId(i)), p);
}

// -----------------------------------------------------------------------------
inline void BoundarySegment::SetPoint(int i, const double p[3])
{
  i = IndexModuloNumberOfPoints(i);
  _Surface->GetPoints()->SetPoint(static_cast<vtkIdType>(PointId(i)), const_cast<double *>(p));
}

// -----------------------------------------------------------------------------
inline Point BoundarySegment::Point(int i) const
{
  double p[3];
  GetPoint(i, p);
  return p;
}

// -----------------------------------------------------------------------------
inline void BoundarySegment::Point(int i, const class Point &p)
{
  SetPoint(i, p);
}

// -----------------------------------------------------------------------------
inline bool BoundarySegment::Contains(int ptId) const
{
  return Find(ptId) != -1;
}

// =============================================================================
// Boundary edges
// =============================================================================

// -----------------------------------------------------------------------------
inline double BoundarySegment::Length() const
{
  if (_Length > 0.) return _Length;
  return ComputeLength();
}

// -----------------------------------------------------------------------------
inline Vector BoundarySegment::EdgeLengths() const
{
  if (_EdgeLengths) return _EdgeLengths;
  return ComputeEdgeLengths();
}

// -----------------------------------------------------------------------------
inline double BoundarySegment::EdgeLength(int i, int di) const
{
  if (_EdgeLengths) {
    return _EdgeLengths(IndexModuloNumberOfPoints(di < 0 ? i-1 : i));
  } else {
    return Point(i).Distance(Point(i + di));
  }
}

// =============================================================================
// Selected points
// =============================================================================

// -----------------------------------------------------------------------------
inline int BoundarySegment::NumberOfSelectedPoints() const
{
  return static_cast<int>(_Selection.size());
}

// -----------------------------------------------------------------------------
inline int BoundarySegment::SelectedPointIndex(int i) const
{
  return _Selection[i];
}

// -----------------------------------------------------------------------------
inline int BoundarySegment::SelectedPointId(int i) const
{
  return _PointIds[_Selection[i]];
}

// -----------------------------------------------------------------------------
inline void BoundarySegment::GetSelectedPoint(int i, double p[3]) const
{
  GetPoint(_Selection[i], p);
}

// -----------------------------------------------------------------------------
inline Point BoundarySegment::SelectedPoint(int i) const
{
  return Point(_Selection[i]);
}

// -----------------------------------------------------------------------------
inline bool BoundarySegment::IsSelected(int i) const
{
  for (auto it = _Selection.begin(); it != _Selection.end(); ++it) {
    if (*it == i) return true;
  }
  return false;
}


} // namespace mirtk

#endif // MIRTK_BoundarySegment_H
