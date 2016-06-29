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

#include "mirtk/BoundarySegment.h"

#include "mirtk/Assert.h"
#include "mirtk/Algorithm.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void BoundarySegment::CopyAttributes(const BoundarySegment &other)
{
  _Surface     = other._Surface;
  _PointIds    = other._PointIds;
  _Index       = other._Index;
  _Selection   = other._Selection;
  _EdgeLengths = other._EdgeLengths;
  _Length      = other._Length;
}

// -----------------------------------------------------------------------------
BoundarySegment::BoundarySegment()
:
  _Length(0.)
{
}

// -----------------------------------------------------------------------------
BoundarySegment::BoundarySegment(vtkPolyData *surface, const Array<int> &ptIds)
:
  _Surface(surface),
  _PointIds(ptIds),
  _Length(0.)
{
}

// -----------------------------------------------------------------------------
BoundarySegment::BoundarySegment(const BoundarySegment &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
BoundarySegment &BoundarySegment::operator =(const BoundarySegment &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
BoundarySegment::~BoundarySegment()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void BoundarySegment::InitializeIndex()
{
  if (_Index.empty()) {
    int i = 0;
    _Index.reserve(_PointIds.size());
    for (auto it = _PointIds.begin(); it != _PointIds.end(); ++it, ++i) {
      _Index[*it] = i;
    }
  }
}

// -----------------------------------------------------------------------------
void BoundarySegment::InitializeLengths()
{
  if (!_EdgeLengths) {
    _EdgeLengths = ComputeEdgeLengths();
    _Length      = _EdgeLengths.Sum();
  } else if (_Length == 0.) {
    _Length = _EdgeLengths.Sum();
  }
}

// -----------------------------------------------------------------------------
Vector BoundarySegment::ComputeEdgeLengths() const
{
  const int npoints = NumberOfPoints();
  if (npoints == 0) return Vector();

  const class Point p0 = Point(0);

  class Point p1, p2;
  Vector      l(npoints);

  p1 = p0;
  for (int i = 1; i < npoints; ++i) {
    p2 = Point(i);
    l(i-1) = p1.Distance(p2);
    p1 = p2;
  }
  l(npoints-1) = p1.Distance(p0);

  return l;
}

// =============================================================================
// Boundary points
// =============================================================================

// -----------------------------------------------------------------------------
int BoundarySegment::Find(int ptId) const
{
  if (_Index.empty()) {
    auto it = find(_PointIds.begin(), _PointIds.end(), ptId);
    if (it == _PointIds.end()) return -1;
    return distance(_PointIds.begin(), it);
  } else {
    auto it = _Index.find(ptId);
    if (it == _Index.end()) return -1;
    return it->second;
  }
}

// -----------------------------------------------------------------------------
int BoundarySegment::FindClosestPoint(const class Point &x, double *dist2) const
{
  int    min_index = -1;
  double min_dist2 = inf;
  double d;

  for (int i = 0; i < NumberOfPoints(); ++i) {
    d = Point(i).SquaredDistance(x);
    if (d < min_dist2) {
      min_dist2 = d;
      min_index = i;
    }
  }

  if (dist2 != nullptr) *dist2 = min_dist2;
  return min_index;
}

// =============================================================================
// Selected points
// =============================================================================

// -----------------------------------------------------------------------------
void BoundarySegment::ReserveSelection(int n)
{
  _Selection.reserve(n);
}

// -----------------------------------------------------------------------------
void BoundarySegment::RemoveSelection(int i)
{
  mirtkAssert(i >= 0 && i < static_cast<int>(_Selection.size()), "valid selection index");
  _Selection.erase(_Selection.begin() + i);
}

// -----------------------------------------------------------------------------
void BoundarySegment::ClearSelection()
{
  _Selection.clear();
}

// -----------------------------------------------------------------------------
void BoundarySegment::SelectPoint(int i)
{
  mirtkAssert(i >= 0 && i < NumberOfPoints(), "valid boundary segment point index");
  if (find(_Selection.begin(), _Selection.end(), i) == _Selection.end()) {
    _Selection.push_back(i);
  }
}

// -----------------------------------------------------------------------------
void BoundarySegment::DeselectPoint(int i)
{
  auto it = find(_Selection.begin(), _Selection.end(), i);
  if (it != _Selection.end()) _Selection.erase(it);
}


} // namespace mirtk
