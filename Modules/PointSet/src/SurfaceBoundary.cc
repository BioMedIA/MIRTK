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

#include "mirtk/SurfaceBoundary.h"

#include "mirtk/List.h"
#include "mirtk/Pair.h"
#include "mirtk/UnorderedSet.h"
#include "mirtk/Algorithm.h"

#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceBoundary::CopyAttributes(const SurfaceBoundary &other)
{
  _Surface   = other._Surface;
  _EdgeTable = other._EdgeTable;
  _PointIds  = other._PointIds;
  _Index     = other._Index;
  _Segments  = other._Segments;
}

// -----------------------------------------------------------------------------
SurfaceBoundary::SurfaceBoundary(vtkPolyData *surface, EdgeTablePointer edgeTable)
:
  _Surface(surface),
  _EdgeTable(edgeTable)
{
  // Get list of boundary edges
  if (!edgeTable) edgeTable = NewShared<class EdgeTable>(_Surface);
  EdgeList boundaryEdges = BoundaryEdges(_Surface, *edgeTable);

  // Get set of boundary point IDs
  UnorderedSet<int> boundaryPtIds;
  for (auto edge = boundaryEdges.begin(); edge != boundaryEdges.end(); ++edge) {
    boundaryPtIds.insert(edge->first);
    boundaryPtIds.insert(edge->second);
  }

  // Store ordered set in container with random access
  _PointIds.clear();
  _PointIds.reserve(boundaryPtIds.size());
  copy(boundaryPtIds.begin(), boundaryPtIds.end(), back_inserter(_PointIds));
  sort(_PointIds.begin(), _PointIds.end());

  // Extract boundary segments
  _Segments.clear();
  int ptId1, ptId2;
  Array<int> ptIds;
  ptIds.reserve(boundaryPtIds.size());
  while (!boundaryPtIds.empty()) {
    ptId1 = *boundaryPtIds.begin();
    ptIds.clear();
    do {
      ptIds.push_back(ptId1);
      boundaryPtIds.erase(ptId1);
      ptId2 = -1;
      for (auto edge = boundaryEdges.begin(); edge != boundaryEdges.end(); ++edge) {
        if (edge->first == ptId1) {
          if (boundaryPtIds.find(edge->second) != boundaryPtIds.end()) {
            ptId2 = edge->second;
            boundaryEdges.erase(edge);
            break;
          }
        } else if (edge->second == ptId1) {
          if (boundaryPtIds.find(edge->first) != boundaryPtIds.end()) {
            ptId2 = edge->first;
            boundaryEdges.erase(edge);
            break;
          }
        }
      }
      ptId1 = ptId2;
    } while (ptId1 != -1);
    _Segments.reserve(_Segments.size() + 1);
    _Segments.push_back(BoundarySegment(_Surface, ptIds));
  }
}

// -----------------------------------------------------------------------------
SurfaceBoundary::SurfaceBoundary(const SurfaceBoundary &other)
:
  Object(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SurfaceBoundary &SurfaceBoundary::operator =(const SurfaceBoundary &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SurfaceBoundary::~SurfaceBoundary()
{
}

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceBoundary::InitializeIndex()
{
  if (_Index.empty()) {
    int i = 0;
    _Index.reserve(_PointIds.size());
    for (auto ptId = _PointIds.begin(); ptId != _PointIds.end(); ++ptId, ++i) {
      _Index[*ptId] = i;
    }
  }
  for (auto segment = _Segments.begin(); segment != _Segments.end(); ++segment) {
    segment->InitializeIndex();
  }
}

// -----------------------------------------------------------------------------
void SurfaceBoundary::InitializeLengths()
{
  for (auto segment = _Segments.begin(); segment != _Segments.end(); ++segment) {
    segment->InitializeLengths();
  }
}

// =============================================================================
// Boundary points
// =============================================================================

// -----------------------------------------------------------------------------
int SurfaceBoundary::Find(int ptId) const
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

// =============================================================================
// Boundary segments
// =============================================================================

// -----------------------------------------------------------------------------
int SurfaceBoundary::FindSegment(int ptId, int *i) const
{
  int idx = -1, n = 0;
  for (const auto &segment : _Segments) {
    idx = segment.Find(ptId);
    if (idx != -1) {
      if (i != nullptr) *i = idx;
      return n;
    }
    ++n;
  }
  if (i != nullptr) *i = -1;
  return -1;
}

// -----------------------------------------------------------------------------
int SurfaceBoundary::FindLongestSegment() const
{
  if (_Segments.size() == 1) return 0;
  int    N = 0, n = 0;
  double L = 0., l;
  for (const auto &segment : _Segments) {
    l = segment.Length();
    if (l > L) {
      L = l;
      N = n;
    }
    ++n;
  }
  return N;
}

// -----------------------------------------------------------------------------
int SurfaceBoundary::FindLargestSegment() const
{
  int N = 0, n = 0;
  int L = 0, l;
  for (const auto &segment : _Segments) {
    l = segment.NumberOfPoints();
    if (l > L) {
      L = l;
      N = n;
    }
    ++n;
  }
  return N;
}

// =============================================================================
// Selected points
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceBoundary::SelectPoint(int i)
{
  const int ptId = PointId(i);
  for (auto && segment : _Segments) {
    auto pos = segment.Find(ptId);
    if (pos >= 0) segment.SelectPoint(pos);
  }
}

// -----------------------------------------------------------------------------
void SurfaceBoundary::DeselectPoint(int i)
{
  const int ptId = PointId(i);
  for (auto && segment : _Segments) {
    auto pos = segment.Find(ptId);
    if (pos >= 0) segment.DeselectPoint(pos);
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
bool SurfaceBoundary::Write(const char *fname) const
{
  const int npts = this->NumberOfPoints();
  vtkSmartPointer<vtkPoints> points;
  points = vtkSmartPointer<vtkPoints>::New();
  points->Allocate(npts);
  vtkSmartPointer<vtkCellArray> lines;
  lines = vtkSmartPointer<vtkCellArray>::New();
  lines->Allocate(lines->EstimateSize(npts, 2));
  for (int s = 0; s < NumberOfSegments(); ++s) {
    const auto &seg = Segment(s);
    if (seg.NumberOfPoints() > 0) {
      lines->InsertNextCell(seg.NumberOfPoints()+1);
      for (int i = 0; i < seg.NumberOfPoints(); ++i) {
        lines->InsertCellPoint(points->InsertNextPoint(Point(i)));
      }
    }
  }
  vtkSmartPointer<vtkPolyData> output;
  output = vtkSmartPointer<vtkPolyData>::New();
  output->SetPoints(points);
  output->SetLines(lines);
  return WritePolyData(fname, output);
}


} // namespace mirtk
