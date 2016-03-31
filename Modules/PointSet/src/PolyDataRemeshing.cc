/*
 * Medical Image Registration ToolKit (MMMIRTK)
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

#include "mirtk/PolyDataRemeshing.h"

#include "mirtk/Config.h" // WINDOWS
#include "mirtk/Vtk.h"
#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Profiling.h"
#include "mirtk/PolyDataSmoothing.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/DataStatistics.h"
#include "mirtk/Transformation.h"
#include "mirtk/VtkMath.h"

#include "vtkIdList.h"
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkMergePoints.h"
#include "vtkPriorityQueue.h"
#include "vtkFloatArray.h"
#include "vtkCellDataToPointData.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void PolyDataRemeshing::CopyAttributes(const PolyDataRemeshing &other)
{
  _Transformation          = other._Transformation;
  _MinFeatureAngle         = other._MinFeatureAngle;
  _MinFeatureAngleCos      = other._MinFeatureAngleCos;
  _MaxFeatureAngle         = other._MaxFeatureAngle;
  _MaxFeatureAngleCos      = other._MaxFeatureAngleCos;
  _MinEdgeLength           = other._MinEdgeLength;
  _MinEdgeLengthSquared    = other._MinEdgeLengthSquared;
  _MaxEdgeLength           = other._MaxEdgeLength;
  _MaxEdgeLengthSquared    = other._MaxEdgeLengthSquared;
  _AdaptiveEdgeLengthArray = other._AdaptiveEdgeLengthArray;
  _MeltingOrder            = other._MeltingOrder;
  _MeltNodes               = other._MeltNodes;
  _MeltTriangles           = other._MeltTriangles;
  _NumberOfMeltedNodes     = other._NumberOfMeltedNodes;
  _NumberOfMeltedEdges     = other._NumberOfMeltedEdges;
  _NumberOfMeltedCells     = other._NumberOfMeltedCells;
  _NumberOfInversions      = other._NumberOfInversions;
  _NumberOfBisections      = other._NumberOfBisections;
  _NumberOfTrisections     = other._NumberOfTrisections;
  _NumberOfQuadsections    = other._NumberOfQuadsections;

  _InvertTrianglesSharingOneLongEdge  = other._InvertTrianglesSharingOneLongEdge;
  _InvertTrianglesToIncreaseMinHeight = other._InvertTrianglesToIncreaseMinHeight;

  // Get output point labels from _Output (copied by PolyDataFilter::Copy)
  if (_Output) {
    _OutputPointLabels  = GetArrayByCaseInsensitiveName(_Output->GetPointData(), "labels");
    _MinEdgeLengthArray = _Output->GetPointData()->GetArray("MinEdgeLength");
    _MaxEdgeLengthArray = _Output->GetPointData()->GetArray("MaxEdgeLength");
  } else {
    _OutputPointLabels  = nullptr;
    _MinEdgeLengthArray = nullptr;
    _MaxEdgeLengthArray = nullptr;
  }
}

// -----------------------------------------------------------------------------
PolyDataRemeshing::PolyDataRemeshing()
:
  _Transformation(nullptr),
  _MinFeatureAngle(180.0),
  _MinFeatureAngleCos(cos(_MinFeatureAngle * rad_per_deg)),
  _MaxFeatureAngle(180.0),
  _MaxFeatureAngleCos(cos(_MaxFeatureAngle * rad_per_deg)),
  _MinEdgeLength(.0),
  _MinEdgeLengthSquared(.0),
  _MaxEdgeLength(numeric_limits<double>::infinity()),
  _MaxEdgeLengthSquared(numeric_limits<double>::infinity()),
  _MeltingOrder(AREA),
  _MeltNodes(true),
  _MeltTriangles(false),
  _InvertTrianglesSharingOneLongEdge(false),
  _InvertTrianglesToIncreaseMinHeight(true),
  _NumberOfMeltedNodes(0),
  _NumberOfMeltedEdges(0),
  _NumberOfMeltedCells(0),
  _NumberOfInversions(0),
  _NumberOfBisections(0),
  _NumberOfTrisections(0),
  _NumberOfQuadsections(0)
{
}

// -----------------------------------------------------------------------------
PolyDataRemeshing::PolyDataRemeshing(const PolyDataRemeshing &other)
:
  PolyDataFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
PolyDataRemeshing &PolyDataRemeshing::operator =(const PolyDataRemeshing &other)
{
  if (this != &other) {
    PolyDataFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
PolyDataRemeshing::~PolyDataRemeshing()
{
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
inline void PolyDataRemeshing::GetPoint(vtkIdType ptId, double p[3]) const
{
  _Output->GetPoint(ptId, p);
  if (_Transformation) _Transformation->Transform(p[0], p[1], p[2]);
}

// -----------------------------------------------------------------------------
inline void PolyDataRemeshing::GetNormal(vtkIdType ptId, double n[3]) const
{
  // TODO: Need to compute normal of transformed surface if _Transformation != nullptr.
  //       Can be computed from normals of adjacent triangles which in turn
  //       must be computed from the transformed triangle corners (cf. GetPoint)
  _Output->GetPointData()->GetNormals()->GetTuple(ptId, n);
}

// -----------------------------------------------------------------------------
inline double PolyDataRemeshing::ComputeArea(vtkIdType cellId) const
{
  vtkIdType npts, *pts;
  _Output->GetCellPoints(cellId, npts, pts);
  if (npts != 3) return numeric_limits<double>::infinity();
  double p1[3], p2[3], p3[3];
  GetPoint(pts[0], p1);
  GetPoint(pts[1], p2);
  GetPoint(pts[2], p3);
  return vtkTriangle::TriangleArea(p1, p2, p3);
}

// -----------------------------------------------------------------------------
inline void PolyDataRemeshing::MiddlePoint(vtkIdType ptId1, vtkIdType ptId2, double p[3]) const
{
  double p2[3];
  _Output->GetPoint(ptId1, p);
  _Output->GetPoint(ptId2, p2);
  p[0] += p2[0], p[1] += p2[1], p[2] += p2[2];
  p[0] *= .5, p[1] *= .5, p[2] *= .5;
}

// -----------------------------------------------------------------------------
inline Point PolyDataRemeshing::MiddlePoint(vtkIdType ptId1, vtkIdType ptId2) const
{
  double p[3];
  MiddlePoint(ptId1, ptId2, p);
  return p;
}

// -----------------------------------------------------------------------------
inline int PolyDataRemeshing::NodeConnectivity(vtkIdType ptId) const
{
  unsigned short ncells;
  vtkIdType      *cells;
  _Output->GetPointCells(ptId, ncells, cells);
  return ncells;
}

// -----------------------------------------------------------------------------
inline void PolyDataRemeshing::DeleteCell(vtkIdType cellId)
{
  _Output->RemoveCellReference(cellId); // before marking it as deleted!
  _Output->DeleteCell(cellId);
}

// -----------------------------------------------------------------------------
inline void PolyDataRemeshing::DeleteCells(vtkIdList *cellIds)
{
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) DeleteCell(cellIds->GetId(i));
}

// -----------------------------------------------------------------------------
inline void PolyDataRemeshing::ReplaceCellPoint(vtkIdType cellId, vtkIdType oldPtId, vtkIdType newPtId)
{
  _Output->RemoveReferenceToCell(oldPtId, cellId); // NOT THREAD SAFE!
  _Output->ReplaceCellPoint(cellId, oldPtId, newPtId);
  _Output->ResizeCellList(newPtId, 1);
  _Output->AddReferenceToCell(newPtId, cellId);
}

// -----------------------------------------------------------------------------
inline vtkIdType PolyDataRemeshing::GetCellEdgeNeighbor(vtkIdType cellId, vtkIdType ptId1, vtkIdType ptId2) const
{
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  _Output->GetCellEdgeNeighbors(cellId, ptId1, ptId2, cellIds);
  vtkIdType neighborCellId = -1;
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
    if (_Output->GetCellType(cellIds->GetId(i)) != VTK_EMPTY_CELL) {
      if (neighborCellId != -1) return -1; // should not happen
      neighborCellId = cellIds->GetId(i);
    }
  }
  return neighborCellId;
}

// -----------------------------------------------------------------------------
inline void PolyDataRemeshing
::GetCellPointNeighbors(vtkIdType cellId, vtkIdType ptId, vtkIdList *ptIds) const
{
  ptIds->Reset();
  unsigned short ncells;
  vtkIdType *cells, *pts, npts;
  _Output->GetPointCells(ptId, ncells, cells);
  for (unsigned short i = 0; i < ncells; ++i) {
    if (cells[i] != cellId) {
      _Output->GetCellPoints(cells[i], npts, pts);
      for (vtkIdType j = 0; j < npts; ++j) {
        if (pts[j] != ptId) ptIds->InsertUniqueId(pts[j]);
      }
    }
  }
}

// -----------------------------------------------------------------------------
inline vtkIdType PolyDataRemeshing
::GetCellEdgeNeighborPoint(vtkIdType cellId, vtkIdType ptId1, vtkIdType ptId2, bool mergeTriples)
{
  unsigned short ncells;
  vtkIdType ptId3, npts, *pts, *cells;

  // Get other cell adjacent to this edge
  const vtkIdType neighborCellId = GetCellEdgeNeighbor(cellId, ptId1, ptId2);
  if (neighborCellId == -1) return -1;

  // Get points which are adjacent to both edge points
  vtkSmartPointer<vtkIdList> ptIds1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> ptIds2 = vtkSmartPointer<vtkIdList>::New();

  GetCellPointNeighbors(cellId, ptId1, ptIds1);
  GetCellPointNeighbors(cellId, ptId2, ptIds2);

  ptIds1->IntersectWith(ptIds2);

  // Intersection contains all points connected to both edge points.
  // This is the third point of cell with ID cellId, and all points on
  // the opposite side of the edge, i.e., the third point defining the
  // neighboring cell. Other points included belong to nodes that are
  // connected to triangles that are fully contained within the triangle
  // defined by the edge points and one of the adjacent points (furthest
  // away from the edge when surface has no self-intersections). These can
  // be removed starting at the third point defining the neighboring cell
  // which in this case has node connectivity three. The next third point
  // defining the new (merged) neighboring cell then has node connectivity
  // three and it's adjacent triangles can be merged as well. This is
  // continued until no more point with node connectivity three remains and
  // the ptIds1 list has only two points left.
  if (ptIds1->GetNumberOfIds() < 2) {
    return -1; // Boundary point, cannot happen in case of closed surface mesh
  }

  if (mergeTriples) {
    while (ptIds1->GetNumberOfIds() > 2) {
      int idx = -1; // index of cell that becomes union of node adjacent cells
      for (vtkIdType i = 0; i < ptIds1->GetNumberOfIds(); ++i) {
        ptId3 = ptIds1->GetId(i);
        _Output->GetPointCells(ptId3, ncells, cells);
        if (ncells != 3) continue; // node connectivity must be three
        for (unsigned short j = 0; j < ncells; ++j) {
          _Output->GetCellPoints(cells[j], npts, pts);
          for (vtkIdType k = 0; k < npts; ++k) {
            if (pts[k] != ptId1 && pts[k] != ptId2 && pts[k] != ptId3) {
              for (unsigned short l = 0; l < ncells; ++l) {
                if (cells[l] == cellId || cells[l] == neighborCellId) {
                  ReplaceCellPoint(cells[l], ptId3, pts[k]);
                  idx = l;
                  break;
                }
              }
              if (idx != -1) {
                for (unsigned short l = 0; l < ncells; ++l) {
                  if (l != idx) DeleteCell(cells[l]);
                }
                ptIds1->DeleteId(ptId3);
                ptIds1->InsertUniqueId(pts[k]); // should be in list already
                if (_MeltingQueue) {
                  for (unsigned short l = 0; l < ncells; ++l) {
                    _MeltingQueue->DeleteId(cells[l]);
                    if (l == idx && cells[l] != cellId) {
                      double priority = MeltingPriority(cells[idx]);
                      if (!IsInf(priority)) {
                        _MeltingQueue->Insert(priority, cells[idx]);
                      }
                    }
                  }
                }
                ++_NumberOfMeltedNodes;
              }
              break;
            }
          }
          break;
        }
        break;
      }
      if (idx == -1) break;
    };
  }

  // Intersection should now contain the other point of this cell
  // and the third point belonging to the cell edge neighbor
  if (ptIds1->GetNumberOfIds() != 2) return -1;

  // Get other cell edge neighbor point
  _Output->GetCellPoints(neighborCellId, npts, pts);
  if (npts == 3) {
    for (vtkIdType i = 0; i < npts; ++i) {
      if (pts[i] != ptId1 && pts[i] != ptId2) {
        return pts[i];
      }
    }
  }
  return -1;
}

// -----------------------------------------------------------------------------
inline double PolyDataRemeshing::MeltingPriority(vtkIdType cellId) const
{
  double priority = numeric_limits<double>::infinity();
  switch (_MeltingOrder) {
    case INDEX: {
      priority = double(cellId);
    } break;
    case AREA: {
      priority = ComputeArea(cellId);
    } break;
    case SHORTEST_EDGE: {
      vtkIdType npts, *pts;
      _Output->GetCellPoints(cellId, npts, pts);
      if (npts == 3) {
        double p1[3], p2[3], p3[3];
        GetPoint(pts[0], p1);
        GetPoint(pts[1], p2);
        GetPoint(pts[2], p3);
        priority = min(min(vtkMath::Distance2BetweenPoints(p1, p2),
                           vtkMath::Distance2BetweenPoints(p1, p3)),
                           vtkMath::Distance2BetweenPoints(p2, p3));
      }
    } break;
  }
  return priority;
}

// -----------------------------------------------------------------------------
inline void PolyDataRemeshing
::InterpolatePointData(vtkIdType newId, vtkIdType ptId1, vtkIdType ptId2)
{
  vtkPointData *pd = _Output->GetPointData();
  if (pd == nullptr) return;
  if (_OutputPointLabels) {
    double label = _OutputPointLabels->GetComponent(ptId1, 0);
    pd->InterpolateEdge(pd, newId, ptId1, ptId2, .5);
    _OutputPointLabels->SetComponent(newId, 0, label);
  } else {
    pd->InterpolateEdge(pd, newId, ptId1, ptId2, .5);
  }
}

// -----------------------------------------------------------------------------
inline void PolyDataRemeshing
::InterpolatePointData(vtkIdType newId, vtkIdList *ptIds, double *weights)
{
  vtkPointData *pd = _Output->GetPointData();
  if (pd == nullptr) return;
  if (_OutputPointLabels) {
    double label = _OutputPointLabels->GetComponent(ptIds->GetId(0), 0);
    pd->InterpolatePoint(pd, newId, ptIds, weights);
    _OutputPointLabels->SetComponent(newId, 0, label);
  } else {
    pd->InterpolatePoint(pd, newId, ptIds, weights);
  }
}

// -----------------------------------------------------------------------------
inline double PolyDataRemeshing
::SquaredMinEdgeLength(vtkIdType ptId1, vtkIdType ptId2) const
{
  if (_MinEdgeLengthArray) {
    double l = .5 * (_MinEdgeLengthArray->GetComponent(ptId1, 0) +
                     _MinEdgeLengthArray->GetComponent(ptId2, 0));
    return l * l;
  }
  return _MinEdgeLengthSquared;
}

// -----------------------------------------------------------------------------
inline double PolyDataRemeshing
::SquaredMaxEdgeLength(vtkIdType ptId1, vtkIdType ptId2) const
{
  if (_MaxEdgeLengthArray) {
    double l = .5 * (_MaxEdgeLengthArray->GetComponent(ptId1, 0) +
                     _MaxEdgeLengthArray->GetComponent(ptId2, 0));
    return l * l;
  }
  return _MaxEdgeLengthSquared;
}

// =============================================================================
// Local remeshing operations
// =============================================================================

// -----------------------------------------------------------------------------
bool PolyDataRemeshing
::MeltEdge(vtkIdType cellId, vtkIdType ptId1, vtkIdType ptId2, vtkIdList *cellIds)
{
  // Check/resolve node connectivity of adjacent points
  vtkIdType neighborPtId = GetCellEdgeNeighborPoint(cellId, ptId1, ptId2, _MeltNodes);
  if (neighborPtId == -1) return false;

  // Get edge neighbor cell
  vtkIdType neighborCellId = GetCellEdgeNeighbor(cellId, ptId1, ptId2);
  if (neighborCellId == -1) return false;

  // Interpolate point data
  InterpolatePointData(ptId1, ptId1, ptId2);

  // Move first point to edge middlepoint
  Point m = MiddlePoint(ptId1, ptId2);
  _Output->GetPoints()->SetPoint(ptId1, m._x, m._y, m._z);
  _Output->GetPoints()->Modified();

  // Replace second point in (remaining) adjacent cells by first point
  _Output->GetPointCells(ptId2, cellIds);
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
    ReplaceCellPoint(cellIds->GetId(i), ptId2, ptId1);
  }

  // Mark cell and its edge neighbor as deleted
  DeleteCell(cellId);
  DeleteCell(neighborCellId);

  // Update priority queue
  _Output->GetPointCells(ptId1, cellIds);
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
    _MeltingQueue->DeleteId(cellIds->GetId(i));
    double priority = MeltingPriority(cellIds->GetId(i));
    if (!IsInf(priority)) _MeltingQueue->Insert(priority, cellIds->GetId(i));
  }
  _MeltingQueue->DeleteId(neighborCellId);

  ++_NumberOfMeltedEdges;
  return true;
}

// -----------------------------------------------------------------------------
bool PolyDataRemeshing::MeltTriangle(vtkIdType cellId, vtkIdList *cellIds)
{
  // Get triangle corners
  vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
  _Output->GetCellPoints(cellId, ptIds);

  // Get points opposite to cell edges and resolve connectivity of 3 if possible
  vtkIdType ptId1 = GetCellEdgeNeighborPoint(cellId, ptIds->GetId(0), ptIds->GetId(1), _MeltNodes);
  if (ptId1 == -1) return false;
  vtkIdType ptId2 = GetCellEdgeNeighborPoint(cellId, ptIds->GetId(1), ptIds->GetId(2), _MeltNodes);
  if (ptId2 == -1) return false;
  vtkIdType ptId3 = GetCellEdgeNeighborPoint(cellId, ptIds->GetId(2), ptIds->GetId(0), _MeltNodes);
  if (ptId3 == -1) return false;

  // Get adjacent triangles
  vtkIdType neighborCellId1 = GetCellEdgeNeighbor(cellId, ptIds->GetId(0), ptIds->GetId(1));
  if (neighborCellId1 == -1) return false;
  vtkIdType neighborCellId2 = GetCellEdgeNeighbor(cellId, ptIds->GetId(1), ptIds->GetId(2));
  if (neighborCellId2 == -1) return false;
  vtkIdType neighborCellId3 = GetCellEdgeNeighbor(cellId, ptIds->GetId(2), ptIds->GetId(0));
  if (neighborCellId3 == -1) return false;

  // Get triangle center point and interpolation weights
  double pcoords[3], c[3], weights[3];
  vtkCell *cell = _Output->GetCell(cellId);
  int subId = cell->GetParametricCenter(pcoords);
  cell->EvaluateLocation(subId, pcoords, c, weights);

  // Interpolate point data
  InterpolatePointData(ptIds->GetId(0), ptIds, weights);

  // Move first point of this triangle to its center
  _Output->GetPoints()->SetPoint(ptIds->GetId(0), c);

  // Link (non-deleted) cells sharing triangle corners to center point
  // instead of any of the other two deleted points
  _Output->GetPointCells(ptIds->GetId(1), cellIds);
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
    ReplaceCellPoint(cellIds->GetId(i), ptIds->GetId(1), ptIds->GetId(0));
  }

  _Output->GetPointCells(ptIds->GetId(2), cellIds);
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
    ReplaceCellPoint(cellIds->GetId(i), ptIds->GetId(2), ptIds->GetId(0));
  }

  // Delete this and adjacent triangles
  DeleteCell(cellId);
  DeleteCell(neighborCellId1);
  DeleteCell(neighborCellId2);
  DeleteCell(neighborCellId3);

  // Update priority queue
  _Output->GetPointCells(ptIds->GetId(0), cellIds);
  for (vtkIdType i = 0, cellId; i < cellIds->GetNumberOfIds(); ++i) {
    cellId = cellIds->GetId(i);
    _MeltingQueue->DeleteId(cellId);
    double priority = MeltingPriority(cellId);
    if (!IsInf(priority)) _MeltingQueue->Insert(priority, cellId);
  }
  _MeltingQueue->DeleteId(neighborCellId1);
  _MeltingQueue->DeleteId(neighborCellId2);
  _MeltingQueue->DeleteId(neighborCellId3);

  ++_NumberOfMeltedCells;
  return true;
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::MeltingOfCells()
{
  MIRTK_START_TIMING();

  int       melt[3], i, j;
  double    p1[3], p2[3], p3[3], length2[3], n1[3], n2[3], n3[3];
  vtkIdType cellId, npts, *pts;

  // Cell ID list shared by melting operation functions to save reallocation
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  cellIds->Allocate(20);

  // Initialize priority queue of cells to process
  _MeltingQueue = vtkSmartPointer<vtkPriorityQueue>::New();
  for (cellId = 0; cellId < _Output->GetNumberOfCells(); ++cellId) {
    double priority = MeltingPriority(cellId);
    if (!IsInf(priority)) _MeltingQueue->Insert(priority, cellId);
  }

  int num_triangle_melting_attempts = 0;

  // Process cells in order determined by priority (e.g., smallest triangle first)
  while (_MeltingQueue->GetNumberOfItems() > 0) {
    cellId = _MeltingQueue->Pop();

    _Output->GetCellPoints(cellId, npts, pts);
    if (npts == 0) continue; // cell marked as deleted (i.e., VTK_EMPTY_CELL)
    mirtkAssert(npts == 3, "surface is triangulated");

    // Get (transformed) point coordinates
    GetPoint(pts[0], p1);
    GetPoint(pts[1], p2);
    GetPoint(pts[2], p3);

    length2[0] = vtkMath::Distance2BetweenPoints(p1, p2);
    length2[1] = vtkMath::Distance2BetweenPoints(p2, p3);
    length2[2] = vtkMath::Distance2BetweenPoints(p3, p1);

    // If triangle-melting is allowed
    if (_MeltTriangles) {

      // Determine which edges are too short
      melt[0] = int(length2[0] < SquaredMinEdgeLength(pts[0], pts[1]));
      melt[1] = int(length2[1] < SquaredMinEdgeLength(pts[1], pts[2]));
      melt[2] = int(length2[2] < SquaredMinEdgeLength(pts[2], pts[0]));

      // Do not melt feature edges (otherwise remeshing may smooth surface too much)
      if ((melt[0] || melt[1] || melt[2]) && _MinFeatureAngle < 180.0) {
        GetNormal(pts[0], n1);
        GetNormal(pts[1], n2);
        GetNormal(pts[2], n3);
        if (melt[0]) melt[0] = int(1.0 - vtkMath::Dot(n1, n2) < _MinFeatureAngleCos);
        if (melt[1]) melt[1] = int(1.0 - vtkMath::Dot(n2, n3) < _MinFeatureAngleCos);
        if (melt[2]) melt[2] = int(1.0 - vtkMath::Dot(n3, n1) < _MinFeatureAngleCos);
      }

      // Perform either edge-melting, triangle-melting, or no operation
      switch (melt[0] + melt[1] + melt[2]) {
        case 1: {
          if      (melt[0]) MeltEdge(cellId, pts[0], pts[1], cellIds);
          else if (melt[1]) MeltEdge(cellId, pts[1], pts[2], cellIds);
          else              MeltEdge(cellId, pts[2], pts[0], cellIds);
        } break;
        case 3: {
          ++num_triangle_melting_attempts;
          if (!MeltTriangle(cellId, cellIds)) {
            i = 0;
            if (length2[1] < length2[0]) i = 1;
            if (length2[2] < length2[i]) i = 2;
            j = (i + 1) % 3;
            MeltEdge(cellId, pts[i], pts[j], cellIds);
          }
        } break;
      }

    // If only individual edges are allowed to be melted
    } else {

      // Determine shortest edge
      i = 0;
      if (length2[1] < length2[0]) i = 1;
      if (length2[2] < length2[i]) i = 2;
      j = (i + 1) % 3;

      // Melt edge if it is too short and not a feature edge
      if (length2[i] < SquaredMinEdgeLength(pts[i], pts[j])) {
        if (_MinFeatureAngle < 180.0) {
          GetNormal(pts[i], n1);
          GetNormal(pts[j], n2);
          if (1.0 - vtkMath::Dot(n1, n2) < _MinFeatureAngleCos) {
            MeltEdge(cellId, pts[i], pts[j], cellIds);
          }
        } else {
          MeltEdge(cellId, pts[i], pts[j], cellIds);
        }
      }

    }
  }

  _MeltingQueue = nullptr;
  MIRTK_DEBUG_TIMING(3, "melting edges and triangles");
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::MeltingOfNodes()
{
  MIRTK_START_TIMING();

  unsigned short             ncells;
  vtkIdType                  ptIdx, npts, *pts, *cells;
  vtkSmartPointer<vtkIdList> ptIds1, ptIds2;
  ptIds1 = vtkSmartPointer<vtkIdList>::New();
  ptIds2 = vtkSmartPointer<vtkIdList>::New();

  bool changed = true;
  for (int iter = 0; iter < 10 && changed; ++iter) {
    changed = false;
    for (vtkIdType ptId = 0; ptId < _Output->GetNumberOfPoints(); ++ptId) {
      _Output->GetPointCells(ptId, ncells, cells);
      switch (ncells) {
        case 1: case 2: {
          for (unsigned short i = 0; i < ncells; ++i) {
            DeleteCell(cells[i]);
          }
          ++_NumberOfMeltedNodes;
          changed = true;
        } break;
        case 3: {
          ptIds1->Reset();
          ptIds2->Reset();
          for (unsigned short i = 0; i < ncells; ++i) {
            _Output->GetCellPoints(cells[i], npts, pts);
            if (i == 0) {
              ptIds1->Allocate(npts);
              for (vtkIdType j = 0; j < npts; ++j) {
                ptIds1->InsertNextId(pts[j]);
              }
            }
            for (vtkIdType j = 0; j < npts; ++j) {
              ptIds2->InsertUniqueId(pts[j]);
            }
          }
          ptIds2->DeleteId(ptId);
          if (ptIds2->GetNumberOfIds() != 3) continue;
          for (ptIdx = 0; ptIdx < ptIds2->GetNumberOfIds(); ++ptIdx) {
            if (ptIds1->IsId(ptIds2->GetId(ptIdx)) == -1) break;
          }
          if (ptIdx == ptIds2->GetNumberOfIds()) continue;
          ReplaceCellPoint(cells[0], ptId, ptIds2->GetId(ptIdx));
          DeleteCell(cells[1]);
          DeleteCell(cells[2]);
          ++_NumberOfMeltedNodes;
          changed = true;
        } break;
      }
    }
  }

  MIRTK_DEBUG_TIMING(3, "melting nodes with connectivity less equal 3");
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing
::Bisect(vtkIdType cellId, vtkIdType ptId1, vtkIdType ptId2, vtkIdType ptId3,
         vtkPoints *points, vtkCellArray *polys)
{
  const Point midPoint = MiddlePoint(ptId1, ptId2);
  const vtkIdType midPtId  = points->InsertNextPoint(midPoint._x, midPoint._y, midPoint._z);
  vtkIdType       pts[3];

  InterpolatePointData(midPtId, ptId1, ptId2);

  pts[0] = ptId1;
  pts[1] = midPtId;
  pts[2] = ptId3;
  polys->InsertNextCell(3, pts);

  pts[0] = midPtId;
  pts[1] = ptId2;
  pts[2] = ptId3;
  polys->InsertNextCell(3, pts);

  ++_NumberOfBisections;
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing
::Trisect(vtkIdType cellId, vtkIdType ptId1, vtkIdType ptId2, vtkIdType ptId3,
          vtkPoints *points, vtkCellArray *polys)
{
  const Point midPoint1 = MiddlePoint(ptId1, ptId2);
  const Point midPoint2 = MiddlePoint(ptId2, ptId3);
  const vtkIdType midPtId1  = points->InsertNextPoint(midPoint1._x, midPoint1._y, midPoint1._z);
  const vtkIdType midPtId2  = points->InsertNextPoint(midPoint2._x, midPoint2._y, midPoint2._z);
  vtkIdType       pts[3];

  InterpolatePointData(midPtId1, ptId1, ptId2);
  InterpolatePointData(midPtId2, ptId2, ptId3);

  pts[0] = ptId1;
  pts[1] = midPtId1;
  pts[2] = ptId3;
  polys->InsertNextCell(3, pts);

  pts[0] = midPtId1;
  pts[1] = ptId2;
  pts[2] = midPtId2;
  polys->InsertNextCell(3, pts);

  pts[0] = midPtId2;
  pts[1] = ptId3;
  pts[2] = midPtId1;
  polys->InsertNextCell(3, pts);

  ++_NumberOfTrisections;
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing
::Quadsect(vtkIdType cellId, vtkIdType ptId1, vtkIdType ptId2, vtkIdType ptId3,
           vtkPoints *points, vtkCellArray *polys)
{
  const Point midPoint1 = MiddlePoint(ptId1, ptId2);
  const Point midPoint2 = MiddlePoint(ptId2, ptId3);
  const Point midPoint3 = MiddlePoint(ptId3, ptId1);
  const vtkIdType midPtId1  = points->InsertNextPoint(midPoint1._x, midPoint1._y, midPoint1._z);
  const vtkIdType midPtId2  = points->InsertNextPoint(midPoint2._x, midPoint2._y, midPoint2._z);
  const vtkIdType midPtId3  = points->InsertNextPoint(midPoint3._x, midPoint3._y, midPoint3._z);
  vtkIdType       pts[3];

  InterpolatePointData(midPtId1, ptId1, ptId2);
  InterpolatePointData(midPtId2, ptId2, ptId3);
  InterpolatePointData(midPtId3, ptId3, ptId1);

  pts[0] = ptId1;
  pts[1] = midPtId1;
  pts[2] = midPtId3;
  polys->InsertNextCell(3, pts);

  pts[0] = midPtId1;
  pts[1] = ptId2;
  pts[2] = midPtId2;
  polys->InsertNextCell(3, pts);

  pts[0] = midPtId1;
  pts[1] = midPtId2;
  pts[2] = midPtId3;
  polys->InsertNextCell(3, pts);

  pts[0] = midPtId2;
  pts[1] = ptId3;
  pts[2] = midPtId3;
  polys->InsertNextCell(3, pts);

  ++_NumberOfQuadsections;
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void PolyDataRemeshing::Initialize()
{
  // Initialize base class
  PolyDataFilter::Initialize();

  // Convert input cell edge length arrays to point data arrays
  vtkSmartPointer<vtkPolyData> input = _Input;
  if (!_AdaptiveEdgeLengthArray && (_MinCellEdgeLengthArray || _MaxCellEdgeLengthArray)) {
    input = vtkSmartPointer<vtkPolyData>::New();
    input->ShallowCopy(_Input);
    input->GetCellData()->Initialize();
    if (_MinCellEdgeLengthArray) input->GetCellData()->AddArray(_MinCellEdgeLengthArray);
    if (_MaxCellEdgeLengthArray) input->GetCellData()->AddArray(_MaxCellEdgeLengthArray);
    vtkNew<vtkCellDataToPointData> c2p;
    SetVTKInput(c2p, input);
    c2p->PassCellDataOff();
    c2p->Update();
    input = c2p->GetPolyDataOutput();

    static int iter = 0; ++iter;
    char fname[64];
    #ifdef WINDOWS
      sprintf_s(fname, 64, "remesh_input_%03d.vtp", iter);
    #else
      sprintf(fname, "remesh_input_%03d.vtp", iter);
    #endif
    WritePolyData(fname, input);
  }

  // Triangulate input surface if necessary and/or make shallow copy
  // as we will modify the point data attribute set of the input
  _TriangulatedInput = Triangulate(input);

  // Compute point normals if needed
  _MinFeatureAngle    = max(.0, min(_MinFeatureAngle, 180.0));
  _MaxFeatureAngle    = max(.0, min(_MaxFeatureAngle, 180.0));
  _MinFeatureAngleCos = 1.0 - cos(_MinFeatureAngle * rad_per_deg);
  _MaxFeatureAngleCos = 1.0 - cos(_MaxFeatureAngle * rad_per_deg);

  if (_MinFeatureAngle < 180.0 || _MaxFeatureAngle < 180.0) {
    vtkNew<vtkPolyDataNormals> calc_normals;
    SetVTKInput(calc_normals, _TriangulatedInput);
    calc_normals->ComputeCellNormalsOff();
    calc_normals->ComputePointNormalsOn();
    calc_normals->AutoOrientNormalsOn();
    calc_normals->SplittingOff();
    calc_normals->Update();
    vtkDataArray *normals = calc_normals->GetOutput()->GetPointData()->GetNormals();
    _TriangulatedInput->GetPointData()->SetNormals(normals);
  }

  // Initialize adaptive edge length range
  _MinEdgeLengthSquared = _MinEdgeLength * _MinEdgeLength;
  _MaxEdgeLengthSquared = _MaxEdgeLength * _MaxEdgeLength;
  this->InitializeEdgeLengthRange();

  // Initialize output
  _Output = vtkSmartPointer<vtkPolyData>::New();
  _Output->DeepCopy(_TriangulatedInput);
  vtkPointData *inputPD  = _TriangulatedInput->GetPointData();
  vtkPointData *outputPD = _Output->GetPointData();
  vtkCellData  *outputCD = _Output->GetCellData();
  outputCD->Initialize(); // TODO: Interpolate cell data during remeshing
  outputPD->InterpolateAllocate(inputPD, _TriangulatedInput->GetNumberOfPoints());
  for (vtkIdType ptId = 0; ptId < _TriangulatedInput->GetNumberOfPoints(); ++ptId) {
    outputPD->CopyData(inputPD, ptId, ptId);
  }

  // Set interpolation of discrete labels to nearest neighbor
  int idx;
  _OutputPointLabels = GetArrayByCaseInsensitiveName(_Output->GetPointData(), "labels", &idx);
  if (_OutputPointLabels) _Output->GetPointData()->SetCopyAttribute(idx, 2);

  // Build links
  _Output->BuildLinks();

  // Reset counters
  _NumberOfMeltedNodes  = 0;
  _NumberOfMeltedEdges  = 0;
  _NumberOfMeltedCells  = 0;
  _NumberOfInversions   = 0;
  _NumberOfBisections   = 0;
  _NumberOfTrisections  = 0;
  _NumberOfQuadsections = 0;
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::InitializeEdgeLengthRange()
{
  // Either derive per-node edge length range from scalar attribute (e.g., curvedness)
  if (_AdaptiveEdgeLengthArray) {

    // Remove previous edge length range arrays
    _TriangulatedInput->GetPointData()->RemoveArray("MinEdgeLength");
    _TriangulatedInput->GetPointData()->RemoveArray("MaxEdgeLength");

    // Global edge length range
    double min_length = _MinEdgeLength;
    double max_length = _MaxEdgeLength;
    if (min_length < .0)   min_length = .0;
    if (IsInf(max_length)) max_length = max(1.0, 10.0 * min_length);

    // Parameters of logistic functions
    const double k  = 8.0;
    const double x1 = 0.4;
    const double x2 = 0.6;

    // Get robust range of scalar values
    using data::statistic::Percentile;
    const double min    = Percentile::Calculate( 5, _AdaptiveEdgeLengthArray);
    const double max    = Percentile::Calculate(95, _AdaptiveEdgeLengthArray);
    const double scale  = 1.0 / (max - min);
    const double offset = - scale * min;

    // Allocate edge length arrays
    if (!_MinEdgeLengthArray) {
      _MinEdgeLengthArray = vtkSmartPointer<vtkFloatArray>::New();
      _MinEdgeLengthArray->SetName("MinEdgeLength");
      _MinEdgeLengthArray->SetNumberOfComponents(1);
    }
    _MinEdgeLengthArray->SetNumberOfTuples(_TriangulatedInput->GetNumberOfPoints());

    if (!_MaxEdgeLengthArray) {
      _MaxEdgeLengthArray = vtkSmartPointer<vtkFloatArray>::New();
      _MaxEdgeLengthArray->SetName("MaxEdgeLength");
      _MaxEdgeLengthArray->SetNumberOfComponents(1);
    }
    _MaxEdgeLengthArray->SetNumberOfTuples(_TriangulatedInput->GetNumberOfPoints());

    // Initialize edge length range
    double x, a, b;
    for (vtkIdType ptId = 0; ptId < _TriangulatedInput->GetNumberOfPoints(); ++ptId) {
      // Rescale adaptive factor to [0, 1]
      x = clamp(scale * _AdaptiveEdgeLengthArray->GetComponent(ptId, 0) + offset, .0, 1.0);
      // Compute lower logistic function interpolation coefficients
      a = 1.0 / (1.0 + exp(- k * (x - x1))), b = 1.0 - a;
      _MinEdgeLengthArray->SetComponent(ptId, 0, a * min_length + b * max_length);
      // Compute upper logistic function interpolation coefficients
      a = 1.0 / (1.0 + exp(- k * (x - x2))), b = 1.0 - a;
      _MaxEdgeLengthArray->SetComponent(ptId, 0, a * min_length + b * max_length);
    }

    // Add edge length point data arrays such that these are
    // interpolated during the local remeshing operations
    _TriangulatedInput->GetPointData()->AddArray(_MinEdgeLengthArray);
    _TriangulatedInput->GetPointData()->AddArray(_MaxEdgeLengthArray);

  // ...or use input cell/point data arrays
  } else {
    _MinEdgeLengthArray = _TriangulatedInput->GetPointData()->GetArray("MinEdgeLength");
    _MaxEdgeLengthArray = _TriangulatedInput->GetPointData()->GetArray("MaxEdgeLength");
  }

  // Smooth edge length ranges
  if (_MinEdgeLengthArray || _MaxEdgeLengthArray) {
    PolyDataSmoothing smoother;
    smoother.Input(_TriangulatedInput);
    if (_MinEdgeLengthArray) smoother.SmoothArray(_MinEdgeLengthArray->GetName());
    if (_MaxEdgeLengthArray) smoother.SmoothArray(_MaxEdgeLengthArray->GetName());
    smoother.NumberOfIterations(_AdaptiveEdgeLengthArray ? 10 : 2);
    smoother.Run();
    vtkPointData * const smoothPD = smoother.Output()->GetPointData();
    if (_MinEdgeLengthArray) {
      _TriangulatedInput->GetPointData()->RemoveArray(_MinEdgeLengthArray->GetName());
      _MinEdgeLengthArray = smoothPD->GetArray(_MinEdgeLengthArray->GetName());
      _TriangulatedInput->GetPointData()->AddArray(_MinEdgeLengthArray);
    }
    if (_MaxEdgeLengthArray) {
      _TriangulatedInput->GetPointData()->RemoveArray(_MaxEdgeLengthArray->GetName());
      _MaxEdgeLengthArray = smoothPD->GetArray(_MaxEdgeLengthArray->GetName());
      _TriangulatedInput->GetPointData()->AddArray(_MaxEdgeLengthArray);
    }
  }
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::Execute()
{
  // Melting pass
  Melting();

  // Inversion pass
  Inversion();

  // Subdivision pass
  Subdivision();
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::Melting()
{
  MIRTK_START_TIMING();

  // Melt triplets of triangles adjacent to nodes with connectivity 3
  // and remove any triangles connected to possibly resulting nodes with
  // connectivity less than 3
  if (_MeltNodes) MeltingOfNodes();

  // Process triangles in chosen melting order
  MeltingOfCells();

  // Melt triplets of triangles adjacent to nodes with connectivity 3
  // and remove any triangles connected to possibly resulting nodes with
  // connectivity less than 3
  if (_MeltNodes) MeltingOfNodes();

  // Remove cells which are marked as deleted
  if (NumberOfMeltings() > 0) {
    _Output->RemoveDeletedCells();
    _Output->BuildLinks();
  }

  MIRTK_DEBUG_TIMING(2, "melting pass");
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::InversionOfTrianglesSharingOneLongEdge()
{
  MIRTK_START_TIMING();

  int       i, j, k;
  double    p[4][3], length2[3], min2[3], max2[3];
  vtkIdType adjCellId, adjPtId, npts, *pts;

  for (vtkIdType cellId = 0; cellId < _Output->GetNumberOfCells(); ++cellId) {

    _Output->GetCellPoints(cellId, npts, pts);
    if (npts == 0) continue; // cell marked as deleted (i.e., VTK_EMPTY_CELL)
    mirtkAssert(npts == 3, "surface is triangulated");

    // Get (transformed) point coordinates
    GetPoint(pts[0], p[0]);
    GetPoint(pts[1], p[1]);
    GetPoint(pts[2], p[2]);

    // Calculate lengths of triangle edges
    length2[0] = vtkMath::Distance2BetweenPoints(p[0], p[1]);
    length2[1] = vtkMath::Distance2BetweenPoints(p[1], p[2]);
    length2[2] = vtkMath::Distance2BetweenPoints(p[2], p[0]);

    min2[0] = SquaredMinEdgeLength(pts[0], pts[1]);
    min2[1] = SquaredMinEdgeLength(pts[1], pts[2]);
    min2[2] = SquaredMinEdgeLength(pts[2], pts[0]);

    max2[0] = SquaredMaxEdgeLength(pts[0], pts[1]);
    max2[1] = SquaredMaxEdgeLength(pts[1], pts[2]);
    max2[2] = SquaredMaxEdgeLength(pts[2], pts[0]);

    // Determine if triangle is elongated
    for (i = 0; i < 3; ++i) {
      if (length2[i] > max2[i]) {
        for (j = 0; j < 3; ++j) {
          if (j != i && (length2[j] < min2[j] || length2[j] > max2[j])) break;
        }
        if (j == 3) break;
      }
    }
    // When long edge found...
    if (i < 3) {                // 1st long edge point index
      j = (i == 2 ? 0 : i + 1); // 2nd long edge point index
      k = (j == 2 ? 0 : j + 1); // 3rd point of this triangle
      // Check connectivity of long edge end points
      if (NodeConnectivity(pts[i]) > 3 && NodeConnectivity(pts[j]) > 3) {
        // Get other vertex of triangle sharing long edge
        adjPtId = GetCellEdgeNeighborPoint(cellId, pts[i], pts[j]);
        if (adjPtId == -1 || _Output->IsEdge(pts[k], adjPtId)) continue;
        // Check if length of other edges are in range
        GetPoint(adjPtId, p[3]);
        length2[0] = vtkMath::Distance2BetweenPoints(p[i], p[3]);
        if (length2[0] < SquaredMinEdgeLength(pts[i], adjPtId) ||
            length2[0] > SquaredMaxEdgeLength(pts[i], adjPtId)) continue;
        length2[1] = vtkMath::Distance2BetweenPoints(p[j], p[3]);
        if (length2[1] < SquaredMinEdgeLength(pts[j], adjPtId) ||
            length2[1] > SquaredMaxEdgeLength(pts[j], adjPtId)) continue;
        // Perform inversion operation
        adjCellId = GetCellEdgeNeighbor(cellId, pts[i], pts[j]);
        if (adjCellId == -1) continue;
        ReplaceCellPoint(cellId,    pts[j], adjPtId);
        ReplaceCellPoint(adjCellId, pts[i], pts[k]);
        ++_NumberOfInversions;
      }
    }
  }

  MIRTK_DEBUG_TIMING(3, "inversion of triangles shared one long edge");
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::InversionOfTrianglesToIncreaseMinHeight()
{
  MIRTK_START_TIMING();

  int       i, j, k;
  double    p[4][3], length2[3], l1, l2, a1, a2, a3, a4;
  vtkIdType adjCellId, adjPtId, npts, *pts;

  for (vtkIdType cellId = 0; cellId < _Output->GetNumberOfCells(); ++cellId) {

    _Output->GetCellPoints(cellId, npts, pts);
    if (npts == 0) continue; // cell marked as deleted (i.e., VTK_EMPTY_CELL)
    mirtkAssert(npts == 3, "surface is triangulated");

    // Get (transformed) point coordinates
    GetPoint(pts[0], p[0]);
    GetPoint(pts[1], p[1]);
    GetPoint(pts[2], p[2]);

    // Calculate lengths of triangle edges
    length2[0] = vtkMath::Distance2BetweenPoints(p[0], p[1]);
    length2[1] = vtkMath::Distance2BetweenPoints(p[1], p[2]);
    length2[2] = vtkMath::Distance2BetweenPoints(p[2], p[0]);

    // Get longest edge in triangle
    i = 0;
    if (length2[1] > length2[0]) i = 1;
    if (length2[2] > length2[i]) i = 2; // 1st long edge point index
    j = (i == 2 ? 0 : i + 1);           // 2nd long edge point index
    k = (j == 2 ? 0 : j + 1);           // 3rd point of this triangle

    a1 = vtkTriangle::TriangleArea(p[0], p[1], p[2]);
    if (a1 > .25 * length2[i]) continue; // ratio of height over l1

    // Check connectivity of long edge end points
    if (NodeConnectivity(pts[i]) > 3 && NodeConnectivity(pts[j]) > 3) {

      // Get other vertex of triangle sharing long edge
      adjCellId = GetCellEdgeNeighbor(cellId, pts[i], pts[j]);
      if (adjCellId == -1) continue;

      adjPtId = GetCellEdgeNeighborPoint(cellId, pts[i], pts[j]);
      if (adjPtId == -1 || _Output->IsEdge(pts[k], adjPtId)) continue;
      GetPoint(adjPtId, p[3]);

      // Do not invert when other edge of neighboring triangle is longer
      if (vtkMath::Distance2BetweenPoints(p[i], p[3]) > length2[i] ||
          vtkMath::Distance2BetweenPoints(p[j], p[3]) > length2[i]) continue;

      // Length of edge before and after inversion
      l1 = sqrt(length2[i]);
      l2 = sqrt(vtkMath::Distance2BetweenPoints(p[k], p[3]));

      // Areas of adjacent triangles after inversion
      a2 = vtkTriangle::TriangleArea(p[i], p[j], p[3]);
      a3 = vtkTriangle::TriangleArea(p[i], p[k], p[3]);
      a4 = vtkTriangle::TriangleArea(p[j], p[k], p[3]);

      // Perform inversion operation when minimum height over edge before
      // and after inversion is greater after the operation
      if ((a1 + a2) * l2 < (a3 + a4) * l1) {
        const vtkIdType ptId[4] = {pts[i], pts[j], pts[k], adjPtId};
        ReplaceCellPoint(adjCellId, ptId[0], ptId[2]);
        ReplaceCellPoint(cellId,    ptId[1], ptId[3]);
        ++_NumberOfInversions;
      }
    }
  }

  MIRTK_DEBUG_TIMING(3, "inversion of triangles to increase height");
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::Inversion()
{
  MIRTK_START_TIMING();

  if (_InvertTrianglesSharingOneLongEdge) {
    InversionOfTrianglesSharingOneLongEdge();
  }
  if (_InvertTrianglesToIncreaseMinHeight) {
    InversionOfTrianglesToIncreaseMinHeight();
  }

  if (_InvertTrianglesSharingOneLongEdge || _InvertTrianglesToIncreaseMinHeight) {
    MIRTK_DEBUG_TIMING(2, "inversion pass");
  }
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::Subdivision()
{
  MIRTK_START_TIMING();

  int       bisect[3];
  vtkIdType npts, *pts;
  double    p1[3], p2[3], p3[3], n1[3], n2[3], n3[3], length2[3], min2[3], max2[3];

  vtkSmartPointer<vtkCellArray> newPolys = vtkSmartPointer<vtkCellArray>::New();
  newPolys->Allocate(_Output->GetNumberOfCells());
  vtkPoints *points = _Output->GetPoints();

  // Get current number of cells (excl. newly added cells)
  const vtkIdType ncells = _Output->GetNumberOfCells();
  for (vtkIdType cellId = 0; cellId < ncells; ++cellId) {

    _Output->GetCellPoints(cellId, npts, pts);
    if (npts == 0) continue; // cell marked as deleted (i.e., VTK_EMPTY_CELL)
    mirtkAssert(npts == 3, "surface is triangulated");

    // Get (transformed) point coordinates
    GetPoint(pts[0], p1);
    GetPoint(pts[1], p2);
    GetPoint(pts[2], p3);

    // Compute squared edge lengths
    length2[0] = vtkMath::Distance2BetweenPoints(p1, p2);
    length2[1] = vtkMath::Distance2BetweenPoints(p2, p3);
    length2[2] = vtkMath::Distance2BetweenPoints(p3, p1);

    // Get desired range of squared edge lengths
    min2[0] = SquaredMinEdgeLength(pts[0], pts[1]);
    min2[1] = SquaredMinEdgeLength(pts[1], pts[2]);
    min2[2] = SquaredMinEdgeLength(pts[2], pts[0]);

    max2[0] = SquaredMaxEdgeLength(pts[0], pts[1]);
    max2[1] = SquaredMaxEdgeLength(pts[1], pts[2]);
    max2[2] = SquaredMaxEdgeLength(pts[2], pts[0]);

    // Determine which edges to bisect
    bisect[0] = int(length2[0] > max2[0]);
    bisect[1] = int(length2[1] > max2[1]);
    bisect[2] = int(length2[2] > max2[2]);

    if (_MaxFeatureAngle < 180.0 && (!bisect[0] || !bisect[1] || !bisect[2])) {
      GetNormal(pts[0], n1);
      GetNormal(pts[1], n2);
      GetNormal(pts[2], n3);
      if (!bisect[0] && length2[0] >= 2.0 * min2[0]) {
        bisect[0] = int(1.0 - vtkMath::Dot(n1, n2) > _MaxFeatureAngleCos);
      }
      if (!bisect[1] && length2[1] >= 2.0 * min2[1]) {
        bisect[1] = int(1.0 - vtkMath::Dot(n2, n3) > _MaxFeatureAngleCos);
      }
      if (!bisect[2] && length2[2] >= 2.0 * min2[2]) {
        bisect[2] = int(1.0 - vtkMath::Dot(n3, n1) > _MaxFeatureAngleCos);
      }
    }

    // Perform subdivision
    switch (bisect[0] + bisect[1] + bisect[2]) {
      case 1:
        if      (bisect[0]) Bisect(cellId, pts[0], pts[1], pts[2], points, newPolys);
        else if (bisect[1]) Bisect(cellId, pts[1], pts[2], pts[0], points, newPolys);
        else                Bisect(cellId, pts[2], pts[0], pts[1], points, newPolys);
        break;
      case 2:
        if      (bisect[0] && bisect[1]) Trisect(cellId, pts[0], pts[1], pts[2], points, newPolys);
        else if (bisect[1] && bisect[2]) Trisect(cellId, pts[1], pts[2], pts[0], points, newPolys);
        else                             Trisect(cellId, pts[2], pts[0], pts[1], points, newPolys);
        break;
      case 3:
        Quadsect(cellId, pts[0], pts[1], pts[2], points, newPolys);
        break;
      default:
        newPolys->InsertNextCell(npts, pts);
    }
  }

  newPolys->Squeeze();
  _Output->DeleteCells();
  _Output->SetPolys(newPolys);
  _Output->BuildLinks();

  MIRTK_DEBUG_TIMING(2, "subdivison pass");
}

// -----------------------------------------------------------------------------
void PolyDataRemeshing::Finalize()
{
  // Finalize base class
  PolyDataFilter::Finalize();

  // If input surface mesh unchanged, set output equal to input
  // Users can check if this->Output() == this->Input() to see if something changed
  if (NumberOfChanges() == 0) {
    _Output = _Input;
    return;
  }

  // Release memory pre-allocated for point data as we may have fewer points now
  _Output->GetPointData()->Squeeze();

  // Remove unused/deleted points and merge middle points of bisected edges
  //
  // FIXME: Under certain circumstances, the cleaner Update causes a segfault.
  //        This is avoided when running the vtkPolyDataNormals filter before.
  //        An inconsistent cell vertex ordering may in particular be caused
  //        by the Quadsect function (maybe also other Subdivision routines).
  vtkNew<vtkPolyDataNormals> prenormals;
  SetVTKInput(prenormals, _Output);
  prenormals->SplittingOff();
  prenormals->AutoOrientNormalsOn();
  prenormals->ConsistencyOn();
  prenormals->ComputeCellNormalsOff();
  prenormals->ComputePointNormalsOn();
  prenormals->NonManifoldTraversalOn();
  prenormals->Update();
  _Output = prenormals->GetOutput();

  vtkNew<vtkCleanPolyData> cleaner;
  SetVTKInput(cleaner, _Output);
  cleaner->PointMergingOn();
  cleaner->SetAbsoluteTolerance(.0);
  cleaner->ConvertPolysToLinesOn();
  cleaner->ConvertLinesToPointsOn();
  cleaner->ConvertStripsToPolysOn();
  cleaner->Update();
  _Output = cleaner->GetOutput();
  _Output->SetVerts(nullptr);
  _Output->SetLines(nullptr);

  // Recompute normals and ensure consistent order of vertices
  vtkNew<vtkPolyDataNormals> normals;
  SetVTKInput(normals, _Output);
  normals->SplittingOff();
  normals->AutoOrientNormalsOn();
  normals->ConsistencyOn();
  normals->ComputeCellNormalsOff();
  normals->ComputePointNormalsOn();
  normals->NonManifoldTraversalOn();
  normals->Update();
  _Output = normals->GetOutput();
}


} // namespace mirtk
