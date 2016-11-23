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

#include "mirtk/PointSetUtils.h"

#include "mirtk/Vtk.h"
#include "mirtk/Utils.h"
#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Stream.h"
#include "mirtk/Parallel.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/DataStatistics.h"
#include "mirtk/PointSet.h"
#include "mirtk/Voxel.h"
#include "mirtk/BaseImage.h"
#include "mirtk/Vector.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/Vector3.h"
#include "mirtk/Algorithm.h"
#include "mirtk/UnorderedSet.h"
#include "mirtk/List.h"

#include "vtkPolyData.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStructuredGrid.h"
#include "vtkDataArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkFloatArray.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataSetAttributes.h"
#include "vtkCell.h"
#include "vtkTriangle.h"
#include "vtkTetra.h"
#include "vtkGenericCell.h"
#include "vtkImageData.h"
#include "vtkImageStencilData.h"
#include "vtkConnectivityFilter.h"
#include "vtkPolyDataConnectivityFilter.h"

#include "vtkImageStencilIterator.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkCleanPolyData.h"
#include "vtkTriangleFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkHull.h"
#include "vtkMaskPoints.h"
#include "vtkDelaunay3D.h"
#include "vtkMassProperties.h"


namespace mirtk {


// =============================================================================
// VTK / MIRTK type conversion
// =============================================================================

// -----------------------------------------------------------------------------
void AddPoints(PointSet &oset, vtkPointSet *iset)
{
  double p[3];
  oset.Reserve(oset.Size() + static_cast<int>(iset->GetNumberOfPoints()));
  for (vtkIdType ptId = 0; ptId < iset->GetNumberOfPoints(); ++ptId) {
    iset->GetPoint(ptId, p);
    oset.Add(Point(p));
  }
}

// =============================================================================
// Point set domain
// =============================================================================

// -----------------------------------------------------------------------------
int Dimension(vtkDataSet *dataset)
{
  double bounds[6];
  dataset->GetBounds(bounds);
  double sx = bounds[1] - bounds[0];
  double sy = bounds[3] - bounds[2];
  double sz = bounds[5] - bounds[4];
  return ((sx > .0 ? 1 : 0) + (sy > .0 ? 1 : 0) + (sz > .0 ? 1 : 0));
}

// -----------------------------------------------------------------------------
ImageAttributes PointSetDomain(vtkPointSet *data, double dx, double dy, double dz)
{
  ImageAttributes attr;
  attr._dx = dx;
  attr._dy = dy >= .0 ? dy : dx;
  attr._dz = dz >= .0 ? dz : dx;
  attr._dt = .0;
  // Compute eigenvectors of covariance matrix
  double c[3], p[3];
  data->GetCenter(c);
  Matrix3x3 covar(.0);
  for (vtkIdType i = 0; i < data->GetNumberOfPoints(); ++i) {
    data->GetPoint(i, p);
    for (int d = 0; d < 3; ++d) p[d] -= c[d];
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        covar[r][c] += p[r] * p[c];
      }
    }
  }
  double  eigen[3];
  Vector3 axis [3];
  covar.EigenSolveSymmetric(eigen, axis);
  Vector3::GenerateOrthonormalBasis(axis[0], axis[1], axis[2]);
  // Set output origin and orientation
  attr._xorigin = c[0];
  attr._yorigin = c[1];
  attr._zorigin = c[2];
  for (int d = 0; d < 3; ++d) {
    attr._xaxis[d] = axis[0][d];
    attr._yaxis[d] = axis[1][d];
    attr._zaxis[d] = axis[2][d];
  }
  // Determine bounds of reoriented data set
  double x, y, z;
  double xmin = numeric_limits<double>::max(), ymin = xmin, zmin = xmin;
  double xmax = -xmin, ymax = -ymin, zmax = -zmin;
  OrderedSet<double> xs, ys, zs;
  for (vtkIdType i = 0; i < data->GetNumberOfPoints(); ++i) {
    data->GetPoint(i, p);
    x = attr._xaxis[0] * p[0] + attr._xaxis[1] * p[1] + attr._xaxis[2] * p[2];
    y = attr._yaxis[0] * p[0] + attr._yaxis[1] * p[1] + attr._yaxis[2] * p[2];
    z = attr._zaxis[0] * p[0] + attr._zaxis[1] * p[1] + attr._zaxis[2] * p[2];
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (z < zmin) zmin = z;
    if (z > zmax) zmax = z;
    if (attr._dx <= .0) xs.insert(x);
    if (attr._dy <= .0) ys.insert(y);
    if (attr._dz <= .0) zs.insert(z);
  }
  // Set output resolution and size
  const double extent[3]  = {xmax - xmin, ymax - ymin, zmax - zmin};
  const double avg_extent = (extent[0] + extent[1] + extent[2]) / 3.0;
  if (attr._dx <= .0) attr._dx = ((extent[0] / avg_extent > 1e-3) ? AverageInterval(xs) : .0);
  if (attr._dy <= .0) attr._dy = ((extent[1] / avg_extent > 1e-3) ? AverageInterval(ys) : .0);
  if (attr._dz <= .0) attr._dz = ((extent[2] / avg_extent > 1e-3) ? AverageInterval(zs) : .0);
  attr._x  = (attr._dx > .0 ? iround(extent[0] / attr._dx) : 0) + 1;
  attr._y  = (attr._dy > .0 ? iround(extent[1] / attr._dy) : 0) + 1;
  attr._z  = (attr._dz > .0 ? iround(extent[2] / attr._dz) : 0) + 1;
  attr._dx = (attr._x  >  1 ? extent[0] / (attr._x - 1) : .0);
  attr._dy = (attr._y  >  1 ? extent[1] / (attr._y - 1) : .0);
  attr._dz = (attr._z  >  1 ? extent[2] / (attr._z - 1) : .0);
  return attr;
}

// -----------------------------------------------------------------------------
ImageAttributes PointSetDomain(vtkPointSet *data, const Vector3D<double> &ds)
{
  return PointSetDomain(data, ds._x, ds._y, ds._z);
}

// -----------------------------------------------------------------------------
ImageAttributes PolyDataDomain(vtkPolyData *data, double dx, double dy, double dz)
{
  return PointSetDomain(data, dx, dy, dz);
}

// -----------------------------------------------------------------------------
ImageAttributes PolyDataDomain(vtkPolyData *data, const Vector3D<double> &ds)
{
  return PointSetDomain(data, ds);
}

// =============================================================================
// Point/cell data
// =============================================================================

// -----------------------------------------------------------------------------
int PolyDataAttributeType(const char *type)
{
  string ltype(type);
  transform(ltype.begin(), ltype.end(), ltype.begin(), ::tolower);
  if (ltype == "scalars")     return vtkDataSetAttributes::SCALARS;
  if (ltype == "vectors")     return vtkDataSetAttributes::VECTORS;
  if (ltype == "normals")     return vtkDataSetAttributes::NORMALS;
  if (ltype == "tcoords")     return vtkDataSetAttributes::TCOORDS;
  if (ltype == "tensors")     return vtkDataSetAttributes::TENSORS;
  if (ltype == "globalids")   return vtkDataSetAttributes::GLOBALIDS;
  if (ltype == "pedigreeids") return vtkDataSetAttributes::PEDIGREEIDS;
  if (ltype == "edgeflag")    return vtkDataSetAttributes::EDGEFLAG;
  return -1;
}

// -----------------------------------------------------------------------------
vtkDataArray *GetArrayByCaseInsensitiveName(vtkDataSetAttributes *data, const char *name, int *loc)
{
  int i;
  vtkDataArray *arr = data->GetArray(name, i);
  if (arr == nullptr) {
    const string lname = ToLower(name);
    for (i = 0; i < data->GetNumberOfArrays(); ++i) {
      const char * const arr_name = data->GetArrayName(i);
      if (arr_name && ToLower(arr_name) == lname) break;
    }
    if (i == data->GetNumberOfArrays()) i = -1;
  }
  if (loc) *loc = i;
  return arr;
}

// -----------------------------------------------------------------------------
int DeepCopyArrayUsingCaseInsensitiveName(vtkDataSetAttributes *dst, vtkDataSetAttributes *src, const char *name, bool deep)
{
  int loc = -1;
  vtkDataArray *src_array = GetArrayByCaseInsensitiveName(src, name);
  if (src_array) {
    vtkDataArray *dst_array = GetArrayByCaseInsensitiveName(dst, name, &loc);
    if (dst_array) {
      dst_array->DeepCopy(src_array);
    } else {
      vtkSmartPointer<vtkDataArray> copy;
      copy.TakeReference(src_array->NewInstance());
      copy->DeepCopy(src_array);
      loc = dst->AddArray(copy);
    }
  }
  return loc;
}

// -----------------------------------------------------------------------------
bool IsCategoricalArrayName(const string &name)
{
  if (name.length() > 2) {
    const string suffix = name.substr(name.length() - 2);
    if (suffix == "Id") return true;
  }
  if (name.length() > 4) {
    const string suffix = name.substr(name.length() - 4);
    if (suffix == "Label" || suffix == "Mask") return true;
  }
  return false;
}

// =============================================================================
// Cells
// =============================================================================

// -----------------------------------------------------------------------------
// Attention: vtkGenericCell is not derived from the respective cell types
double ComputeArea(vtkCell *cell)
{
  if (cell->GetCellDimension() < 2) return .0;
  vtkPoints *points = cell->GetPoints();
  switch (cell->GetCellType()) {
    case VTK_TRIANGLE: {
      double p1[3], p2[3], p3[3];
      points->GetPoint(0, p1);
      points->GetPoint(1, p2);
      points->GetPoint(2, p3);
      return vtkTriangle::TriangleArea(p1, p2, p3);
    }
    default: return numeric_limits<double>::quiet_NaN();
  }
}

// -----------------------------------------------------------------------------
// Attention: vtkGenericCell is not derived from the respective cell types
double ComputeVolume(vtkCell *cell)
{
  if (cell->GetCellDimension() < 3) return .0;
  vtkPoints *points = cell->GetPoints();
  switch (cell->GetCellType()) {
    case VTK_TETRA: {
      double p1[3], p2[3], p3[3], p4[3];
      points->GetPoint(0, p1);
      points->GetPoint(1, p2);
      points->GetPoint(2, p3);
      points->GetPoint(3, p4);
      return vtkTetra::ComputeVolume(p1, p2, p3, p4);
    }
    default: return numeric_limits<double>::quiet_NaN();
  }
}

// =============================================================================
// Basic point set manipulation
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> DataSetSurface(vtkSmartPointer<vtkDataSet> dataset, bool passPtIds, bool passCellIds)
{
  if (IsSurfaceMesh(dataset)) {
    vtkSmartPointer<vtkPolyData> surface = vtkSmartPointer<vtkPolyData>::New();
    surface->ShallowCopy(dataset);
    return surface;
  } else {
    vtkSmartPointer<vtkDataSetSurfaceFilter> surface;
    surface = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    SetVTKInput(surface, dataset);
    surface->SetPassThroughPointIds(passPtIds);
    surface->SetPassThroughCellIds(passCellIds);
    surface->Update();
    return surface->GetOutput();
  }
}

// -----------------------------------------------------------------------------
void Center(vtkSmartPointer<vtkPointSet> pointset)
{
  double c[3], p[3];
  pointset->GetCenter(c);
  vtkPoints *points = pointset->GetPoints();
  for (vtkIdType ptId = 0; ptId < pointset->GetNumberOfPoints(); ++ptId) {
    points->GetPoint(ptId, p);
    p[0] -= c[0];
    p[1] -= c[1];
    p[2] -= c[2];
    points->SetPoint(ptId, p);
  }
  pointset->Modified();
}

// -----------------------------------------------------------------------------
void Scale(vtkSmartPointer<vtkPointSet> pointset, double scale)
{
  double c[3], p[3];
  pointset->GetCenter(c);
  vtkPoints *points = pointset->GetPoints();
  for (vtkIdType ptId = 0; ptId < pointset->GetNumberOfPoints(); ++ptId) {
    points->GetPoint(ptId, p);
    p[0] = c[0] + scale * (p[0] - c[0]);
    p[1] = c[1] + scale * (p[1] - c[1]);
    p[2] = c[2] + scale * (p[2] - c[2]);
    points->SetPoint(ptId, p);
  }
  pointset->Modified();
}

// =============================================================================
// Surface meshes
// =============================================================================

// -----------------------------------------------------------------------------
UnorderedSet<int> BoundaryPoints(vtkDataSet *dataset, const EdgeTable *edgeTable)
{
  UniquePtr<EdgeTable> tmpEdgeTable;
  if (edgeTable == nullptr) {
    tmpEdgeTable.reset(new EdgeTable(dataset));
    edgeTable = tmpEdgeTable.get();
  }
  UnorderedSet<int> boundaryPtIds;
  vtkSmartPointer<vtkIdList> cellIds1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellIds2 = vtkSmartPointer<vtkIdList>::New();
  vtkIdType ptId1, ptId2;
  EdgeIterator edgeIt(*edgeTable);
  for (edgeIt.InitTraversal(); edgeIt.GetNextEdge(ptId1, ptId2) != -1;) {
    dataset->GetPointCells(ptId1, cellIds1);
    dataset->GetPointCells(ptId2, cellIds2);
    cellIds1->IntersectWith(cellIds2);
    if (cellIds1->GetNumberOfIds() < 2) {
      boundaryPtIds.insert(static_cast<int>(ptId1));
      boundaryPtIds.insert(static_cast<int>(ptId2));
    }
  }
  return boundaryPtIds;
}

// -----------------------------------------------------------------------------
EdgeList BoundaryEdges(vtkDataSet *dataset, const EdgeTable &edgeTable)
{
  EdgeList boundaryEdges;
  vtkSmartPointer<vtkIdList> cellIds1 = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkIdList> cellIds2 = vtkSmartPointer<vtkIdList>::New();
  vtkIdType ptId1, ptId2;
  EdgeIterator edgeIt(edgeTable);
  for (edgeIt.InitTraversal(); edgeIt.GetNextEdge(ptId1, ptId2) != -1;) {
    dataset->GetPointCells(ptId1, cellIds1);
    dataset->GetPointCells(ptId2, cellIds2);
    cellIds1->IntersectWith(cellIds2);
    if (cellIds1->GetNumberOfIds() < 2) {
      boundaryEdges.push_back(MakePair(static_cast<int>(ptId1), static_cast<int>(ptId2)));
    }
  }
  return boundaryEdges;
}

// -----------------------------------------------------------------------------
EdgeList BoundaryEdges(vtkDataSet *dataset)
{
  EdgeTable edgeTable(dataset);
  return BoundaryEdges(dataset, edgeTable);
}

// -----------------------------------------------------------------------------
EdgeList GetPointEdges(const EdgeList &edges, int ptId)
{
  EdgeList result;
  for (const auto &edge : edges) {
    if (edge.first == ptId) {
      result.push_back(edge);
    } else if (edge.second == ptId) {
      result.push_back(MakePair(edge.second, edge.first));
    }
  }
  return result;
}

// -----------------------------------------------------------------------------
EdgeList PopPointEdges(EdgeList &edges, int ptId)
{
  EdgeList result;
  for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
    if (edge->first == ptId) {
      result.push_back(*edge);
      edge = edges.erase(edge);
    } else if (edge->second == ptId) {
      result.push_back(MakePair(edge->second, edge->first));
      edge = edges.erase(edge);
    }
  }
  return result;
}

// -----------------------------------------------------------------------------
Array<Array<int> > BoundarySegments(vtkDataSet *dataset, const EdgeTable *edgeTable)
{
  Array<Array<int> > boundaries;
  UniquePtr<EdgeTable> tmpEdgeTable;
  if (edgeTable == nullptr) {
    tmpEdgeTable.reset(new EdgeTable(dataset));
    edgeTable = tmpEdgeTable.get();
  }
  int ptId1, ptId2;
  UnorderedSet<int> boundaryPtIds;
  EdgeList boundaryEdges = BoundaryEdges(dataset, *edgeTable);
  for (auto edge = boundaryEdges.begin(); edge != boundaryEdges.end(); ++edge) {
    boundaryPtIds.insert(edge->first);
    boundaryPtIds.insert(edge->second);
  }
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
    boundaries.reserve(boundaries.size() + 1);
    boundaries.push_back(ptIds);
  }
  return boundaries;
}

// -----------------------------------------------------------------------------
int NumberOfPoints(vtkDataSet *dataset)
{
  int n = 0;
  bool deleted;
  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  for (vtkIdType ptId = 0; ptId < dataset->GetNumberOfPoints(); ++ptId) {
    dataset->GetPointCells(ptId, cellIds);
    deleted = true;
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
      if (dataset->GetCellType(cellIds->GetId(i)) != VTK_EMPTY_CELL) {
        deleted = false;
        break;
      }
    }
    if (!deleted) {
      ++n;
    }
  }
  return n;
}

// -----------------------------------------------------------------------------
int NumberOfEdges(vtkDataSet *dataset, const EdgeTable *edgeTable)
{
  if (edgeTable == nullptr) {
    EdgeTable edgeTable(dataset);
    return edgeTable.NumberOfEdges();
  } else {
    return edgeTable->NumberOfEdges();
  }
}

// -----------------------------------------------------------------------------
int NumberOfFaces(vtkDataSet *dataset)
{
  int nfaces = 0;
  for (vtkIdType cellId = 0; cellId < dataset->GetNumberOfCells(); ++cellId) {
    switch (dataset->GetCellType(cellId)) {
      case VTK_TRIANGLE:
      case VTK_QUAD:
      case VTK_POLYGON: {
        nfaces += 1;
      } break;
      case VTK_TRIANGLE_STRIP: {
        nfaces += dataset->GetCell(cellId)->GetNumberOfFaces();
      } break;
    }
  }
  return nfaces;
}

// -----------------------------------------------------------------------------
int NumberOfEmptyCells(vtkDataSet *dataset)
{
  int n = 0;
  for (vtkIdType cellId = 0; cellId < dataset->GetNumberOfCells(); ++cellId) {
    if (dataset->GetCellType(cellId) == VTK_EMPTY_CELL) {
      ++n;
    }
  }
  return n;
}

// -----------------------------------------------------------------------------
int NumberOfConnectedComponents(vtkDataSet *dataset)
{
  vtkPolyData *polydata = vtkPolyData::SafeDownCast(dataset);
  if (polydata != nullptr) {
    vtkNew<vtkPolyDataConnectivityFilter> connectivity;
    connectivity->SetExtractionModeToAllRegions();
    connectivity->ScalarConnectivityOff();
    SetVTKInput(connectivity, polydata);
    connectivity->Update();
    return connectivity->GetNumberOfExtractedRegions();
  } else {
    vtkNew<vtkConnectivityFilter> connectivity;
    connectivity->SetExtractionModeToAllRegions();
    connectivity->ScalarConnectivityOff();
    SetVTKInput(connectivity, dataset);
    connectivity->Update();
    return connectivity->GetNumberOfExtractedRegions();
  }
}

// -----------------------------------------------------------------------------
int NumberOfBoundarySegments(vtkDataSet *dataset, const EdgeTable *edgeTable)
{
  return static_cast<int>(BoundarySegments(dataset, edgeTable).size());
}

// -----------------------------------------------------------------------------
int EulerCharacteristic(vtkDataSet *dataset, const EdgeTable &edgeTable,
                        int *npoints, int *nedges, int *nfaces)
{
  int num_points = NumberOfPoints(dataset);
  int num_edges  = NumberOfEdges(dataset, &edgeTable);
  int num_faces  = NumberOfFaces(dataset);
  if (npoints != nullptr) *npoints = num_points;
  if (nedges  != nullptr) *nedges  = num_edges;
  if (nfaces  != nullptr) *nfaces  = num_faces;
  return num_points - num_edges + num_faces;
}

// -----------------------------------------------------------------------------
int EulerCharacteristic(vtkDataSet *dataset, int *npoints, int *nedges, int *nfaces)
{
  EdgeTable edgeTable(dataset);
  return EulerCharacteristic(dataset, edgeTable, npoints, nedges, nfaces);
}

// -----------------------------------------------------------------------------
double Genus(vtkDataSet *dataset, const EdgeTable &edgeTable,
             int *npoints, int *nedges, int *nfaces, int *nbounds, int *ncomps, int *euler)
{
  int num_points = NumberOfPoints(dataset);
  int num_edges  = NumberOfEdges(dataset);
  int num_faces  = NumberOfFaces(dataset);
  int num_comps  = NumberOfConnectedComponents(dataset);
  int num_bounds = NumberOfBoundarySegments(dataset, &edgeTable);
  int chi        = num_points - num_edges + num_faces;
  if (npoints != nullptr) *npoints = num_points;
  if (nedges  != nullptr) *nedges  = num_edges;
  if (nfaces  != nullptr) *nfaces  = num_faces;
  if (nbounds != nullptr) *nbounds = num_bounds;
  if (ncomps  != nullptr) *ncomps  = num_comps;
  if (euler   != nullptr) *euler   = chi;
  return double(2 * num_comps - num_bounds - chi) / 2.0;
}

// -----------------------------------------------------------------------------
double Genus(vtkDataSet *dataset, int *npoints, int *nedges, int *nfaces, int *nbounds, int *ncomps, int *euler)
{
  EdgeTable edgeTable(dataset);
  return Genus(dataset, edgeTable, npoints, nedges, nfaces, nbounds, ncomps, euler);
}

// -----------------------------------------------------------------------------
namespace AreaUtils {

struct ComputeSurfaceArea
{
  vtkPolyData  *_Surface;
  vtkDataArray *_CellArea;
  double        _SurfaceArea;

  ComputeSurfaceArea(vtkPolyData *surface)
  :
    _Surface(surface),
    _CellArea(surface->GetCellData()->GetArray("Area")),
    _SurfaceArea(.0)
  {}

  ComputeSurfaceArea(const ComputeSurfaceArea &other, split)
  :
    _Surface(other._Surface), _CellArea(other._CellArea), _SurfaceArea(.0)
  {}

  void join(const ComputeSurfaceArea &other)
  {
    _SurfaceArea += other._SurfaceArea;
  }

  void operator ()(const blocked_range<vtkIdType> cellIds)
  {
    vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
    vtkIdType npts, *pts;
    double p1[3], p2[3], p3[3], area;
    for (vtkIdType cellId = cellIds.begin(); cellId != cellIds.end(); ++cellId)
    {
      _Surface->GetCellPoints(cellId, npts, pts);
      if (npts == 3) {
        _Surface->GetPoint(pts[0], p1);
        _Surface->GetPoint(pts[1], p2);
        _Surface->GetPoint(pts[2], p3);
        area = vtkTriangle::TriangleArea(p1, p2, p3);
      } else {
        _Surface->GetCell(cellId, cell);
        area = ComputeArea(cell);
      }
      if (_CellArea) _CellArea->SetComponent(cellId, 0, area);
      _SurfaceArea += area;
    }
  }
};

} // namespace AreaUtils

// -----------------------------------------------------------------------------
double Area(vtkPolyData *polydata, bool per_cell)
{
  vtkSmartPointer<vtkDataArray> area;
  if (per_cell) {
    area = vtkSmartPointer<vtkFloatArray>::New();
    area->SetName("Area");
    area->SetNumberOfComponents(1);
    area->SetNumberOfTuples(polydata->GetNumberOfCells());
    polydata->GetCellData()->AddArray(area);
  }
  AreaUtils::ComputeSurfaceArea eval(polydata);
  blocked_range<vtkIdType> cellIds(0, polydata->GetNumberOfCells());
  parallel_reduce(cellIds, eval);
  return eval._SurfaceArea;
}

// -----------------------------------------------------------------------------
double Area(vtkSmartPointer<vtkPointSet> pointset)
{
  return Area(DataSetSurface(pointset), false);
}


namespace EdgeLengthUtils {

// -----------------------------------------------------------------------------
/// Compute lengths of all edges and store it in C array
struct ComputeEdgeLength
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double          *_EdgeLength;

  void operator ()(const blocked_range<int> &re) const
  {
    int    ptId1, ptId2, edgeId;
    double p1[3], p2[3];

    EdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); (edgeId = it.GetNextEdge(ptId1, ptId2)) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      _EdgeLength[edgeId] = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute squared lengths of all edges and store it in C array
struct ComputeSquaredEdgeLength
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double          *_EdgeLength;

  void operator ()(const blocked_range<int> &re) const
  {
    int    ptId1, ptId2, edgeId;
    double p1[3], p2[3];

    EdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); (edgeId = it.GetNextEdge(ptId1, ptId2)) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      _EdgeLength[edgeId] = vtkMath::Distance2BetweenPoints(p1, p2);
    }
  }
};

// -----------------------------------------------------------------------------
/// Determine minimum/maximum edge length
struct MinMaxEdgeLength
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _Min;
  double           _Max;

  MinMaxEdgeLength() : _Min(+inf), _Max(-inf) {}

  MinMaxEdgeLength(const MinMaxEdgeLength &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Min(other._Min),
    _Max(other._Max)
  {}

  void join(const MinMaxEdgeLength &other)
  {
    if (other._Min < _Min) _Min = other._Min;
    if (other._Max > _Max) _Max = other._Max;
  }

  void operator ()(const blocked_range<int> &re)
  {
    int    ptId1, ptId2;
    double p1[3], p2[3], d;

    EdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); it.GetNextEdge(ptId1, ptId2) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      d = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      if (d < _Min) _Min = d;
      if (d > _Max) _Max = d;
    }
  }
};

// -----------------------------------------------------------------------------
/// Calculate sum of edge lengths
struct SumEdgeLengths
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _Sum;

  SumEdgeLengths() : _Sum(.0) {}

  SumEdgeLengths(const SumEdgeLengths &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Sum(.0)
  {}

  void join(const SumEdgeLengths &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &re)
  {
    int    ptId1, ptId2;
    double p1[3], p2[3];

    EdgeIterator it(*_EdgeTable);
    for (it.InitTraversal(re); it.GetNextEdge(ptId1, ptId2) != -1;) {
      _Points->GetPoint(ptId1, p1);
      _Points->GetPoint(ptId2, p2);
      _Sum += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    }
  }
};

// -----------------------------------------------------------------------------
/// Calculate sum of edge lengths
struct SumEdgeLengthsWithDuplicates
{
  vtkPoints       *_Points;
  const EdgeTable *_EdgeTable;
  double           _Sum;

  SumEdgeLengthsWithDuplicates() : _Sum(.0) {}

  SumEdgeLengthsWithDuplicates(const SumEdgeLengthsWithDuplicates &other, split)
  :
    _Points(other._Points),
    _EdgeTable(other._EdgeTable),
    _Sum(.0)
  {}

  void join(const SumEdgeLengthsWithDuplicates &other)
  {
    _Sum += other._Sum;
  }

  void operator ()(const blocked_range<int> &ptIds)
  {
    const int *adjPtIds;
    int        numAdjPts;
    double     p1[3], p2[3];

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      _Points->GetPoint(ptId, p1);
      _EdgeTable->GetAdjacentPoints(ptId, numAdjPts, adjPtIds);
      for (int i = 0; i < numAdjPts; ++i) {
        _Points->GetPoint(adjPtIds[i], p2);
        _Sum += sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      }
    }
  }
};


} // namespace EdgeLengthUtils

// -----------------------------------------------------------------------------
Vector EdgeLengths(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  const int n = edgeTable.NumberOfEdges();
  Vector edgeLengths(n);
  if (n == 0) return edgeLengths;
  EdgeLengthUtils::ComputeEdgeLength eval;
  eval._Points     = points;
  eval._EdgeTable  = &edgeTable;
  eval._EdgeLength = edgeLengths.RawPointer();
  parallel_for(blocked_range<int>(0, n), eval);
  return edgeLengths;
}

// -----------------------------------------------------------------------------
Vector SquaredEdgeLengths(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  const int n = edgeTable.NumberOfEdges();
  Vector edgeLengths(n);
  if (n == 0) return edgeLengths;
  EdgeLengthUtils::ComputeSquaredEdgeLength eval;
  eval._Points     = points;
  eval._EdgeTable  = &edgeTable;
  eval._EdgeLength = edgeLengths.RawPointer();
  parallel_for(blocked_range<int>(0, n), eval);
  return edgeLengths;
}

// -----------------------------------------------------------------------------
double AverageEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  if (edgeTable.NumberOfEdges() == 0) return .0;
  // FIXME: Parallel reduction of sum of edge lengths yields different result
  //        from the code below or a single-threaded execution. This
  //        could indicate a bug in irtkEdgeIterator::InitTraversal.
//  EdgeLengthUtils::SumEdgeLengths eval;
//  eval._Points    = points;
//  eval._EdgeTable = &edgeTable;
//  parallel_reduce(blocked_range<int>(0, edgeTable.NumberOfEdges()), eval);
//  eval(blocked_range<int>(0, edgeTable.NumberOfEdges()));
//  return eval._Sum / edgeTable.NumberOfEdges();
  EdgeLengthUtils::SumEdgeLengthsWithDuplicates eval;
  eval._Points    = points;
  eval._EdgeTable = &edgeTable;
  parallel_reduce(blocked_range<int>(0, static_cast<int>(points->GetNumberOfPoints())), eval);
  return eval._Sum / (2 * edgeTable.NumberOfEdges());
}

// -----------------------------------------------------------------------------
double AverageEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return AverageEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
double RobustAverageEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  const int n = edgeTable.NumberOfEdges();
  if (n == 0) return .0;
  double *edgeLength = new double[n];
  EdgeLengthUtils::ComputeEdgeLength eval;
  eval._Points     = points;
  eval._EdgeTable  = &edgeTable;
  eval._EdgeLength = edgeLength;
  parallel_for(blocked_range<int>(0, n), eval);
  double mean = data::statistic::RobustMean::Calculate(5, n, edgeLength);
  delete[] edgeLength;
  return mean;
}

// -----------------------------------------------------------------------------
double RobustAverageEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return RobustAverageEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
double MedianEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  const int n = edgeTable.NumberOfEdges();
  if (n == 0) return .0;
  double *edgeLength = new double[n];
  EdgeLengthUtils::ComputeEdgeLength eval;
  eval._Points     = points;
  eval._EdgeTable  = &edgeTable;
  eval._EdgeLength = edgeLength;
  parallel_for(blocked_range<int>(0, n), eval);
  sort(edgeLength, edgeLength + n);
  double median = edgeLength[n / 2];
  delete[] edgeLength;
  return median;
}

// -----------------------------------------------------------------------------
double MedianEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return MedianEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
void GetMinMaxEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable, double &min, double &max)
{
  min = max = .0;
  if (edgeTable.NumberOfEdges() == 0) return;
  EdgeLengthUtils::MinMaxEdgeLength eval;
  eval._Points    = points;
  eval._EdgeTable = &edgeTable;
  parallel_reduce(blocked_range<int>(0, edgeTable.NumberOfEdges()), eval);
  min = eval._Min;
  max = eval._Max;
}

// -----------------------------------------------------------------------------
void GetMinMaxEdgeLength(vtkSmartPointer<vtkPointSet> pointset, double &min, double &max)
{
  EdgeTable edgeTable(pointset);
  GetMinMaxEdgeLength(pointset->GetPoints(), edgeTable, min, max);
}

// -----------------------------------------------------------------------------
double MinEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  double min, max;
  GetMinMaxEdgeLength(points, edgeTable, min, max);
  return min;
}

// -----------------------------------------------------------------------------
double MinEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return MinEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
double MaxEdgeLength(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable)
{
  double min, max;
  GetMinMaxEdgeLength(points, edgeTable, min, max);
  return max;
}

// -----------------------------------------------------------------------------
double MaxEdgeLength(vtkSmartPointer<vtkPointSet> pointset)
{
  EdgeTable edgeTable(pointset);
  return MaxEdgeLength(pointset->GetPoints(), edgeTable);
}

// -----------------------------------------------------------------------------
void EdgeLengthNormalDistribution(vtkSmartPointer<vtkPoints> points, const EdgeTable &edgeTable, double &mean, double &sigma)
{
  mean = sigma = .0;
  if (edgeTable.NumberOfEdges() == 0) return;

  int    ptId1, ptId2, n = 0;
  double p1[3], p2[3], d, delta;

  EdgeIterator it(edgeTable);
  for (it.InitTraversal(); it.GetNextEdge(ptId1, ptId2) != -1;) {
    points->GetPoint(ptId1, p1);
    points->GetPoint(ptId2, p2);
    d = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
    ++n;
    delta = d - mean;
    mean  += delta / n;
    sigma += delta * (d - mean);
  }

  if (n > 1) sigma /= n - 1;
  sigma = sqrt(sigma);
}

// -----------------------------------------------------------------------------
void EdgeLengthNormalDistribution(vtkSmartPointer<vtkPointSet> pointset, double &mean, double &sigma)
{
  EdgeTable edgeTable(pointset);
  EdgeLengthNormalDistribution(pointset->GetPoints(), edgeTable, mean, sigma);
}

// -----------------------------------------------------------------------------
double Volume(vtkSmartPointer<vtkPolyData> surface)
{
  vtkSmartPointer<vtkMassProperties> mp = vtkMassProperties::New();
  SetVTKInput(mp, surface);
  return mp->GetVolume();
}

// -----------------------------------------------------------------------------
bool IsSurfaceMesh(vtkDataSet *dataset)
{
  vtkPolyData *polydata = vtkPolyData::SafeDownCast(dataset);
  if (polydata == nullptr) return false;
  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();
  for (vtkIdType cellId = 0; cellId < polydata->GetNumberOfCells(); ++cellId) {
    polydata->GetCell(cellId, cell);
    if (cell->GetCellDimension() != 2) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> Triangulate(vtkSmartPointer<vtkPolyData> mesh)
{
  vtkSmartPointer<vtkPolyData> output;
  if (IsTriangularMesh(mesh)) {
    output.TakeReference(mesh->NewInstance());
    output->ShallowCopy(mesh);
  } else {
    vtkSmartPointer<vtkTriangleFilter> filter;
    filter = vtkSmartPointer<vtkTriangleFilter>::New();
    SetVTKInput(filter, mesh);
    filter->PassVertsOff();
    filter->PassLinesOff();
    vtkSmartPointer<vtkCleanPolyData> cleaner;
    cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    cleaner->PointMergingOff();
    cleaner->ConvertLinesToPointsOn();
    cleaner->ConvertPolysToLinesOn();
    cleaner->ConvertStripsToPolysOn();
    SetVTKConnection(cleaner, filter);
    cleaner->Update();
    cleaner->GetOutput()->SetLines(NULL);
    cleaner->GetOutput()->SetVerts(NULL);
    output = cleaner->GetOutput();
  }
  return output;
}

// -----------------------------------------------------------------------------
// Implements two different filter pipelines to extract the convex hull
vtkSmartPointer<vtkPolyData> ConvexHull(vtkSmartPointer<vtkPointSet> pointset, int levels)
{
  // Spatially stratify points to prevent "Unable to factor linear system"
  // warning of vtkDelaunay3D filter due to numerical imprecisions
  vtkNew<vtkMaskPoints> stratify;
  stratify->RandomModeOn();
  stratify->SetRandomModeType(2);
  stratify->SetMaximumNumberOfPoints(.75 * pointset->GetNumberOfPoints());
  SetVTKInput(stratify, pointset);

  // Get convex hull of largest component
  vtkNew<vtkHull> hull;
  hull->AddRecursiveSpherePlanes(levels);
  //SetVTKConnection(hull, stratify);
  SetVTKInput(hull, pointset);

  // Compute Delaunay triangulation
  vtkNew<vtkDelaunay3D> delaunay;
  SetVTKConnection(delaunay, hull);

  // Construct surface mesh
  vtkNew<vtkDataSetSurfaceFilter> mesher;
  SetVTKConnection(mesher, delaunay);

  mesher->Update();
  return mesher->GetOutput();
}

// -----------------------------------------------------------------------------
bool IsTriangularMesh(vtkDataSet *input)
{
  if (vtkPointSet::SafeDownCast(input) == nullptr) return false;
  for (vtkIdType cellId = 0; cellId < input->GetNumberOfCells(); ++cellId) {
    int type = input->GetCellType(cellId);
    if (type != VTK_EMPTY_CELL && type != VTK_TRIANGLE) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
bool IsTetrahedralMesh(vtkDataSet *input)
{
  if (vtkPointSet::SafeDownCast(input) == nullptr) return false;
  for (vtkIdType cellId = 0; cellId < input->GetNumberOfCells(); ++cellId) {
    int type = input->GetCellType(cellId);
    if (type != VTK_EMPTY_CELL && type != VTK_TETRA) return false;
  }
  return true;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> Tetrahedralize(vtkSmartPointer<vtkPointSet> input)
{
  vtkSmartPointer<vtkPointSet> mesh;
  if (IsTetrahedralMesh(input)) {
    mesh.TakeReference(input->NewInstance());
    mesh->ShallowCopy(input);
  } else {
    // TODO: Use TetGen library to tetrahedralize interior of input PLC
    cerr << "irtkPolyDataUtils::Tetrahedralize: Not implemented, use TetGen command instead" << endl;
    exit(1);
  }
  return mesh;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkImageData> NewVtkMask(int nx, int ny, int nz)
{
  vtkSmartPointer<vtkImageData> imagedata = vtkSmartPointer<vtkImageData>::New();
  imagedata->SetOrigin(.0, .0, .0);
  imagedata->SetDimensions(nx, ny, nz);
  imagedata->SetSpacing(1.0, 1.0, 1.0);
#if VTK_MAJOR_VERSION >= 6
  imagedata->AllocateScalars(ToVTKDataType(MIRTK_VOXEL_BINARY), 1);
#else
  imagedata->SetScalarType(ToVTKDataType(MIRTK_VOXEL_BINARY));
  imagedata->AllocateScalars();
#endif
  return imagedata;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> WorldToImage(vtkSmartPointer<vtkPointSet> pointset,
                                          const BaseImage             *image)
{
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(pointset->GetNumberOfPoints());
  double p[3];
  for (vtkIdType ptId = 0; ptId < points->GetNumberOfPoints(); ++ptId) {
    pointset->GetPoint(ptId, p);
    image->WorldToImage(p[0], p[1], p[2]);
    points->SetPoint(ptId, p);
  }
  vtkSmartPointer<vtkPointSet> output;
  output.TakeReference(pointset->NewInstance());
  output->ShallowCopy(pointset);
  output->SetPoints(points);
  return output;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkImageStencilData> ImageStencil(vtkSmartPointer<vtkImageData> image,
                                                  vtkSmartPointer<vtkPointSet>  pointset)
{
  vtkSmartPointer<vtkPolyData> surface = DataSetSurface(pointset);
  vtkSmartPointer<vtkPolyDataToImageStencil> filter;
  filter = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
  SetVTKInput(filter, surface);
  filter->SetOutputOrigin(image->GetOrigin());
  filter->SetOutputSpacing(image->GetSpacing());
  filter->SetOutputWholeExtent(image->GetExtent());
  filter->Update();
  return filter->GetOutput();
}

// -----------------------------------------------------------------------------
void ImageStencilToMask(vtkSmartPointer<vtkImageStencilData> stencil,
                        vtkSmartPointer<vtkImageData>        image)
{
  if (image->GetScalarType() != ToVTKDataType(MIRTK_VOXEL_BINARY)) {
    cerr << "ImageStencilToMask: vtkImageData must have scalar type MIRTK_VOXEL_BINARY" << endl;
    exit(1);
  }
  const vtkIdType nvox = image->GetNumberOfPoints();
  memset(image->GetScalarPointer(), 0, nvox * sizeof(BinaryPixel));
  vtkImageStencilIterator<BinaryPixel> it;
  it.Initialize(image, stencil, image->GetExtent());
  while (!it.IsAtEnd()) {
    if (it.IsInStencil()) {
      for (BinaryPixel *cur = it.BeginSpan(); cur != it.EndSpan(); ++cur) {
        *cur = static_cast<BinaryPixel>(true);
      }
    }
    it.NextSpan();
  }
}


} // namespace mirtk
