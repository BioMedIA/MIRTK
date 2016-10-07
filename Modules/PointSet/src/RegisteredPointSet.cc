/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "mirtk/RegisteredPointSet.h"

#include "mirtk/Vtk.h"
#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/PointSet.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"

#include "mirtk/LinearInterpolateImageFunction.hxx"

#include "vtkSmartPointer.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkIdTypeArray.h"
#include "vtkFloatArray.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkStructuredGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkPolyDataNormals.h"
#include "vtkPointLocator.h"
#include "vtkCellLocator.h"


namespace mirtk {


// =============================================================================
// Copy points to PointSet structure
// =============================================================================

namespace RegisteredPointSetUtils {

/// Type of interpolator used to interpolate dense displacement field
typedef GenericLinearInterpolateImageFunction<GenericImage<double> > DisplacementInterpolator;

// -----------------------------------------------------------------------------
/// Copy VTK points to MIRTK point set
class CopyVtkPointsToIrtkPointSet
{
  vtkPoints *_Input;
  PointSet  *_Output;

public:

  CopyVtkPointsToIrtkPointSet(vtkPoints *in, PointSet &out)
  :
    _Input(in), _Output(&out)
  {}

  void operator ()(const blocked_range<int> &re) const
  {
    double p[3];
    for (int i = re.begin(); i != re.end(); ++i) {
      _Input ->GetPoint(i, p);
      _Output->SetPoint(i, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Update points of extracted point set surface
class UpdateSurfacePoints
{
  vtkPoints      *_Points;
  vtkPolyData    *_Surface;
  vtkIdTypeArray *_PtIds;

public:

  UpdateSurfacePoints(vtkPoints *points, vtkPolyData *surface, vtkIdTypeArray *ptIds)
  :
    _Points(points), _Surface(surface), _PtIds(ptIds)
  {}

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double p[3];
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      _Points->GetPoint(_PtIds->GetComponent(ptId, 0), p);
      _Surface->GetPoints()->SetPoint(ptId, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Update point data of extracted point set surface
class UpdateSurfacePointData
{
  vtkPointData *_PointData;
  vtkPolyData  *_Surface;
  vtkDataArray *_PtIds;
  int           _MaxNumberOfComponents;

public:

  UpdateSurfacePointData(vtkPointData *pd, vtkPolyData *surface, vtkIdTypeArray *ptIds)
  :
    _PointData(pd), _Surface(surface), _PtIds(ptIds)
  {
    mirtkAssert(_PointData->GetNumberOfArrays() == _Surface->GetPointData()->GetNumberOfArrays(),
                "surface has expected number of point data arrays");
    _MaxNumberOfComponents = 0;
    for (vtkIdType i = 0; i < _PointData->GetNumberOfArrays(); ++i) {
      int n = _PointData->GetArray(i)->GetNumberOfComponents();
      mirtkAssert(_Surface->GetPointData()->GetArray(i)->GetNumberOfComponents() == n,
                  "surface arrays have identical number of components");
      if (n > _MaxNumberOfComponents) _MaxNumberOfComponents = n;
    }
  }

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    vtkIdType origPtId;
    double *tuple = new double[_MaxNumberOfComponents];
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      origPtId = static_cast<vtkIdType>(_PtIds->GetComponent(ptId, 0));
      for (vtkIdType i = 0; i < _PointData->GetNumberOfArrays(); ++i) {
        _PointData->GetArray(i)->GetTuple(origPtId, tuple);
        _Surface->GetPointData()->GetArray(i)->SetTuple(ptId, tuple);
      }
    }
    delete[] tuple;
  }
};

// -----------------------------------------------------------------------------
/// Transform each point by applying the transformation evaluated at this point
struct TransformPoint
{
  vtkPoints            *_InputPoints;
  vtkPoints            *_OutputPoints;
  const Transformation *_Transformation;
  double                _t, _t0;

  void operator ()(const blocked_range<vtkIdType> &ids) const
  {
    double p[3];
    for (vtkIdType id = ids.begin(); id != ids.end(); ++id) {
      _InputPoints->GetPoint(id, p);
      _Transformation->Transform(p[0], p[1], p[2], _t, _t0);
      _OutputPoints->SetPoint(id, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Transform each point by applying the displacement interpolated at this point
struct TransformPointUsingInterpolatedDisplacements
{
  vtkPoints                      *_InputPoints;
  vtkPoints                      *_OutputPoints;
  const DisplacementInterpolator *_Displacement;

  void operator ()(const blocked_range<vtkIdType> &ids) const
  {
    double p[3], d[3];
    for (vtkIdType id = ids.begin(); id != ids.end(); ++id) {
      _InputPoints->GetPoint(id, p);
      memcpy(d, p, 3 * sizeof(double));
      _Displacement->WorldToImage(d[0], d[1], d[2]);
      _Displacement->Evaluate(d, d[0], d[1], d[2]);
      p[0] += d[0], p[1] += d[1], p[2] += d[2];
      _OutputPoints->SetPoint(id, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Rescale points
struct RescalePoints
{
  vtkPoints *_Points;
  double     _Slope;
  double     _Intercept;

  RescalePoints(vtkPoints *points, double m, double t)
  :
    _Points(points), _Slope(m), _Intercept(t)
  {}

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double p[3];
    for (vtkIdType i = re.begin(); i != re.end(); ++i) {
      _Points->GetPoint(i, p);
      p[0] = p[0] * _Slope + _Intercept;
      p[1] = p[1] * _Slope + _Intercept;
      p[2] = p[2] * _Slope + _Intercept;
      _Points->SetPoint(i, p);
    }
  }
};

// -----------------------------------------------------------------------------
/// Rescale tuples of given data array
struct RescaleData
{
  vtkDataArray *_Data;
  double        _Slope;
  double        _Intercept;

  RescaleData(vtkDataArray *data, double m, double t)
  :
    _Data(data), _Slope(m), _Intercept(t)
  {}

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double *v = new double[_Data->GetNumberOfComponents()];
    for (vtkIdType i = re.begin(); i != re.end(); ++i) {
      _Data->GetTuple(i, v);
      for (int c = 0; c < _Data->GetNumberOfComponents(); ++c) {
        v[c] = v[c] * _Slope + _Intercept;
      }
      _Data->SetTuple(i, v);
    }
    delete[] v;
  }
};


} // namespace RegisteredPointSetUtils
using namespace RegisteredPointSetUtils;

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
RegisteredPointSet::RegisteredPointSet(vtkPointSet                *data,
                                       const class Transformation *t)
:
  _InputPointSet           (data),
  _InputSurfacePoints      (NULL),
  _InputTime               (1.0),
  _Time                    (.0),
  _Transformation          (t),
  _CopyAll                 (true),
  _NeighborhoodRadius      (0),
  _SelfUpdate              (true),
  _UpdateSurfaceNormals    (false),
  _UpdateSurfaceFaceNormals(false),
  _ExternalDisplacement    (NULL),
  _Displacement            (NULL)
{
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::CopyAttributes(const RegisteredPointSet &other)
{
  if (_InputSurfacePoints != &_InputPoints) Delete(_InputSurfacePoints);
  Delete(_Displacement);

  _InputPointSet                 = other._InputPointSet;
  _InputSurface                  = other._InputSurface;
  _InputPoints                   = other._InputPoints;
  _InputTime                     = other._InputTime;
  _Time                          = other._Time;
  _Transformation                = other._Transformation;
  _CopyAll                       = other._CopyAll;
  _NeighborhoodRadius            = other._NeighborhoodRadius;
  _AverageInputEdgeLength        = other._AverageInputEdgeLength;
  _AverageInputSurfaceEdgeLength = other._AverageInputSurfaceEdgeLength;
  _InputDiameter                 = other._InputDiameter;
  _Diameter                      = other._Diameter;
  _EdgeTable                     = other._EdgeTable;
  _SurfaceEdgeTable              = other._SurfaceEdgeTable;
  _NodeNeighbors                 = other._NodeNeighbors;
  _SurfaceNodeNeighbors          = other._SurfaceNodeNeighbors;
  _SelfUpdate                    = other._SelfUpdate;
  _UpdateSurfaceNormals          = other._UpdateSurfaceNormals;
  _UpdateSurfaceFaceNormals      = other._UpdateSurfaceFaceNormals;
  _Domain                        = other._Domain;
  _ExternalDisplacement          = other._ExternalDisplacement;

  if (other._InputSurfacePoints) {
    if (other._InputSurfacePoints == &other._InputPoints) {
      _InputSurfacePoints = &_InputPoints;
    } else {
      _InputSurfacePoints = new class PointSet(*other._InputSurfacePoints);
    }
  }
  if (other._OutputPointSet) {
    _OutputPointSet.TakeReference(other._OutputPointSet->NewInstance());
    _OutputPointSet->DeepCopy(other._OutputPointSet);
  } else {
    _OutputPointSet = NULL;
  }
  if (other._OutputSurface) {
    _OutputSurface.TakeReference(other._OutputSurface->NewInstance());
    _OutputSurface->DeepCopy(other._OutputSurface);
  } else {
    _OutputSurface = NULL;
  }
}

// -----------------------------------------------------------------------------
RegisteredPointSet::RegisteredPointSet(const RegisteredPointSet &other)
:
  Object(other),
  _InputSurfacePoints(NULL),
  _Displacement(NULL)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
RegisteredPointSet &RegisteredPointSet::operator =(const RegisteredPointSet &other)
{
  if (this != &other) {
    Object::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
RegisteredPointSet::~RegisteredPointSet()
{
  if (_InputSurfacePoints != &_InputPoints) Delete(_InputSurfacePoints);
  Delete(_Displacement);
}

// -----------------------------------------------------------------------------
void RegisteredPointSet
::Initialize(bool deep_copy_points)
{
  // Clear previous tables
  _EdgeTable = nullptr;
  _SurfaceEdgeTable = nullptr;
  _NodeNeighbors.Clear();
  _SurfaceNodeNeighbors.Clear();

  // Check input
  if (!_InputPointSet) {
    cerr << "RegisteredPointSet::Initialize: Missing input point set!" << endl;
    exit(1);
  }

  // Extract input point set surface
  _InputSurface  = vtkPolyData::SafeDownCast(_InputPointSet);
  _IsSurfaceMesh = (_InputSurface != NULL);
  if (!_IsSurfaceMesh) {
    _InputSurface = DataSetSurface(_InputPointSet, true, true);
  }
  _InputSurface->BuildLinks();

  // Make shallow copy of input point set and its surface
  //
  // Note that the vtkPoints of the point set are replaced by a new instance
  // upon the first Update of a transformed point set. If no transformation is
  // set or if the point set is never updated, no copy of input points is made.
  // This is overridden by the deep_copy_points argument.
  _OutputPointSet.TakeReference(_InputPointSet->NewInstance());
  _OutputPointSet->ShallowCopy(_InputPointSet);
  if (deep_copy_points) {
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->DeepCopy(_InputPointSet->GetPoints());
    _OutputPointSet->SetPoints(points);
  }

  _OutputSurface = vtkPolyData::SafeDownCast(_OutputPointSet);
  if (!_OutputSurface) {
    _OutputSurface.TakeReference(_InputSurface->NewInstance());
    _OutputSurface->ShallowCopy(_InputSurface);
    // ShallowCopy also copies the previously built cells and links
    if (deep_copy_points) {
      vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      points->DeepCopy(_InputSurface->GetPoints());
      _OutputSurface->SetPoints(points);
    }
  }

  // Reset point and cell data
  if (!_CopyAll || !_PointDataToCopy.empty()) {
    _OutputPointSet->GetPointData()->Initialize();
    _OutputSurface ->GetPointData()->Initialize();
  }
  if (!_CopyAll) {
    _OutputPointSet->GetCellData()->Initialize();
    _OutputSurface ->GetCellData()->Initialize();
  }

  // Initialize structures depending on input points
  InputPointsChanged();

  // Mark on demand output properties as invalid
  PointsChanged();

  // Initialize locators
  _PointLocator = vtkSmartPointer<vtkPointLocator>::New();
  _CellLocator  = vtkSmartPointer<vtkCellLocator>::New();
  _PointLocator->SetDataSet(_OutputPointSet);
  _CellLocator ->SetDataSet(_OutputPointSet);

  if (_IsSurfaceMesh) {
    _SurfacePointLocator = _PointLocator;
    _SurfaceCellLocator  = _CellLocator;
  } else {
    _SurfacePointLocator = vtkSmartPointer<vtkPointLocator>::New();
    _SurfaceCellLocator  = vtkSmartPointer<vtkCellLocator>::New();
    _SurfacePointLocator->SetDataSet(_OutputSurface);
    _SurfaceCellLocator ->SetDataSet(_OutputSurface);
  }
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::PointsChanged()
{
  _UpdateSurfaceNormals     = true;
  _UpdateSurfaceFaceNormals = true;
  _SurfaceArea              = .0;

  _OutputPointSet->Modified();
  _OutputSurface ->Modified();

  double bounds[6];
  _OutputPointSet->GetBounds(bounds);
  _Diameter = sqrt(pow(bounds[1] - bounds[0], 2) +
                        pow(bounds[3] - bounds[2], 2) +
                        pow(bounds[5] - bounds[4], 2));
  if (!_IsSurfaceMesh) _OutputSurface->ComputeBounds();
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::InputPointsChanged()
{
  // Copy input points to MIRTK structure needed by the energy terms to compute
  // the gradient w.r.t. the transformation parameters
  // Transformation::ParametricGradient and thus used by energy terms
  this->GetInputPoints(_InputPoints);
  if (_InputSurface == _InputPointSet) {
    if (_InputSurfacePoints != &_InputPoints) {
      delete _InputSurfacePoints;
      _InputSurfacePoints = &_InputPoints;
    }
  } else {
    if (!_InputSurfacePoints || _InputSurfacePoints == &_InputPoints) {
      _InputSurfacePoints = new class PointSet;
    }
    this->GetInputSurfacePoints(*_InputSurfacePoints);
  }

  // Precompute average input point set edge length
  {
    SharedPtr<const EdgeTable> edgeTable = SharedEdgeTable();
    _AverageInputEdgeLength = AverageEdgeLength(_InputPointSet->GetPoints(), *edgeTable);
  }
  if (_InputSurface == _InputPointSet) {
    _AverageInputSurfaceEdgeLength = _AverageInputEdgeLength;
  } else {
    SharedPtr<const EdgeTable> edgeTable = SharedSurfaceEdgeTable();
    _AverageInputSurfaceEdgeLength = AverageEdgeLength(_InputSurface->GetPoints(), *edgeTable);
  }

  // Invalidate input surface area
  _InputSurfaceArea = .0;

  // Recompute input bounds
  double bounds[6];
  _InputPointSet->GetBounds(bounds);
  _InputDiameter = sqrt(pow(bounds[1] - bounds[0], 2) +
                             pow(bounds[3] - bounds[2], 2) +
                             pow(bounds[5] - bounds[4], 2));
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::BuildEdgeTables()
{
  if (!_EdgeTable) {
    _EdgeTable = NewShared<EdgeTable>(_InputPointSet);
  }
  if (_InputSurface == _InputPointSet) {
    _SurfaceEdgeTable = _EdgeTable;
  } else if (!_SurfaceEdgeTable) {
    _SurfaceEdgeTable = NewShared<EdgeTable>(_InputSurface);
  }
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::BuildNeighborhoodTables(int n)
{
  if (n < _NeighborhoodRadius) n = _NeighborhoodRadius;
  if (_NodeNeighbors.Rows() == 0 || n > _NodeNeighbors.Maximum()) {
    SharedPtr<const EdgeTable> edgeTable = SharedEdgeTable();
//    if (n <= 0) {
//      double r = _NeighborhoodRadius * _AverageInputEdgeLength;
//      _NodeNeighbors.Initialize(_InputPointSet, r, edgeTable());
//    } else {
//      _NodeNeighbors.Initialize(_InputPointSet, n, edgeTable());
//    }
    _NodeNeighbors.Initialize(_InputPointSet, n, edgeTable.get());
  }
  if (_InputSurface != _InputPointSet) {
    if (_SurfaceNodeNeighbors.Rows() == 0 || n > _SurfaceNodeNeighbors.Maximum()) {
      SharedPtr<const EdgeTable> edgeTable = SharedSurfaceEdgeTable();
//      if (n <= 0) {
//        double r = _NeighborhoodRadius * _AverageInputSurfaceEdgeLength;
//        _SurfaceNodeNeighbors.Initialize(_InputSurface, r, edgeTable.get());
//      } else {
//        _SurfaceNodeNeighbors.Initialize(_InputSurface, n, edgeTable.get());
//      }
      _SurfaceNodeNeighbors.Initialize(_InputSurface, n, edgeTable.get());
    }
  }
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::BuildLocators()
{
  _PointLocator->BuildLocator();
  _CellLocator ->BuildLocator();
  if (!_IsSurfaceMesh) {
    _SurfacePointLocator->BuildLocator();
    _SurfaceCellLocator ->BuildLocator();
  }
}

// =============================================================================
// Copy points to PointSet structure
// =============================================================================

// -----------------------------------------------------------------------------
vtkIdTypeArray *RegisteredPointSet::OriginalSurfacePointIds() const
{
  return vtkIdTypeArray::SafeDownCast(_InputSurface->GetPointData()->GetArray("vtkOriginalPointIds"));
}

// -----------------------------------------------------------------------------
vtkIdTypeArray *RegisteredPointSet::OriginalSurfaceCellIds() const
{
  return vtkIdTypeArray::SafeDownCast(_InputSurface->GetCellData()->GetArray("vtkOriginalCellIds"));
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::GetInputPoints(class PointSet &pset) const
{
  const int n = this->NumberOfPoints();
  pset.Reserve(n);
  pset.Resize(n);
  pset.ShrinkToFit();
  CopyVtkPointsToIrtkPointSet copy(_InputPointSet->GetPoints(), pset);
  parallel_for(blocked_range<int>(0, pset.Size()), copy);
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::GetInputSurfacePoints(class PointSet &pset) const
{
  const int n = this->NumberOfSurfacePoints();
  pset.Reserve(n);
  pset.Resize (n);
  pset.ShrinkToFit();
  CopyVtkPointsToIrtkPointSet copy(_InputSurface->GetPoints(), pset);
  parallel_for(blocked_range<int>(0, pset.Size()), copy);
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::GetPoints(class PointSet &pset) const
{
  const int n = this->NumberOfPoints();
  pset.Reserve(n);
  pset.Resize(n);
  pset.ShrinkToFit();
  CopyVtkPointsToIrtkPointSet copy(_OutputPointSet->GetPoints(), pset);
  parallel_for(blocked_range<int>(0, pset.Size()), copy);
}

// -----------------------------------------------------------------------------
vtkDataArray *RegisteredPointSet::InitialStatus() const
{
  return _OutputPointSet->GetPointData()->GetArray("InitialStatus");
}

// -----------------------------------------------------------------------------
vtkDataArray *RegisteredPointSet::Status() const
{
  return _OutputPointSet->GetPointData()->GetArray("Status");
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::GetSurfacePoints(class PointSet &pset) const
{
  const int n = this->NumberOfSurfacePoints();
  pset.Reserve(n);
  pset.Resize(n);
  pset.ShrinkToFit();
  CopyVtkPointsToIrtkPointSet copy(_OutputSurface->GetPoints(), pset);
  parallel_for(blocked_range<int>(0, pset.Size()), copy);
}

// -----------------------------------------------------------------------------
vtkDataArray *RegisteredPointSet::SurfaceNormals() const
{
  if (_OutputSurface->GetPointData()->GetNormals() == NULL || _UpdateSurfaceNormals) {
    vtkSmartPointer<vtkPolyDataNormals> filter = vtkSmartPointer<vtkPolyDataNormals>::New();
    SetVTKInput(filter, _OutputSurface);
    filter->SplittingOff();
    filter->ComputePointNormalsOn();
    filter->ComputeCellNormalsOff();
    filter->ConsistencyOff();
    filter->AutoOrientNormalsOff();
    filter->FlipNormalsOff();
    filter->Update();
    vtkDataArray *normals = filter->GetOutput()->GetPointData()->GetNormals();
    _OutputSurface->GetPointData()->SetNormals(normals);
  }
  return _OutputSurface->GetPointData()->GetNormals();
}

// -----------------------------------------------------------------------------
vtkDataArray *RegisteredPointSet::SurfaceFaceNormals() const
{
  if (_OutputSurface->GetCellData()->GetNormals() == NULL || _UpdateSurfaceFaceNormals) {
    vtkSmartPointer<vtkPolyDataNormals> filter = vtkSmartPointer<vtkPolyDataNormals>::New();
    SetVTKInput(filter, _OutputSurface);
    filter->SplittingOff();
    filter->ComputePointNormalsOff();
    filter->ComputeCellNormalsOn();
    filter->ConsistencyOff();
    filter->AutoOrientNormalsOff();
    filter->FlipNormalsOff();
    filter->Update();
    vtkDataArray *normals = filter->GetOutput()->GetCellData()->GetNormals();
    _OutputSurface->GetCellData()->SetNormals(normals);
  }
  return _OutputSurface->GetCellData()->GetNormals();
}

// -----------------------------------------------------------------------------
SharedPtr<const RegisteredPointSet::EdgeTable> RegisteredPointSet::SharedEdgeTable() const
{
  if (!_EdgeTable) {
    _EdgeTable = NewShared<EdgeTable>(_InputPointSet);
  }
  return _EdgeTable;
}

// -----------------------------------------------------------------------------
const RegisteredPointSet::EdgeTable *RegisteredPointSet::Edges() const
{
  return SharedEdgeTable().get();
}

// -----------------------------------------------------------------------------
SharedPtr<const RegisteredPointSet::EdgeTable> RegisteredPointSet::SharedSurfaceEdgeTable() const
{
  if (!_SurfaceEdgeTable) {
    if (_InputSurface == _InputPointSet) {
      _SurfaceEdgeTable = SharedEdgeTable();
    } else {
      _SurfaceEdgeTable = NewShared<EdgeTable>(_InputSurface);
    }
  }
  return _SurfaceEdgeTable;
}

// -----------------------------------------------------------------------------
const RegisteredPointSet::EdgeTable *RegisteredPointSet::SurfaceEdges() const
{
  return SharedSurfaceEdgeTable().get();
}

// -----------------------------------------------------------------------------
const RegisteredPointSet::NodeNeighbors *RegisteredPointSet::Neighbors(int n) const
{
  if (_NodeNeighbors.Rows() == 0 || n > _NodeNeighbors.Maximum()) {
    SharedPtr<const EdgeTable> edgeTable = SharedEdgeTable();
//    if (n <= 0) {
//      double r = _NeighborhoodRadius * _AverageInputEdgeLength;
//      _NodeNeighbors.Initialize(_InputPointSet, r, edgeTable.get());
//    } else {
//      _NodeNeighbors.Initialize(_InputPointSet, n, edgeTable.get());
//    }
    _NodeNeighbors.Initialize(_InputPointSet, max(n, _NeighborhoodRadius), edgeTable.get());
  }
  return &_NodeNeighbors;
}

// -----------------------------------------------------------------------------
const RegisteredPointSet::NodeNeighbors *RegisteredPointSet::SurfaceNeighbors(int n) const
{
  if (_InputSurface == _InputPointSet) return Neighbors(n);
  if (_SurfaceNodeNeighbors.Rows() == 0 || n > _SurfaceNodeNeighbors.Maximum()) {
    SharedPtr<const EdgeTable> edgeTable = SharedSurfaceEdgeTable();
//    if (n <= 0) {
//      double r = _NeighborhoodRadius * _AverageInputSurfaceEdgeLength;
//      _SurfaceNodeNeighbors.Initialize(_InputSurface, r, edgeTable.get());
//    } else {
//      _SurfaceNodeNeighbors.Initialize(_InputSurface, n, edgeTable.get());
//    }
    _SurfaceNodeNeighbors.Initialize(_InputSurface, max(n, _NeighborhoodRadius), edgeTable.get());
  }
  return &_SurfaceNodeNeighbors;
}

// -----------------------------------------------------------------------------
vtkDataArray *RegisteredPointSet::InitialSurfaceStatus() const
{
  return _OutputSurface->GetPointData()->GetArray("InitialStatus");
}

// -----------------------------------------------------------------------------
vtkDataArray *RegisteredPointSet::SurfaceStatus() const
{
  return _OutputSurface->GetPointData()->GetArray("Status");
}

// -----------------------------------------------------------------------------
double RegisteredPointSet::InputSurfaceArea() const
{
  if (_InputSurfaceArea <= .0) _InputSurfaceArea = Area(_InputSurface);
  return _InputSurfaceArea;
}

// -----------------------------------------------------------------------------
double RegisteredPointSet::SurfaceArea() const
{
  if (_SurfaceArea <= .0) _SurfaceArea = Area(_OutputSurface);
  return _SurfaceArea;
}

// -----------------------------------------------------------------------------
vtkAbstractPointLocator *RegisteredPointSet::PointLocator() const
{
  _PointLocator->BuildLocator();
  return _PointLocator;
}

// -----------------------------------------------------------------------------
vtkAbstractCellLocator *RegisteredPointSet::CellLocator() const
{
  _CellLocator->BuildLocator();
  return _CellLocator;
}

// -----------------------------------------------------------------------------
vtkAbstractPointLocator *RegisteredPointSet::SurfacePointLocator() const
{
  _SurfacePointLocator->BuildLocator();
  return _SurfacePointLocator;
}

// -----------------------------------------------------------------------------
vtkAbstractCellLocator *RegisteredPointSet::SurfaceCellLocator() const
{
  _SurfaceCellLocator->BuildLocator();
  return _SurfaceCellLocator;
}

// =============================================================================
// Update
// =============================================================================

// -----------------------------------------------------------------------------
void RegisteredPointSet::Update(bool force)
{
  // Do nothing if self-update is disabled (i.e., external process is
  // responsible for update of registered data)
  if (!force && !_SelfUpdate) return;
  MIRTK_START_TIMING();

  // ---------------------------------------------------------------------------
  // Transform points
  if (_Transformation) {
    // Update points of output point set
    vtkPoints *inputPoints  = _InputPointSet ->GetPoints();
    vtkPoints *outputPoints = _OutputPointSet->GetPoints();
    const vtkIdType npoints = inputPoints->GetNumberOfPoints();
    if (outputPoints == inputPoints) {
      vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      _OutputPointSet->SetPoints(points);
      outputPoints = points;
    }
    outputPoints->SetNumberOfPoints(npoints);
    if (_ExternalDisplacement) {
      MIRTK_START_TIMING();
      DisplacementInterpolator disp;
      disp.Input(_ExternalDisplacement);
      disp.Initialize();
      TransformPointUsingInterpolatedDisplacements transform;
      transform._InputPoints  = inputPoints;
      transform._OutputPoints = outputPoints;
      transform._Displacement = &disp;
      parallel_for(blocked_range<vtkIdType>(0, npoints), transform);
      MIRTK_DEBUG_TIMING(7, "transforming points");
    } else if (_Transformation->RequiresCachingOfDisplacements() && _Domain) {
      MIRTK_START_TIMING();
      if (!_Displacement) _Displacement = new GenericImage<double>();
      _Displacement->Initialize(_Domain, 3);
      _Transformation->Displacement(*_Displacement, _Time, _InputTime);
      DisplacementInterpolator disp;
      disp.Input(_Displacement);
      disp.Initialize();
      MIRTK_DEBUG_TIMING(7, "caching displacements");
      MIRTK_RESET_TIMING();
      TransformPointUsingInterpolatedDisplacements transform;
      transform._InputPoints  = inputPoints;
      transform._OutputPoints = outputPoints;
      transform._Displacement = &disp;
      parallel_for(blocked_range<vtkIdType>(0, npoints), transform);
      MIRTK_DEBUG_TIMING(7, "transforming points");
    } else {
      MIRTK_START_TIMING();
      TransformPoint transform;
      transform._InputPoints    = inputPoints;
      transform._OutputPoints   = outputPoints;
      transform._Transformation = _Transformation;
      transform._t              = _Time;
      transform._t0             = _InputTime;
      parallel_for(blocked_range<vtkIdType>(0, npoints), transform);
      MIRTK_DEBUG_TIMING(7, "transforming points");
    }
    // Update points of output point set surface
    if (_OutputSurface->GetPoints() != outputPoints) {
      if (_OutputSurface->GetPoints() == _InputSurface->GetPoints()) {
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        _OutputSurface->SetPoints(points);
      }
      const vtkIdType nspoints = _InputSurface->GetNumberOfPoints();
      _OutputSurface->GetPoints()->SetNumberOfPoints(nspoints);
      UpdateSurfacePoints update(outputPoints, _OutputSurface, OriginalSurfacePointIds());
      parallel_for(blocked_range<vtkIdType>(0, nspoints), update);
    }
    // Invalidate cached attributes which are recomputed on demand only
    PointsChanged();
  }

  // ---------------------------------------------------------------------------
  // Copy and optionally normalize point data
  if (!_PointDataToCopy.empty()) {
    MIRTK_START_TIMING();
    // Update point data of output point set
    vtkPointData * const inputPD  = _InputPointSet ->GetPointData();
    vtkPointData * const outputPD = _OutputPointSet->GetPointData();
    vtkDataArray *input_array;
    vtkSmartPointer<vtkDataArray> output_array;
    Array<ScalingFunction>::iterator it;
    for (it = _PointDataToCopy.begin(); it != _PointDataToCopy.end(); ++it) {
      if (it->_Slope == .0) continue;
      if (it->_InputIndex < 0 || it->_InputIndex >= inputPD->GetNumberOfArrays()) {
        cerr << "RegisteredPointSet::Update: Invalid point data index: " << it->_InputIndex << endl;
        exit(1);
      }
      int attr = inputPD->IsArrayAnAttribute(it->_InputIndex);
      input_array = inputPD->GetArray(it->_InputIndex);
      if (0 <= it->_OutputIndex && it->_OutputIndex < outputPD->GetNumberOfArrays()) {
        output_array = outputPD->GetArray(it->_OutputIndex);
      } else {
        if (it->_Slope == 1.0 && it->_Intercept == .0) {
          output_array.TakeReference(input_array->NewInstance());
        } else {
          output_array = vtkSmartPointer<vtkFloatArray>::New();
        }
        output_array->SetName(input_array->GetName());
        it->_OutputIndex = -1;
      }
      output_array->DeepCopy(input_array);
      if (attr < 0) {
        parallel_for(blocked_range<vtkIdType>(0, outputPD->GetNumberOfTuples()),
                     RescaleData(output_array, it->_Slope, it->_Intercept));
      }
      if (it->_OutputIndex == -1) {
        it->_OutputIndex = outputPD->AddArray(output_array);
        if (attr >= 0) outputPD->SetActiveAttribute(it->_OutputIndex, attr);
      }
    }
    // Update point data of output point set surface
    if (_OutputSurface->GetPointData() != outputPD) {
      UpdateSurfacePointData update(outputPD, _OutputSurface, OriginalSurfacePointIds());
      parallel_for(blocked_range<vtkIdType>(0, _OutputSurface->GetNumberOfPoints()), update);
    }
    MIRTK_DEBUG_TIMING(7, "copying and rescaling point data");
  }

  MIRTK_DEBUG_TIMING(2, "update of" << (_Transformation ? " moving " : " fixed ")
                                    << (_IsSurfaceMesh  ? "surface"  : "point set"));
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
const char *RegisteredPointSet::DefaultExtension() const
{
  return mirtk::DefaultExtension(_InputPointSet);
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::Write(const char       *fname,
                               vtkAbstractArray *pointdata,
                               vtkAbstractArray *celldata) const
{
  Write(fname, &pointdata, 1, &celldata, 1);
}

// -----------------------------------------------------------------------------
void RegisteredPointSet::Write(const char        *fname,
                               vtkAbstractArray **pointdata, int npointdata,
                               vtkAbstractArray **celldata,  int ncelldata) const
{
  // Add point data
  int *pointdataidx = NULL;
  if (npointdata > 0) {
    pointdataidx = new int[npointdata];
    for (int i = 0; i < npointdata; ++i) {
      pointdataidx[i] = _OutputPointSet->GetPointData()->AddArray(pointdata[i]);
    }
  }
  // Add cell data
  int *celldataidx = NULL;
  if (ncelldata > 0) {
    celldataidx = new int[ncelldata];
    for (int i = 0; i < ncelldata; ++i) {
      celldataidx[i] = _OutputPointSet->GetCellData()->AddArray(celldata[i]);
    }
  }
  // Write data set with additional point and cell data
  WritePointSet(fname, _OutputPointSet);
  // Remove point data
  for (int i = 0; i < npointdata; ++i) {
    _OutputPointSet->GetPointData()->RemoveArray(pointdataidx[i]);
  }
  delete[] pointdataidx;
  // Remove cell data
  for (int i = 0; i < ncelldata; ++i) {
    _OutputPointSet->GetCellData()->RemoveArray(celldataidx[i]);
  }
  delete[] celldataidx;
}


} // namespace mirtk
