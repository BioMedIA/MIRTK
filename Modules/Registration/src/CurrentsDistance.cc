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

#include "mirtk/CurrentsDistance.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Vector3D.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/ObjectFactory.h"
#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkCellArray.h"
#include "vtkFloatArray.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
// Due to a bug in vtkKdTreePointLocator, calling BuildLocator
// is not sufficient to make FindClosestPoint thread-safe as it does
// not call vtkBSPIntersections::BuildRegionsList
// (cf. http://www.vtk.org/Bug/view.php?id=15206 ).
// Therefore, use Octree instead of Kd-tree in the meantime.
#include "vtkOctreePointLocator.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(CurrentsDistance);


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
CurrentsDistance::CurrentsDistance(const char *name, double weight)
:
  PointSetDistance(name, weight),
  _Sigma(-0.05),
  _Symmetric(true),
  _TargetNormSquared(.0)
{
}

// -----------------------------------------------------------------------------
CurrentsDistance::CurrentsDistance(const CurrentsDistance &other)
:
  PointSetDistance(other),
  _TargetNormSquared(other._TargetNormSquared)
{
  if (other._TargetCurrent) {
    _TargetCurrent = vtkSmartPointer<vtkPolyData>::New();
    _TargetCurrent->DeepCopy(other._TargetCurrent);
  }
  if (other._SourceCurrent) {
    _SourceCurrent = vtkSmartPointer<vtkPolyData>::New();
    _SourceCurrent->DeepCopy(other._SourceCurrent);
  }
}

// -----------------------------------------------------------------------------
CurrentsDistance &CurrentsDistance::operator =(const CurrentsDistance &other)
{
  PointSetDistance::operator =(other);
  _TargetNormSquared = other._TargetNormSquared;
  if (other._TargetCurrent) {
    _TargetCurrent = vtkSmartPointer<vtkPolyData>::New();
    _TargetCurrent->DeepCopy(other._TargetCurrent);
  } else {
    _TargetCurrent = NULL;
  }
  if (other._SourceCurrent) {
    _SourceCurrent = vtkSmartPointer<vtkPolyData>::New();
    _SourceCurrent->DeepCopy(other._SourceCurrent);
  } else {
    _SourceCurrent = NULL;
  }
  return *this;
}

// -----------------------------------------------------------------------------
CurrentsDistance::~CurrentsDistance()
{
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> CurrentsDistance::SurfaceToCurrent(vtkPolyData *data)
{
  MIRTK_START_TIMING();

  vtkCellArray *faces = data ->GetPolys();
  vtkIdType num_faces = faces->GetNumberOfCells();

  vtkSmartPointer<vtkFloatArray> normals = vtkSmartPointer<vtkFloatArray>::New();
  vtkSmartPointer<vtkPoints>     centers = vtkSmartPointer<vtkPoints>    ::New();

  normals->SetName("normals");
  normals->SetNumberOfComponents(3);
  normals->SetNumberOfTuples(num_faces);
  centers->SetNumberOfPoints(num_faces);

  faces->InitTraversal();
  vtkIdType npts, *id;
  double    v1[3], v2[3], v3[3], e2[3], e3[3], c[3], n[3];
  for (vtkIdType i = 0; i < num_faces; ++i) {
    faces->GetNextCell(npts, id);
    if (npts != 3) {
      cerr << "CurrentsDistance::SurfaceToCurrent: Surface cells must have three points each!" << endl;
      exit(1);
    }
    data->GetPoint(id[0], v1);
    data->GetPoint(id[1], v2);
    data->GetPoint(id[2], v3);
    for (int d = 0; d < 3; ++d) {
      e2[d] = v3[d] - v1[d];
      e3[d] = v2[d] - v1[d];
      c [d] = (v1[d] + v2[d] + v3[d]) / 3.0;
    }
    vtkMath::Cross(e2, e3, n);
    centers->SetPoint (i, c);
    normals->SetTuple3(i, 0.5 * n[0], 0.5 * n[1], 0.5 * n[2]);
  }

  vtkSmartPointer<vtkPolyData> current = vtkSmartPointer<vtkPolyData>::New();
  current->SetPoints(centers);
  current->GetPointData()->AddArray(normals);

  MIRTK_DEBUG_TIMING(3, "conversion of surface mesh to current");
  return current;
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPolyData> CurrentsDistance::ToCurrent(vtkPointSet *data)
{
  // TODO: Support also points (0D) and curves (1D)
  vtkPolyData *surface = vtkPolyData::SafeDownCast(data);
  if (surface->GetPolys()->GetNumberOfCells() == 0) {
    cerr << "CurrentsDistance::ToCurrent: Only surface meshes supported at the moment" << endl;
    exit(1);
  }
  return SurfaceToCurrent(surface);
}

// -----------------------------------------------------------------------------
class CurrentsDistanceDotProduct
{
private:

  vtkPoints               *_CentersA;
  vtkFloatArray           *_WeightsA;
  vtkPoints               *_CentersB;
  vtkFloatArray           *_WeightsB;
  vtkAbstractPointLocator *_LocatorB;
  double                   _Variance;
  double                   _Radius;
  vtkDataArray            *_Value;
  double                   _Sum;

public:

  CurrentsDistanceDotProduct(double sigma)
  :
    _CentersA(NULL), _WeightsA(NULL),
    _CentersB(NULL), _WeightsB(NULL), _LocatorB(NULL),
    _Variance(sigma * sigma), _Radius(2.5 * sigma),
    _Value(NULL), _Sum(.0)
  {}

  inline double Evaluate(vtkPolyData *a, vtkPolyData *b, vtkDataArray *value = NULL)
  {
    MIRTK_START_TIMING();
    // Get centers of n-currents
    _CentersA = a->GetPoints();
    _CentersB = b->GetPoints();
    // Get normals of 2-currents (i.e., surfaces)
    _WeightsA = vtkFloatArray::SafeDownCast(a->GetPointData()->GetArray("normals"));
    _WeightsB = vtkFloatArray::SafeDownCast(b->GetPointData()->GetArray("normals"));
    // Get segments of 1-currents (i.e., curves)
    if (!_WeightsA && !_WeightsB) {
      _WeightsA = vtkFloatArray::SafeDownCast(a->GetPointData()->GetArray("segments"));
      _WeightsB = vtkFloatArray::SafeDownCast(b->GetPointData()->GetArray("segments"));
    } else if ((_WeightsA && !_WeightsB) || (!_WeightsA && _WeightsB)) {
      cerr << "Cannot compute inner product between different types of currents" << endl;
      exit(1);
    }
    // Get weights of 0-currents (i.e., points)
    if (!_WeightsA && !_WeightsB) {
      _WeightsA = vtkFloatArray::SafeDownCast(a->GetPointData()->GetArray("weights"));
      _WeightsB = vtkFloatArray::SafeDownCast(b->GetPointData()->GetArray("weights"));
    } else if ((_WeightsA && !_WeightsB) || (!_WeightsA && _WeightsB)) {
      cerr << "Cannot compute inner product between different types of currents" << endl;
      exit(1);
    }
    // Initialize point locator
    _LocatorB = vtkOctreePointLocator::New();
    _LocatorB->SetDataSet(b);
    _LocatorB->BuildLocator();
    // Evaluate inner product
    _Value = value;
    _Sum   = .0;
    blocked_range<vtkIdType> cellsA(0, _CentersA->GetNumberOfPoints());
    parallel_reduce(cellsA, *this);
    // Free point locator
    _LocatorB->Delete(), _LocatorB = NULL;
    MIRTK_DEBUG_TIMING(3, "evaluation of dot product of currents");
    return _Sum;
  }

/* private: */

  CurrentsDistanceDotProduct(CurrentsDistanceDotProduct &other, split)
  :
    _CentersA(other._CentersA),
    _WeightsA(other._WeightsA),
    _CentersB(other._CentersB),
    _WeightsB(other._WeightsB),
    _LocatorB(other._LocatorB),
    _Variance(other._Variance),
    _Radius  (other._Radius),
    _Value   (other._Value),
    _Sum     (.0)
  {}

  void join(CurrentsDistanceDotProduct &other)
  {
    _Sum += other._Sum;
  }

  inline double EvaluateKernel(double ca[3], double cb[3])
  {
    return exp(- vtkMath::Distance2BetweenPoints(ca, cb) / _Variance);
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    // In case of point clouds, the _Weights(A|B) arrays contain
    // scalar tuples only, i.e., GetTuple does not change the second
    // and third component of da and db, respectively. The dot product
    // <da, db> is thus equal to the scalar product of the point weights.
    double ca[3], cb[3], da[3] = {1.0, .0, .0}, db[3] = {1.0, .0, .0};
    for (vtkIdType i = re.begin(); i != re.end(); ++i) {
      _CentersA->GetPoint(i, ca);
      _WeightsA->GetTuple(i, da);
      _LocatorB->FindPointsWithinRadius(_Radius, ca, ids);
      double value = .0;
      for (vtkIdType k = 0; k < ids->GetNumberOfIds(); ++k) {
        vtkIdType j = ids->GetId(k);
        _CentersB->GetPoint(j, cb);
        _WeightsB->GetTuple(j, db);
        value += EvaluateKernel(ca, cb) * vtkMath::Dot(da, db);
      }
      if (_Value) _Value->SetTuple1(i, _Value->GetTuple1(i) + value);
      _Sum += value;
    }
  }
};

// -----------------------------------------------------------------------------
class CurrentsDistanceGradient
{
private:

  vtkPointSet             *_SurfaceA;
  vtkPoints               *_CentersA;
  vtkFloatArray           *_WeightsA;
  vtkAbstractPointLocator *_LocatorA;
  vtkPointSet             *_SurfaceB;
  vtkPoints               *_CentersB;
  vtkFloatArray           *_WeightsB;
  vtkAbstractPointLocator *_LocatorB;
  double                   _Variance;
  double                   _Radius;
  Vector3D<double>        *_Gradient;

public:

  CurrentsDistanceGradient(double sigma)
  :
    _SurfaceA(NULL), _CentersA(NULL), _WeightsA(NULL), _LocatorA(NULL),
    _SurfaceB(NULL), _CentersB(NULL), _WeightsB(NULL), _LocatorB(NULL),
    _Variance(sigma * sigma), _Radius(2.5 * sigma), _Gradient(NULL)
  {}

  inline void EvaluateGradient(vtkPointSet *sa, vtkPolyData *ca,
                               vtkPointSet *sb, vtkPolyData *cb,
                               Vector3D<double> *g)
  {
    MIRTK_START_TIMING();
    _SurfaceA = sa;
    _SurfaceB = sb;
    // Get centers of n-currents
    _CentersA = ca->GetPoints();
    _CentersB = cb->GetPoints();
    // Get normals of 2-currents (i.e., surfaces)
    _WeightsA = vtkFloatArray::SafeDownCast(ca->GetPointData()->GetArray("normals"));
    _WeightsB = vtkFloatArray::SafeDownCast(cb->GetPointData()->GetArray("normals"));
    // Get segments of 1-currents (i.e., curves)
    if (!_WeightsA && !_WeightsB) {
      _WeightsA = vtkFloatArray::SafeDownCast(ca->GetPointData()->GetArray("segments"));
      _WeightsB = vtkFloatArray::SafeDownCast(cb->GetPointData()->GetArray("segments"));
    } else if ((_WeightsA && !_WeightsB) || (!_WeightsA && _WeightsB)) {
      cerr << "Cannot compute inner product between different types of currents" << endl;
      exit(1);
    }
    // Get weights of 0-currents (i.e., points)
    if (!_WeightsA && !_WeightsB) {
      _WeightsA = vtkFloatArray::SafeDownCast(ca->GetPointData()->GetArray("weights"));
      _WeightsB = vtkFloatArray::SafeDownCast(cb->GetPointData()->GetArray("weights"));
    } else if ((_WeightsA && !_WeightsB) || (!_WeightsA && _WeightsB)) {
      cerr << "Cannot compute inner product between different types of currents" << endl;
      exit(1);
    }
    // Initialize point locators
    _LocatorA = vtkOctreePointLocator::New();
    _LocatorA->SetDataSet(ca);
    _LocatorA->BuildLocator();
    _LocatorB = vtkOctreePointLocator::New();
    _LocatorB->SetDataSet(cb);
    _LocatorB->BuildLocator();
    // Evaluate gradient of currents distance measure
    _Gradient = g;
    for (int i = 0; i < _SurfaceA->GetNumberOfPoints(); ++i) {
      _Gradient[i]._x = _Gradient[i]._y = _Gradient[i]._z = .0;
    }
    blocked_range<vtkIdType> cellsA(0, _SurfaceA->GetNumberOfCells());
    parallel_reduce(cellsA, *this);
    // Free point locators
    _LocatorA->Delete(), _LocatorA = NULL;
    _LocatorB->Delete(), _LocatorB = NULL;
    MIRTK_DEBUG_TIMING(3, "evaluation of gradient of currents distance");
  }

/* private: */

  CurrentsDistanceGradient(CurrentsDistanceGradient &other)
  :
    _SurfaceA(other._SurfaceA),
    _CentersA(other._CentersA),
    _WeightsA(other._WeightsA),
    _LocatorA(other._LocatorA),
    _SurfaceB(other._SurfaceB),
    _CentersB(other._CentersB),
    _WeightsB(other._WeightsB),
    _LocatorB(other._LocatorB),
    _Variance(other._Variance),
    _Radius  (other._Radius),
    _Gradient(other._Gradient)
  {}

  CurrentsDistanceGradient(CurrentsDistanceGradient &lhs, split)
  :
    _SurfaceA(lhs._SurfaceA),
    _CentersA(lhs._CentersA),
    _WeightsA(lhs._WeightsA),
    _LocatorA(lhs._LocatorA),
    _SurfaceB(lhs._SurfaceB),
    _CentersB(lhs._CentersB),
    _WeightsB(lhs._WeightsB),
    _LocatorB(lhs._LocatorB),
    _Variance(lhs._Variance),
    _Radius  (lhs._Radius)
  {
    CAllocate(_Gradient, _SurfaceA->GetNumberOfPoints());
  }

  void join(CurrentsDistanceGradient &rhs)
  {
    for (int i = 0; i < _SurfaceA->GetNumberOfPoints(); ++i) {
      _Gradient[i] += rhs._Gradient[i];
    }
    Deallocate(rhs._Gradient);
  }

  inline double EvaluateKernel(double ca[3], double cb[3])
  {
    return exp(- vtkMath::Distance2BetweenPoints(ca, cb) / _Variance);
  }

  inline void EvaluateKernelGradient(double g[3], double ca[3], double cb[3])
  {
    const double w = -2.0 * EvaluateKernel(ca, cb) / _Variance;
    g[0] = w * (ca[0] - cb[0]);
    g[1] = w * (ca[1] - cb[1]);
    g[2] = w * (ca[2] - cb[2]);
  }

  inline void Add(Vector3D<double> &v, double g[3])
  {
    v._x += g[0], v._y += g[1], v._z += g[2];
  }

  void operator ()(const blocked_range<vtkIdType> &re)
  {
    vtkIdType j, i1, i2, i3;       // vertex indices
    double    v1[3], v2[3], v3[3]; // vertex coordinates
    double    e1[3], e2[3], e3[3]; // edge vectors
    double    c1[3], c2[3];        // center coordinates
    double    n1[3], n2[3];        // surface normals
    double    kws[3], dks[3][3], kwt[3], dkt[3][3], kw[3];
    double    w, g[3];

    const double _2over3 = 2.0 / 3.0;

    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();

    // Loop over transformed triangles
    for (vtkIdType i = re.begin(); i != re.end(); ++i) {
      // Get vertex indices
      _SurfaceA->GetCellPoints(i, ids);
      i1 = ids->GetId(0);
      i2 = ids->GetId(1);
      i3 = ids->GetId(2);
      Vector3D<double> &g1 = _Gradient[i1];
      Vector3D<double> &g2 = _Gradient[i2];
      Vector3D<double> &g3 = _Gradient[i3];
      // Get coordinates of vertices
      _SurfaceA->GetPoint(i1, v1);
      _SurfaceA->GetPoint(i2, v2);
      _SurfaceA->GetPoint(i3, v3);
      // Get center and normal
      _CentersA->GetPoint(i, c1);
      _WeightsA->GetTuple(i, n1);
      // Compute kds = KtauS and dks = gradKtauS.transpose()
      // (cf. Deformetrica 2.0 OrientedSurfaceMesh::ComputeMatchGradient)
      _LocatorA->FindPointsWithinRadius(_Radius, c1, ids);
      memset(kws, 0, 3 * sizeof(double));
      memset(dks, 0, 9 * sizeof(double));
      for (vtkIdType k = 0; k < ids->GetNumberOfIds(); ++k) {
        j = ids->GetId(k);
        _CentersA->GetPoint(j, c2);
        _WeightsA->GetTuple(j, n2);
        w = EvaluateKernel(c1, c2);
        EvaluateKernelGradient(g, c1, c2);
        for (int d = 0; d < 3; ++d) {
          kws   [d] += w    * n2[d];
          dks[0][d] += g[0] * n2[d];
          dks[1][d] += g[1] * n2[d];
          dks[2][d] += g[2] * n2[d];
        }
      }
      // Compute kwt = KtauT and dkt = gradKtauT.transpose()
      // (cf. Deformetrica 2.0 OrientedSurfaceMesh::ComputeMatchGradient)
      _LocatorB->FindPointsWithinRadius(_Radius, c1, ids);
      memset(kwt, 0, 3 * sizeof(double));
      memset(dkt, 0, 9 * sizeof(double));
      for (vtkIdType k = 0; k < ids->GetNumberOfIds(); ++k) {
        j = ids->GetId(k);
        _CentersB->GetPoint(j, c2);
        _WeightsB->GetTuple(j, n2);
        w = EvaluateKernel(c1, c2);
        EvaluateKernelGradient(g, c1, c2);
        for (int d = 0; d < 3; ++d) {
          kwt   [d] += w    * n2[d];
          dkt[0][d] += g[0] * n2[d];
          dkt[1][d] += g[1] * n2[d];
          dkt[2][d] += g[2] * n2[d];
        }
      }
      // Add gradient
      for (int d = 0; d < 3; ++d) {
        e1[d] =  v3[d] -  v2[d];
        e2[d] =  v1[d] -  v3[d];
        e3[d] =  v2[d] -  v1[d];
        kw[d] = kws[d] - kwt[d];
        g [d] = ((dks[d][0] - dkt[d][0]) * n1[0] +
                 (dks[d][1] - dkt[d][1]) * n1[1] +
                 (dks[d][2] - dkt[d][2]) * n1[2]) * _2over3;
      }
      Add(g1, g), Add(g2, g), Add(g3, g);
      vtkMath::Cross(e1, kw, g), Add(g1, g);
      vtkMath::Cross(e2, kw, g), Add(g2, g);
      vtkMath::Cross(e3, kw, g), Add(g3, g);
    }
  }
};

// =============================================================================
// Initialization
// =============================================================================

// -----------------------------------------------------------------------------
void CurrentsDistance::Init()
{
  // Derive kernel sigma value, assuming points distributed in a spherical arrangement
  if (_Sigma == .0) {
    cerr << "CurrentsDistance::Initialize: Kernel width may not be zero!" << endl;
    exit(1);
  }
  if (_Sigma < .0) {
    double ba[6], bb[6];
    _Target->InputPointSet()->GetBounds(ba);
    _Source->InputPointSet()->GetBounds(bb);
    double va = (ba[1] - ba[0]) * (ba[3] - ba[2]) * (ba[5] - ba[4]);
    double vb = (bb[1] - bb[0]) * (bb[3] - bb[2]) * (bb[5] - bb[4]);
    double ra = pow(va * 3.0 / 4.0 / pi, (1.0/3.0));
    double rb = pow(vb * 3.0 / 4.0 / pi, (1.0/3.0));
    _Sigma = 0.5 * (ra + rb) * abs(_Sigma);
  }

  // Get currents representation of input data sets
  _TargetCurrent = ToCurrent(_Target->InputPointSet());
  _SourceCurrent = ToCurrent(_Source->InputPointSet());

  // Compute squared norm of fixed current(s)
  _TargetNormSquared = .0;
  CurrentsDistanceDotProduct dot_product(_Sigma);
  if (!_Target->Transformation()) {
    _TargetNormSquared += dot_product.Evaluate(_TargetCurrent, _TargetCurrent);
  }
  if (!_Source->Transformation()) {
    _TargetNormSquared += dot_product.Evaluate(_SourceCurrent, _SourceCurrent);
  }
}

// -----------------------------------------------------------------------------
void CurrentsDistance::Initialize()
{
  // Initialize base class
  PointSetDistance::Initialize();
  // Initialize this class
  CurrentsDistance::Init();
}

// -----------------------------------------------------------------------------
void CurrentsDistance::Reinitialize()
{
  // Reinitialize base class
  PointSetDistance::Reinitialize();
  // Reinitialize this class
  CurrentsDistance::Init();
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool CurrentsDistance::SetWithPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Currents kernel width") == 0) {
    return FromString(value, _Sigma) && _Sigma != .0;
  }
  if (strcmp(param, "Symmetric currents distance") == 0) {
    return FromString(value, _Symmetric);
  }
  return PointSetDistance::SetWithPrefix(param, value);
}

// -----------------------------------------------------------------------------
bool CurrentsDistance::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Kernel width") == 0) {
    return FromString(value, _Sigma) && _Sigma != .0;
  }
  if (strcmp(param, "Symmetric") == 0) {
    return FromString(value, _Symmetric);
  }
  return PointSetDistance::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList CurrentsDistance::Parameter() const
{
  ParameterList params = PointSetDistance::Parameter();
  InsertWithPrefix(params, "Kernel width", _Sigma);
  InsertWithPrefix(params, "Symmetric",    _Symmetric);
  return params;
}

// =============================================================================
// Evaluation
// =============================================================================

// -----------------------------------------------------------------------------
void CurrentsDistance::Update(bool)
{
  if (_Target->Transformation()) {
    _Target->Update();
    _TargetCurrent = ToCurrent(_Target->PointSet());
  }
  if (_Source->Transformation()) {
    _Source->Update();
    _SourceCurrent = ToCurrent(_Source->PointSet());
  }
}

// -----------------------------------------------------------------------------
double CurrentsDistance::Evaluate()
{
  MIRTK_START_TIMING();
  double d = _TargetNormSquared;
  CurrentsDistanceDotProduct dot_product(_Sigma);
  if (_Target->Transformation()) {
    d += dot_product.Evaluate(_TargetCurrent, _TargetCurrent);
  }
  if (_Source->Transformation()) {
    d += dot_product.Evaluate(_SourceCurrent, _SourceCurrent);
  }
  if (_Symmetric) {
    d -= dot_product.Evaluate(_TargetCurrent, _SourceCurrent);
    d -= dot_product.Evaluate(_SourceCurrent, _TargetCurrent);
  } else {
    d -= 2.0 * dot_product.Evaluate(_TargetCurrent, _SourceCurrent);
  }
  MIRTK_DEBUG_TIMING(2, "evaluation of currents distance");
  return d / ((_TargetCurrent->GetNumberOfPoints() + _SourceCurrent->GetNumberOfPoints()) / 2);
}

// -----------------------------------------------------------------------------
void CurrentsDistance::NonParametricGradient(const RegisteredPointSet *target,
                                             GradientType             *gradient)
{
  vtkPointSet *sa = _Target->PointSet();
  vtkPolyData *ca = _TargetCurrent;
  vtkPointSet *sb = _Source->PointSet();
  vtkPolyData *cb = _SourceCurrent;
  if (target == _Source) swap(sa, sb), swap(ca, cb);
  CurrentsDistanceGradient d(_Sigma);
  d.EvaluateGradient(sa, ca, sb, cb, gradient);
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
static inline void NegateTuples1(vtkDataArray *array)
{
  for (vtkIdType i = 0; i < array->GetNumberOfTuples(); ++i) {
    array->SetTuple1(i, -array->GetTuple1(i));
  }
}

// -----------------------------------------------------------------------------
void CurrentsDistance::WriteDataSets(const char *p, const char *suffix, bool all) const
{
  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  if (_Target->Transformation() || all) {
    vtkSmartPointer<vtkFloatArray> dist;
    if (_Target->Transformation()) {
      CurrentsDistanceDotProduct dot_product(_Sigma);
      dist = vtkSmartPointer<vtkFloatArray>::New();
      dist->SetName("distance");
      dist->SetNumberOfComponents(1);
      dist->SetNumberOfTuples(_TargetCurrent->GetNumberOfPoints());
      dist->FillComponent(0, .0);
      dot_product.Evaluate(_TargetCurrent, _SourceCurrent, dist);
      NegateTuples1(dist);
      dot_product.Evaluate(_TargetCurrent, _TargetCurrent, dist);
    }
    snprintf(fname, sz, "%starget%s%s", prefix, suffix, _Target->DefaultExtension());
    _Target->Write(fname, NULL, dist);
  }

  if (_Source->Transformation() || all) {
    vtkSmartPointer<vtkFloatArray> dist;
    if (_Source->Transformation()) {
      CurrentsDistanceDotProduct dot_product(_Sigma);
      dist = vtkSmartPointer<vtkFloatArray>::New();
      dist->SetName("distance");
      dist->SetNumberOfComponents(1);
      dist->SetNumberOfTuples(_SourceCurrent->GetNumberOfPoints());
      dist->FillComponent(0, .0);
      dot_product.Evaluate(_SourceCurrent, _TargetCurrent, dist);
      NegateTuples1(dist);
      dot_product.Evaluate(_SourceCurrent, _SourceCurrent, dist);
    }
    snprintf(fname, sz, "%ssource%s%s", prefix, suffix, _Source->DefaultExtension());
    _Source->Write(fname, NULL, dist);
  }
}


} // namespace mirtk
