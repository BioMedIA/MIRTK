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

#include "mirtk/SurfaceCurvature.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"

#include "mirtk/Matrix3x3.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/MeshSmoothing.h"
#include "mirtk/PointSetUtils.h"

#include "mirtk/Vtk.h"
#include "mirtk/VtkMath.h"

#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyDataNormals.h"
#include "vtkCurvatures.h"
#include "vtkIdList.h"

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"


namespace mirtk {


// =============================================================================
// Names of data arrays
// =============================================================================

MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::MINIMUM            = "Minimum_Curvature";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::MAXIMUM            = "Maximum_Curvature";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::PRINCIPAL          = "Principal_Curvatures";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::MEAN               = "Mean_Curvature";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::GAUSS              = "Gauss_Curvature";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::CURVEDNESS         = "Curvedness";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::MINIMUM_DIRECTION  = "Minimum_Curvature_Direction";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::MAXIMUM_DIRECTION  = "Maximum_Curvature_Direction";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::TENSOR             = "Curvature_Tensor";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::INVERSE_TENSOR     = "Inverse_Curvature_Tensor";
MIRTK_PointSet_EXPORT const char * const SurfaceCurvature::NORMALS            = "Normals";

// =============================================================================
// Auxiliary functors
// =============================================================================

namespace SurfaceCurvatureUtils {


// -----------------------------------------------------------------------------
/// Compute curvature tensors for each undirected edge
class ComputeEdgeTensors
{
  vtkPolyData     *_Surface;
  vtkDataArray    *_Normals;
  const EdgeTable *_EdgeTable;
  double           _DistNorm;
  vtkDataArray    *_Tensors;

public:

  /// Compute edge tensors for specified range of edges
  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    vtkIdType npts, *pts, ptId1, ptId2, cellId1, cellId2;
    double    d, p1[3], p2[3], e[3], n1[3], n2[3], dp, cp[3], beta, s, T[6];

    vtkSmartPointer<vtkIdList> adjCellIds = vtkSmartPointer<vtkIdList>::New();

    // Iterate over faces
    for (cellId1 = re.begin(); cellId1 != re.end(); ++cellId1) {
      _Surface->GetCellPoints(cellId1, npts, pts);
      for (vtkIdType i = 0; i < npts; ++i) {
        // Get edge points
        ptId1 = pts[i];
        ptId2 = pts[(i + 1) % npts];
        // Process edge only if there is really just one other neighbour, i.e.,
        // the edge has two adjacent faces with IDs cellId1 and cellId2, and
        // when the edge was not visited before (ensured by cellId1 < cellId2)
        _Surface->GetCellEdgeNeighbors(cellId1, ptId1, ptId2, adjCellIds);
        if (adjCellIds->GetNumberOfIds() == 1 && cellId1 < (cellId2 = adjCellIds->GetId(0))) {
          // Compute normalized edge vector
          _Surface->GetPoint(ptId1, p1);
          _Surface->GetPoint(ptId2, p2);
          e[0] = p2[0] - p1[0];
          e[1] = p2[1] - p1[1];
          e[2] = p2[2] - p1[2];
          d = vtkMath::Normalize(e);
          // Compute signed angle of face normals
          _Normals->GetTuple(cellId1, n1);
          _Normals->GetTuple(cellId2, n2);
          dp = max(-1.0, min(vtkMath::Dot(n1, n2), 1.0));
          vtkMath::Cross(n1, n2, cp);
          beta = sgn(vtkMath::Dot(cp, e)) * acos(dp);
          // Compute upper triangular part of symmetric curvature tensor
          s  = beta * d;
          s /= _DistNorm; // Avoid too large numerics
          T[0] = s * e[0] * e[0]; // ParaView compatible order:
          T[1] = s * e[1] * e[1]; // XX, YY, ZZ, XY, YZ, XZ
          T[2] = s * e[2] * e[2];
          T[3] = s * e[0] * e[1];
          T[4] = s * e[1] * e[2];
          T[5] = s * e[0] * e[2];
          _Tensors->SetTuple(_EdgeTable->EdgeId(ptId1, ptId2), T);
        }
      }
    }
  }

  /// Compute curvature tensor for each edge of surface mesh
  static vtkSmartPointer<vtkDataArray> Run(vtkPolyData *surface, const EdgeTable &edgeTable)
  {
    vtkSmartPointer<vtkDataArray> tensors = vtkSmartPointer<vtkDoubleArray>::New();
    tensors->SetName(SurfaceCurvature::TENSOR);
    tensors->SetNumberOfComponents(6);
    tensors->SetNumberOfTuples(edgeTable.NumberOfEdges());
    ComputeEdgeTensors body;
    body._EdgeTable = &edgeTable;
    body._Surface   = surface;
    body._Normals   = surface->GetCellData()->GetNormals();
    body._DistNorm  = AverageEdgeLength(surface->GetPoints(), edgeTable);
    body._Tensors   = tensors;
    parallel_for(blocked_range<vtkIdType>(0, surface->GetNumberOfCells()), body);
    return tensors;
  }
};

// -----------------------------------------------------------------------------
/// Average curvature tensors of adjacent edges to obtain curvature tensor at nodes
struct AverageEdgeTensors
{
  const EdgeTable *_EdgeTable;
  vtkDataArray    *_EdgeTensors;
  vtkDataArray    *_NodeTensors;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double    Te[6], Tp[6];
    int       npts, edgeId, i, c;
    const int *pts;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      _EdgeTable->GetAdjacentPoints(ptId, npts, pts);
      memset(Tp, 0, 6 * sizeof(double));
      for (i = 0; i < npts; ++i) {
        edgeId = _EdgeTable->EdgeId(ptId, pts[i]);
        _EdgeTensors->GetTuple(edgeId, Te);
        for (c = 0; c < 6; ++c) Tp[c] += Te[c];
      }
      for (c = 0; c < 6; ++c) Tp[c] /= npts;
      _NodeTensors->SetTuple(ptId, Tp);
    }
  }

  static void Run(vtkPolyData *surface, const EdgeTable &edgeTable, vtkDataArray *edgeTensors)
  {
    vtkSmartPointer<vtkDataArray> tensors = vtkSmartPointer<vtkDoubleArray>::New();
    tensors->SetName(SurfaceCurvature::TENSOR);
    tensors->SetNumberOfComponents(6);
    tensors->SetNumberOfTuples(surface->GetNumberOfPoints());
    AverageEdgeTensors body;
    body._EdgeTable   = &edgeTable;
    body._EdgeTensors = edgeTensors;
    body._NodeTensors = tensors;
    parallel_for(blocked_range<vtkIdType>(0, surface->GetNumberOfPoints()), body);
    surface->GetPointData()->AddArray(tensors);
  }
};

// -----------------------------------------------------------------------------
/// Perform eigenanalysis of curvature tensors
struct DecomposeTensors
{
  typedef Eigen::Matrix3d                            EigenMatrix;
  typedef Eigen::Vector3d                            EigenVector;
  typedef Eigen::SelfAdjointEigenSolver<EigenMatrix> EigenSolver;
  typedef EigenSolver::RealVectorType                EigenValues;
  typedef EigenSolver::MatrixType                    EigenVectors;

  vtkDataArray *_InputTensors;
  vtkDataArray *_OutputTensors;
  vtkDataArray *_InverseTensors;
  vtkDataArray *_Minimum;
  vtkDataArray *_Maximum;
  vtkDataArray *_Normals;          // Input normals
  vtkDataArray *_NormalDirection;  // Output normals
  vtkDataArray *_MinimumDirection;
  vtkDataArray *_MaximumDirection;

  static void ComposeTensor(double T[6], EigenValues lambda, EigenVectors e, int i, int j, int k)
  {
    T[SurfaceCurvature::XX] = lambda(i) * e(i, i) * e(i, i)
                            + lambda(j) * e(i, j) * e(i, j)
                            + lambda(k) * e(i, k) * e(i, k);
    T[SurfaceCurvature::YY] = lambda(i) * e(j, i) * e(j, i)
                            + lambda(j) * e(j, j) * e(j, j)
                            + lambda(k) * e(j, k) * e(j, k);
    T[SurfaceCurvature::ZZ] = lambda(i) * e(k, i) * e(k, i)
                            + lambda(j) * e(k, j) * e(k, j)
                            + lambda(k) * e(k, k) * e(k, k);
    T[SurfaceCurvature::XY] = lambda(i) * e(i, i) * e(j, i)
                            + lambda(j) * e(i, j) * e(j, j)
                            + lambda(k) * e(i, k) * e(j, k);
    T[SurfaceCurvature::YZ] = lambda(i) * e(j, i) * e(k, i)
                            + lambda(j) * e(j, j) * e(k, j)
                            + lambda(k) * e(j, k) * e(k, k);
    T[SurfaceCurvature::XZ] = lambda(i) * e(i, i) * e(k, i)
                            + lambda(j) * e(i, j) * e(k, j)
                            + lambda(k) * e(i, k) * e(k, k);
  }

  static void ComposeTensor(double T[6], EigenValues lambda, double a[3], double b[3], double c[3], int i, int j, int k)
  {
    T[SurfaceCurvature::XX] = lambda(i) * a[0] * a[0]
                            + lambda(j) * b[0] * b[0]
                            + lambda(k) * c[0] * c[0];
    T[SurfaceCurvature::YY] = lambda(i) * a[1] * a[1]
                            + lambda(j) * b[1] * b[1]
                            + lambda(k) * c[1] * c[1];
    T[SurfaceCurvature::ZZ] = lambda(i) * a[2] * a[2]
                            + lambda(j) * b[2] * b[2]
                            + lambda(k) * c[2] * c[2];
    T[SurfaceCurvature::XY] = lambda(i) * a[0] * a[1]
                            + lambda(j) * b[0] * b[1]
                            + lambda(k) * c[0] * c[1];
    T[SurfaceCurvature::YZ] = lambda(i) * a[1] * a[2]
                            + lambda(j) * b[1] * b[2]
                            + lambda(k) * c[1] * c[2];
    T[SurfaceCurvature::XZ] = lambda(i) * a[0] * a[2]
                            + lambda(j) * b[0] * b[2]
                            + lambda(k) * c[0] * c[2];
  }

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    int    perm[3], i, j, k;
    double T[6], n[3], e1[3], e2[3], dp[3], v[3], sign;

    EigenVector normal;
    EigenMatrix tensor, inverse;
    EigenSolver eigensolver;

    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      // Get 3x3 symmetric real curvature tensor
      _InputTensors->GetTuple(ptId, T);
      tensor << T[SurfaceCurvature::XX],
                T[SurfaceCurvature::XY],
                T[SurfaceCurvature::XZ],
                T[SurfaceCurvature::YX],
                T[SurfaceCurvature::YY],
                T[SurfaceCurvature::YZ],
                T[SurfaceCurvature::ZX],
                T[SurfaceCurvature::ZY],
                T[SurfaceCurvature::ZZ];
      // Decompose symmetric 3x3 curvature tensor
      eigensolver.compute(tensor);
      const EigenValues  &lambda = eigensolver.eigenvalues();
      const EigenVectors &V      = eigensolver.eigenvectors();
      // Find permutation of eigenvalues so they are sorted by increasing magnitude
      perm[0] = 0, perm[1] = 1, perm[2] = 2;
      for (i = 0; i < 3; ++i) {
        for (j = i, k = i + 1; k < 3; ++k) {
          if (fabs(lambda[perm[k]]) < fabs(lambda[perm[j]])) j = k;
        }
        swap(perm[i], perm[j]);
      }
      i = perm[0]; // normal
      j = perm[1]; // k_min
      k = perm[2]; // k_max
      // If input normals provided...
      if (_Normals) {
        // Ensure that normal direction is the one which is most similar
        // to the input surface normal because in flat areas the magnitude
        // of all eigenvalues is very similar and can lead to errors
        _Normals->GetTuple(ptId, n);
        normal << n[0], n[1], n[2];
        dp[i] = fabs(normal.dot(V.col(i)));
        dp[j] = fabs(normal.dot(V.col(j)));
        dp[k] = fabs(normal.dot(V.col(k)));
        if (dp[i] > dp[j]) {
          if (dp[i] < dp[k]) swap(i, k);
        } else {
          if (dp[j] > dp[k]) swap(i, j);
          else               swap(i, k);
        }
        // Preserve input surface normal and find new curvature directions
        // which are orthogonal to this normal vector instead, keeping the
        // found direction of maximum curvature which is most reliable.
        //
        // If output normals based on the curvature tensor where requested,
        // replace the input normal, however, and ensure that inward/outward
        // relationship is unchanged.
        if (_NormalDirection) {
          sign = sgn(n[0] * V(0, i) + n[1] * V(1, i) + n[2] * V(2, i));
          n [0] = V(0, i), n [1] = V(1, i), n [2] = V(2, i);
          e1[0] = V(0, j), e1[1] = V(1, j), e1[2] = V(2, j);
          e2[0] = V(0, k), e2[1] = V(1, k), e2[2] = V(2, k);
          vtkMath::MultiplyScalar(n, sign);
          vtkMath::Cross(n, e1, v);
          if (vtkMath::Dot(v, e2) < .0) vtkMath::MultiplyScalar(e2, -1.0);
        } else {
          e2[0] = V(0, k), e2[1] = V(1, k), e2[2] = V(2, k);
          vtkMath::Cross(e2, n, e1);
        }
      } else {
        e1[0] = V(0, j), e1[1] = V(1, j), e1[2] = V(2, j);
        e2[0] = V(0, k), e2[1] = V(1, k), e2[2] = V(2, k);
        vtkMath::Cross(n, e1, v);
        if (vtkMath::Dot(v, e2) < .0) vtkMath::MultiplyScalar(e2, -1.0);
      }
      // Set minimum and maximum principle curvature
      if (_Minimum) _Minimum->SetComponent(ptId, 0, lambda[j]);
      if (_Maximum) _Maximum->SetComponent(ptId, 0, lambda[k]);
      // Set direction vectors
      if (_NormalDirection ) _NormalDirection ->SetTuple(ptId, n);
      if (_MinimumDirection) _MinimumDirection->SetTuple(ptId, e1);
      if (_MaximumDirection) _MaximumDirection->SetTuple(ptId, e2);
      // Compose tensor from reordered and reoriented eigenvalues/-vectors
      if (_OutputTensors || _InverseTensors) {
        ComposeTensor(T, lambda, n, e1, e2, i, j, k);
      }
      // Set (reordered and reoriented) curvature tensor
      if (_OutputTensors) _OutputTensors->SetTuple(ptId, T);
      // Compute inverse of (reordered and reoriented) curvature tensor
      if (_InverseTensors) {
        tensor << T[SurfaceCurvature::XX],
                  T[SurfaceCurvature::XY],
                  T[SurfaceCurvature::XZ],
                  T[SurfaceCurvature::YX],
                  T[SurfaceCurvature::YY],
                  T[SurfaceCurvature::YZ],
                  T[SurfaceCurvature::ZX],
                  T[SurfaceCurvature::ZY],
                  T[SurfaceCurvature::ZZ];
        inverse = tensor.inverse();
        T[SurfaceCurvature::XX] = inverse(0, 0);
        T[SurfaceCurvature::YY] = inverse(1, 1);
        T[SurfaceCurvature::ZZ] = inverse(2, 2);
        T[SurfaceCurvature::XY] = inverse(0, 1);
        T[SurfaceCurvature::YZ] = inverse(1, 2);
        T[SurfaceCurvature::XZ] = inverse(0, 2);
        _InverseTensors->SetTuple(ptId, T);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute minimum and maximum curvature from mean and Gauss curvatures
struct CalculateMinMaxCurvature
{
  vtkDataArray *_Mean;
  vtkDataArray *_Gauss;
  vtkDataArray *_Minimum;
  vtkDataArray *_Maximum;
  int           _MinimumComponentIndex;
  int           _MaximumComponentIndex;

  CalculateMinMaxCurvature(vtkDataArray *min, vtkDataArray *max)
  :
    _Minimum(min),
    _Maximum(max),
    _MinimumComponentIndex(0),
    _MaximumComponentIndex(0)
  {}

  CalculateMinMaxCurvature(vtkDataArray *min, int j1, vtkDataArray *max, int j2)
  :
    _Minimum(min),
    _Maximum(max),
    _MinimumComponentIndex(j1),
    _MaximumComponentIndex(j2)
  {}

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double h, k, h2_k;
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      h    = _Mean ->GetComponent(ptId, 0);
      k    = _Gauss->GetComponent(ptId, 0);
      h2_k = h * h - k;
      h2_k = (h2_k > .0 ? sqrt(h2_k) : .0);
      if (_Minimum) {
        _Minimum->SetComponent(ptId, _MinimumComponentIndex, h - h2_k);
      }
      if (_Maximum) {
        _Maximum->SetComponent(ptId, _MaximumComponentIndex, h + h2_k);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute mean curvature from minium and maximum eigenvalues
struct CalculateMeanCurvature
{
  vtkDataArray *_Minimum;
  vtkDataArray *_Maximum;
  vtkDataArray *_Mean;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double k_min, k_max;
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      k_min = _Minimum->GetComponent(ptId, 0);
      k_max = _Maximum->GetComponent(ptId, 0);
      _Mean->SetComponent(ptId, 0, .5 * (k_min + k_max));
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute Gauss curvature from minium and maximum curvature
struct CalculateGaussCurvature
{
  vtkDataArray *_Minimum;
  vtkDataArray *_Maximum;
  vtkDataArray *_Gauss;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double k_min, k_max;
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      k_min = _Minimum->GetComponent(ptId, 0);
      k_max = _Maximum->GetComponent(ptId, 0);
      _Gauss->SetComponent(ptId, 0, k_min * k_max);
    }
  }
};

// -----------------------------------------------------------------------------
/// Compute curvendess from minium and maximum curvature
struct CalculateCurvedness
{
  vtkDataArray *_Minimum;
  vtkDataArray *_Maximum;
  vtkDataArray *_Curvedness;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    double k_min, k_max;
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      k_min = _Minimum->GetComponent(ptId, 0);
      k_max = _Maximum->GetComponent(ptId, 0);
      _Curvedness->SetComponent(ptId, 0, sqrt(.5 * (k_min * k_min + k_max * k_max)));
    }
  }
};

// -----------------------------------------------------------------------------
/// Multiply tuples by scalar
struct MultiplyScalar
{
  vtkDataArray *_Input;
  vtkDataArray *_Output;
  double        _Scalar;

  void operator ()(const blocked_range<vtkIdType> &re) const
  {
    for (vtkIdType ptId = re.begin(); ptId != re.end(); ++ptId) {
      for (int c = 0; c < _Input->GetNumberOfComponents(); ++c) {
        _Output->SetComponent(ptId, c, _Scalar * _Input->GetComponent(ptId, c));
      }
    }
  }

  static void Run(vtkDataArray *output, vtkDataArray *input, double s)
  {
    output->SetNumberOfComponents(input->GetNumberOfComponents());
    output->SetNumberOfTuples(input->GetNumberOfTuples());
    MultiplyScalar body;
    body._Input  = input;
    body._Output = output;
    body._Scalar = s;
    parallel_for(blocked_range<vtkIdType>(0, input->GetNumberOfTuples()), body);
  }
};


} // namespace SurfaceCurvatureUtils
using namespace SurfaceCurvatureUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceCurvature::CopyAttributes(const SurfaceCurvature &other)
{
  _CurvatureType   = other._CurvatureType;
  _VtkCurvatures   = other._VtkCurvatures;
  _TensorAveraging = other._TensorAveraging;
  _Normalize       = other._Normalize;
  _Volume          = other._Volume;
  _Radius          = other._Radius;
}

// -----------------------------------------------------------------------------
SurfaceCurvature::SurfaceCurvature(int type)
:
  _CurvatureType(type),
  _VtkCurvatures(false),
  _TensorAveraging(3),
  _Normalize(false),
  _Volume(.0),
  _Radius(.0)
{
}

// -----------------------------------------------------------------------------
SurfaceCurvature::SurfaceCurvature(const SurfaceCurvature &other)
:
  SurfaceFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
SurfaceCurvature &SurfaceCurvature::operator =(const SurfaceCurvature &other)
{
  if (this != &other) {
    SurfaceFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
SurfaceCurvature::~SurfaceCurvature()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void SurfaceCurvature::Initialize()
{
  // Initialize base class
  SurfaceFilter::Initialize();

  // Check if vtkCurvatures can be used
  if (_CurvatureType > Curvedness) {
    _VtkCurvatures = false;
  }

  // Compute face normals
  if (!_VtkCurvatures && _Output->GetCellData()->GetNormals() == nullptr) {
    vtkNew<vtkPolyDataNormals> filter;
    SetVTKInput(filter, _Output);
    filter->ComputePointNormalsOn();
    filter->ComputeCellNormalsOn();
    filter->SplittingOff();
    filter->AutoOrientNormalsOff();
    filter->ConsistencyOn();
    filter->Update();
    _Output = filter->GetOutput();
    _Output->BuildLinks();
  }

  // Compute volume and approximate radius of convex hull
  if (_Normalize && (_CurvatureType & (Mean | Gauss)) != 0) {
    _Volume = mirtk::Volume(ConvexHull(_Input));
    _Radius = pow( 3 * _Volume / (4.0 * pi), 1.0/3.0);
  }

  // Remove requested outputs from input to force recomputation
  vtkPointData *pd = _Output->GetPointData();
  if ((_CurvatureType & Minimum         ) != 0) pd->RemoveArray(MINIMUM);
  if ((_CurvatureType & Maximum         ) != 0) pd->RemoveArray(MAXIMUM);
  if ((_CurvatureType & Mean            ) != 0) pd->RemoveArray(MEAN);
  if ((_CurvatureType & Gauss           ) != 0) pd->RemoveArray(GAUSS);
  if ((_CurvatureType & MinimumDirection) != 0) pd->RemoveArray(MINIMUM_DIRECTION);
  if ((_CurvatureType & MaximumDirection) != 0) pd->RemoveArray(MAXIMUM_DIRECTION);
  if ((_CurvatureType & Tensor          ) != 0) pd->RemoveArray(TENSOR);
  if ((_CurvatureType & InverseTensor   ) != 0) pd->RemoveArray(TENSOR);
}

// -----------------------------------------------------------------------------
vtkDataArray *SurfaceCurvature::GetMinimumCurvature()
{
  ComputeMinMaxCurvature(true, false);
  return _Output->GetPointData()->GetArray(MINIMUM);
}

// -----------------------------------------------------------------------------
vtkDataArray *SurfaceCurvature::GetMaximumCurvature()
{
  ComputeMinMaxCurvature(false, true);
  return _Output->GetPointData()->GetArray(MAXIMUM);
}

// -----------------------------------------------------------------------------
vtkDataArray *SurfaceCurvature::GetPrincipalCurvatures()
{
  ComputePrincipalCurvatures();
  return _Output->GetPointData()->GetArray(PRINCIPAL);
}

// -----------------------------------------------------------------------------
vtkDataArray *SurfaceCurvature::GetMeanCurvature()
{
  ComputeMeanCurvature();
  return _Output->GetPointData()->GetArray(MEAN);
}

// -----------------------------------------------------------------------------
vtkDataArray *SurfaceCurvature::GetGaussCurvature()
{
  ComputeGaussCurvature();
  return _Output->GetPointData()->GetArray(GAUSS);
}

// -----------------------------------------------------------------------------
vtkDataArray *SurfaceCurvature::GetCurvedness()
{
  ComputeCurvedness();
  return _Output->GetPointData()->GetArray(CURVEDNESS);
}

// -----------------------------------------------------------------------------
void SurfaceCurvature::Execute()
{
  // Compute curvature tensor field and its eigenvalues if needed
  if (!_VtkCurvatures) {
    ComputeTensorField();
    if ((_CurvatureType - (_CurvatureType & Tensor)) != 0) {
      DecomposeTensorField();
    }
  }

  // Compute mean and/or Gauss curvature either using vtkCurvatures
  // or from minimum and maximum curvature values
  if ((_CurvatureType & Mean ) != 0) ComputeMeanCurvature();
  if ((_CurvatureType & Gauss) != 0) ComputeGaussCurvature();

  // Compute minimum and/or maximum curvature from mean and Gauss curvature
  // if not obtained from curvature tensor before
  if (_VtkCurvatures) {
    ComputeMinMaxCurvature((_CurvatureType & Minimum) != 0, (_CurvatureType & Maximum) != 0);
  }

  // Compute curvedness from minimum and maximum curvature
  if ((_CurvatureType & Curvedness) != 0) ComputeCurvedness();
}

// -----------------------------------------------------------------------------
void SurfaceCurvature::ComputeTensorField()
{
  MIRTK_START_TIMING();

  // Construct edge table if none provided
  InitializeEdgeTable();

  // Compute curvature tensor for each edge
  vtkSmartPointer<vtkDataArray> tensors;
  tensors = ComputeEdgeTensors::Run(_Output, *_EdgeTable);

  // Locally integrate edge tensors
  AverageEdgeTensors::Run(_Output, *_EdgeTable, tensors);

  MeshSmoothing smoother;
  smoother.Input(_Output);
  smoother.EdgeTable(_EdgeTable);
  smoother.SmoothPointsOff();
  smoother.SmoothArray(TENSOR);
  smoother.Weighting(MeshSmoothing::Combinatorial);
  smoother.AdjacentValuesOnly(false);
  smoother.NumberOfIterations(_TensorAveraging);
  smoother.Sigma(.0);
  smoother.Lambda(1.0);
  smoother.Run();
  _Output = smoother.Output();

  MIRTK_DEBUG_TIMING(3, "computation of curvature tensors");
}

// -----------------------------------------------------------------------------
void SurfaceCurvature::DecomposeTensorField()
{
  MIRTK_START_TIMING();

  const vtkIdType n = _Output->GetNumberOfPoints();

  vtkSmartPointer<vtkDataArray> minimum, maximum;
  vtkSmartPointer<vtkDataArray> minimum_direction, maximum_direction;
  vtkSmartPointer<vtkDataArray> tensors, inverse_tensors;
  vtkSmartPointer<vtkDataArray> input_normals, output_normals;

  tensors = _Output->GetPointData()->GetArray(TENSOR);

  if ((_CurvatureType & (Minimum | Mean | Gauss | Curvedness)) != 0) {
    minimum = NewArray(MINIMUM, n, 1);
    _Output->GetPointData()->AddArray(minimum);
  }
  if ((_CurvatureType & (Maximum | Mean | Gauss | Curvedness)) != 0) {
    maximum = NewArray(MAXIMUM, n, 1);
    _Output->GetPointData()->AddArray(maximum);
  }
  if ((_CurvatureType & MinimumDirection) != 0) {
    minimum_direction = NewArray(MINIMUM_DIRECTION, n, 3);
    _Output->GetPointData()->AddArray(minimum_direction);
  }
  if ((_CurvatureType & MaximumDirection) != 0) {
    maximum_direction = NewArray(MAXIMUM_DIRECTION, n, 3);
    _Output->GetPointData()->AddArray(maximum_direction);
  }
  input_normals = _Output->GetPointData()->GetNormals();
  if ((_CurvatureType & Normal) != 0) output_normals = input_normals;
  if ((_CurvatureType & InverseTensor) != 0) {
    if ((_CurvatureType & Tensor) != 0) {
      inverse_tensors = NewArray(INVERSE_TENSOR, n, 6);
      _Output->GetPointData()->AddArray(inverse_tensors);
    } else {
      inverse_tensors = tensors;
      inverse_tensors->SetName(INVERSE_TENSOR);
    }
  }

  DecomposeTensors eval;
  eval._InputTensors     = tensors;
  eval._OutputTensors    = (_CurvatureType & Tensor) ? tensors : NULL;
  eval._InverseTensors   = inverse_tensors;
  eval._Minimum          = minimum;
  eval._Maximum          = maximum;
  eval._Normals          = input_normals;
  eval._NormalDirection  = output_normals;
  eval._MinimumDirection = minimum_direction;
  eval._MaximumDirection = maximum_direction;
  parallel_for(blocked_range<vtkIdType>(0, n), eval);

  if (_Normalize) {
    if (minimum) MultiplyScalar::Run(minimum, minimum, _Radius);
    if (maximum) MultiplyScalar::Run(maximum, maximum, _Radius);
  }

  MIRTK_DEBUG_TIMING(3, "decomposition of curvature tensors");
}

// -----------------------------------------------------------------------------
void SurfaceCurvature::ComputeMeanCurvature()
{
  const vtkIdType n = _Output->GetNumberOfPoints();

  vtkSmartPointer<vtkDataArray> output;
  output = _Output->GetPointData()->GetArray(MEAN);
  if (output) return;

  MIRTK_START_TIMING();

  vtkDataArray *minimum = _Output->GetPointData()->GetArray(MINIMUM);
  vtkDataArray *maximum = _Output->GetPointData()->GetArray(MAXIMUM);

  if (_VtkCurvatures || !minimum || !maximum) {
    vtkSmartPointer<vtkCurvatures> curvatures;
    curvatures = vtkSmartPointer<vtkCurvatures>::New();
    SetVTKInput(curvatures, _Output);
    curvatures->SetCurvatureTypeToMean();
    curvatures->Update();
    output = curvatures->GetOutput()->GetPointData()->GetArray("Mean_Curvature");
    output->SetName(MEAN);
    if (_Normalize) MultiplyScalar::Run(output, output, _Radius);
  } else {
    output = NewArray(MEAN, n, 1);
    CalculateMeanCurvature eval;
    eval._Minimum = minimum;
    eval._Maximum = maximum;
    eval._Mean    = output;
    parallel_for(blocked_range<vtkIdType>(0, n), eval);
  }

  _Output->GetPointData()->AddArray(output);

  MIRTK_DEBUG_TIMING(3, "computation of mean curvature");
}

// -----------------------------------------------------------------------------
void SurfaceCurvature::ComputeGaussCurvature()
{
  const vtkIdType n = _Output->GetNumberOfPoints();

  vtkSmartPointer<vtkDataArray> output;
  output = _Output->GetPointData()->GetArray(GAUSS);
  if (output) return;

  MIRTK_START_TIMING();

  vtkDataArray *minimum = _Output->GetPointData()->GetArray(MINIMUM);
  vtkDataArray *maximum = _Output->GetPointData()->GetArray(MAXIMUM);

  if (_VtkCurvatures || !minimum || !maximum) {
    vtkSmartPointer<vtkCurvatures> curvatures;
    curvatures = vtkSmartPointer<vtkCurvatures>::New();
    SetVTKInput(curvatures, _Output);
    curvatures->SetCurvatureTypeToGaussian();
    curvatures->Update();
    output = curvatures->GetOutput()->GetPointData()->GetArray("Gauss_Curvature");
    output->SetName(GAUSS);
    if (_Normalize) MultiplyScalar::Run(output, output, _Radius * _Radius);
  } else {
    output = NewArray(GAUSS, n, 1);
    CalculateGaussCurvature eval;
    eval._Minimum = minimum;
    eval._Maximum = maximum;
    eval._Gauss   = output;
    parallel_for(blocked_range<vtkIdType>(0, n), eval);
  }

  _Output->GetPointData()->AddArray(output);

  MIRTK_DEBUG_TIMING(3, "computation of Gauss curvature");
}

// -----------------------------------------------------------------------------
void SurfaceCurvature::ComputePrincipalCurvatures()
{
  const vtkIdType n = _Output->GetNumberOfPoints();

  // Skip if requested output arrays already exist
  vtkSmartPointer<vtkDataArray> principal;
  principal = _Output->GetPointData()->GetArray(PRINCIPAL);
  if (principal) return;

  vtkDataArray *minimum = _Output->GetPointData()->GetArray(MINIMUM);
  vtkDataArray *maximum = _Output->GetPointData()->GetArray(MAXIMUM);

  if (minimum && maximum) {

    // Allocate output arrays
    principal = NewArray(PRINCIPAL, n, 2);
    _Output->GetPointData()->AddArray(principal);

    // Copy principal curvatures
    for (vtkIdType ptId = 0; ptId < n; ++ptId) {
      principal->SetComponent(ptId, 0, minimum->GetComponent(ptId, 0));
      principal->SetComponent(ptId, 1, minimum->GetComponent(ptId, 1));
    }

  } else {

    // Get mean and Gauss curvature
    vtkDataArray *mean  = GetMeanCurvature();
    vtkDataArray *gauss = GetGaussCurvature();

    MIRTK_START_TIMING();

    // Allocate output arrays
    principal = NewArray(PRINCIPAL, n, 2);
    _Output->GetPointData()->AddArray(principal);

    // Compute principal curvatures
    CalculateMinMaxCurvature eval(principal, 0, principal, 1);
    eval._Mean  = mean;
    eval._Gauss = gauss;
    parallel_for(blocked_range<vtkIdType>(0, n), eval);

    MIRTK_DEBUG_TIMING(3, "computation of principal curvatures");
  }
}

// -----------------------------------------------------------------------------
void SurfaceCurvature::ComputeMinMaxCurvature(bool min, bool max)
{
  const vtkIdType n = _Output->GetNumberOfPoints();

  // Skip if requested output arrays already exist
  vtkSmartPointer<vtkDataArray> minimum, maximum;
  minimum = _Output->GetPointData()->GetArray(MINIMUM);
  maximum = _Output->GetPointData()->GetArray(MAXIMUM);
  if ((!min || minimum) && (!max || maximum)) return;

  // Get mean and Gauss curvature
  vtkDataArray *mean  = GetMeanCurvature();
  vtkDataArray *gauss = GetGaussCurvature();

  MIRTK_START_TIMING();

  // Allocate output arrays
  if (min) {
    minimum = NewArray(MINIMUM, n, 1);
    _Output->GetPointData()->AddArray(minimum);
  }
  if (max) {
    maximum = NewArray(MAXIMUM, n, 1);
    _Output->GetPointData()->AddArray(maximum);
  }

  // Compute minimum/maximum curvature
  CalculateMinMaxCurvature eval(minimum, maximum);
  eval._Mean  = mean;
  eval._Gauss = gauss;
  parallel_for(blocked_range<vtkIdType>(0, n), eval);

  MIRTK_DEBUG_TIMING(3, "computation of min/max curvature");
}

// -----------------------------------------------------------------------------
void SurfaceCurvature::ComputeCurvedness()
{
  const vtkIdType n = _Output->GetNumberOfPoints();

  vtkSmartPointer<vtkDataArray> output;
  output = _Output->GetPointData()->GetArray(CURVEDNESS);
  if (output) return;

  ComputeMinMaxCurvature(true, true);
  vtkDataArray *minimum = GetMinimumCurvature();
  vtkDataArray *maximum = GetMaximumCurvature();

  MIRTK_START_TIMING();

  output = NewArray(CURVEDNESS, n, 1);

  CalculateCurvedness eval;
  eval._Minimum    = minimum;
  eval._Maximum    = maximum;
  eval._Curvedness = output;
  parallel_for(blocked_range<vtkIdType>(0, n), eval);

  _Output->GetPointData()->AddArray(output);

  MIRTK_DEBUG_TIMING(3, "computation of curvedness");
}

// -----------------------------------------------------------------------------
void SurfaceCurvature::Finalize()
{
  // Remove not requested output which was needed to derive the requested output
  vtkPointData *pd = _Output->GetPointData();
  if ((_CurvatureType & Minimum         ) == 0) pd->RemoveArray(MINIMUM);
  if ((_CurvatureType & Maximum         ) == 0) pd->RemoveArray(MAXIMUM);
  if ((_CurvatureType & Mean            ) == 0) pd->RemoveArray(MEAN);
  if ((_CurvatureType & Gauss           ) == 0) pd->RemoveArray(GAUSS);
  if ((_CurvatureType & Curvedness      ) == 0) pd->RemoveArray(CURVEDNESS);
  if ((_CurvatureType & MinimumDirection) == 0) pd->RemoveArray(MINIMUM_DIRECTION);
  if ((_CurvatureType & MaximumDirection) == 0) pd->RemoveArray(MAXIMUM_DIRECTION);
  if ((_CurvatureType & Tensor          ) == 0) pd->RemoveArray(TENSOR);
  if ((_CurvatureType & InverseTensor   ) == 0) pd->RemoveArray(INVERSE_TENSOR);

  // Finalize base class
  SurfaceFilter::Finalize();
}


} // namespace mirtk
