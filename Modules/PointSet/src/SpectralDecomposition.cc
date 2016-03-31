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

#include "mirtk/SpectralDecomposition.h"

#include "mirtk/Assert.h"
#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Algorithm.h"
#include "mirtk/Memory.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/Matrix.h"
#include "mirtk/SparseMatrix.h"
#include "mirtk/Vector.h"
#include "mirtk/PointSet.h"
#include "mirtk/EdgeTable.h"

#include "vtkSmartPointer.h"
#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkOctreePointLocator.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"

#include <random>

#define USE_mlxHungarian 0


namespace mirtk { namespace SpectralDecomposition {


// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
inline int NumberOfPoints(vtkPolyData *data, const Array<int> *ids)
{
  return ((ids && !ids->empty()) ? static_cast<int>(ids->size()) : data->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline int PointIndex(vtkPolyData *data, const Array<int> *ids, int i)
{
  return (ids && !ids->empty()) ? (*ids)[i] : i;
}

// -----------------------------------------------------------------------------
/// Generate random permutation of indices [first, last)
Array<int> RandPerm(int first, int last, int n = -1, int seed = 0)
{
  Array<int> idx(last - first);
  for (int i = first; i < last; ++i) idx[i - first] = i;
  shuffle(idx.begin(), idx.end(), std::default_random_engine(seed));
  if (n >= 0) idx.resize(n);
  return idx;
}

// -----------------------------------------------------------------------------
/// Get points at specified indices
void GetPoints(PointSet &lp, const PointSet &p, const Array<int> &idx)
{
  const int n = static_cast<int>(idx.size());
  lp.Resize(n);
  for (int i = 0; i < n; ++i) {
    lp.SetPoint(i, p.GetPoint(idx[i]));
  }
}

// -----------------------------------------------------------------------------
/// Get points at specified indices
void GetPoints(PointSet &lp, vtkPoints *p, const Array<int> &idx)
{
  const int n = static_cast<int>(idx.size());
  lp.Resize(n);
  double pt[3];
  for (int i = 0; i < n; ++i) {
    p->GetPoint(idx[i], pt);
    lp.SetPoint(i,      pt);
  }
}

// -----------------------------------------------------------------------------
/// Get values of eigenmodes at specified nodes
void GetEigenmodes(Matrix &lm, const Matrix &m, const Array<int> &idx)
{
  const int n = static_cast<int>(idx.size());
  lm.Initialize(n, m.Cols());
  for (int r = 0; r < n; ++r)
  for (int c = 0; c < m.Cols(); ++c) {
    lm(r, c) = m(idx[r], c);
  }
}

// -----------------------------------------------------------------------------
/// Flig sign of specified eigenmodes
void FlipSign(Matrix &m, const Array<int> &idx)
{
  for (Array<int>::const_iterator c = idx.begin(); c != idx.end(); ++c) {
    for (int r = 0; r < m.Rows(); ++r) m(r, *c) = -m(r, *c);
  }
}

// -----------------------------------------------------------------------------
/// Subtract minimum and divide each dimension by the range of its point coordinates
class NormalizePoints
{
  PointSet &_Point;
  Point     _Min;
  Point     _Range;

  NormalizePoints(PointSet &p) : _Point(p) {}

public:

  void operator ()(const blocked_range<int> &re) const
  {
    for (int i = re.begin(); i != re.end(); ++i) {
      (_Point(i) -= _Min) /= _Range;
    }
  }

  static void Run(PointSet &p, Point *min = NULL, Point *max = NULL)
  {
    if (p.Size() == 0) return;
    NormalizePoints norm(p);
    if (!min || !max) {
      p.BoundingBox(norm._Min, norm._Range);
      norm._Range -= norm._Min;
    }
    if (min) norm._Min   = (*min);
    if (max) norm._Range = (*max) - norm._Min;
    blocked_range<int> pts(0, p.Size());
    parallel_for(pts, norm);
  }
};

// -----------------------------------------------------------------------------
/// Subtract minimum and divide each eigenmode by the range of its values
class RescaleEigenmodes
{
  Matrix &_Matrix;

  RescaleEigenmodes(Matrix &m) : _Matrix(m) {}

public:

  void operator ()(const blocked_range<int> &re) const
  {
    double cmin, cmax, crange;
    for (int c = re.begin(); c != re.end(); ++c) {
      _Matrix.ColRange(c, cmin, cmax);
      if (cmax > cmin) {
        crange = cmax - cmin;
        for (int r = 0; r < _Matrix.Rows(); ++r) {
          double &m = _Matrix(r, c);
          m = (m - cmin) / crange - .5;
        }
      } else {
        for (int r = 0; r < _Matrix.Rows(); ++r) {
          _Matrix(r, c) -= cmin;
        }
      }
    }
  }

  static void Run(Matrix &m)
  {
    if (m.NumberOfElements() == 0) return;
    RescaleEigenmodes norm(m);
    blocked_range<int> cols(0, m.Cols());
    parallel_for(cols, norm);
  }
};

// -----------------------------------------------------------------------------
/// Scale columns of matrix by specified weights
class ScaleCols
{
  Matrix       &_Matrix;
  const Vector &_Weight;

  ScaleCols(Matrix &m, const Vector &w) : _Matrix(m), _Weight(w) {}

public:

  void operator ()(const blocked_range<int> &re) const
  {
    for (int c = re.begin(); c != re.end(); ++c) {
      _Matrix.ScaleCol(c, _Weight(c));
    }
  }

  static void Run(Matrix &m, const Vector &w)
  {
    ScaleCols scale(m, w);
    blocked_range<int> cols(0, m.Cols());
    parallel_for(cols, scale);
  }
};

// =============================================================================
// Bipartite graph matching / Optimal assignment problem
// =============================================================================

// -----------------------------------------------------------------------------
/// Compute for each pair of eigenmodes and sign flip the assignment cost
class ComputeAssignmentCosts
{
  const Matrix                         &_E1; /// First eigenmodes
  const Matrix                         &_E2; /// Second eigenmodes
  const Array<Array<int> > &_J;  /// Map first row index to second row index
  Matrix                               &_C;  /// Assignment costs
  Array<int>                      _S;  /// Indices of pairs with sign flip

  ComputeAssignmentCosts(const Matrix  &m1, const Matrix &m2,
                         const Array<Array<int> > &r12, Matrix &c12)
  :
    _E1(m1), _E2(m2), _J(r12), _C(c12)
  {}

public:

  ComputeAssignmentCosts(const ComputeAssignmentCosts &lhs, split)
  :
    _E1(lhs._E1), _E2(lhs._E2), _J(lhs._J), _C(lhs._C)
  {}

  void join(const ComputeAssignmentCosts &rhs)
  {
    _S.insert(_S.end(), rhs._S.begin(), rhs._S.end());
  }

  void operator ()(const blocked_range2d<int> &re)
  {
    const int n = static_cast<int>(_J   .size());
    const int k = static_cast<int>(_J[0].size());
    double c1, c2, sum1, sum2;
    for (int c = re.rows().begin(); c != re.rows().end(); ++c) {
    for (int r = re.cols().begin(); r != re.cols().end(); ++r) {
      c1 = c2 = .0;
      for (int i = 0; i < n; ++i) {
        sum1 = sum2 = .0;
        for (int j = 0; j < k; ++j) {
          sum1 += pow(_E1(i, r) - _E2(_J[i][j], c), 2);
          sum2 += pow(_E1(i, r) + _E2(_J[i][j], c), 2);
        }
        c1 += sum1 / k;
        c2 += sum2 / k;
      }
      if (c1 <= c2) {
        _C(r, c) = sqrt(c1) / n;
      } else {
        _C(r, c) = sqrt(c2) / n;
        _S.push_back(_C.Index(r, c));
      }
    }
    }
  }

  static Matrix Run(const PointSet &p1, const Matrix &m1,
                    const PointSet &p2, const Matrix &m2,
                    const Array<Array<int> > &r12, Array<int> *flip = NULL)
  {
    Matrix w12(m1.Cols(), m2.Cols());
    ComputeAssignmentCosts eval(m1, m2, r12, w12);
    blocked_range2d<int> idx(0, m1.Cols(), 0, m2.Cols());
    parallel_reduce(idx, eval);
    if (flip) {
      sort(eval._S.begin(), eval._S.end());
      (*flip) = eval._S;
    }
    return w12;
  }
};

// -----------------------------------------------------------------------------
// FIXME: Implement hungarian algorithm in C++
Array<int> WeightedMatching(const Matrix &cost)
{
  Array<int> match(cost.Rows(), -1); // Maps row index to matching column index
  MIRTK_START_TIMING();
#if USE_mlxHungarian
  const int nlhs = 1, nrhs = 1; // M = hungarian(cost)
  mxArray *lhs[nlhs] = { NULL };
  mxArray *rhs[nrhs] = { cost.MxArray() };
  bool ok = mlxHungarian(nlhs, lhs, nrhs, rhs);
  for (int i = 0; i < nrhs; ++i) mxDestroyArray(rhs[i]);
  if (!ok) {
    cerr << "WeightedMatching: Failed to call mlxHungarian, ensure that libirtkmwutilsInitialize was called" << endl;
    exit(1);
  }
  mirtkAssert(mxIsDouble(lhs[0]), "mlxHungarian returns match in matrix of type DOUBLE");
  mirtkAssert(static_cast<int>(mxGetM(lhs[0])) == cost.Rows(), "same size as cost matrix");
  mirtkAssert(static_cast<int>(mxGetN(lhs[0])) == cost.Cols(), "same size as cost matrix");
  double *m = mxGetPr(lhs[0]);
  for (int c = 0; c < cost.Cols(); ++c) {
    for (int r = 0; r < cost.Rows(); ++r) {
      if (*m++ != .0) match[r] = c;
    }
  }
  for (int i = 0; i < nlhs; ++i) mxDestroyArray(lhs[i]);
#else
  cerr << "No C/C++ hungarian implementation selected during build configuration" << endl;
  exit(1);
#endif
  MIRTK_DEBUG_TIMING(10, "bipartite matching");
  return match;
}

// =============================================================================
// Find closest points
// =============================================================================

// -----------------------------------------------------------------------------
FindClosestPoints::FindClosestPoints(const PointSet          &p1,
                                     vtkAbstractPointLocator *p2,
                                     Array<int>        &idx)
:
  _P1(p1), _P2(p2), _Idx(idx)
{
}

// -----------------------------------------------------------------------------
FindClosestNPoints::FindClosestNPoints(const PointSet                 &p1,
                                       vtkAbstractPointLocator        *p2,
                                       Array<Array<int> > &idx,
                                       int                             num)
:
  _P1(p1), _P2(p2), _Idx(idx), _Num(num)
{
}

// -----------------------------------------------------------------------------
Array<int> FindClosestPoints::Run(const PointSet &p1, const PointSet &p2)
{
  Array<int> idx(p1.Size());
  if (!idx.empty()) {
    vtkSmartPointer<vtkPoints> points;
    points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(p2.Size());
    for (int i = 0; i < p2.Size(); ++i) {
      const Point &p = p2(i);
      points->SetPoint(i, p._x, p._y, p._z);
    }
    vtkSmartPointer<vtkPolyData> dataset;
    dataset = vtkSmartPointer<vtkPolyData>::New();
    dataset->SetPoints(points);
    vtkSmartPointer<vtkOctreePointLocator> locator;
    locator = vtkSmartPointer<vtkOctreePointLocator>::New();
    locator->SetDataSet(dataset);
    locator->BuildLocator();
    FindClosestPoints find_closest_point(p1, locator, idx);
    blocked_range<int> i(0, p1.Size());
    parallel_for(i, find_closest_point);
  }
  return idx;
}

// -----------------------------------------------------------------------------
Array<Array<int> > FindClosestNPoints::Run(const PointSet &p1, const PointSet &p2, int num)
{
  Array<Array<int> > idx(p1.Size(), Array<int>(num, -1));
  if (!idx.empty()) {
    vtkSmartPointer<vtkPoints> points;
    points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(p2.Size());
    for (int i = 0; i < p2.Size(); ++i) {
      const Point &p = p2(i);
      points->SetPoint(i, p._x, p._y, p._z);
    }
    vtkSmartPointer<vtkPolyData> dataset;
    dataset = vtkSmartPointer<vtkPolyData>::New();
    dataset->SetPoints(points);
    vtkSmartPointer<vtkOctreePointLocator> locator;
    locator = vtkSmartPointer<vtkOctreePointLocator>::New();
    locator->SetDataSet(dataset);
    locator->BuildLocator();
    FindClosestNPoints find_closest_points(p1, locator, idx, num);
    blocked_range<int> i(0, p1.Size());
    parallel_for(i, find_closest_points);
  }
  return idx;
}

// =============================================================================
// Node/edge feature weights
// =============================================================================

// -----------------------------------------------------------------------------
/// Function \rho used to enforce positive node weights
inline double PositiveNodeWeight(double w)
{
  return exp(w);
}

// -----------------------------------------------------------------------------
FeatureWeights EdgeWeights(vtkPolyData                    *dataset,
                           const Array<string> &names,
                           const Array<double>      *weights)
{
  mirtkAssert(weights == NULL || weights->size() == 1 || names.size() == weights->size(),
              "number of weights is either 0, 1, or equal number of feature names");

  vtkPointData * const data = dataset->GetPointData();

  // Get indices of feature arrays
  FeatureWeights edge_weights(names.size());
  for (size_t i = 0; i < names.size(); ++i) {
    if (data->GetArray(names[i].c_str(), edge_weights[i].first) == NULL) {
      cerr << "Missing node feature named: " << names[i] << endl;
      exit(1);
    }
    if      (weights == NULL)      edge_weights[i].second = .2;
    else if (weights->size() == 1) edge_weights[i].second = weights->at(0);
    else                           edge_weights[i].second = weights->at(i);
  }

  // Adjusts weights of edge features
  if (!edge_weights.empty()) {

    double range[2];
    FeatureWeights::iterator ew;
    for (ew = edge_weights.begin(); ew != edge_weights.end(); ++ew) {
      data->GetArray(ew->first)->GetRange(range, -1);
      ew->second /= (range[1] - range[0]);
    }
  }

  return edge_weights;
}

// -----------------------------------------------------------------------------
FeatureWeights NodeWeights(vtkPolyData                    *dataset,
                           const Array<string> &names,
                           const Array<double>      *weights)
{
  mirtkAssert(weights == NULL || weights->size() == 1 || names.size() == weights->size(),
              "number of weights is either 0, 1, or equal number of feature names");

  vtkPointData * const data = dataset->GetPointData();

  // Get indices of feature arrays
  FeatureWeights node_weights(names.size());
  for (size_t i = 0; i < names.size(); ++i) {
    if (data->GetArray(names[i].c_str(), node_weights[i].first) == NULL) {
      cerr << "Missing node feature named: " << names[i] << endl;
      exit(1);
    }
    if      (weights == NULL)      node_weights[i].second = 1.0;
    else if (weights->size() == 1) node_weights[i].second = weights->at(0);
    else                           node_weights[i].second = weights->at(i);
  }

  // Adjusts weights of node features
  if (!node_weights.empty()) {
    double range[2];
    FeatureWeights::iterator nw;
    for (nw = node_weights.begin(); nw != node_weights.end(); ++nw) {
      data->GetArray(nw->first)->GetRange(range, -1);
      nw->second /= (PositiveNodeWeight(range[1]) - PositiveNodeWeight(range[0]));
    }
  }

  return node_weights;
}

// =============================================================================
// Graph connectivity
// =============================================================================

// -----------------------------------------------------------------------------
void AdjacencyMatrix(SparseMatrix::Entries        adjw[],
                     SparseMatrix::StorageLayout  layout,
                     int r1, int c1, vtkPolyData *dataset,
                     FeatureWeights               weights)
{
  MIRTK_START_TIMING();
  if (layout == SparseMatrix::CCS) swap(r1, c1);

  // Get number of nodes and number of features per node
  const int n = dataset->GetNumberOfPoints();
  const int d = NumberOfFeatures(dataset, &weights);
  if (n == 0 || d == 0) return;

  // Get edge table
  EdgeTable edgeTable(dataset);
  if (edgeTable.NumberOfEdges() == 0) return;

  // Allocate memory for node feature vectors
  double *p1 = new double[d];
  double *p2 = new double[d];

  // Iterate over edges
  vtkIdType ptId1, ptId2;
  SparseMatrix::EntryType w;

  EdgeIterator it(edgeTable);
  for (it.InitTraversal(); it.GetNextEdge(ptId1, ptId2) != -1;) {
    GetFeatures(dataset, ptId1, p1, &weights);
    GetFeatures(dataset, ptId2, p2, &weights);
    w = static_cast<SparseMatrix::EntryType>(1.0 / (sqrt(Distance2BetweenPoints(p1, p2, d)) + EPSILON));
    adjw[r1 + ptId1].push_back(MakePair(c1 + static_cast<int>(ptId2), w));
    adjw[r1 + ptId2].push_back(MakePair(c1 + static_cast<int>(ptId1), w));
  }

  // Free memory
  delete[] p1;
  delete[] p2;

  MIRTK_DEBUG_TIMING(7, "calculating adjacency weights");
}

// -----------------------------------------------------------------------------
void AdjacencyMatrix(SparseMatrix  &adjw,
                     vtkPolyData   *dataset,
                     FeatureWeights weights)
{
  const int n = dataset->GetNumberOfPoints();
  if (n == 0) {
    adjw.Initialize(0);
    return;
  }
  SparseMatrix::Entries *entries = new SparseMatrix::Entries[n];
  AdjacencyMatrix(entries, adjw.Layout(), 0, 0, dataset, weights);
  adjw.Initialize(n, n, entries);
  adjw.Index();
  delete[] entries;
}

// -----------------------------------------------------------------------------
SparseMatrix AdjacencyMatrix(vtkPolyData                *dataset,
                             SparseMatrix::StorageLayout layout,
                             FeatureWeights              weights)
{
  SparseMatrix adjw(layout);
  AdjacencyMatrix(adjw, dataset, weights);
  return adjw;
}

// -----------------------------------------------------------------------------
/// Calculate inter-mesh affinity weights (with links at subsampled points only)
template <typename TCorWeight>
struct CalculateConnectivityWeights
{
  typedef GenericSparseMatrix<TCorWeight>                WeightMatrix;
  typedef typename WeightMatrix::Entries                 WeightEntries;
  typedef typename WeightMatrix::Entries::const_iterator WeightIterator;

  vtkPolyData                *_Target;
  const Array<int>     *_TargetSample;
  vtkPolyData                *_Source;
  const Array<int>     *_SourceSample;
  const WeightMatrix         *_CorWeight;
  bool                        _CorSamples;
  bool                        _CorTranspose;
  SparseMatrix::Entries      *_LnkWeight;
  SparseMatrix::StorageLayout _Layout;
  int                         _M, _N;
  int                         _RowOffset;
  int                         _ColOffset;

  void operator ()(const blocked_range<int> &re) const
  {
    int            i1, i2, r, c;
    double         p1[3], p2[3], d;
    WeightEntries  weight;
    WeightIterator wend;

    for (int i = re.begin(); i != re.end(); ++i) {
      i1 = PointIndex(_Target, _TargetSample, i);
      if (_CorTranspose) {
        if (_Layout == SparseMatrix::CRS) _CorWeight->GetCol(_CorSamples ? i : i1, weight);
        else                              _CorWeight->GetRow(_CorSamples ? i : i1, weight);
      } else {
        if (_Layout == SparseMatrix::CRS) _CorWeight->GetRow(_CorSamples ? i : i1, weight);
        else                              _CorWeight->GetCol(_CorSamples ? i : i1, weight);
      }
      if (!weight.empty()) {
        wend = weight.end();
        do { --wend; } while (wend->first >= _N && wend != weight.begin());
        if (wend->first < _N) ++wend;
        r = i1 + _RowOffset;
        double wsum = .0;
        for (WeightIterator w = weight.begin(); w != wend; ++w) {
          wsum += w->second;
        }
        if (wsum > .0) {
          for (WeightIterator w = weight.begin(); w != wend; ++w) {
            i2 = (_CorSamples ? PointIndex(_Source, _SourceSample, w->first) : w->first);
            c  = i2 + _ColOffset;
            _Source->GetPoint(i2, p2);
            d = sqrt(Distance2BetweenPoints(p1, p2));
            _LnkWeight[r].push_back(MakePair(c, (w->second / wsum) / (d + EPSILON)));
          }
        }
      }
    }
  }
};

// -----------------------------------------------------------------------------
template <typename TCorWeight>
void ConnectivityMatrix(SparseMatrix::Entries                  lnkw[],
                        SparseMatrix::StorageLayout            layout,
                        int                                    r1,
                        int                                    c1,
                        vtkPolyData                           *target,
                        const Array<int>                *target_sample,
                        vtkPolyData                           *source,
                        const Array<int>                *source_sample,
                        const GenericSparseMatrix<TCorWeight> *corw,
                        bool                                   transpose)
{
  MIRTK_START_TIMING();
  CalculateConnectivityWeights<TCorWeight> body;
  body._M = NumberOfPoints(target, target_sample);
  body._N = NumberOfPoints(source, source_sample);
  if (body._M == 0 || body._N == 0) return;
  body._Target       = target;
  body._TargetSample = target_sample;
  body._Source       = source;
  body._SourceSample = source_sample;
  body._CorWeight    = corw;
  body._CorTranspose = transpose;
  body._LnkWeight    = lnkw;
  body._Layout       = layout;
  body._RowOffset    = r1;
  body._ColOffset    = c1;
  if (layout == SparseMatrix::CCS) {
    swap(body._Target,       body._Source);
    swap(body._TargetSample, body._SourceSample);
    swap(body._M,            body._N);
    swap(body._RowOffset,    body._ColOffset);
  }
  bool all_points = false;
  if (body._CorTranspose) {
    if (body._Layout == SparseMatrix::CRS) {
      all_points = (body._CorWeight->Rows() >= body._Source->GetNumberOfPoints());
    } else {
      all_points = (body._CorWeight->Cols() >= body._Source->GetNumberOfPoints());
    }
  } else {
    if (body._Layout == SparseMatrix::CRS) {
      all_points = (body._CorWeight->Cols() >= body._Source->GetNumberOfPoints());
    } else {
      all_points = (body._CorWeight->Rows() >= body._Source->GetNumberOfPoints());
    }
  }
  if (all_points) {
    body._N          = body._Source->GetNumberOfPoints();
    body._CorSamples = false;
  } else {
    body._CorSamples = true;
  }
  blocked_range<int> range(0, body._M);
  parallel_for(range, body);
  MIRTK_DEBUG_TIMING(7, "calculating connectivity weights");
}

template void ConnectivityMatrix<float>(
  SparseMatrix::Entries [], SparseMatrix::StorageLayout, int, int,
  vtkPolyData *, const Array<int> *, vtkPolyData *, const Array<int> *,
  const GenericSparseMatrix<float> *, bool
);
template void ConnectivityMatrix<double>(
  SparseMatrix::Entries [], SparseMatrix::StorageLayout, int, int,
  vtkPolyData *, const Array<int> *, vtkPolyData *, const Array<int> *,
  const GenericSparseMatrix<double> *, bool
);

// -----------------------------------------------------------------------------
template <typename TCorWeight>
void ConnectivityMatrix(SparseMatrix                          &lnkw,
                        vtkPolyData                           *target,
                        const Array<int>                *target_sample,
                        vtkPolyData                           *source,
                        const Array<int>                *source_sample,
                        const GenericSparseMatrix<TCorWeight> *weight,
                        bool                                   transpose)
{
  const int m = target->GetNumberOfPoints();
  const int n = source->GetNumberOfPoints();
  if (m == 0 || n == 0) return;
  SparseMatrix::Entries *entries;
  entries = new SparseMatrix::Entries[lnkw.Layout() == SparseMatrix::CRS ? m : n];
  ConnectivityMatrix(entries, lnkw.Layout(), 0, 0,
                     target, target_sample, source, source_sample, weight, transpose);
  lnkw.Initialize(m, n, entries);
  delete[] entries;
}

template void ConnectivityMatrix<float >(
  SparseMatrix &, vtkPolyData *, const Array<int> *,
                  vtkPolyData *, const Array<int> *, const GenericSparseMatrix<float> *, bool
);
template void ConnectivityMatrix<double>(
  SparseMatrix &, vtkPolyData *, const Array<int> *,
                  vtkPolyData *, const Array<int> *, const GenericSparseMatrix<double> *, bool
);

// -----------------------------------------------------------------------------
template <typename TCorWeight>
SparseMatrix ConnectivityMatrix(vtkPolyData                           *target,
                                const Array<int>                *target_sample,
                                vtkPolyData                           *source,
                                const Array<int>                *source_sample,
                                const GenericSparseMatrix<TCorWeight> *weight,
                                bool                                   transpose,
                                SparseMatrix::StorageLayout            layout)
{
  SparseMatrix lnkw(layout);
  ConnectivityMatrix(lnkw, target, target_sample, source, source_sample, weight, transpose);
  return lnkw;
}

template SparseMatrix ConnectivityMatrix<float>(
  vtkPolyData *, const Array<int> *, vtkPolyData *, const Array<int> *,
  const GenericSparseMatrix<float> *, bool, SparseMatrix::StorageLayout
);
template SparseMatrix ConnectivityMatrix<double>(
  vtkPolyData *, const Array<int> *, vtkPolyData *, const Array<int> *,
  const GenericSparseMatrix<double> *, bool, SparseMatrix::StorageLayout
);

// -----------------------------------------------------------------------------
void Degree(Vector &D, const SparseMatrix &A)
{
  mirtkAssert(A.Rows() == A.Cols(), "adjacency matrix must be square");
  const int n = A.Rows();
  D.Initialize(n);
  if (A.Layout() == SparseMatrix::CRS) {
    for (int i = 0; i < n; ++i) D(i) = A.RowSum(i);
  } else {
    for (int i = 0; i < n; ++i) D(i) = A.ColSum(i);
  }
}

// -----------------------------------------------------------------------------
Vector Degree(const SparseMatrix &A)
{
  Vector D;
  Degree(D, A);
  return D;
}

// -----------------------------------------------------------------------------
void DegreeMatrix(SparseMatrix &D, const SparseMatrix &A)
{
  mirtkAssert(A.Rows() == A.Cols(), "adjacency matrix must be square");
  const int n = A.Rows();
  D.Initialize(n, n, n);
  if (A.Layout() == SparseMatrix::CRS) {
    for (int i = 0; i < n; ++i) D(i, i) = A.RowSum(i);
  } else {
    for (int i = 0; i < n; ++i) D(i, i) = A.ColSum(i);
  }
}

// -----------------------------------------------------------------------------
SparseMatrix DegreeMatrix(const SparseMatrix &A)
{
  SparseMatrix D(A.Layout());
  DegreeMatrix(D, A);
  return D;
}

// -----------------------------------------------------------------------------
/// Auxiliary functor used by GeneralLaplacian
struct NegateAffinitiesAndDivideByNodeDegree
{
  SparseMatrix     *_A;
  const Vector *_D;
  const Vector *_G;
  bool              _DOwner;

  NegateAffinitiesAndDivideByNodeDegree(SparseMatrix &A, const Vector *D = NULL, const Vector *G = NULL)
  :
    _A(&A), _D(D), _G(G && G->Rows() > 0 ? G : NULL), _DOwner(false)
  {
    if (_A->Layout() == SparseMatrix::CCS) _A->Index();
    if (_D == NULL || _D->Rows() == 0) {
      _D      = new Vector(_A->Rows());
      _DOwner = true;
      Degree(*const_cast<Vector *>(_D), A);
    }
  }

  NegateAffinitiesAndDivideByNodeDegree(const NegateAffinitiesAndDivideByNodeDegree &other)
  :
    _A(other._A), _D(other._D), _G(other._G), _DOwner(false)
  {}

  ~NegateAffinitiesAndDivideByNodeDegree()
  {
    if (_DOwner) delete _D;
  }

  void operator ()(const blocked_range<int> &row) const
  {
    if (_G) {
      for (int r = row.begin(); r != row.end(); ++r) {
        _A->ScaleRow(r, -_G->Get(r) / _D->Get(r));
      }
    } else {
      for (int r = row.begin(); r != row.end(); ++r) {
        _A->ScaleRow(r, -1.0 / _D->Get(r));
      }
    }
  }
};

// -----------------------------------------------------------------------------
void NormalizedLaplacian(SparseMatrix &L, const SparseMatrix &A,
                         const Vector *G, const Vector *D)
{
  const int n = A.Rows();

  mirtkAssert(A.Cols() == A.Rows(), "adjacency matrix must be square");
  mirtkAssert(G == NULL || G->Rows() == 0 || G->Rows() == n,
              "weight either given for every node or none at all");
  mirtkAssert(D == NULL || D->Rows() == 0 || D->Rows() == n,
              "degree either given for every node or none at all");

  if (G && G->Rows() == 0) G = NULL;

  MIRTK_START_TIMING();

  // - G D^-1 A
  L = A;
  NegateAffinitiesAndDivideByNodeDegree mul(L, D, G);
  blocked_range<int> rows(0, n);
  parallel_for(rows, mul);

  // + G
  if (G) L.Diag( *G);
  else   L.Diag(1.0);

  MIRTK_DEBUG_TIMING(7, "calculating the normalized Laplacian");
}

// -----------------------------------------------------------------------------
void Laplacian(SparseMatrix &L, vtkPolyData *dataset,
               FeatureWeights edge_weights, FeatureWeights node_weights)
{
  mirtkAssert(dataset->GetNumberOfPoints() > 1, "dataset must have at least two points");

  MIRTK_START_TIMING();

  // Scale edge weights by range of spatial coordinates
  if (!edge_weights.empty()) {
    double p[3], pmin[3], pmax[3];
    dataset->GetPoint(0, pmin);
    dataset->GetPoint(0, pmax);
    for (vtkIdType i = 1; i < dataset->GetNumberOfPoints(); ++i) {
      dataset->GetPoint(i, p);
      for (int d = 0; d < 3; ++d) {
        if (p[d] < pmin[d]) pmin[d] = p[d];
        if (p[d] > pmax[d]) pmax[d] = p[d];
      }
    }
    double max_xyz_range = max(max(pmax[0] - pmin[0],
                                   pmax[1] - pmin[1]),
                                   pmax[2] - pmin[2]);
    FeatureWeights::iterator ew;
    for (ew = edge_weights.begin(); ew != edge_weights.end(); ++ew) {
      ew->second *= max_xyz_range;
    }
  }

  // Weighted adjacency matrix (with edge weights)
  AdjacencyMatrix(L, dataset, edge_weights);

  // Node weighting
  Vector D, G;
  Degree(D, L);
  if (!node_weights.empty()) {
    const int n = dataset->GetNumberOfPoints();
    vtkPointData * const data = dataset->GetPointData();

    // Scale node weights by range of node degrees
    double min_degree = D(0), max_degree = D(0);
    for (int i = 1; i < n; ++i) {
      if (D(i) < min_degree) min_degree = D(i);
      if (D(i) > max_degree) max_degree = D(i);
    }
    double degree_range = max_degree - min_degree;
    FeatureWeights::iterator nw;
    for (nw = node_weights.begin(); nw != node_weights.end(); ++nw) {
      nw->second *= degree_range;
    }

    // Sum weighted node weights
    G.Initialize(n, .0);
    for (nw = node_weights.begin(); nw != node_weights.end(); ++nw) {
      vtkDataArray * const feature = data->GetArray(nw->first);
      for (int i = 0; i < n; ++i) {
        for (int j = 0; j < feature->GetNumberOfComponents(); ++j) {
          G(i) += nw->second * PositiveNodeWeight(feature->GetComponent(i, j));
        }
      }
    }

    // Divide by number of weights in summation
    int num = 0;
    for (nw = node_weights.begin(); nw != node_weights.end(); ++nw) {
      num += data->GetArray(nw->first)->GetNumberOfComponents();
    }
    if (num > 0) {
      for (int i = 0; i < n; ++i) G(i) /= num;
    }
  }

  // Laplacian matrix
  NormalizedLaplacian(L, L, &G, &D);

  MIRTK_DEBUG_TIMING(7, "calculating the graph Laplacian");
}

// =============================================================================
// Spectral decomposition
// =============================================================================

// -----------------------------------------------------------------------------
int ComputeEigenmodes(const SparseMatrix &L, int k, Matrix &m, Vector &v)
{
  MIRTK_START_TIMING();

  // Repeat eigen decomposition until enough eigenmodes are found
  const double     eps     = 1e-10; // Discard eigenvectors with eigenvalue below
  int              niter   = 0;     // Number of iterations
  int              nev     = 0;     // Number of eigenvectors (incl. null vector)
  int              nconv   = 0;     // Number of converged eigenvalues
  int              fiedler = 1;     // Index of Fiedler vector
  Array<int> idx;             // Ascending order of eigenvalues
  Matrix           vec;             // Eigenvectors
  Vector           val;             // Eigenvalues

  do {

    // Request k plus expected number of null vectors
    nev = k + fiedler;

    // Perform iterative eigen decomposition
    nconv = L.Eigenvectors(vec, val, nev, "sm");

    // Determine ascending order
    idx.resize(nconv);
    Array<bool> done(nconv, false);
    for (int i = 0; i < nconv; ++i) {
      double min = numeric_limits<double>::max();
      for (int j = 0; j < nconv; ++j) {
        if (!done[j] && val(j) < min) min = val(j), idx[i] = j;
      }
      done[idx[i]] = true;
    }

    // Determine Fiedler vector
    fiedler = 1;
    while (fiedler < nconv && val(idx[fiedler]) <= eps) ++fiedler;

  } while ((++niter) < 3 && nconv == nev && nconv - fiedler < k);

  if (nconv - fiedler < min(3, k)) {
    cerr << "Eigen decomposition failed to find sufficient number of eigenmodes" << endl;
    exit(1);
  }
  k = nconv - fiedler;

  // Remove indices of duplicate eigenvalues
  for (int i = fiedler+1; i < nconv; ++i) {
    if (fequal(val(idx[i]), val(idx[i-1]))) {
      idx.erase(idx.begin() + i);
      --i, --nconv, --k;
    }
  }

  // Get eigenmodes
  m.Initialize(vec.Rows(), k);
  v.Initialize(k);
  for (int i = fiedler, j = 0; j < k; ++i, ++j) {
    int c = idx[i];
    for (int r = 0; r < vec.Rows(); ++r) m(r, j) = vec(r, c);
    v(j) = val(c);
  }

  MIRTK_DEBUG_TIMING(7, "calculating the eigenmodes (#iter=" << niter << ")");
  return k;
}

// -----------------------------------------------------------------------------
int ComputeEigenmodes(vtkPolyData *d, int k, Matrix &m, Vector &v,
                      FeatureWeights ew, FeatureWeights nw)
{
  SparseMatrix L(SparseMatrix::CCS);
  Laplacian(L, d, ew, nw);
  return ComputeEigenmodes(L, k, m, v);
}

// -----------------------------------------------------------------------------
int ComputeEigenmodes(vtkPolyData *d, int k, FeatureWeights ew, FeatureWeights nw)
{
  SparseMatrix L(SparseMatrix::CCS);
  Laplacian(L, d, ew, nw);
  Matrix m;
  Vector v;
  k = ComputeEigenmodes(L, k, m, v);
  if (k > 0) SetEigenmodes(d, m);
  return k;
}

// -----------------------------------------------------------------------------
/// Reorder eigenmodes to correct for sign ambiguity and multiplicity
///
/// \param[in]     lp1 First set of downsampled points.
/// \param[in,out] lm1 First set of downsampled eigenmodes which are rescaled to [-.5, .5].
/// \param[in]      m1 First set of original eigenmodes.
/// \param[in]      v1 Eigenvalues corresponding to first set of eigenmodes.
/// \param[in]     lp2 Second set of downsampled points.
/// \param[in,out] lm2 Second set of downsampled eigenmodes which are rescaled to [-.5, .5].
/// \param[in,out]  m2 Second set of original eigenmodes which are reordered.
/// \param[in]      v2 Eigenvalues corresponding to second set of eigenmodes.
/// \param[in]     knn Number of spatially nearest neighbors to consider.
///
/// \returns Confidence/cost of each match.
Vector MatchEigenmodes(const PointSet &lp1, Matrix       &lm1,
                       const Matrix   &m1,  const Vector &v1,
                       const PointSet &lp2, Matrix       &lm2,
                       Matrix         &m2,  Vector       &v2,
                       int knn = 5)
{
  mirtkAssert(m1.Cols() <= m2.Cols(), "second dataset must have at least as many eigenmodes as first dataset");
  const int k = m1.Cols();
  // FOCUSR MATLAB implementation normalizes point coordinates of the two sets
  // independent of each other to [0, 1]. Does this not change the Euclidean
  // distance relationship of the points? All we need the spatial coordinates
  // for is to determine the closest points in space to set up the cost matrix.
  // Therefore, no normalization of the point coordinates is being done here.
#if 0
  NormalizePoints::Run(lp1);
  NormalizePoints::Run(lp2);
#endif
  // Normalize (downsampled) spectral coordinates to [-.5, .5]
  RescaleEigenmodes::Run(lm1);
  RescaleEigenmodes::Run(lm2);
  // Compute for each pair of eigenmodes and sign flip the assignment cost
  Array<int> flip;
  Array<Array<int> > r12 = FindClosestNPoints::Run(lp1, lp2, 1);
  Matrix w12 = ComputeAssignmentCosts::Run(lp1, lm1, lp2, lm2, r12, &flip);
  // Find optimal assignment
  Array<int> c12 = WeightedMatching(w12);
  // Flip signs
  Array<int> idx(k);
  for (int i = 0; i < k; ++i) idx[i] = w12.Index(i, c12[i]);
  sort(idx.begin(), idx.end());
  flip.resize(set_intersection(flip.begin(), flip.end(),
                               idx .begin(), idx .end(),
                               flip.begin()) - flip.begin());
  for (size_t i = 0; i < flip.size(); ++i) flip[i] = w12.ColIndex(flip[i]);
  FlipSign(m2, flip);
  // Reorder eigenmodes
  m2.PermuteCols(c12);
  v2.PermuteRows(c12);
  // Return confidence of matches
  Vector w(k);
  for (int i = 0; i < k; ++i) w(i) = w12(i, c12[i]);
  return w;
}

// -----------------------------------------------------------------------------
Vector MatchEigenmodes(const PointSet &p1, const Matrix &m1, const Vector &v1,
                       const PointSet &p2,       Matrix &m2,       Vector &v2,
                       int ratio, int knn)
{
  // Downsample points
  //
  // FIXME: The downsampling must ensure a uniform distribution of the point
  //        samples. Otherwise, the random samples may be to unbalanced and
  //        result in a misleading bipartite match, i.e., sign flip.
#if 0
  if (ratio < 1) ratio = 1;
  Array<int> perm;
  PointSet lp1, lp2;
  Matrix   lm1, lm2;
  perm = RandPerm(0, p1.Size(), p1.Size() / ratio);
  GetPoints    (lp1, p1, perm);
  GetEigenmodes(lm1, m1, perm);
  perm = RandPerm(0, p2.Size(), p2.Size() / ratio);
  GetPoints    (lp2, p2, perm);
  GetEigenmodes(lm2, m2, perm);
#else
  PointSet lp1(p1), lp2(p2);
  Matrix   lm1(m1), lm2(m2);
#endif
  // Match eigenmodes
  return MatchEigenmodes(lp1, lm1, m1, v1, lp2, lm2, m2, v2, knn);
}

// -----------------------------------------------------------------------------
Vector MatchEigenmodes(vtkPoints *p1, const Matrix &m1, const Vector &v1,
                       vtkPoints *p2,       Matrix &m2,       Vector &v2,
                       int ratio, int knn)
{
  // Downsample points
  //
  // FIXME: The downsampling must ensure a uniform distribution of the point
  //        samples. Otherwise, the random samples may be to unbalanced and
  //        result in a misleading bipartite match, i.e., sign flip.
#if 0
  if (ratio < 1) ratio = 1;
  Array<int> perm;
  PointSet lp1, lp2;
  Matrix   lm1, lm2;
  perm = RandPerm(0, p1->GetNumberOfPoints(), p1->GetNumberOfPoints() / ratio);
  GetPoints    (lp1, p1, perm);
  GetEigenmodes(lm1, m1, perm);
  perm = RandPerm(0, p2->GetNumberOfPoints(), p2->GetNumberOfPoints() / ratio);
  GetPoints    (lp2, p2, perm);
  GetEigenmodes(lm2, m2, perm);
#else
  Matrix lm1(m1), lm2(m2);
  PointSet lp1(static_cast<int>(p1->GetNumberOfPoints()));
  PointSet lp2(static_cast<int>(p2->GetNumberOfPoints()));
  double p[3];
  for (vtkIdType i = 0; i < p1->GetNumberOfPoints(); ++i) {
    p1->GetPoint(i, p);
    lp1.SetPoint(i, p);
  }
  for (vtkIdType i = 0; i < p2->GetNumberOfPoints(); ++i) {
    p2->GetPoint(i, p);
    lp2.SetPoint(i, p);
  }
#endif
  // Match eigenmodes
  return MatchEigenmodes(lp1, lm1, m1, v1, lp2, lm2, m2, v2, knn);
}

// -----------------------------------------------------------------------------
void WeightEigenmodes(const Matrix &m1, const Vector &v1,
                      const Matrix &m2, const Vector &v2, Vector &w)
{
  mirtkAssert(m1.Cols() <= m2.Cols(), "second dataset must have at least as many eigenmodes as first dataset");
  const int k = w.Rows();
  double wsum = .0;
  for (int i = 0; i < k; ++i) {
    // Note: FOCUSR MATLAB implementation uses sum(M' * Q')', which is
    //       equal to sum(Q, 2) no matter what the permutation matrix M is.
    //       This is presumably a mistake in the MATLAB code.
    wsum += (w(i) *= max(v1(i), v2(i)));
  }
  (wsum /= k) *= 2.0;
  for (int i = 0; i < k; ++i) {
    w(i) = exp(- w(i) * w(i) / wsum);
  }
}

// -----------------------------------------------------------------------------
void NormalizeEigenmodes(Matrix &m)
{
  RescaleEigenmodes::Run(m);
}

// -----------------------------------------------------------------------------
void ScaleEigenmodes(Matrix &m, const Vector &w)
{
  ScaleCols::Run(m, w);
}

// -----------------------------------------------------------------------------
void ScaleEigenmodes(Matrix &m1, Matrix &m2, const Vector &w)
{
  ScaleCols::Run(m1, w);
  ScaleCols::Run(m2, w);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkDataArray> ToPointDataArray(const Matrix &m, int r1, int n, int k, const char *name)
{
  if (n <  0) n = m.Rows();
  if (k <= 0 || k > m.Cols()) k = m.Cols();
  mirtkAssert(r1 + n <= m.Rows(), "spectral components given for each point");
  vtkSmartPointer<vtkDataArray> s = vtkSmartPointer<vtkFloatArray>::New();
  if (name) s->SetName(name);
  s->SetNumberOfComponents(k);
  s->SetNumberOfTuples(n);
  for (int r = 0; r < n; ++r) {
    for (int c = 0; c < k; ++c) s->SetComponent(r, c, m(r1 + r, c));
  }
  return s;
}

// -----------------------------------------------------------------------------
int SetEigenmodes(vtkPolyData *d, const Matrix &m, const char *name)
{
  return SetEigenmodes(d, m, 0, -1, name);
}

// -----------------------------------------------------------------------------
int SetEigenmodes(vtkPolyData *d, const Matrix &m, int r1, int k, const char *name)
{
  const int n = d->GetNumberOfPoints();
  if (k <= 0 || k > m.Cols()) k = m.Cols();
  mirtkAssert(r1 + n <= m.Rows(), "spectral components given for each point");
  vtkSmartPointer<vtkDataArray> s;
  int                           i;
  if (name) s = d->GetPointData()->GetArray(name, i);
  if (!s) {
    s = vtkSmartPointer<vtkFloatArray>::New();
    if (name) s->SetName(name);
    i = d->GetPointData()->AddArray(s);
  }
  s->SetNumberOfComponents(k);
  s->SetNumberOfTuples(n);
  for (int r = 0; r < n; ++r) {
    for (int c = 0; c < k; ++c) s->SetComponent(r, c, m(r1 + r, c));
  }
  return i;
}

// -----------------------------------------------------------------------------
Matrix GetEigenmodes(vtkPolyData *d, const char *name)
{
  mirtkAssert(name != NULL, "array name must be specified");
  vtkDataArray *s = d->GetPointData()->GetArray(name);
  if (s == NULL) return Matrix();
  const int n = static_cast<int>(s->GetNumberOfTuples());
  const int k = static_cast<int>(s->GetNumberOfComponents());
  Matrix m(n, k);
  double *row = new double[k];
  for (vtkIdType r = 0; r < n; ++r) {
    s->GetTuple(r, row);
    for (int c = 0; c < k; ++c) m(r, c) = row[c];
  }
  delete[] row;
  return m;
}

// -----------------------------------------------------------------------------
/// Compute spectral components and account for sign ambiguity and multiplicity
///
/// The signs of the computed eigenmodes of the second dataset are flipped
/// and the eigenmodes reordered in order to minimize the bipartite assignment
/// cost between the eigenmodes of two datasets. This accounts for the sign
/// ambiguity and multiplicity of the eigenvalue problem. The assignment costs
/// are then used to scale the eigenmodes. The resulting spectral components
/// are added to the point data of the respective vtkPolyData instance with
/// name "eigenmodes". An existing point data array with this name will be
/// overwritten.
///
/// \param[in,out] d1  Surface mesh or contour with optional extra
///                    features used as edge and/or node weights.
/// \param[in,out] d2  Surface mesh or contour with optional extra
///                    features used as edge and/or node weights.
/// \oaram[in]     k   Number of desired spectral components.
/// \param[in]     ew  Weights of edge features if any (cf. EdgeWeights).
/// \param[in]     nw  Weights of node features if any (cf. NodeWeights).
///
/// \returns Actual number of spectral components found for both datasets.
int ComputeEigenmodes(vtkPolyData *d1, vtkPolyData *d2, int k,
                      FeatureWeights ew, FeatureWeights nw)
{
  Matrix m1, m2;
  Vector v1, v2, w;
  // Compute eigenmodes of each dataset
  int k1 = ComputeEigenmodes(d1, k, m1, v1, ew, nw);
  int k2 = ComputeEigenmodes(d2, k, m2, v2, ew, nw);
  k = min(k1, k2);
  // Match sign and order of eigenmodes
  if (k1 <= k2) {
    w = MatchEigenmodes(d1->GetPoints(), m1, v1, d2->GetPoints(), m2, v2);
  } else {
    w = MatchEigenmodes(d2->GetPoints(), m2, v2, d1->GetPoints(), m1, v1);
  }
  // Scale eigenmodes
  WeightEigenmodes(m1, v1, m2, v2, w);
  ScaleEigenmodes (m1,     m2,     w);
  // Add eigenmodes to point data of datasets
  SetEigenmodes(d1, m1, 0, k);
  SetEigenmodes(d2, m2, 0, k);
  return k;
}


} } // namespace mirtk::SpectralDecomposition
