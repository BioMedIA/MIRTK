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

#ifndef MIRTK_SpectralDecomposition_H
#define MIRTK_SpectralDecomposition_H

#include "mirtk/Pair.h"
#include "mirtk/Array.h"
#include "mirtk/PointSet.h"
#include "mirtk/Vector.h"
#include "mirtk/Matrix.h"
#include "mirtk/SparseMatrix.h"
#include "mirtk/Parallel.h"

#include "mirtk/PointSetExport.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"
#include "vtkIdList.h"
#include "vtkAbstractPointLocator.h"


namespace mirtk {
namespace SpectralDecomposition {


// =============================================================================
// Types and constants
// =============================================================================

typedef GenericSparseMatrix<double>  SparseMatrix;
typedef Array<Pair<int, double> >    FeatureWeights;

const double EPSILON = 1e-6;

// =============================================================================
// Bipartite graph matching / Optimal assignment problem
// =============================================================================

/// Find optimal assignment given cost matrix using the Hungarian method
Array<int> WeightedMatching(const Matrix &cost);

// =============================================================================
// Find closest points
// =============================================================================

// -----------------------------------------------------------------------------
/// For each point in the first set, find closest point in the second set
class FindClosestPoints
{
  // Attention: vtkKdTree/vtkKdTreePointLocator is not thread-safe
  //            (cf. http://www.vtk.org/Bug/view.php?id=15206 )!
  const PointSet          &_P1;
  vtkAbstractPointLocator *_P2;
  Array<int>              &_Idx;

  FindClosestPoints(const PointSet &, vtkAbstractPointLocator *, Array<int> &);

public:

  /// Internal function executed by worker threads
  void operator ()(const blocked_range<int> &re) const
  {
    for (int i = re.begin(); i != re.end(); ++i) {
      const Point &p = _P1(i);
      _Idx[i] = static_cast<int>(_P2->FindClosestPoint(p._x, p._y, p._z));
    }
  }

  /// For each point in the first set, find closest point in the second set
  static Array<int> Run(const PointSet &p1, const PointSet &p2);
};

// -----------------------------------------------------------------------------
/// For each point in the first set, find closest point in the second set
class FindClosestNPoints
{
  // Attention: vtkKdTree/vtkKdTreePointLocator is not thread-safe
  //            (cf. http://www.vtk.org/Bug/view.php?id=15206 )!
  const PointSet          &_P1;
  vtkAbstractPointLocator *_P2;
  Array<Array<int> >      &_Idx;
  int                      _Num;

  FindClosestNPoints(const PointSet &, vtkAbstractPointLocator *, Array<Array<int> > &, int);

public:

  /// Internal function executed by worker threads
  void operator ()(const blocked_range<int> &re) const
  {
    double p[3];
    vtkIdList *ids = vtkIdList::New();
    for (int i = re.begin(); i != re.end(); ++i) {
      _P1.GetPoint(i, p);
      _P2->FindClosestNPoints(_Num, p, ids);
      for (vtkIdType j = 0; j < ids->GetNumberOfIds(); ++j) {
        _Idx[i][j] = static_cast<int>(ids->GetId(j));
      }
    }
    ids->Delete();
  }

  /// For each point in the first set, find closest point in the second set
  static Array<Array<int> > Run(const PointSet &p1, const PointSet &p2, int);
};

// =============================================================================
// Node/edge feature weights
// =============================================================================

/// Get absolute edge weights
///
/// \param[in] dataset  Surface mesh.
/// \param[in] names    Names of feature arrays in \p dataset, e.g., "curvature",
///                     "sulcal depth", "cortical thickness", ..., to use.
/// \param[in] weights  Relative weights of edge features. If \c NULL, the
///                     default weight is used for all features. If only one
///                     weight value given, it is used for all features.
///                     Otherwise, a separate weight must be given for each
///                     named feature array.
///
/// \returns Pairs of feature array indices and corresponding absolute weights.
FeatureWeights EdgeWeights(vtkPolyData         *dataset,
                           const Array<string> &names,
                           const Array<double> *weights = NULL);

/// Get absolute node weights
///
/// \param[in] dataset  Surface mesh.
/// \param[in] names    Names of feature arrays in \p dataset, e.g., "curvature",
///                     "sulcal depth", "cortical thickness", ..., to use.
/// \param[in] weights  Relative weights of node features. If \c NULL, the
///                     default weight is used for all features. If only one
///                     weight value given, it is used for all features.
///                     Otherwise, a separate weight must be given for each
///                     named feature array.
///
/// \returns Pairs of feature array indices and corresponding absolute weights.
FeatureWeights NodeWeights(vtkPolyData         *dataset,
                           const Array<string> &names,
                           const Array<double> *weights = NULL);

// -----------------------------------------------------------------------------
/// Size of node feature vectors
inline int NumberOfFeatures(vtkPolyData *dataset, const FeatureWeights *weights = NULL)
{
  int d = 3;
  if (weights) {
    vtkPointData * const data = dataset->GetPointData();
    FeatureWeights::const_iterator weight = weights->begin();
    while (weight != weights->end()) {
      vtkDataArray * const feature = data->GetArray(weight->first);
      d += feature->GetNumberOfComponents();
      ++weight;
    }
  }
  return d;
}

// -----------------------------------------------------------------------------
/// Get node features
inline void GetFeatures(vtkPolyData *dataset, vtkIdType i, double *p, const FeatureWeights *weights = NULL)
{
  dataset->GetPoint(i, p);
  p += 3;
  if (weights) {
    vtkPointData * const data = dataset->GetPointData();
    FeatureWeights::const_iterator weight = weights->begin();
    while (weight != weights->end()) {
      vtkDataArray * const feature = data->GetArray(weight->first);
      feature->GetTuple(i, p);
      for (int c = 0; c < feature->GetNumberOfComponents(); ++c) {
        p[c] *= weight->second;
      }
      p += feature->GetNumberOfComponents();
      ++weight;
    }
  }
}

// -----------------------------------------------------------------------------
/// Calculate squared distance between 3-dimensional points
inline double Distance2BetweenPoints(const double *p1, const double *p2)
{
  double dx, dist2;
  dx = p2[0] - p1[0], dist2  = dx * dx;
  dx = p2[1] - p1[1], dist2 += dx * dx;
  dx = p2[2] - p1[2], dist2 += dx * dx;
  return dist2;
}

// -----------------------------------------------------------------------------
/// Calculate squared distance between d-dimensional points
inline double Distance2BetweenPoints(const double *p1, const double *p2, int d)
{
  double dx, dist2 = .0;
  for (int i = 0; i < d; ++i) {
    dx = p2[i] - p1[i];
    dist2 += dx * dx;
  }
  return dist2;
}

// =============================================================================
// Graph connectivity
// =============================================================================

/// Calculate intra-mesh affinity weights
void AdjacencyMatrix(SparseMatrix::Entries       adjw[],
                     SparseMatrix::StorageLayout layout,
                     int r1, int c1, vtkPolyData *dataset,
                     FeatureWeights weights = FeatureWeights());

/// Calculate intra-mesh affinity weights
void AdjacencyMatrix(SparseMatrix  &adjw,
                     vtkPolyData   *dataset,
                     FeatureWeights weights = FeatureWeights());

/// Calculate intra-mesh affinity weights
SparseMatrix AdjacencyMatrix(vtkPolyData                *dataset,
                             SparseMatrix::StorageLayout layout  = SparseMatrix::CCS,
                             FeatureWeights              weights = FeatureWeights());

/// Calculate inter-mesh affinity weights (with links at subsampled points only)
template <typename TCorWeight>
void ConnectivityMatrix(SparseMatrix::Entries                  lnkw[],
                        SparseMatrix::StorageLayout            layout,
                        int                                    r1,
                        int                                    c1,
                        vtkPolyData                           *target,
                        const Array<int>                      *target_sample,
                        vtkPolyData                           *source,
                        const Array<int>                      *source_sample,
                        const GenericSparseMatrix<TCorWeight> *weight,
                        bool                                   weight_transpose = false);

/// Calculate inter-mesh affinity weights (with links at subsampled points only)
template <typename TCorWeight>
void ConnectivityMatrix(SparseMatrix                          &lnkw,
                        vtkPolyData                           *target,
                        const Array<int>                      *target_sample,
                        vtkPolyData                           *source,
                        const Array<int>                      *source_sample,
                        const GenericSparseMatrix<TCorWeight> *weight,
                        bool                                   transpose = false);

/// Calculate inter-mesh affinity weights (with links at subsampled points only)
template <typename TCorWeight>
SparseMatrix ConnectivityMatrix(vtkPolyData                           *target,
                                const Array<int>                      *target_sample,
                                vtkPolyData                           *source,
                                const Array<int>                      *source_sample,
                                const GenericSparseMatrix<TCorWeight> *weight,
                                bool                                   transpose = false,
                                SparseMatrix::StorageLayout            layout = SparseMatrix::CCS);

/// Compute node degrees
void Degree(Vector &D, const SparseMatrix &A);

/// Compute node degrees
Vector Degree(const SparseMatrix &A);

/// Compute degree matrix
void DegreeMatrix(SparseMatrix &D, const SparseMatrix &A);

/// Compute degree matrix
SparseMatrix DegreeMatrix(const SparseMatrix &A);

/// Compute normalized Laplacian
///
/// \param[out] L Normalized Laplacian matrix.
/// \param[in]  A (Weighted) Adjacency matrix.
/// \param[in]  G Optional node weights. If \c NULL, G = Id.
/// \param[in]  D Main diagonal of (precomputed) degree matrix.
///
/// \returns L = G D^-1 (D - A) = G - D^-1 A
void NormalizedLaplacian(SparseMatrix &L, const SparseMatrix &A,
                         const Vector *G = NULL, const Vector *D = NULL);

/// Compute general graph Laplacian
///
/// \param[out] L            General graph Laplacian matrix.
/// \param[in]  dataset      Surface mesh or contour with optional extra
///                          feature arrays used as edge and/or node weights.
/// \param[in]  edge_weights Weights of edge features if any (cf. EdgeWeights).
/// \param[in]  node_weights Weights of node features if any (cf. NodeWeights).
void Laplacian(SparseMatrix &L, vtkPolyData *dataset,
               FeatureWeights edge_weights = FeatureWeights(),
               FeatureWeights node_weights = FeatureWeights());

// =============================================================================
// Spectral decomposition
// =============================================================================

/// Compute spectral components
///
/// This function performs a spectral analysis of the given graph Laplacian matrix.
///
/// \param[in]  L General graph Laplacian matrix.
/// \param[in]  k Number of spectral components.
/// \param[out] m Eigenmodes of \c d, i.e., spectral coordinates.
/// \param[out] v Eigenvalues of \c d, i.e., resonance frequencies.
///
/// \returns Actual number of spectral components found.
int ComputeEigenmodes(const SparseMatrix &L, int k, Matrix &m, Vector &v);

/// Compute spectral components
///
/// This function performs a spectral analysis of the general graph Laplacian
/// matrix of the specified surface mesh.
///
/// \param[in]  d   Surface mesh or contour with optional extra
///                 features used as edge and/or node weights.
/// \param[in]  k   Number of spectral components.
/// \param[out] m   Eigenmodes of \c d, i.e., spectral coordinates.
/// \param[out] v   Eigenvalues of \c d, i.e., resonance frequencies.
/// \param[in]  ew  Weights of edge features if any (cf. EdgeWeights).
/// \param[in]  nw  Weights of node features if any (cf. NodeWeights).
///
/// \returns Actual number of spectral components found.
int ComputeEigenmodes(vtkPolyData *d, int k, Matrix &m, Vector &v,
                      FeatureWeights ew = FeatureWeights(),
                      FeatureWeights nw = FeatureWeights());

/// Compute spectral components
///
/// This function performs a spectral analysis of the general graph Laplacian
/// matrix of the specified surface mesh. The components of the resulting \p k
/// eigenmodes are stored as point data of the input mesh in the array named
/// \"eigenmodes\".
///
/// \param[in,out] d   Surface mesh or contour with optional extra
///                    features used as edge and/or node weights.
/// \param[in]     k   Number of spectral components.
/// \param[in]     ew  Weights of edge features if any (cf. EdgeWeights).
/// \param[in]     nw  Weights of node features if any (cf. NodeWeights).
///
/// \returns Actual number of spectral components found.
int ComputeEigenmodes(vtkPolyData *d, int k,
                      FeatureWeights ew = FeatureWeights(),
                      FeatureWeights nw = FeatureWeights());

/// Reorder eigenmodes to correct for sign ambiguity and multiplicity
///
/// \param[in]      p1 First set of spatial coordinates.
/// \param[in]      m1 First set of original eigenmodes.
/// \param[in]      v1 Eigenvalues corresponding to first set of eigenmodes.
/// \param[in]      p2 Second set of spatial coordinates.
/// \param[in,out]  m2 Second set of original eigenmodes which are reordered.
/// \param[in]      v2 Eigenvalues corresponding to second set of eigenmodes.
/// \param[in]   ratio Downsample ratio, e.g., 2 to use only half of the points.
/// \param[in]     knn Number of spatially nearest neighbors to consider.
///
/// \returns Confidence/cost of each match.
Vector MatchEigenmodes(const PointSet &p1, const Matrix &m1, const Vector &v1,
                       const PointSet &p2, Matrix       &m2, Vector       &v2,
                       int ratio = 10, int knn = 1);

/// Reorder eigenmodes to correct for sign ambiguity and multiplicity
///
/// \param[in]      p1 First set of spatial coordinates.
/// \param[in]      m1 First set of original eigenmodes.
/// \param[in]      v1 Eigenvalues corresponding to first set of eigenmodes.
/// \param[in]      p2 Second set of spatial coordinates.
/// \param[in,out]  m2 Second set of original eigenmodes which are reordered.
/// \param[in]      v2 Eigenvalues corresponding to second set of eigenmodes.
/// \param[in]   ratio Downsample ratio, e.g., 2 to use only half of the points.
/// \param[in]     knn Number of spatially nearest neighbors to consider.
///
/// \returns Confidence/cost of each match.
Vector MatchEigenmodes(vtkPoints *p1, const Matrix &m1, const Vector &v1,
                       vtkPoints *p2, Matrix       &m2, Vector       &v2,
                       int ratio = 10, int knn = 1);

/// Compute weights of eigenmodes given the match confidence
///
/// Reduce importance of eigenmode by lowering weight for
/// - Eigenmodes corresponding to high frequencies.
/// - Pairs of eigenmodes with low confidence in the match.
///
/// \param[in]     m1 First set of eigenmodes.
/// \param[in]     v1 Eigenvalues corresponding to first set of eigenmodes.
/// \param[in]     m2 Second set of eigenmodes.
/// \param[in]     v2 Eigenvalues corresponding to second set of eigenmodes.
/// \param[in,out] w  Input: Match confidence/cost. Output: Weights of eigenmodes.
void WeightEigenmodes(const Matrix &m1, const Vector &v1,
                      const Matrix &m2, const Vector &v2, Vector &w);

/// Scale eigenmodes given their individual weights
///
/// \param[in,out] m  (Joint) Eigenmodes.
/// \param[in]     w  Weights of eigenmodes.
void ScaleEigenmodes(Matrix &m, const Vector &w);

/// Normalize eigenmodes individually to the range [-1, 1]
void NormalizeEigenmodes(Matrix &m);

/// Scale eigenmodes given their individual weights
///
/// \param[in,out] m1 First set of eigenmodes.
/// \param[in,out] m2 Second set of eigenmodes.
/// \param[in]     w  Weights of eigenmodes (cf. WeightEigenmodes).
void ScaleEigenmodes(Matrix &m1, Matrix &m2, const Vector &w);

/// Convert eigenmodes to dataset point data array
vtkSmartPointer<vtkDataArray>
ToPointDataArray(const Matrix &m, int r1 = 0,
                 int n = -1, int k = -1, const char *name = "eigenmodes");

/// Set spectral components as point data of vtkPolyData instance
///
/// \returns Index of eigenmodes point data array.
int SetEigenmodes(vtkPolyData *d, const Matrix &m, const char *name = "eigenmodes");

/// Set spectral components as point data of vtkPolyData instance
///
/// \returns Index of eigenmodes point data array.
int SetEigenmodes(vtkPolyData *d, const Matrix &m,
                  int r1, int k = -1, const char *name = "eigenmodes");

/// Get spectral components from point data of vtkPolyData instance
Matrix GetEigenmodes(vtkPolyData *d, const char *name = "eigenmodes");

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
/// \param[in]     k   Number of desired spectral components.
/// \param[in]     ew  Weights of edge features if any (cf. EdgeWeights).
/// \param[in]     nw  Weights of node features if any (cf. NodeWeights).
///
/// \returns Actual number of spectral components found for both datasets.
int ComputeEigenmodes(vtkPolyData *d1, vtkPolyData *d2, int k,
                      FeatureWeights ew = FeatureWeights(),
                      FeatureWeights nw = FeatureWeights());


} } // namespace mirtk::SpectralDecomposition

#endif // MIRTK_SpectralDecomposition_H
