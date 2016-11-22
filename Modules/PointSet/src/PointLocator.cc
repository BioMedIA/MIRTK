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

#if defined(HAVE_FLANN) && defined(_MSC_VER)
  #pragma warning(disable: 4267) // conversion from 'size_t'
  #ifndef _SCL_SECURE_NO_WARNINGS
    #define _SCL_SECURE_NO_WARNINGS
  #endif
#endif

#include "mirtk/PointLocator.h"

#include "mirtk/Assert.h"
#include "mirtk/Array.h"
#include "mirtk/ArrayHeap.h"
#include "mirtk/Pair.h"
#include "mirtk/Allocate.h"
#include "mirtk/Deallocate.h"
#include "mirtk/Parallel.h"

#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkIdList.h"
#include "vtkOctreePointLocator.h"


#ifdef HAVE_FLANN
  // Compiler warning "Using integer absolute value function 'abs' when
  // argument is of floating point type" caused by use of abs function
  // in kdtree_index.h of FLANN. Using std::abs fixes the issue.
  using std::abs;
  // Disable MSVC warnings
  #if defined(_MSC_VER)
    #pragma warning(push)
    #pragma warning(disable: 4267) // conversion from 'size_t'
    #pragma warning(disable: 4291) // no matching operator delete found
    #pragma warning(disable: 4996) // _CRT_SECURE_NO_WARNINGS
  // Disable "Unused typedef 'ElementType'" in flann/ground_truth.h
  #elif defined(__clang__) // *also* defines __GNUG__!
    #pragma clang diagnostic push
    #if !defined(__has_warning) || __has_warning("-Wunused-local-typedefs")
      #pragma clang diagnostic ignored "-Wunused-local-typedefs"
    #endif
  #elif defined(__GNUG__)
    #pragma GCC diagnostic push
    #if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ > 7)
      #pragma GCC diagnostic ignored "-Wunused-local-typedefs"
    #endif
  #endif
  // Include FLANN
  #include "flann/flann.hpp"
  // Enable warnings again
  #if defined(_MSC_VER)
    #pragma warning(pop)
  #elif defined(__clang__)
    #pragma clang diagnostic pop
  #elif defined(__GNUG__)
    #pragma GCC diagnostic pop
  #endif
#endif // HAVE_FLANN


namespace mirtk {


////////////////////////////////////////////////////////////////////////////////
// class: FlannPointLocator (using FLANN)
////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_FLANN


/**
 * Auxiliary class to wrap Kd-tree structure of third-party library
 */
class FlannPointLocator
{
public:

  /// Datatype used for internal FLANN index
  typedef float FlannType;

  /// List of point features to use for nearest neighbor search
  typedef PointLocator::FeatureList FeatureList;

protected:

  /// Dataset for which search structure is build
  mirtkPublicAggregateMacro(vtkPointSet, DataSet);

  /// Indices of points to consider only or NULL
  mirtkPublicAggregateMacro(const Array<int>, Sample);

  /// Indices/names and rescaling parameters of point data arrays
  mirtkPublicAttributeMacro(FeatureList, Features);

  /// Dimension of feature Arrays/points
  mirtkPublicAttributeMacro(int, PointDimension);

public:

  /// FLANN Kd tree used for n-dimensional feature spaces
  flann::Index<flann::L2_Simple<FlannType> > _FlannTree;

  /// Constructor
  FlannPointLocator();

  /// Destructor
  virtual ~FlannPointLocator();

  /// Initialize FLANN index
  void Initialize();

  /// Find nearest neighbor
  int FindClosestPoint(double *, double *);

  /// Find nearest neighbor
  Array<int> FindClosestPoint(vtkPointSet *, const Array<int> *,
                              const FeatureList *, Array<double> *);

  /// Find nearest neighbors
  Array<int> FindClosestNPoints(int, double *, Array<double> *);

  /// Find nearest neighbors
  Array<Array<int> >
  FindClosestNPoints(int, vtkPointSet *, const Array<int> *,
                     const FeatureList *, Array<Array<double> > *);

  /// Find points within radius
  Array<int> FindPointsWithinRadius(double, double *, Array<double> *);

  /// Find points within radius
  Array<Array<int> > FindPointsWithinRadius(double, vtkPointSet *,
                                                    const Array<int> *,
                                                    const FeatureList *,
                                                    Array<Array<double> > *);

protected:

  /// Create FLANN matrix containing feature Arrays of points of given dataset
  flann::Matrix<FlannType> FlannMatrix(vtkPointSet *,
                                       const Array<int> * = NULL,
                                       const FeatureList * = NULL);

  /// Create FLANN matrix from single feature Array
  flann::Matrix<FlannType> FlannMatrix(double *);

};

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
FlannPointLocator::FlannPointLocator()
:
  _PointDimension(0),
  _FlannTree(flann::KDTreeSingleIndexParams(10))
{
}

// ---------------------------------------------------------------------------
FlannPointLocator::~FlannPointLocator()
{
}

// -----------------------------------------------------------------------------
void FlannPointLocator::Initialize()
{
  mirtkAssert(_PointDimension > 0, "_PointDimension attribute set");
  flann::Matrix<FlannType> points = FlannMatrix(_DataSet, _Sample, &_Features);
  _FlannTree.buildIndex(points);
  delete[] points.ptr();
}

// =============================================================================
// Conversion helpers
// =============================================================================

// ---------------------------------------------------------------------------
flann::Matrix<FlannPointLocator::FlannType>
FlannPointLocator::FlannMatrix(vtkPointSet       *dataset,
                               const Array<int>  *sample,
                               const FeatureList *features)
{
  const int n = PointLocator::GetNumberOfPoints(dataset, sample);
  double *point = Allocate<double>(_PointDimension);
  flann::Matrix<FlannType> matrix(Allocate<FlannType>(n * _PointDimension), n, _PointDimension);
  for (int i = 0; i < n; ++i) {
    PointLocator::GetPoint(point, dataset, sample, i, features);
    FlannType *v = matrix[i];
    for (int j = 0; j < _PointDimension; ++j, ++v) {
      (*v) = static_cast<FlannType>(point[j]);
    }
  }
  Deallocate(point);
  return matrix;
}

// ---------------------------------------------------------------------------
flann::Matrix<FlannPointLocator::FlannType> FlannPointLocator::FlannMatrix(double *point)
{
  float *m = Allocate<FlannType>(_PointDimension);
  flann::Matrix<FlannType> matrix(m, 1, _PointDimension);
  float *v = matrix[0];
  for (int j = 0; j < _PointDimension; ++j, ++v) {
    (*v) = static_cast<FlannType>(point[j]);
  }
  return matrix;
}

// =============================================================================
// Closest point
// =============================================================================

// -----------------------------------------------------------------------------
int FlannPointLocator::FindClosestPoint(double *point, double *dist2)
{
  Array<Array<int> > indices;
  Array<Array<FlannType> > dists;
  flann::Matrix<FlannType> queries = FlannMatrix(point);
  _FlannTree.knnSearch(queries, indices, dists, 1,
                       flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
  delete[] queries.ptr();
  if (dist2) (*dist2) = dists[0][0];
  return indices[0][0];
}

// -----------------------------------------------------------------------------
Array<int> FlannPointLocator
::FindClosestPoint(vtkPointSet *dataset, const Array<int> *sample,
                   const FeatureList *features, Array<double> *dist2)
{
  Array<Array<int   > >  indices;
  Array<Array<double> > *dists;
  dists   = (dist2 ? new Array<Array<double> >() : NULL);
  indices = FindClosestNPoints(1, dataset, sample, features, dists);
  Array<int> index(indices.size());
  for (size_t i = 0; i < indices.size(); ++i) {
    index[i] = indices[i][0];
  }
  if (dist2) {
    dist2->resize(dists->size());
    for (size_t i = 0; i < dists->size(); ++i) {
      (*dist2)[i] = (*dists)[i][0];
    }
  }
  delete dists;
  return index;
}

// =============================================================================
// Nearest neighbors
// =============================================================================

// -----------------------------------------------------------------------------
Array<int> FlannPointLocator
::FindClosestNPoints(int k, double *point, Array<double> *dist2)
{
  Array<Array<int> > indices;
  Array<Array<FlannType> > dists;
  flann::Matrix<FlannType> queries = FlannMatrix(point);
  _FlannTree.knnSearch(queries, indices, dists, k,
                       flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
  delete[] queries.ptr();
  if (dist2) {
    dist2->resize(dists[0].size());
    for (size_t i = 0; i < dists[0].size(); ++i) {
      (*dist2)[i] = dists[0][i];
    }
  }
  return indices[0];
}

// -----------------------------------------------------------------------------
Array<Array<int> > FlannPointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const Array<int> *sample,
                     const FeatureList *features, Array<Array<double> > *dist2)
{
  Array<Array<int> > indices;
  Array<Array<FlannType> > dists;
  flann::Matrix<FlannType> queries = FlannMatrix(dataset, sample, features);
  _FlannTree.knnSearch(queries, indices, dists, k,
                       flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
  delete[] queries.ptr();
  if (dist2) {
    dist2->resize(dists.size());
    for (size_t i = 0; i < dists.size(); ++i) {
      Array<double> &row = (*dist2)[i];
      row.resize(dists[i].size());
      for (size_t j = 0; j < dists[i].size(); ++j) {
        row[j] = static_cast<double>(dists[i][j]);
      }
    }
  }
  return indices;
}

// =============================================================================
// Radius search
// =============================================================================

// -----------------------------------------------------------------------------
Array<int> FlannPointLocator
::FindPointsWithinRadius(double radius, double *point, Array<double> *dist2)
{
  Array<Array<int> > indices;
  Array<Array<FlannType> > dists;
  flann::Matrix<FlannType> queries = FlannMatrix(point);
  _FlannTree.radiusSearch(queries, indices, dists, radius * radius,
                          flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
  delete[] queries.ptr();
  if (dist2) {
    dist2->resize(dists[0].size());
    for (size_t i = 0; i < dist2->size(); ++i) {
      (*dist2)[i] = dists[0][i];
    }
  }
  return indices[0];
}

// -----------------------------------------------------------------------------
Array<Array<int> > FlannPointLocator
::FindPointsWithinRadius(double radius, vtkPointSet           *dataset,
                                        const Array<int>      *sample,
                                        const FeatureList     *features,
                                        Array<Array<double> > *dist2)
{
  Array<Array<int> > indices;
  Array<Array<FlannType> > dists;
  flann::Matrix<FlannType> queries = FlannMatrix(dataset, sample, features);
  _FlannTree.radiusSearch(queries, indices, dists, radius,
                          flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));
  delete[] queries.ptr();
  if (dist2) {
    dist2->resize(dists.size());
    for (size_t i = 0; i < dists.size(); ++i) {
      Array<double> &row = (*dist2)[i];
      row.resize(dists[0].size());
      for (size_t j = 0; j < row.size(); ++j) {
        row[j] = static_cast<double>(dists[i][j]);
      }
    }
  }
  return indices;
}

#endif // HAVE_FLANN
////////////////////////////////////////////////////////////////////////////////
// class: PointLocator
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
PointLocator::PointLocator()
:
  _DataSet(NULL),
  _Sample(NULL),
  _GlobalIndices(false),
  _NumberOfPoints(0),
  _PointDimension(0)
{
}

// -----------------------------------------------------------------------------
PointLocator::~PointLocator()
{
}

// -----------------------------------------------------------------------------
void PointLocator::Initialize()
{
  // Destruct previous internal locator(s)
  _VtkLocator = NULL;
#ifdef HAVE_FLANN
  _FlannLocator = nullptr;
#endif
  // Check inputs
  if (!_DataSet) {
    cerr << "PointLocator: Missing dataset!" << endl;
    exit(1);
  }
  _NumberOfPoints = GetNumberOfPoints(_DataSet, _Sample);
  if (_NumberOfPoints == 0) {
    cerr << "PointLocator: No points in search tree!" << endl;
    exit(1);
  }
  _PointDimension = GetPointDimension(_DataSet, &_Features);
  if (_PointDimension == 0) {
    cerr << "PointLocator: Point feature vector size is zero!" << endl;
    exit(1);
  }
  // Build VTK locator for 3-D feature vectors
  if (_PointDimension <= 3) {
    // Dataset of 2/3-D sample feature points
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(_NumberOfPoints);
    double point[3] = {0};
    for (int i = 0; i < _NumberOfPoints; ++i) {
      GetPoint(point, _DataSet, _Sample, i, &_Features);
      points->SetPoint(i, point);
    }
    vtkSmartPointer<vtkPolyData> dataset = vtkSmartPointer<vtkPolyData>::New();
    dataset->SetPoints(points);
    // Note: vtkOctreeLocator preferred over vtkKdTree because the
    //       latter is not thread safe even after BuildLocator was called!
    _VtkLocator = vtkSmartPointer<vtkOctreePointLocator>::New();
    _VtkLocator->SetDataSet(dataset);
    _VtkLocator->BuildLocator();
  } else {
#ifdef HAVE_FLANN
    // Build FLANN tree for N-D feature vectors
    _FlannLocator = NewShared<FlannPointLocator>();
    _FlannLocator->DataSet(_DataSet);
    _FlannLocator->Sample(_Sample);
    _FlannLocator->Features(_Features);
    _FlannLocator->PointDimension(_PointDimension);
    _FlannLocator->Initialize();
#endif
  }
}

// -----------------------------------------------------------------------------
PointLocator *PointLocator::New(vtkPointSet       *dataset,
                                const Array<int>  *sample,
                                const FeatureList *features)
{
  UniquePtr<PointLocator> locator(new PointLocator());
  locator->DataSet(dataset);
  locator->Sample(sample);
  if (features) locator->Features(*features);
  locator->Initialize();
  return locator.release();
}

// =============================================================================
// Closest point
// =============================================================================

// -----------------------------------------------------------------------------
int PointLocator::FindClosestPoint(double *point, double *dist2)
{
  int index;

  // ---------------------------------------------------------------------------
  // Using VTK
  if (_VtkLocator) {
    vtkIdType j = _VtkLocator->FindClosestPoint(point);
    if (dist2) {
      double p[3] = {0};
      _VtkLocator->GetDataSet()->GetPoint(j, p);
      *dist2 = Distance2BetweenPoints(p, point);
    }
    index = static_cast<int>(j);
  }

  // ---------------------------------------------------------------------------
  // Using FLANN
#ifdef HAVE_FLANN
  else if (_FlannLocator) {
    index = _FlannLocator->FindClosestPoint(point, dist2);
  }
#endif

  // ---------------------------------------------------------------------------
  // Brute force
  else {
    index = -1;
    Array<double> p(_PointDimension);
    double d2, mind2 = inf;
    for (int i = 0, idx; i < _NumberOfPoints; ++i) {
      idx = GetPointIndex(_DataSet, _Sample, i);
      GetPoint(p.data(), _DataSet, idx, &_Features);
      d2 = Distance2BetweenPoints(p.data(), point, _PointDimension);
      if (d2 < mind2) {
        index = i;
        mind2 = d2;
      }
    }
    if (dist2) (*dist2) = mind2;
  }

  // ---------------------------------------------------------------------------
  // Map to global point ID
  if (_Sample && _GlobalIndices) {
    index = (*_Sample)[index];
  }
  return index;
}

// -----------------------------------------------------------------------------
namespace PointLocatorUtils {
struct FindClosestPoint
{
  vtkPointSet                     *_DataSet;
  const Array<int>                *_Sample;
  const PointLocator::FeatureList *_Features;
  PointLocator                    *_Locator;
  Array<int>                      *_Index;
  Array<double>                   *_Dist2;

  void operator ()(const blocked_range<int> &idx) const
  {
    if (_Dist2) {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Index)[i] = _Locator->FindClosestPoint(_DataSet, _Sample, i, _Features, &((*_Dist2)[i]));
      }
    } else {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Index)[i] = _Locator->FindClosestPoint(_DataSet, _Sample, i, _Features);
      }
    }
  }
};
} // namespace PointLocatorUtils

// -----------------------------------------------------------------------------
Array<int> PointLocator
::FindClosestPoint(vtkPointSet *dataset, const Array<int> *sample,
                   const FeatureList *features, Array<double> *dist2)
{
#ifdef HAVE_FLANN
  if (!_VtkLocator && _FlannLocator) {
    Array<int> index;
    index = _FlannLocator->FindClosestPoint(dataset, sample, features, dist2);
    if (_Sample && _GlobalIndices) {
      for (auto &&i : index) {
        i = (*_Sample)[i];
      }
    }
    return index;
  }
#endif

  Array<int> index(GetNumberOfPoints(dataset, sample));
  if (dist2) dist2->resize(index.size());
  PointLocatorUtils::FindClosestPoint query;
  query._DataSet  = dataset;
  query._Sample   = sample;
  query._Features = features;
  query._Locator  = this;
  query._Index    = &index;
  query._Dist2    = dist2;
  parallel_for(blocked_range<int>(0, static_cast<int>(index.size())), query);
  return index;
}

// =============================================================================
// Nearest neighbors
// =============================================================================

// -----------------------------------------------------------------------------
Array<int> PointLocator::FindClosestNPoints(int k, double *point, Array<double> *dist2)
{
  Array<int> indices;

  // ---------------------------------------------------------------------------
  // Using VTK
  if (_VtkLocator) {
    double p[3];
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    _VtkLocator->FindClosestNPoints(k, point, ids);
    indices.resize(ids->GetNumberOfIds());
    if (dist2) dist2->resize(ids->GetNumberOfIds());
    for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
      indices[i] = static_cast<int>(ids->GetId(i));
      if (dist2) {
        _VtkLocator->GetDataSet()->GetPoint(ids->GetId(i), p);
        (*dist2)[i] = Distance2BetweenPoints(point, p);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Using FLANN
#ifdef HAVE_FLANN
  else if (_FlannLocator) {
    indices = _FlannLocator->FindClosestNPoints(k, point, dist2);
  }
#endif

  // ---------------------------------------------------------------------------
  // Brute force
  else {
    struct Comp
    {
      bool operator()(const Pair<double, int> &a, const Pair<double, int> &b)
      {
        return a.first > b.first;
      }
    } comp;

    Array<double> p(_PointDimension);
    Array<Pair<double, int> > dists(_NumberOfPoints);
    for (int i = 0; i < _NumberOfPoints; ++i) {
      int idx = GetPointIndex(_DataSet, _Sample, i);
      GetPoint(p.data(), _DataSet, idx, &_Features);
      dists[i] = MakePair(Distance2BetweenPoints(point, p.data(), _PointDimension), i);
    }
    make_heap(dists.begin(), dists.end(), comp);
    indices.resize(k);
    for (int i = 0; i < k; ++i) {
      Pair<double, int> &min = dists.front();
      if (dist2) (*dist2)[i] = min.first;
      indices[i] = min.second;
      pop_heap(dists.begin(), dists.end(), comp);
      dists.pop_back();
    }
  }

  // ---------------------------------------------------------------------------
  // Map to global point IDs
  if (_Sample && _GlobalIndices) {
    for (auto &&i : indices) {
      i = (*_Sample)[i];
    }
  }
  return indices;
}

// -----------------------------------------------------------------------------
namespace PointLocatorUtils {
struct FindClosestNPoints
{
  vtkPointSet                     *_DataSet;
  const Array<int>                *_Sample;
  const PointLocator::FeatureList *_Features;
  PointLocator                    *_Locator;
  Array<Array<int> >              *_Indices;
  Array<Array<double> >           *_Dist2;
  int                              _K;

  void operator ()(const blocked_range<int> &idx) const
  {
    if (_Dist2) {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Indices)[i] = _Locator->FindClosestNPoints(_K, _DataSet, _Sample, i, _Features, &((*_Dist2)[i]));
      }
    } else {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Indices)[i] = _Locator->FindClosestNPoints(_K, _DataSet, _Sample, i, _Features);
      }
    }
  }
};
} // namespace PointLocatorUtils

// -----------------------------------------------------------------------------
Array<Array<int> > PointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const Array<int> *sample,
                     const FeatureList *features, Array<Array<double> > *dist2)
{
  mirtkAssert(GetPointDimension(dataset, features) == _PointDimension,
              "Query points must have same dimension as feature points");

  if (k > _NumberOfPoints) {
    cerr << "PointLocator::FindClosestNPoints: Cannot find more points than there are in total" << endl;
    exit(1);
  }

#ifdef HAVE_FLANN
  if (!_VtkLocator && _FlannLocator) {
    Array<Array<int> > indices;
    indices = _FlannLocator->FindClosestNPoints(k, dataset, sample, features, dist2);
    if (_Sample && _GlobalIndices) {
      for (auto &&index : indices) {
        for (auto &&i : index) {
          i = (*_Sample)[i];
        }
      }
    }
    return indices;
  }
#endif

  Array<Array<int> > indices(GetNumberOfPoints(dataset, sample));
  if (dist2) dist2->resize(indices.size());
  PointLocatorUtils::FindClosestNPoints query;
  query._DataSet  = dataset;
  query._Sample   = sample;
  query._Features = features;
  query._Locator  = this;
  query._Indices  = &indices;
  query._Dist2    = dist2;
  query._K        = k;
  parallel_for(blocked_range<int>(0, static_cast<int>(indices.size())), query);
  return indices;
}

// =============================================================================
// Radius search
// =============================================================================

// -----------------------------------------------------------------------------
Array<int> PointLocator
::FindPointsWithinRadius(double radius, double *point, Array<double> *dist2)
{
  Array<int> indices;

  // ---------------------------------------------------------------------------
  // Using VTK
  if (_VtkLocator) {
    double p[3] = {0};
    vtkNew<vtkIdList> ids;
    _VtkLocator->FindPointsWithinRadius(radius, point, ids.GetPointer());
    indices.resize(ids->GetNumberOfIds());
    if (dist2) dist2->resize(ids->GetNumberOfIds());
    for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
      indices[i] = static_cast<int>(ids->GetId(i));
      if (dist2) {
        _VtkLocator->GetDataSet()->GetPoint(ids->GetId(i), p);
        (*dist2)[i] = Distance2BetweenPoints(point, p);
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Using FLANN
#ifdef HAVE_FLANN
  else if (_FlannLocator) {
    indices = _FlannLocator->FindPointsWithinRadius(radius, point, dist2);
  }
#endif

  // ---------------------------------------------------------------------------
  // Brute force
  else {
    struct Comp
    {
      bool operator()(const Pair<double, int> &a, const Pair<double, int> &b)
      {
        return a.first > b.first;
      }
    } comp;

    const double maxdist2 = radius * radius;
    Array<double> p(_PointDimension);
    Array<Pair<double, int> > dists(_NumberOfPoints);
    int idx, k = 0;
    for (int i = 0; i < _NumberOfPoints; ++i) {
      idx = GetPointIndex(_DataSet, _Sample, i);
      GetPoint(p.data(), _DataSet, idx, &_Features);
      dists[i] = MakePair(Distance2BetweenPoints(point, p.data(), _PointDimension), i);
      if (dists[i].first <= maxdist2) ++k;
    }
    make_heap(dists.begin(), dists.end(), comp);
    indices.resize(k);
    if (dist2) dist2->resize(k);
    for (int i = 0; i < k; ++i) {
      Pair<double, int> &min = dists.front();
      if (dist2) (*dist2)[i] = min.first;
      indices[i] = min.second;
      pop_heap(dists.begin(), dists.end(), comp);
      dists.pop_back();
    }
  }

  // ---------------------------------------------------------------------------
  // Map to global point IDs
  if (_Sample && _GlobalIndices) {
    for (auto &&i : indices) i = (*_Sample)[i];
  }
  return indices;
}

// -----------------------------------------------------------------------------
namespace PointLocatorUtils {
struct FindPointsWithinRadius
{
  vtkPointSet                     *_DataSet;
  const Array<int>                *_Sample;
  const PointLocator::FeatureList *_Features;
  PointLocator                    *_Locator;
  Array<Array<int> >              *_Indices;
  Array<Array<double> >           *_Dist2;
  double                           _Radius;

  void operator ()(const blocked_range<int> &idx) const
  {
    if (_Dist2) {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Indices)[i] = _Locator->FindPointsWithinRadius(_Radius, _DataSet, _Sample, i, _Features, &((*_Dist2)[i]));
      }
    } else {
      for (int i = idx.begin(); i != idx.end(); ++i) {
        (*_Indices)[i] = _Locator->FindPointsWithinRadius(_Radius, _DataSet, _Sample, i, _Features);
      }
    }
  }
};
} // namespace PointLocatorUtils

// -----------------------------------------------------------------------------
Array<Array<int> > PointLocator
::FindPointsWithinRadius(double radius, vtkPointSet           *dataset,
                                        const Array<int>      *sample,
                                        const FeatureList     *features,
                                        Array<Array<double> > *dist2)
{
  mirtkAssert(GetPointDimension(dataset, features) == _PointDimension,
              "Query points must have same dimension as feature points");

#ifdef HAVE_FLANN
  if (!_VtkLocator && _FlannLocator) {
    Array<Array<int> > indices;
    indices = _FlannLocator->FindPointsWithinRadius(radius, dataset, sample, features, dist2);
    if (_Sample && _GlobalIndices) {
      for (auto &&index : indices) {
        for (auto &&i : index) {
          i = (*_Sample)[i];
        }
      }
    }
    return indices;
  }
#endif

  Array<Array<int> > indices(GetNumberOfPoints(dataset, sample));
  if (dist2) dist2->resize(indices.size());
  PointLocatorUtils::FindPointsWithinRadius query;
  query._DataSet  = dataset;
  query._Sample   = sample;
  query._Features = features;
  query._Locator  = this;
  query._Indices  = &indices;
  query._Dist2    = dist2;
  query._Radius   = radius;
  parallel_for(blocked_range<int>(0, static_cast<int>(indices.size())), query);
  return indices;
}


} // namespace mirtk
