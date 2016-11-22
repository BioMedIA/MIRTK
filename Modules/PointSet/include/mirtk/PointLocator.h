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

#ifndef MIRTK_PointLocator_H
#define MIRTK_PointLocator_H

#include "mirtk/Object.h"

#include "mirtk/Array.h"
#include "mirtk/Point.h"
#include "mirtk/Memory.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkDataArray.h"


// Forward declaration of implementation using VTK point locator
class vtkOctreePointLocator;


namespace mirtk {


// Forward declaration of implementation using FLANN (if available)
class FlannPointLocator;


/**
 * Point search structure for establishing point correspondences
 *
 * This point locator implementation is specialized for use by
 * PointCorrespondence subclass implementations. It allows the search of
 * nearest neighbors within the n-dimensional feature space spanned by the
 * feature arrays used to establish point correspondences.
 *
 * The implementation uses either VTK for point search in three dimensions or
 * FLANN for higher dimensional feature spaces if available. Alternatively,
 * an ITK Kd tree is used if the library is available. As last resort, a brute
 * force search without actual Kd tree search structure is performed.
 *
 * \todo Consider actual implementation of own thread-safe N-dimensional Kd tree?
 *
 * \attention When only a subset of the points is used, i.e., an index array
 *            of point set samples is given (_Sample attribute), the indices
 *            returned by this locator are the respective closest sample
 *            indices into this sample point index array. Set the _GlobalIndices
 *            flag to enable the internal mapping of closest sample indices to
 *            global point indices.
 */
class PointLocator : public Object
{
  mirtkObjectMacro(PointLocator);

  // ---------------------------------------------------------------------------
  // Feature Arrays
public:

  /// Type storing name and/or index of point feature together with
  /// an optional linear rescaling function
  struct FeatureInfo
  {
    string _Name;      ///< Name of feature/point data array
    int    _Index;     ///< Index of feature/point data array
    double _Weight;    ///< Weight of feature/point data
    double _Slope;     ///< Rescaling slope of feature/point data
    double _Intercept; ///< Rescaling intercept of feature/point data

    /// Constructor
    FeatureInfo(int index = -2, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(""), _Index(index), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}

    /// Constructor
    FeatureInfo(const char *name, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(name), _Index(-2), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}

    /// Constructor
    FeatureInfo(const string &name, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(name), _Index(-2), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}

    /// Constructor
    FeatureInfo(const char *name, int index, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(name), _Index(index), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}

    /// Constructor
    FeatureInfo(const string &name, int index, double weight = 1.0, double slope = 1.0, double intercept = .0)
    :
      _Name(name), _Index(index), _Weight(weight), _Slope(slope), _Intercept(intercept)
    {}
  };

  /// List of point features to use for nearest neighbor search
  typedef Array<FeatureInfo> FeatureList;

  /// Get point data array
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] feature Index/name and weight/rescaling parameters.
  ///
  /// \returns Pointer to point data array of dataset or NULL.
  static vtkDataArray *GetDataArray(vtkPointSet *dataset, const FeatureInfo &feature);

  /// Get number of points
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  ///
  /// \returns Number of (sample) points.
  static int GetNumberOfPoints(vtkPointSet *dataset, const Array<int> *sample = NULL);

  /// Get size of feature Arrays
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] feature Indices/names and weights/rescaling parameters of features.
  static int GetPointDimension(vtkPointSet *dataset, const FeatureList *feature);

  /// Get index of feature Array
  ///
  /// \param[in] dataset Dataset.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in] index   Index of point in \p dataset or \p sample (if not NULL).
  ///
  /// \returns Index of feature Array/point in \p dataset.
  static int GetPointIndex(vtkPointSet *dataset, const Array<int> *sample, int index);

  /// Get spatial (sample) point
  ///
  /// \param[out] point   Spatial point.
  /// \param[in]  dataset Dataset.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  static void GetPoint(Point &point, vtkPointSet *dataset, const Array<int> *sample, int index);

  /// Get feature Array
  ///
  /// \param[out] point   Feature Array.
  /// \param[in]  dataset Dataset.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  feature Indices/names and weights/rescaling parameters of features.
  static void GetPoint(double *point, vtkPointSet *dataset, const Array<int> *sample,
                       int index, const FeatureList *feature = NULL);

  /// Get feature Array
  ///
  /// \param[out] point   Feature Array.
  /// \param[in]  dataset Dataset.
  /// \param[in]  index   Index of point in \p dataset.
  /// \param[in]  feature Indices/names and weights/rescaling parameters of features.
  static void GetPoint(double *point, vtkPointSet *dataset, int index, const FeatureList *feature = NULL);

  /// Calculate squared Euclidean distance between feature Arrays
  ///
  /// \param[in] a First  feature Array.
  /// \param[in] b Second feature Array.
  /// \param[in] d Dimension of feature Arrays.
  ///
  /// \returns Squared Euclidean distance, i.e., dot product of a and b.
  static double Distance2BetweenPoints(const double *a, const double *b, int d = 3);

  // ---------------------------------------------------------------------------
  // Attributes

  /// Dataset for which search structure is build
  mirtkPublicAggregateMacro(vtkPointSet, DataSet);

  /// Indices of points to consider only or NULL
  mirtkPublicAggregateMacro(const Array<int>, Sample);

  /// Indices/names and rescaling parameters of point data arrays
  mirtkPublicAttributeMacro(FeatureList, Features);

  /// Return global _DataSet point indices instead of _Sample indices
  mirtkPublicAttributeMacro(bool, GlobalIndices);

  /// Number of points in search structure
  mirtkReadOnlyAttributeMacro(int, NumberOfPoints);

  /// Dimension of feature Arrays/points
  mirtkReadOnlyAttributeMacro(int, PointDimension);

  /// VTK point locator used for three-dimensional feature spaces
  vtkSmartPointer<vtkOctreePointLocator> _VtkLocator;

  /// FLANN point locator used for higher-dimensional feature spaces when available
  SharedPtr<FlannPointLocator> _FlannLocator;

  // ---------------------------------------------------------------------------
  // Construction/Destruction

protected:

  /// Constructor
  PointLocator();

  /// Initialize point locator
  void Initialize();

private:

  /// Copy constructor
  /// \note Intentionally not implemented
  PointLocator(const PointLocator &);

  /// Assignment operator
  /// \note Intentionally not implemented
  void operator =(const PointLocator &);

public:

  /// Construct new point locator for search in given dataset
  ///
  /// \param[in] dataset Dataset in which points are searched.
  /// \param[in] sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in] feature Indices and weights of point data in \p dataset to use.
  static PointLocator *New(vtkPointSet       *dataset,
                           const Array<int>  *sample  = NULL,
                           const FeatureList *feature = NULL);

  /// Destructor
  virtual ~PointLocator();

  // ---------------------------------------------------------------------------
  // Closest point

  /// Find closest point
  ///
  /// \param[in]  point Query point.
  /// \param[out] dist2 Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the query \p point.
  int FindClosestPoint(double *point, double *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the i-th point in \p dataset.
  int FindClosestPoint(vtkPointSet       *dataset,
                       const Array<int>  *sample,
                       int                index,
                       const FeatureList *features,
                       double            *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset  Dataset whose nearest neighbor is queried.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the i-th point in \p dataset.
  int FindClosestPoint(vtkPointSet       *dataset,
                       int                index,
                       const FeatureList *features,
                       double            *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset Dataset whose nearest neighbor is queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  /// \param[out] dist2   Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the i-th point in \p dataset.
  int FindClosestPoint(vtkPointSet      *dataset,
                       const Array<int> *sample,
                       int               index,
                       double           *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset Dataset whose nearest neighbor is queried.
  /// \param[in]  index   Index of point in \p dataset.
  /// \param[out] dist2   Squared Euclidean distance of closest point.
  ///
  /// \returns Index of point/sample in dataset closest to the i-th point in \p dataset.
  int FindClosestPoint(vtkPointSet *dataset,
                       int          index,
                       double      *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples closest to each point in \p dataset.
  Array<int> FindClosestPoint(vtkPointSet       *dataset,
                              const Array<int>  *sample,
                              const FeatureList *features,
                              Array<double>     *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset  Dataset whose nearest neighbor is queried.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples closest to each point in \p dataset.
  Array<int> FindClosestPoint(vtkPointSet       *dataset,
                              const FeatureList *features,
                              Array<double>     *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset Dataset whose nearest neighbor is queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[out] dist2   Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples closest to each point in \p dataset.
  Array<int> FindClosestPoint(vtkPointSet      *dataset,
                              const Array<int> *sample,
                              Array<double>    *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset Dataset whose nearest neighbor is queried.
  /// \param[out] dist2   Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples closest to each point in \p dataset.
  Array<int> FindClosestPoint(vtkPointSet   *dataset,
                              Array<double> *dist2 = NULL);

  /// Find closest point
  ///
  /// \param[in]  dataset1  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample1   Indices of points in \p dataset1 to consider only or NULL for all.
  /// \param[in]  features1 Indices and weights of point data in \p dataset1 to use.
  /// \param[in]  dataset2  Dataset in which to perform nearest neighbor search.
  /// \param[in]  sample2   Indices of points in \p dataset2 to consider only or NULL for all.
  /// \param[in]  features2 Indices and weights of point data in \p dataset2 to use.
  /// \param[out] dist2     Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points/samples in \p dataset2 which are closest to each point in \p dataset1.
  static Array<int> FindClosestPoint(vtkPointSet       *dataset1,
                                     const Array<int>  *sample1,
                                     const FeatureList *features1,
                                     vtkPointSet       *dataset2,
                                     const Array<int>  *sample2,
                                     const FeatureList *features2,
                                     Array<double>     *dist2 = NULL);

  // ---------------------------------------------------------------------------
  // Nearest neighbors (kNN)

  /// Find k nearest neighbors
  ///
  /// \param[in]  k     Number of nearest neighbors to find.
  /// \param[in]  point Query point.
  /// \param[out] dist2 Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindClosestNPoints(int k, double *point, Array<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindClosestNPoints(int k, vtkPointSet       *dataset,
                                       const Array<int>  *sample,
                                       int                index,
                                       const FeatureList *features,
                                Array<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindClosestNPoints(int k, vtkPointSet      *dataset,
                                       const Array<int> *sample,
                                       int               index,
                                Array<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  index    Index of point in \p dataset.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindClosestNPoints(int k, vtkPointSet       *dataset,
                                       int                index,
                                       const FeatureList *features,
                                Array<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k       Number of nearest neighbors to find.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  index   Index of point in \p dataset.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindClosestNPoints(int k, vtkPointSet *dataset, int index,
                                Array<double> *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of closest \p k points.
  Array<Array<int> > FindClosestNPoints(int k, vtkPointSet       *dataset,
                                               const Array<int>  *sample,
                                               const FeatureList *features,
                                        Array<Array<double> > *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of closest \p k points.
  Array<Array<int> > FindClosestNPoints(int k, vtkPointSet      *dataset,
                                               const Array<int> *sample,
                                        Array<Array<double> > *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k        Number of nearest neighbors to find.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of closest \p k points.
  Array<Array<int> > FindClosestNPoints(int k, vtkPointSet       *dataset,
                                               const FeatureList *features,
                                        Array<Array<double> > *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k       Number of nearest neighbors to find.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of closest \p k points.
  Array<Array<int> > FindClosestNPoints(int k, vtkPointSet *dataset,
                                        Array<Array<double> > *dist2 = NULL);

  /// Find k nearest neighbors
  ///
  /// \param[in]  k         Number of nearest neighbors to find.
  /// \param[in]  dataset1  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample1   Indices of points in \p dataset1 to consider only or NULL for all.
  /// \param[in]  features1 Indices and weights of point data in \p dataset1 to use.
  /// \param[in]  dataset2  Dataset in which to perform nearest neighbor search.
  /// \param[in]  sample2   Indices of points in \p dataset2 to consider only or NULL for all.
  /// \param[in]  features2 Indices and weights of point data in \p dataset2 to use.
  /// \param[out] dist2     Squared Euclidean distance of closest point.
  ///
  /// \returns For each point in \p dataset1, indices of closest \p k points in \p dataset2.
  static Array<Array<int> > FindClosestNPoints(int                    k,
                                               vtkPointSet           *dataset1,
                                               const Array<int>      *sample1,
                                               const FeatureList     *features1,
                                               vtkPointSet           *dataset2,
                                               const Array<int>      *sample2,
                                               const FeatureList     *features2,
                                               Array<Array<double> > *dist2 = NULL);

  // ---------------------------------------------------------------------------
  // Radius search

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  point   Query point.
  /// \param[out] dist2   Squared Euclidean distance of closest point.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindPointsWithinRadius(double radius, double *point, Array<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius   Search radius in N-D point feature space.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index    Index of point in \p dataset or \p sample (if not NULL).
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                   const Array<int>  *sample,
                                                   int                index,
                                                   const FeatureList *features,
                                    Array<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  index   Index of point in \p dataset or \p sample (if not NULL).
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindPointsWithinRadius(double radius, vtkPointSet      *dataset,
                                                   const Array<int> *sample,
                                                   int               index,
                                    Array<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius   Search radius in N-D point feature space.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  index    Index of point in \p dataset.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                   int                index,
                                                   const FeatureList *features,
                                    Array<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  index   Index of point in \p dataset.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns Indices of points in dataset closests to the i-th point in \p d.
  Array<int> FindPointsWithinRadius(double radius, vtkPointSet *dataset, int index,
                                    Array<double> *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius   Search radius in N-D point feature space.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  sample   Indices of points in \p dataset to consider only or NULL for all.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of points within search \p radius.
  Array<Array<int> > FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                           const Array<int>  *sample,
                                                           const FeatureList *features,
                                            Array<Array<double> > *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[in]  sample  Indices of points in \p dataset to consider only or NULL for all.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of points within search \p radius.
  Array<Array<int> > FindPointsWithinRadius(double radius, vtkPointSet      *dataset,
                                                           const Array<int> *sample,
                                            Array<Array<double> > *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius   Search radius in N-D point feature space.
  /// \param[in]  dataset  Dataset whose nearest neighbors are queried.
  /// \param[in]  features Indices and weights of point data in \p dataset to use.
  /// \param[out] dist2    Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of points within search \p radius.
  Array<Array<int> > FindPointsWithinRadius(double radius, vtkPointSet       *dataset,
                                                           const FeatureList *features,
                                            Array<Array<double> > *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius  Search radius in N-D point feature space.
  /// \param[in]  dataset Dataset whose nearest neighbors are queried.
  /// \param[out] dist2   Squared Euclidean distance of found nearest neighbors.
  ///
  /// \returns For each point in \p dataset, indices of points within search \p radius.
  Array<Array<int> > FindPointsWithinRadius(double radius, vtkPointSet *dataset,
                                            Array<Array<double> > *dist2 = NULL);

  /// Find neighbors within search radius
  ///
  /// \param[in]  radius    Search radius in N-D point feature space.
  /// \param[in]  dataset1  Dataset whose nearest neighbor is queried.
  /// \param[in]  sample1   Indices of points in \p dataset1 to consider only or NULL for all.
  /// \param[in]  features1 Indices and weights of point data in \p dataset1 to use.
  /// \param[in]  dataset2  Dataset in which to perform nearest neighbor search.
  /// \param[in]  sample2   Indices of points in \p dataset2 to consider only or NULL for all.
  /// \param[in]  features2 Indices and weights of point data in \p dataset2 to use.
  /// \param[out] dist2     Squared Euclidean distance of closest point.
  ///
  /// \returns For each point in \p dataset1, indices of points in \p dataset2 within search \p radius.
  static Array<Array<int> > FindPointsWithinRadius(double                 radius,
                                                   vtkPointSet           *dataset1,
                                                   const Array<int>      *sample1,
                                                   const FeatureList     *features1,
                                                   vtkPointSet           *dataset2,
                                                   const Array<int>      *sample2,
                                                   const FeatureList     *features2,
                                                   Array<Array<double> > *dist2 = NULL);

};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// =============================================================================
// Static helpers
// =============================================================================

// -----------------------------------------------------------------------------
inline vtkDataArray *PointLocator::GetDataArray(vtkPointSet *dataset, const FeatureInfo &feature)
{
  vtkPointData * const pd = dataset->GetPointData();
  if (feature._Index >= 0) return pd->GetArray(feature._Index);
  return pd->GetArray(feature._Name.c_str());
}

// -----------------------------------------------------------------------------
inline int PointLocator::GetNumberOfPoints(vtkPointSet *dataset, const Array<int> *sample)
{
  return (sample && sample->size() > 0 ? static_cast<int>(sample->size()) : dataset->GetNumberOfPoints());
}

// -----------------------------------------------------------------------------
inline int PointLocator::GetPointDimension(vtkPointSet *dataset, const FeatureList *features)
{
  int dim;
  if (features && features->size() > 0) {
    dim = 0;
    vtkDataArray *feature_array;
    for (FeatureList::const_iterator feature = features->begin(); feature != features->end(); ++feature) {
      if (feature->_Weight != .0) {
        if (feature->_Index == -1) {
          dim += 3;
        } else {
          feature_array = GetDataArray(dataset, *feature);
          dim += feature_array->GetNumberOfComponents();
        }
      }
    }
  } else {
    dim = 3;
  }
  return dim;
}

// -----------------------------------------------------------------------------
inline int PointLocator::GetPointIndex(vtkPointSet *dataset, const Array<int> *sample, int index)
{
  return (sample && sample->size() > 0 ? (*sample)[index] : index);
}

// -----------------------------------------------------------------------------
inline void PointLocator
::GetPoint(Point &point, vtkPointSet *dataset, const Array<int> *sample, int index)
{
  double pt[3];
  const int i = GetPointIndex(dataset, sample, index);
  dataset->GetPoint(i, pt);
  point._x = pt[0], point._y = pt[1], point._z = pt[2];
}

// -----------------------------------------------------------------------------
inline void PointLocator
::GetPoint(double *point, vtkPointSet *dataset, const Array<int> *sample, int index, const FeatureList *features)
{
  const int i = GetPointIndex(dataset, sample, index);
  if (features && features->size() > 0) {
    vtkDataArray *feature_array;
    for (FeatureList::const_iterator feature = features->begin(); feature != features->end(); ++feature) {
      if (feature->_Weight != .0) {
        if (feature->_Index == -1) {
          dataset->GetPoint(i, point);
          for (int d = 0; d < 3; ++d, ++point) {
            (*point) = feature->_Weight * (feature->_Slope * (*point) + feature->_Intercept);
          }
        } else {
          feature_array = GetDataArray(dataset, *feature);
          feature_array->GetTuple(i, point);
          for (int d = 0; d < feature_array->GetNumberOfComponents(); ++d, ++point) {
            (*point) = feature->_Weight * (feature->_Slope * (*point) + feature->_Intercept);
          }
        }
      }
    }
  } else {
    dataset->GetPoint(i, point);
  }
}

// -----------------------------------------------------------------------------
inline void PointLocator
::GetPoint(double *point, vtkPointSet *dataset, int index, const FeatureList *features)
{
  GetPoint(point, dataset, NULL, index, features);
}

// -----------------------------------------------------------------------------
inline double PointLocator
::Distance2BetweenPoints(const double *a, const double *b, int d)
{
  double dx, dist2 = .0;
  for (int i = 0; i < d; ++i) {
    dx = b[i] - a[i];
    dist2 += dx * dx;
  }
  return dist2;
}

// =============================================================================
// Closest point
// =============================================================================

// -----------------------------------------------------------------------------
inline int PointLocator
::FindClosestPoint(vtkPointSet *dataset, const Array<int> *sample, int index,
                   const FeatureList *features, double *dist2)
{
  double *point = new double[_PointDimension];
  GetPoint(point, dataset, sample, index, features);
  int i = this->FindClosestPoint(point, dist2);
  delete[] point;
  return i;
}

// -----------------------------------------------------------------------------
inline int PointLocator
::FindClosestPoint(vtkPointSet *dataset, int index,
                   const FeatureList *features, double *dist2)
{
  return FindClosestPoint(dataset, NULL, index, features, dist2);
}

// -----------------------------------------------------------------------------
inline int PointLocator
::FindClosestPoint(vtkPointSet *dataset, const Array<int> *sample,
                   int index, double *dist2)
{
  return FindClosestPoint(dataset, sample, index, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline int PointLocator
::FindClosestPoint(vtkPointSet *dataset, int index, double *dist2)
{
  return FindClosestPoint(dataset, NULL, index, dist2);
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindClosestPoint(vtkPointSet *dataset, const FeatureList *features, Array<double> *dist2)
{
  return FindClosestPoint(dataset, NULL, features, dist2);
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindClosestPoint(vtkPointSet *dataset, const Array<int> *sample, Array<double> *dist2)
{
  return FindClosestPoint(dataset, sample, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindClosestPoint(vtkPointSet *dataset, Array<double> *dist2)
{
  return FindClosestPoint(dataset, NULL, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindClosestPoint(vtkPointSet *dataset1, const Array<int> *sample1, const FeatureList *features1,
                   vtkPointSet *dataset2, const Array<int> *sample2, const FeatureList *features2,
                   Array<double> *dist2)
{
  UniquePtr<PointLocator> locator(PointLocator::New(dataset2, sample2, features2));
  return locator->FindClosestPoint(dataset1, sample1, features1, dist2);
}

// =============================================================================
// Nearest neighbors (kNN)
// =============================================================================

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const Array<int> *sample, int index,
                     const FeatureList *features, Array<double> *dist2)
{
  double *point = new double[_PointDimension];
  GetPoint(point, dataset, sample, index, features);
  Array<int> indices = this->FindClosestNPoints(k, point, dist2);
  delete[] point;
  return indices;
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, int index,
                     const FeatureList *features, Array<double> *dist2)
{
  return FindClosestNPoints(k, dataset, NULL, index, features, dist2);
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const Array<int> *sample,
                     int index, Array<double> *dist2)
{
  return FindClosestNPoints(k, dataset, sample, index, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, int index, Array<double> *dist2)
{
  return FindClosestNPoints(k, dataset, NULL, index, dist2);
}

// -----------------------------------------------------------------------------
inline Array<Array<int> > PointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const FeatureList *features, Array<Array<double> > *dist2)
{
  return FindClosestNPoints(k, dataset, NULL, features, dist2);
}

// -----------------------------------------------------------------------------
inline Array<Array<int> > PointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, const Array<int> *sample, Array<Array<double> > *dist2)
{
  return FindClosestNPoints(k, dataset, sample, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline Array<Array<int> > PointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset, Array<Array<double> > *dist2)
{
  return FindClosestNPoints(k, dataset, NULL, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline Array<Array<int> > PointLocator
::FindClosestNPoints(int k, vtkPointSet *dataset1, const Array<int> *sample1, const FeatureList *features1,
                            vtkPointSet *dataset2, const Array<int> *sample2, const FeatureList *features2,
                     Array<Array<double> > *dist2)
{
  UniquePtr<PointLocator> locator(PointLocator::New(dataset2, sample2, features2));
  return locator->FindClosestNPoints(k, dataset1, sample1, features1, dist2);
}

// =============================================================================
// Radius search
// =============================================================================

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, const Array<int> *sample,
                         int index, const FeatureList *features, Array<double> *dist2)
{
  double *point = new double[_PointDimension];
  GetPoint(point, dataset, sample, index, features);
  Array<int> indices = this->FindPointsWithinRadius(radius, point, dist2);
  delete[] point;
  return indices;
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, int index,
                         const FeatureList *features, Array<double> *dist2)
{
  return FindPointsWithinRadius(radius, dataset, NULL, index, features, dist2);
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, const Array<int> *sample,
                         int index, Array<double> *dist2)
{
  return FindPointsWithinRadius(radius, dataset, sample, index, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline Array<int> PointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, int index, Array<double> *dist2)
{
  return FindPointsWithinRadius(radius, dataset, NULL, index, dist2);
}

// -----------------------------------------------------------------------------
inline Array<Array<int> > PointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset,
                         const FeatureList *features, Array<Array<double> > *dist2)
{
  return FindPointsWithinRadius(radius, dataset, NULL, features, dist2);
}

// -----------------------------------------------------------------------------
inline Array<Array<int> > PointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, const Array<int> *sample, Array<Array<double> > *dist2)
{
  return FindPointsWithinRadius(radius, dataset, sample, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline Array<Array<int> > PointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset, Array<Array<double> > *dist2)
{
  return FindPointsWithinRadius(radius, dataset, NULL, NULL, dist2);
}

// -----------------------------------------------------------------------------
inline Array<Array<int> > PointLocator
::FindPointsWithinRadius(double radius, vtkPointSet *dataset1, const Array<int> *sample1, const FeatureList *features1,
                                        vtkPointSet *dataset2, const Array<int> *sample2, const FeatureList *features2,
                         Array<Array<double> > *dist2)
{
  UniquePtr<PointLocator> locator(PointLocator::New(dataset2, sample2, features2));
  return locator->FindPointsWithinRadius(radius, dataset1, sample1, features1, dist2);
}


} // namespace mirtk

#endif // MIRTK_PointLocator_H
