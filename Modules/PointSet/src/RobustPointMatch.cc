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

#include "mirtk/RobustPointMatch.h"

#include "mirtk/Math.h"
#include "mirtk/Memory.h"
#include "mirtk/Pair.h"
#include "mirtk/Array.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/PointLocator.h"
#include "mirtk/SparseMatrix.h"
#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkIdList.h"
#include "vtkPointSet.h"
#include "vtkOctreePointLocator.h"


namespace mirtk {


// =============================================================================
// Utilities
// =============================================================================

namespace RobustPointMatchUtils {


// -----------------------------------------------------------------------------
class SquaredDistance
{
  vtkSmartPointer<vtkPointSet>           _Target;
  const Array<int>                      *_Sample;
  vtkSmartPointer<vtkOctreePointLocator> _Source;
  int                                    _Num;
  double                                 _Max;
  double                                 _Sum, _Sum2;
  int                                    _Cnt;

public:

  SquaredDistance()
  :
    _Max(.0), _Sum(.0), _Sum2(.0), _Cnt(0)
  {}

  SquaredDistance(const SquaredDistance &lhs, split)
  :
    _Target(lhs._Target), _Sample(lhs._Sample), _Source(lhs._Source),
    _Num(lhs._Num), _Max(.0), _Sum(.0), _Sum2(.0), _Cnt(0)
  {}

  void join(const SquaredDistance &rhs)
  {
    _Max   = max(_Max, rhs._Max);
    _Sum  += rhs._Sum;
    _Sum2 += rhs._Sum2;
    _Cnt  += rhs._Cnt;
  }

  void operator()(const blocked_range<int> &re)
  {
    double p1[3], p2[3], dist2 = .0;
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    for (int r = re.begin(); r != re.end(); ++r) {
      _Target->GetPoint(PointCorrespondence::GetPointIndex(_Target, _Sample, r), p1);
      _Source->FindClosestNPoints(_Num, p1, ids);
      for (vtkIdType i = 0; i < ids->GetNumberOfIds(); ++i) {
        _Source->GetDataSet()->GetPoint(ids->GetId(i), p2);
        dist2 = vtkMath::Distance2BetweenPoints(p1, p2);
        _Sum  += dist2;
        _Sum2 += dist2 * dist2;
      }
      _Cnt += static_cast<int>(ids->GetNumberOfIds());
      // Note: ids are sorted from closest to farthest
      if (dist2 > _Max) _Max = dist2;
    }
  }

  void Add(vtkPointSet      *target,
           const Array<int> *sample,
           vtkPointSet      *source,
           int               num = 0)
  {
    const int m = PointLocator::GetNumberOfPoints(target, sample);
    const int n = source->GetNumberOfPoints();
    if (m == 0 || n == 0) return;
    if (num <= 0) num = n;
    _Num    = num;
    _Target = target;
    _Sample = sample;
    _Source = vtkSmartPointer<vtkOctreePointLocator>::New();
    _Source->SetDataSet(source);
    _Source->BuildLocator();
    blocked_range<int> range(0, m);
    parallel_reduce(range, *this);
    _Target = NULL;
    _Source = NULL;
  }

  void Add(vtkPointSet *target, const Array<int> *sample, vtkPointSet *source, double source_frac, int min_num  = 1)
  {
    int num = max(min_num, iround(source_frac * source->GetNumberOfPoints()));
    Add(target, sample, source, num);
  }

  double Mean()  { return _Cnt > 0 ? _Sum / _Cnt : .0; }
  double Sigma() { return _Cnt > 0 ? sqrt(_Sum2 / _Cnt - pow(_Sum / _Cnt, 2)) : .0; }
  double Max()   { return _Max; }
};

// -----------------------------------------------------------------------------
/// Calculate fuzzy correspondence weights
class CalculateCorrespondenceWeights
{
  typedef RobustPointMatch::WeightMatrix WeightMatrix;
  typedef RobustPointMatch::FeatureList  FeatureList;

  const RegisteredPointSet *_Target;
  const Array<int>         *_TargetSample;
  const FeatureList        *_TargetFeatures;
  const RegisteredPointSet *_Source;
  const Array<int>         *_SourceSample;
  const FeatureList        *_SourceFeatures;
  int                       _NumberOfPoints;
  int                       _NumberOfFeatures;
  double                    _MaxDist;
  double                    _Temperature;
  double                    _VarianceOfFeatures;
  WeightMatrix::Entries    *_CorrWeights;
  vtkAbstractPointLocator  *_Locator;

  CalculateCorrespondenceWeights() {}

public:

  CalculateCorrespondenceWeights(const CalculateCorrespondenceWeights &lhs)
  :
    _Target            (lhs._Target),
    _TargetSample      (lhs._TargetSample),
    _TargetFeatures    (lhs._TargetFeatures),
    _Source            (lhs._Source),
    _SourceSample      (lhs._SourceSample),
    _SourceFeatures    (lhs._SourceFeatures),
    _NumberOfPoints    (lhs._NumberOfPoints),
    _NumberOfFeatures  (lhs._NumberOfFeatures),
    _MaxDist           (lhs._MaxDist),
    _Temperature       (lhs._Temperature),
    _VarianceOfFeatures(lhs._VarianceOfFeatures),
    _CorrWeights       (lhs._CorrWeights),
    _Locator           (lhs._Locator)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    WeightMatrix::Entries::iterator weight;
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    double *p1 = Allocate<double>(_NumberOfFeatures); // spatial + optional extra features
    double *p2 = Allocate<double>(_NumberOfFeatures);
    for (int i = re.begin(); i != re.end(); ++i) {
      PointCorrespondence::GetPoint(p1, _Target, _TargetSample, i, _TargetFeatures);
      _Locator->FindPointsWithinRadius(_MaxDist, p1, ids);
      if (ids->GetNumberOfIds() > 0) {
        _CorrWeights[i].resize(_CorrWeights[i].size() + ids->GetNumberOfIds());
        weight = _CorrWeights[i].end() - ids->GetNumberOfIds();
        for (vtkIdType j = 0; j < ids->GetNumberOfIds(); ++j, ++weight) {
          PointCorrespondence::GetPoint(p2, _Source, _SourceSample, ids->GetId(j), _SourceFeatures);
          weight->first  = PointCorrespondence::GetPointIndex(_Source, _SourceSample, ids->GetId(j));
          weight->second = exp(- vtkMath::Distance2BetweenPoints(p1, p2) / _Temperature);
          if (_NumberOfFeatures > 3) {
            weight->second *= exp(- PointCorrespondence::Distance2BetweenPoints(p1+3, p2+3, _NumberOfFeatures-3) / _VarianceOfFeatures);
          }
        }
      }
    }
    Deallocate(p1);
    Deallocate(p2);
  }

  static void Run(const RegisteredPointSet    *target,
                  const Array<int>            *target_sample,
                  const FeatureList           *target_features,
                  const RegisteredPointSet    *source,
                  const Array<int>            *source_sample,
                  const FeatureList           *source_features,
                  WeightMatrix::Entries       *corrw,
                  WeightMatrix::StorageLayout  layout,
                  double                       temperature,
                  double                       var_of_features,
                  double                       min_weight)
  {
    if (min_weight >= 1.0) {
      cerr << "RobustPointMatch::CalculateWeights: Minimum weight must be less than 1" << endl;
      exit(1);
    }
    if (min_weight <= .0) {
      cerr << "RobustPointMatch::CalculateWeights: Minimum weight must be positive" << endl;
      exit(1);
    }
    if (layout == WeightMatrix::CCS) {
      swap(target,          source);
      swap(target_sample,   source_sample);
      swap(target_features, source_features);
    }
    const int m = PointCorrespondence::GetNumberOfPoints(target, target_sample);
    const int n = PointCorrespondence::GetNumberOfPoints(source, source_sample);
    if (m == 0 || n == 0) return;
    MIRTK_START_TIMING();
    CalculateCorrespondenceWeights body;
    body._Target             = target;
    body._TargetSample       = target_sample;
    body._TargetFeatures     = target_features;
    body._Source             = source;
    body._SourceSample       = source_sample;
    body._SourceFeatures     = source_features;
    body._NumberOfPoints     = n;
    body._NumberOfFeatures   = PointCorrespondence::GetPointDimension(target, target_features);
    body._MaxDist            = sqrt(temperature * -log(min_weight));
    body._Temperature        = temperature;
    body._VarianceOfFeatures = var_of_features;
    body._CorrWeights        = corrw;
    vtkSmartPointer<vtkOctreePointLocator> locator;
    vtkSmartPointer<vtkPointSet>           dataset;
    dataset = PointCorrespondence::GetPointSet(source, source_sample);
    locator = vtkSmartPointer<vtkOctreePointLocator>::New();
    locator->SetDataSet(dataset);
    locator->BuildLocator();
    body._Locator = locator;
    parallel_for(blocked_range<int>(0, m), body);
    for (int i = 0; i < m; ++i) sort(corrw[i].begin(), corrw[i].end());
    MIRTK_DEBUG_TIMING(7, "calculating weight for each pair of points");
  }
};

// -----------------------------------------------------------------------------
class GetCentroidOfPoints
{
  vtkPointSet      *_DataSet;
  const Array<int> *_Sample;
  Point             _Sum;

  GetCentroidOfPoints() : _Sum(.0,.0,.0) {}

public:

  GetCentroidOfPoints(const GetCentroidOfPoints &lhs, split)
  :
    _DataSet(lhs._DataSet),
    _Sample (lhs._Sample),
    _Sum    (lhs._Sum)
  {}

  void join(const GetCentroidOfPoints &rhs)
  {
    _Sum += rhs._Sum;
  }

  ~GetCentroidOfPoints()
  {
  }

  void operator ()(const blocked_range<int> &re)
  {
    double p[3];
    for (int i = re.begin(); i != re.end(); ++i) {
      _DataSet->GetPoint(PointCorrespondence::GetPointIndex(_DataSet, _Sample, i), p);
      _Sum += Point(p);
    }
  }

  static void Run(vtkPointSet      *dataset,
                  const Array<int> *sample,
                  Point            &centroid)
  {
    const int n = PointCorrespondence::GetNumberOfPoints(dataset, sample);
    if (n == 0) return;
    MIRTK_START_TIMING();
    GetCentroidOfPoints body;
    body._DataSet = dataset;
    body._Sample  = sample;
    blocked_range<int> i(0, n);
    parallel_reduce(i, body);
    centroid = body._Sum / n;
    MIRTK_DEBUG_TIMING(7, "getting centroid of points");
  }
};

// -----------------------------------------------------------------------------
/// Calculate outlier weights
class CalculateOutlierWeights1
{
  typedef RobustPointMatch::WeightMatrix WeightMatrix;

  const RegisteredPointSet *_DataSet;
  const Array<int>         *_Sample;
  Point                     _Cluster;
  double                    _Temperature;
  WeightMatrix::Entries    *_CorrWeights;
  int                       _N;

  CalculateOutlierWeights1() {}

public:

  CalculateOutlierWeights1(const CalculateOutlierWeights1 &lhs)
  :
    _DataSet    (lhs._DataSet),
    _Sample     (lhs._Sample),
    _Cluster    (lhs._Cluster),
    _Temperature(lhs._Temperature),
    _CorrWeights(lhs._CorrWeights),
    _N          (lhs._N)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    double weight;
    Point p;
    for (int i = re.begin(); i != re.end(); ++i) {
      _DataSet->GetPoint(PointCorrespondence::GetPointIndex(_DataSet, _Sample, i), p);
      weight = exp(- p.SquaredDistance(_Cluster) / _Temperature);
      _CorrWeights[i].push_back(MakePair(_N, static_cast<WeightMatrix::EntryType>(weight)));
    }
  }

  static void Run(const RegisteredPointSet *dataset,
                  const Array<int>         *sample,
                  const Point              &cluster,
                  WeightMatrix::Entries    *corrw,
                  int                       n,
                  double                    temperature)
  {
    const int m = PointCorrespondence::GetNumberOfPoints(dataset, sample);
    if (m == 0) return;
    MIRTK_START_TIMING();
    CalculateOutlierWeights1 body;
    body._DataSet     = dataset;
    body._Sample      = sample;
    body._Cluster     = cluster;
    body._Temperature = temperature;
    body._CorrWeights = corrw;
    body._N           = n;
    blocked_range<int> i(0, m);
    parallel_for(i, body);
    MIRTK_DEBUG_TIMING(7, "calculating outlier weights 1");
  }
};

// -----------------------------------------------------------------------------
/// Calculate outlier weights
class CalculateOutlierWeights2
{
  typedef RobustPointMatch::WeightMatrix WeightMatrix;

  const RegisteredPointSet *_DataSet;
  const Array<int>         *_Sample;
  Point                     _Cluster;
  double                    _Temperature;
  WeightMatrix::Entries    *_CorrWeights;

  CalculateOutlierWeights2() {}

public:

  CalculateOutlierWeights2(const CalculateOutlierWeights2 &lhs)
  :
    _DataSet    (lhs._DataSet),
    _Sample     (lhs._Sample),
    _Cluster    (lhs._Cluster),
    _Temperature(lhs._Temperature),
    _CorrWeights(lhs._CorrWeights)
  {}

  void operator()(const blocked_range<int> &re) const
  {
    double weight;
    Point p;
    for (int i = re.begin(); i != re.end(); ++i) {
      _DataSet->GetPoint(PointCorrespondence::GetPointIndex(_DataSet, _Sample, i), p);
      weight = exp(- p.SquaredDistance(_Cluster) / _Temperature);
      (*_CorrWeights)[i] = MakePair(i, static_cast<WeightMatrix::EntryType>(weight));
    }
  }

  static void Run(const RegisteredPointSet *dataset,
                  const Array<int>         *sample,
                  const Point              &cluster,
                  WeightMatrix::Entries    &corrw,
                  double                    temperature)
  {
    const int n = PointCorrespondence::GetNumberOfPoints(dataset, sample);
    if (n == 0) return;
    MIRTK_START_TIMING();
    corrw.resize(n);
    CalculateOutlierWeights2 body;
    body._DataSet     = dataset;
    body._Sample      = sample;
    body._Cluster     = cluster;
    body._Temperature = temperature;
    body._CorrWeights = &corrw;
    blocked_range<int> i(0, n);
    parallel_for(i, body);
    MIRTK_DEBUG_TIMING(7, "calculating outlier weights 2");
  }
};


} // namespace RobustPointMatchUtils
using namespace RobustPointMatchUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
RobustPointMatch::RobustPointMatch()
:
  _InitialTemperature(numeric_limits<double>::quiet_NaN()),
  _AnnealingRate     (numeric_limits<double>::quiet_NaN()),
  _FinalTemperature  (numeric_limits<double>::quiet_NaN()),
  _Temperature       (numeric_limits<double>::quiet_NaN()),
  _VarianceOfFeatures(numeric_limits<double>::quiet_NaN())
{
}

// -----------------------------------------------------------------------------
RobustPointMatch::RobustPointMatch(const RegisteredPointSet *target,
                                           const RegisteredPointSet *source)
:
  _InitialTemperature(numeric_limits<double>::quiet_NaN()),
  _AnnealingRate     (numeric_limits<double>::quiet_NaN()),
  _FinalTemperature  (numeric_limits<double>::quiet_NaN()),
  _Temperature       (numeric_limits<double>::quiet_NaN()),
  _VarianceOfFeatures(numeric_limits<double>::quiet_NaN())
{
  Target(target);
  Source(source);
  Initialize();
}

// -----------------------------------------------------------------------------
RobustPointMatch::RobustPointMatch(const RobustPointMatch &other)
:
  FuzzyCorrespondence(other),
  _InitialTemperature  (other._InitialTemperature),
  _AnnealingRate       (other._AnnealingRate),
  _FinalTemperature    (other._FinalTemperature),
  _Temperature         (other._Temperature),
  _VarianceOfFeatures  (other._VarianceOfFeatures),
  _TargetOutlierCluster(other._TargetOutlierCluster),
  _SourceOutlierCluster(other._SourceOutlierCluster)
{
}

// -----------------------------------------------------------------------------
PointCorrespondence *RobustPointMatch::NewInstance() const
{
  return new RobustPointMatch(*this);
}

// -----------------------------------------------------------------------------
RobustPointMatch::~RobustPointMatch()
{
}

// -----------------------------------------------------------------------------
RobustPointMatch::TypeId RobustPointMatch::Type() const
{
  return TypeId::RobustPointMatch;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool RobustPointMatch::Set(const char *name, const char *value)
{
  // Initial temperature
  if (strcmp(name, "Initial temperature") == 0 || strcmp(name, "Temperature") == 0) {
    return FromString(value, _InitialTemperature);
  }
  // Annealing rate, negative value indicates kNN search
  if (strcmp(name, "Annealing rate") == 0) {
    return FromString(value, _AnnealingRate) && _AnnealingRate < 1.0;
  }
  // Final temperature
  if (strcmp(name, "Final temperature") == 0) {
    return FromString(value, _FinalTemperature);
  }
  // Variance of extra features (resp. of their differences)
  if (strcmp(name, "Variance of features") == 0) {
    return FromString(value, _VarianceOfFeatures);
  }
  return FuzzyCorrespondence::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList RobustPointMatch::Parameter() const
{
  ParameterList params = FuzzyCorrespondence::Parameter();
  Insert(params, "Initial temperature",  _InitialTemperature);
  Insert(params, "Annealing rate",       _AnnealingRate);
  Insert(params, "Final temperature",    _FinalTemperature);
  Insert(params, "Variance of features", _VarianceOfFeatures);
  return params;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void RobustPointMatch::Initialize()
{
  // Initialize base class
  FuzzyCorrespondence::Initialize();

  // Ensure that spatial coordinates are first components of feature vector
  size_t i;
  for (i = 0; i < _TargetFeatures.size(); ++i) {
    if (_TargetFeatures[i]._Index == -1) break;
  }
  if (i == _TargetFeatures.size()) {
    FeatureList::value_type info;
    info._Name      = "spatial coordinates";
    info._Index     = -1;
    info._Slope     = 1.0;
    info._Intercept = .0;
    _TargetFeatures.insert(_TargetFeatures.begin(), info);
  } else if (i != 0) {
    swap(_TargetFeatures[0], _TargetFeatures[i]);
    if (_TargetFeatures[0]._Slope == .0) _TargetFeatures[0]._Slope = 1.0;
  }

  for (i = 0; i < _SourceFeatures.size(); ++i) {
    if (_SourceFeatures[i]._Index == -1) break;
  }
  if (i == _SourceFeatures.size()) {
    FeatureList::value_type info;
    info._Name      = "spatial coordinates";
    info._Index     = -1;
    info._Slope     = 1.0;
    info._Intercept = .0;
    _SourceFeatures.insert(_SourceFeatures.begin(), info);
  } else if (i != 0) {
    swap(_SourceFeatures[0], _SourceFeatures[i]);
    if (_SourceFeatures[0]._Slope == .0) _SourceFeatures[0]._Slope = 1.0;
  }

  // Initialize annealing process
  this->InitializeAnnealing();

  // TODO: Determine variance
  if (IsNaN(_VarianceOfFeatures)) _VarianceOfFeatures = 1.0;
}

// -----------------------------------------------------------------------------
void RobustPointMatch::InitializeAnnealing()
{
  MIRTK_START_TIMING();

  // Annealing rate
  if (IsNaN(_AnnealingRate)) _AnnealingRate = 0.93;
  if (_AnnealingRate >= 1.0) {
    cerr << this->NameOfClass() << "::Initialize: ";
    cerr << "Annealing rate must be less than 1, normally it is in the range [0.9, 0.99]" << endl;
    cerr << "    Alternatively, set it to a negative integer to specify the number" << endl;
    cerr << "    of nearest neighbors to consider when choosing a new temperature." << endl;
    exit(1);
  }

  // Annealing rate used below to adjust temperature range
  double annealing_rate = (_AnnealingRate <= .0 ? .93 : _AnnealingRate);

  // Temperature range
  if (IsNaN(_InitialTemperature) || _InitialTemperature <= .0) {
    MIRTK_LOG_EVENT("Initializing annealing process...\n");

    // Number of nearest neighbors
    int n = 10;
    if (_InitialTemperature < .0) n = ceil(-_InitialTemperature);
    else if (_AnnealingRate < .0) n = ceil(-_AnnealingRate);
    const int k = min(n, min(_M, _N));
    MIRTK_LOG_EVENT("  Considering " << k << " nearest neighors\n");

    // Determine mean and standard deviation of distances to k nearest neighbors
    SquaredDistance dist2;
    vtkSmartPointer<vtkPointSet> source = GetPointSet(_Source, _SourceSample);
    dist2.Add(_Target->PointSet(), _TargetSample, source, k);
    MIRTK_LOG_EVENT("  Mean squared distance = " << dist2.Mean() << " (sigma = " << dist2.Sigma() << ")\n");

    // Set initial temperature of annealing process
    if (IsNaN(_InitialTemperature) || _InitialTemperature <= .0) {
      _InitialTemperature = dist2.Mean() + 1.5 * dist2.Sigma();
    }

    MIRTK_LOG_EVENT("Initializing annealing proces... done\n");
  }

  // Set final temperature of annealing process
  if (IsNaN(_FinalTemperature)) {
    _FinalTemperature = pow(annealing_rate, 50) * _InitialTemperature;
  }

  MIRTK_LOG_EVENT("Initial temperature = " << _InitialTemperature << "\n");
  MIRTK_LOG_EVENT("Final   temperature = " << _FinalTemperature   << "\n");

  // Set initial temperature
  _Temperature = _InitialTemperature;
  MIRTK_LOG_EVENT("Temperature = " << _InitialTemperature << "\n");

  MIRTK_DEBUG_TIMING(6, "initialization of annealing process");
}

// -----------------------------------------------------------------------------
bool RobustPointMatch::Upgrade()
{
  // Have base class check ratio of outliers
  if (!FuzzyCorrespondence::Upgrade()) return false;

  // Reduce temperature
  if (_AnnealingRate < .0) {

    // Number of nearest neighbors
    const int k = min(int(ceil(-_AnnealingRate)), min(_M, _N));

    // Compute statistics of (squared) distances
    SquaredDistance dist2;
    vtkSmartPointer<vtkPointSet> source = GetPointSet(_Source, _SourceSample);
    dist2.Add(_Target->PointSet(), _TargetSample, source, k);

    // Adjust temperature further to speed up annealing process
    _Temperature = min(0.96 * _Temperature, dist2.Mean());
  } else {
    _Temperature *= _AnnealingRate;
  }

  // Stop if final temperature reached
  if (_Temperature < _FinalTemperature) {
    MIRTK_LOG_EVENT("Final temperature reached\n");
    return false;
  }

  // Log current temperature
  MIRTK_LOG_EVENT("Temperature = " << _Temperature << "\n");
  return true;
}

// -----------------------------------------------------------------------------
void RobustPointMatch::CalculateWeights()
{
  MIRTK_START_TIMING();

  // Size of weight matrix
  const int m = _M + 1;
  const int n = _N + 1;

  // Allocate lists for non-zero weight entries
  const int nentries = (_Weight.Layout() == WeightMatrix::CRS ? m : n);
  WeightMatrix::Entries *entries = new WeightMatrix::Entries[nentries];

  // Calculate correspondence weights
  CalculateCorrespondenceWeights::Run(_Target, _TargetSample, &_TargetFeatures,
                                      _Source, _SourceSample, &_SourceFeatures,
                                      entries, _Weight.Layout(),
                                      _Temperature, _VarianceOfFeatures, _MinWeight);

  // Compute cluster to which source outliers are matched
  // (i.e., centroid of target points!)
  GetCentroidOfPoints::Run(_Target->PointSet(), _TargetSample, _SourceOutlierCluster);

  // Compute cluster to which target outliers are matched
  // (i.e., centroid of source points!)
  GetCentroidOfPoints::Run(_Source->PointSet(), _SourceSample, _TargetOutlierCluster);

  // Calculate weights of point assignment to outlier clusters
  if (_Weight.Layout() == WeightMatrix::CRS) {

    CalculateOutlierWeights1::Run(_Target, _TargetSample, _TargetOutlierCluster,
                                  entries, _N, _InitialTemperature);

    CalculateOutlierWeights2::Run(_Source, _SourceSample, _SourceOutlierCluster,
                                  entries[_M], _InitialTemperature);

  } else {

    CalculateOutlierWeights1::Run(_Source, _SourceSample, _SourceOutlierCluster,
                                  entries, _M, _InitialTemperature);

    CalculateOutlierWeights2::Run(_Target, _TargetSample, _TargetOutlierCluster,
                                  entries[_N], _InitialTemperature);

  }

  MIRTK_DEBUG_TIMING(6, "calculating correspondence weights (" << m << "x" << n << ")");

  // Initialize correspondence matrix
  MIRTK_RESET_TIMING();
  _Weight.Initialize(m, n, entries, true);
  delete[] entries;
  MIRTK_DEBUG_TIMING(7, "copying sparse matrix entries (NNZ=" << _Weight.NNZ() << ")");
}


} // namespace mirtk
