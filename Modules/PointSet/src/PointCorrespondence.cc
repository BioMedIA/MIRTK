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

#include "mirtk/PointCorrespondence.h"

#include "mirtk/Vtk.h"
#include "mirtk/Math.h"
#include "mirtk/Pair.h"
#include "mirtk/Array.h"
#include "mirtk/Algorithm.h"
#include "mirtk/Vector.h"
#include "mirtk/Matrix.h"
#include "mirtk/SparseMatrix.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/Transformation.h"
#include "mirtk/SpectralDecomposition.h"
#include "mirtk/PointSetIO.h"

#include "mirtk/CommonExport.h"

#include "vtkSmartPointer.h"
#include "vtkDataArray.h"
#include "vtkPointData.h"
#include "vtkOctreePointLocator.h"
#include "vtkPMaskPoints.h"
#include "vtkFloatArray.h"
#include "vtkCharArray.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkWindowedSincPolyDataFilter.h"

#include "mirtk/FiducialMatch.h"
#include "mirtk/ClosestPoint.h"
#include "mirtk/ClosestPointLabel.h"
#include "mirtk/ClosestCell.h"
#include "mirtk/SpectralMatch.h"
#include "mirtk/RobustClosestPoint.h"
#include "mirtk/RobustPointMatch.h"


namespace mirtk {


// Global debug level (cf. mirtk/Options.h)
MIRTK_Common_EXPORT extern int debug;


// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
PointCorrespondence *PointCorrespondence::New(TypeId type)
{
  switch (type) {
    case FiducialMatch:       return new class FiducialMatch();
    case ClosestPoint:        return new class ClosestPoint();
    case ClosestPointLabel:   return new class ClosestPointLabel();
    case ClosestCell:         return new class ClosestCell();
    case SpectralMatch:       return new class SpectralMatch();
    case RobustClosestPoint:  return new class RobustClosestPoint();
    case RobustPointMatch:    return new class RobustPointMatch();
    default:
      cerr << "PointCorrespondence::New: Unknown type = " << type << endl;
      exit(1);
  }
}

// -----------------------------------------------------------------------------
PointCorrespondence *PointCorrespondence::New(const char *type_name)
{
  TypeId type = Unknown;
  if (!FromString(type_name, type)) {
    cerr << "PointCorrespondence::New: Unknown type = " << type_name << endl;
    exit(1);
  }
  return New(type);
}

// -----------------------------------------------------------------------------
template <> bool FromString(const char *str, PointCorrespondence::TypeId &type)
{
  if (strcmp(str, "Index") == 0 || strcmp(str, "Fiducial") == 0) {
    type = PointCorrespondence::FiducialMatch;
  } else if (strcmp(str, "CP")            == 0 ||
             strcmp(str, "ClosestPoint")  == 0 ||
             strcmp(str, "Closest Point") == 0 ||
             strcmp(str, "Closest point") == 0) {
    type = PointCorrespondence::ClosestPoint;
  } else if (strcmp(str, "ClosestPointLabel")   == 0 ||
             strcmp(str, "Closest Point Label") == 0 ||
             strcmp(str, "Closest point label") == 0) {
    type = PointCorrespondence::ClosestPointLabel;
  } else if (strcmp(str, "CSP")                   == 0 ||
             strcmp(str, "ClosestCell")           == 0 ||
             strcmp(str, "Closest Cell")          == 0 ||
             strcmp(str, "Closest cell")          == 0 ||
             strcmp(str, "Closest Surface Point") == 0 ||
             strcmp(str, "Closest surface point") == 0 ||
             strcmp(str, "Closest")               == 0) {
    type = PointCorrespondence::ClosestCell;
  } else if (strcmp(str, "SM")                   == 0 ||
             strcmp(str, "SpectralMatch")   == 0 ||
             strcmp(str, "Spectral Match") == 0 ||
             strcmp(str, "Spectral match") == 0) {
    type = PointCorrespondence::SpectralMatch;
  } else if (strcmp(str, "RCP")                  == 0 ||
             strcmp(str, "RobustClosestPoint")   == 0 ||
             strcmp(str, "Robust Closest Point") == 0 ||
             strcmp(str, "Robust closest point") == 0) {
    type = PointCorrespondence::RobustClosestPoint;
  } else if (strcmp(str, "RPM")                == 0 ||
             strcmp(str, "RobustPointMatch")   == 0 ||
             strcmp(str, "Robust Point Match") == 0 ||
             strcmp(str, "Robust point match") == 0) {
    type = PointCorrespondence::RobustPointMatch;
  } else {
    type = PointCorrespondence::Unknown;
  }
  return (type != PointCorrespondence::Unknown);
}

// -----------------------------------------------------------------------------
template <> string ToString(const PointCorrespondence::TypeId &type, int w, char c, bool left)
{
  string str;
  switch (type) {
    case PointCorrespondence::FiducialMatch:      str = "Index";              break;
    case PointCorrespondence::ClosestPoint:       str = "ClosestPoint";       break;
    case PointCorrespondence::ClosestPointLabel:  str = "ClosestPointLabel";  break;
    case PointCorrespondence::ClosestCell:        str = "ClosestCell";        break;
    case PointCorrespondence::SpectralMatch:      str = "SpectralMatch";      break;
    case PointCorrespondence::RobustClosestPoint: str = "RobustClosestPoint"; break;
    case PointCorrespondence::RobustPointMatch:   str = "RobustPointMatch";   break;
    default:                                      str = "Unknown";            break;
  }
  return ToString(str, w, c, left);
}

// =============================================================================
// Auxiliary functions
// =============================================================================

namespace PointCorrespondenceUtils {


// -----------------------------------------------------------------------------
/// Auxiliary functor used by SamplePoints
class FindIndicesOfPointsClosestToSamples
{
  vtkPointSet             *_Samples;
  vtkAbstractPointLocator *_Locator;
  Array<int>              *_Index;

public:

  void operator()(const blocked_range<vtkIdType> &re) const
  {
    double pt[3];
    Array<int> &index = (*_Index);
    for (vtkIdType i = re.begin(); i != re.end(); ++i) {
      _Samples->GetPoint(i, pt);
      index[i] = _Locator->FindClosestPoint(pt);
    }
  }

  static void Run(vtkAbstractPointLocator *locator, vtkPointSet *samples, Array<int> &index)
  {
    index.resize(samples->GetNumberOfPoints());
    FindIndicesOfPointsClosestToSamples body;
    body._Samples = samples;
    body._Locator = locator;
    body._Index   = &index;
    blocked_range<vtkIdType> pts(0, samples->GetNumberOfPoints());
    parallel_for(pts, body);
    sort(index.begin(), index.end());
    index.erase(unique(index.begin(), index.end()), index.end());
  }
};

// -----------------------------------------------------------------------------
void SamplePoints(vtkPointSet *pointset, Array<int> &indices,
                  int maxnum, double maxdist, bool stratified)
{
  MIRTK_START_TIMING();

  // Reset set of drawn samples (empty --> use all points)
  indices.clear();

  // Skip if input data set contains no points (failsafe)
  if (pointset->GetNumberOfPoints() == 0) return;

  // Build point locator
  //
  // Due to a bug in vtkKdTreePointLocator, calling BuildLocator
  // is not sufficient to make FindClosestPoint thread-safe as it does
  // not call vtkBSPIntersections::BuildRegionsList
  // (cf. http://www.vtk.org/Bug/view.php?id=15206 ).
  vtkSmartPointer<vtkOctreePointLocator> locator;
  if (maxnum > 0 || maxdist > .0) {
    locator = vtkSmartPointer<vtkOctreePointLocator>::New();
    locator->SetDataSet(pointset);
    locator->BuildLocator();
  }

  // Determine maximum number of samples
  int nsamples;
  if (maxnum == 0 && maxdist > .0) {
    const double r = maxdist / 2.0;
    // Count number of points within radius r of each input point
    double p[3];
    vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
    double sum = .0;
    for (vtkIdType i = 0; i < pointset->GetNumberOfPoints(); ++i) {
      pointset->GetPoint(i, p);
      locator->FindPointsWithinRadius(r, p, ids);
      sum += ids->GetNumberOfIds();
    }
    // Divide number of points by average number of points within radius r
    nsamples = iround(pointset->GetNumberOfPoints() / (sum / pointset->GetNumberOfPoints()));
  } else {
    nsamples = maxnum;
  }

  // Extract specified maximum number of samples
  if (nsamples > 0) {
    vtkSmartPointer<vtkPMaskPoints> sampler;
    sampler = vtkSmartPointer<vtkPMaskPoints>::New();
    sampler->GenerateVerticesOff();
    sampler->SetMaximumNumberOfPoints(nsamples);
    sampler->SetRandomModeType(stratified ? 2 : 1);
    sampler->RandomModeOn();
    sampler->ProportionalMaximumNumberOfPointsOn();
    SetVTKInput(sampler, pointset);
    MIRTK_START_TIMING();
    sampler->Update();
    MIRTK_DEBUG_TIMING(6, "uniformly subsampling point set");
    MIRTK_RESET_TIMING();
    FindIndicesOfPointsClosestToSamples::Run(locator, sampler->GetOutput(), indices);
    MIRTK_DEBUG_TIMING(6, "finding closest sample points");
  }

  MIRTK_DEBUG_TIMING(5, "subsampling of point set (#samples = " << nsamples << ")");
}

// -----------------------------------------------------------------------------
/// Maximum range of spatial coordinates
double MaxSpatialRange(vtkPointSet *target, vtkPointSet *source)
{
  double target_range[6], source_range[6], range[3];
  target->GetBounds(target_range);
  source->GetBounds(source_range);
  range[0] = max(target_range[1] - target_range[0], source_range[1] - source_range[0]);
  range[1] = max(target_range[3] - target_range[2], source_range[3] - source_range[2]);
  range[2] = max(target_range[5] - target_range[4], source_range[5] - source_range[4]);
  return max(max(range[0], range[1]), range[2]);
}

// -----------------------------------------------------------------------------
vtkSmartPointer<vtkPointSet> SpectralPoints(vtkPointSet *input)
{
  vtkDataArray *eigenmodes = input->GetPointData()->GetArray("eigenmodes");
  if (!eigenmodes) return NULL;

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(input->GetNumberOfPoints());

  double *p = new double[max(3, eigenmodes->GetNumberOfComponents())];
  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
    eigenmodes->GetTuple(i, p);
    points->SetPoint(i, p);
  }
  delete[] p;

  vtkSmartPointer<vtkPointSet> output;
  output.TakeReference(input->NewInstance());
  output->DeepCopy (input);
  output->SetPoints(points);
  return output;
}

// -----------------------------------------------------------------------------
int ComputeNormals(vtkPolyData *dataset)
{
  // Smooth polydata
  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother;
  smoother = vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
  smoother->FeatureEdgeSmoothingOff();
  smoother->SetNumberOfIterations(25);
  smoother->SetPassBand(.1);
  smoother->NormalizeCoordinatesOn();
  SetVTKInput(smoother, dataset);
  // Calculate normals
  vtkSmartPointer<vtkPolyDataNormals> calculator;
  calculator = vtkSmartPointer<vtkPolyDataNormals>::New();
  calculator->SplittingOff();
  calculator->ConsistencyOn();
  calculator->AutoOrientNormalsOn();
  SetVTKConnection(calculator, smoother);
  calculator->Update();
  vtkDataArray *normals = calculator->GetOutput()->GetPointData()->GetNormals();
  normals->SetName("Normals");
  dataset->GetPointData()->SetNormals(normals);
  int i = -1;
  dataset->GetPointData()->GetArray("Normals", i);
  return i;
}

// -----------------------------------------------------------------------------
int ComputeSpectralNormals(vtkPolyData *dataset)
{
  vtkSmartPointer<vtkPointSet> spectral_pointset = SpectralPoints(dataset);
  vtkPolyData *spectral_polydata = vtkPolyData::SafeDownCast(spectral_pointset);
  int i = -1;
  if (spectral_polydata) {
    ComputeNormals(spectral_polydata);
    vtkDataArray *normals = dataset->GetPointData()->GetArray("spectral normals", i);
    if (normals) normals->DeepCopy(spectral_polydata->GetPointData()->GetNormals());
    else {
      i = dataset->GetPointData()->AddArray(spectral_polydata->GetPointData()->GetNormals());
    }
  }
  return i;
}

// -----------------------------------------------------------------------------
int GetEigenmodes(const RegisteredPointSet *dataset, int k, Matrix &m, Vector &v)
{
  int           i = -2;
  vtkDataArray *s = dataset->PointSet()->GetPointData()->GetArray("eigenmodes", i);
  if (dataset->Transformation() || s == NULL || s->GetNumberOfComponents() < k || v.Rows() < k) {
    vtkPolyData *polydata = vtkPolyData::SafeDownCast(dataset->PointSet());
    if (!polydata) {
      cerr << "PointCorrespondence: Can only compute eigenmodes of surface mesh" << endl;
      exit(1);
    }
    if (SpectralDecomposition::ComputeEigenmodes(polydata, k, m, v) != k) {
      cerr << "PointCorrespondence: Failed to compute " << k << " eigenmodes" << endl;
      exit(1);
    }
    return -2;
  } else {
    m.Initialize(static_cast<int>(s->GetNumberOfTuples()), k);
    double *row = new double[k];
    for (vtkIdType r = 0; r < s->GetNumberOfTuples(); ++r) {
      s->GetTuple(r, row);
      for (int c = 0; c < k; ++c) m(r, c) = row[c];
    }
    delete[] row;
    return i;
  }
}

// -----------------------------------------------------------------------------
void GetNormalizationParametersOfEigenmodes(double &slope, double &intercept,
                                            vtkDataArray *m1, vtkDataArray *m2,
                                            Vector w = Vector())
{
  const int k = min(m1->GetNumberOfComponents(), m2->GetNumberOfComponents());
  if (w.Rows() < k) w.Resize(k, 1.0);
  double range[2], minval = numeric_limits<double>::infinity(), maxrange = .0;
  for (int c = 0; c < k; ++c) {
    m1->GetRange(range, c);
    minval   = min(minval,   w(c) * range[0]);
    maxrange = max(maxrange, w(c) * (range[1] - range[0]));
    m2->GetRange(range, c);
    minval   = min(minval,   w(c) * range[0]);
    maxrange = max(maxrange, w(c) * (range[1] - range[0]));
  }
  slope     = 2.0 / maxrange;
  intercept = - slope * minval - 1.0;
}

// -----------------------------------------------------------------------------
/// (Re-)compute eigenmodes of datasets independent of each other and
/// correct for sign ambiguity and order afterwards using bipartite matching
void UpdateEigenmodes(const RegisteredPointSet         *target,
                      Vector                           &target_eigenvalues,
                      PointCorrespondence::FeatureList &target_features,
                      const RegisteredPointSet         *source,
                      Vector                           &source_eigenvalues,
                      PointCorrespondence::FeatureList &source_features,
                      int                               k)
{
  using namespace SpectralDecomposition;
  // Get feature entries
  PointCorrespondence::FeatureList::iterator target_info;
  for (target_info = target_features.begin(); target_info != target_features.end(); ++target_info) {
    if (target_info->_Name == "eigenmodes") break;
  }
  if (target_info == target_features.end()) return;
  PointCorrespondence::FeatureList::iterator source_info;
  for (source_info = source_features.begin(); source_info != source_features.end(); ++source_info) {
    if (source_info->_Name == "eigenmodes") break;
  }
  if (source_info == source_features.end()) return;
  // Cast to vtkPolyData
  vtkPolyData *target_polydata = vtkPolyData::SafeDownCast(target->PointSet());
  vtkPolyData *source_polydata = vtkPolyData::SafeDownCast(source->PointSet());
  if (!target_polydata || !source_polydata) {
    cerr << "FiducialRegistrationError::UpdateEigenmodes: Point sets must be surface meshes" << endl;
    exit(1);
  }
  // Get/compute eigenmodes of each dataset
  Matrix m1, m2;
  target_info->_Index = GetEigenmodes(target, k, m1, target_eigenvalues);
  source_info->_Index = GetEigenmodes(source, k, m2, source_eigenvalues);
  // Match sign and order of eigenmodes
  if (source->Transformation() != NULL && target->Transformation() == NULL) {
    MatchEigenmodes(target->Points(), m1, target_eigenvalues,
                    source->Points(), m2, source_eigenvalues);
    source_info->_Index = -2; // indicate that eigenmodes need to be updated
  } else {
    MatchEigenmodes(source->Points(), m2, source_eigenvalues,
                    target->Points(), m1, target_eigenvalues);
    target_info->_Index = -2; // indicate that eigenmodes need to be updated
  }
  // Set eigenmodes as point data of datasets to be registered
  if (target_info->_Index < 0) {
    target_info->_Index = SetEigenmodes(target_polydata, m1, 0, k, "eigenmodes");
  }
  if (source_info->_Index < 0) {
    source_info->_Index = SetEigenmodes(source_polydata, m2, 0, k, "eigenmodes");
  }
  // Set rescaling parameters s.t. eigenmodes are normalized to range [-1 +1]
  double slope, intercept;
  vtkDataArray *target_eigenmodes = target->PointSet()->GetPointData()->GetArray(target_info->_Index);
  vtkDataArray *source_eigenmodes = source->PointSet()->GetPointData()->GetArray(source_info->_Index);
  GetNormalizationParametersOfEigenmodes(slope, intercept, target_eigenmodes, source_eigenmodes);
  target_info->_Slope     = source_info->_Slope     = slope;
  target_info->_Intercept = source_info->_Intercept = intercept;
}

// -----------------------------------------------------------------------------
void UpdateEigenmodes(const RegisteredPointSet         *target,
                      PointCorrespondence::FeatureList &target_features,
                      const RegisteredPointSet         *source,
                      PointCorrespondence::FeatureList &source_features,
                      Vector                           &eigenvalues,
                      int                               k)
{
  using namespace SpectralDecomposition;
  typedef GenericSparseMatrix<double> SparseMatrix;
  vtkDataArray *target_eigenmodes, *source_eigenmodes;
  // Get feature entries
  PointCorrespondence::FeatureList::iterator target_info;
  for (target_info = target_features.begin(); target_info != target_features.end(); ++target_info) {
    if (target_info->_Name == "eigenmodes") break;
  }
  if (target_info == target_features.end()) return;
  PointCorrespondence::FeatureList::iterator source_info;
  for (source_info = source_features.begin(); source_info != source_features.end(); ++source_info) {
    if (source_info->_Name == "eigenmodes") break;
  }
  if (source_info == source_features.end()) return;
  // Cast to vtkPolyData
  vtkPolyData *target_polydata = vtkPolyData::SafeDownCast(target->PointSet());
  vtkPolyData *source_polydata = vtkPolyData::SafeDownCast(source->PointSet());
  if (!target_polydata || !source_polydata) {
    cerr << "FiducialRegistrationError::UpdateEigenmodes: Point sets must be surface meshes" << endl;
    exit(1);
  }
  // Compute initial spectral coordinates
  if (ComputeEigenmodes(target_polydata, source_polydata, k) < k) {
    cerr << "FiducialRegistrationError::UpdateEigenmodes: Failed to find " << k << " initial eigenmodes" << endl;
    exit(1);
  }
  // Sample points for which inter-dataset links will be added
  Array<int> target_sample, source_sample;
  SamplePoints(target, target_sample, max(10, target->NumberOfPoints() / 10));
  SamplePoints(source, source_sample, max(10, source->NumberOfPoints() / 10));
  // Find initial point correspondences
  PointLocator::FeatureList target_match_features;
  PointLocator::FeatureList source_match_features;
  PointLocator::FeatureInfo feature;
  vtkPointData * const targetPD = target_polydata->GetPointData();
  vtkPointData * const sourcePD = source_polydata->GetPointData();
  // Use computed spectral coordinates for initial match
  int target_index, source_index;
  target_eigenmodes = targetPD->GetArray(feature._Name.c_str(), target_index);
  source_eigenmodes = sourcePD->GetArray(feature._Name.c_str(), source_index);
  GetNormalizationParametersOfEigenmodes(feature._Slope, feature._Intercept,
                                         target_eigenmodes, source_eigenmodes);
  feature._Name   = "eigenmodes";
  feature._Weight = target_info->_Weight;
  feature._Index  = target_index;
  target_match_features.push_back(feature);
  feature._Weight = source_info->_Weight;
  feature._Index  = source_index;
  source_match_features.push_back(feature);
  // Optionally also use normals of 3D spectral points
  feature._Name = "spectral normals";
  for (size_t i = 0; i < target_features.size(); ++i) {
    if (target_features[i]._Name == feature._Name) {
      ComputeSpectralNormals(target_polydata);
      feature._Weight    = target_features[i]._Weight;
      feature._Slope     = target_features[i]._Slope;
      feature._Intercept = target_features[i]._Intercept;
      targetPD->GetArray(feature._Name.c_str(), feature._Index);
      target_match_features.push_back(feature);
      break;
    }
  }
  for (size_t i = 0; i < source_features.size(); ++i) {
    if (source_features[i]._Name == feature._Name) {
      ComputeSpectralNormals(source_polydata);
      sourcePD->GetArray(feature._Name.c_str(), feature._Index);
      feature._Weight    = source_features[i]._Weight;
      feature._Slope     = source_features[i]._Slope;
      feature._Intercept = source_features[i]._Intercept;
      source_match_features.push_back(feature);
      break;
    }
  }
  // Find nearest neighbors of point samples
  Array<int>    corr12, corr21;
  Array<double> dist12, dist21;
  corr12 = PointLocator::FindClosestPoint(target->PointSet(), &target_sample, &target_match_features,
                                          source->PointSet(), &source_sample, &source_match_features, &dist12);
  corr21 = PointLocator::FindClosestPoint(source->PointSet(), &source_sample, &source_match_features,
                                          target->PointSet(), &target_sample, &target_match_features, &dist21);
  // Set intra-mesh affinity weights
  const int m = target->NumberOfPoints();
  const int n = source->NumberOfPoints();
  SparseMatrix::Entries *cols = new SparseMatrix::Entries[m + n];
  AdjacencyMatrix(cols, SparseMatrix::CCS, 0, 0, vtkPolyData::SafeDownCast(target->InputPointSet()));
  AdjacencyMatrix(cols, SparseMatrix::CCS, m, m, vtkPolyData::SafeDownCast(source->InputPointSet()));
  // Set inter-mesh affinity weights
  const int ncorr12 = static_cast<int>(corr12.size());
  for (int i = 0; i < ncorr12; ++i) {
    dist12[i] = 1.0 / (sqrt(dist12[i]) + EPSILON);
    const int r = corr12[i] + m;
    const int c = PointCorrespondence::GetPointIndex(target, &target_sample, i);
    cols[c].push_back(MakePair(r, dist12[i]));
    cols[r].push_back(MakePair(c, dist12[i]));
  }
  const int ncorr21 = static_cast<int>(corr21.size());
  for (int i = 0; i < ncorr21; ++i) {
    dist21[i] = 1.0 / (sqrt(dist21[i]) + EPSILON);
    const int r = corr21[i];
    const int c = PointCorrespondence::GetPointIndex(source, &source_sample, i) + m;
    cols[c].push_back(MakePair(r, dist21[i]));
    cols[r].push_back(MakePair(c, dist21[i]));
  }
  // Compute graph Laplacian of joint connectivity graph
  SparseMatrix L(SparseMatrix::CCS);
  L.Initialize(m + n, m + n, cols);
  NormalizedLaplacian(L, L);
  delete[] cols;
  // Compute eigenmodes of joint graph Laplacian
  Matrix eigenmodes;
  if (ComputeEigenmodes(L, k+2, eigenmodes, eigenvalues) < k) {
    cerr << "FiducialRegistrationError::UpdateEigenmodes: Failed to compute " << k << " eigenmodes" << endl;
    exit(1);
  }
  eigenvalues.Resize(k);
  // Set eigenmodes as point data
  target_info->_Index = SetEigenmodes(target_polydata, eigenmodes, 0, k, "eigenmodes");
  source_info->_Index = SetEigenmodes(source_polydata, eigenmodes, m, k, "eigenmodes");
  target_eigenmodes   = targetPD->GetArray(target_info->_Index);
  source_eigenmodes   = sourcePD->GetArray(source_info->_Index);
  // Calculate weight of eigenmodes, more importance to lower frequencies
  double wsum = .0;
  for (int c = 0; c < k; ++c) {
    wsum += 1.0 / sqrt(eigenvalues(c));
  }
  Vector weight(k);
  for (int c = k; c >= 0; --c) {
    weight(c) = (1.0 / sqrt(eigenvalues(c))) / wsum;
  }
  // Set rescaling parameters s.t. eigenmodes are normalized to range [-1 +1]
  double slope, intercept;
  GetNormalizationParametersOfEigenmodes(slope, intercept, target_eigenmodes, source_eigenmodes, weight);
  target_info->_Slope     = source_info->_Slope     = slope;
  target_info->_Intercept = source_info->_Intercept = intercept;
}

// -----------------------------------------------------------------------------
void UpdateEigenmodesIfUsed(const RegisteredPointSet         *target,
                            Vector                           &target_eigenvalues,
                            PointCorrespondence::FeatureList &target_features,
                            const RegisteredPointSet         *source,
                            Vector                           &source_eigenvalues,
                            PointCorrespondence::FeatureList &source_features,
                            int k, bool diffeo)
{
  // Skip if no spectral coordinates used as features
  if (k == 0) {
    target_eigenvalues.Clear();
    source_eigenvalues.Clear();
    return;
  }
  size_t i;
  for (i = 0; i < target_features.size(); ++i) {
    if ((target_features[i]._Name == "eigenmodes" ||
         target_features[i]._Name == "spectral normals")
        && target_features[i]._Weight != .0) break;
  }
  if (i == target_features.size()) {
    target_eigenvalues.Clear();
    source_eigenvalues.Clear();
    return;
  }
  // Compute eigenmodes
  if (diffeo) {
    UpdateEigenmodes(target, target_features,
                     source, source_features, source_eigenvalues, k);
    target_eigenvalues = source_eigenvalues;
  } else {
    UpdateEigenmodes(target, target_eigenvalues, target_features,
                     source, source_eigenvalues, source_features, k);
  }
}


} // namespace PointCorrespondenceUtils
using namespace PointCorrespondenceUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
PointCorrespondence::PointCorrespondence(const RegisteredPointSet *target,
                                         const RegisteredPointSet *source)
:
  _M(0), _N(0), _NumberOfFeatures(0),
  _Target(target), _TargetSample(NULL), _TargetFeatures(),
  _Source(source), _SourceSample(NULL), _SourceFeatures(),
  _DimensionOfSpectralPoints(5),
  _DiffeomorphicSpectralDecomposition(true),
  _UpdateSpectralPoints(false),
  _FromTargetToSource(true),
  _FromSourceToTarget(true),
  _DefaultDirection(TargetToSource)
{
}

// -----------------------------------------------------------------------------
PointCorrespondence::PointCorrespondence(const PointCorrespondence &other)
:
  _M(other._M), _N(other._N), _NumberOfFeatures(other._NumberOfFeatures),
  _Target(other._Target), _TargetSample(other._TargetSample), _TargetFeatures(other._TargetFeatures),
  _Source(other._Source), _SourceSample(other._SourceSample), _SourceFeatures(other._SourceFeatures),
  _DimensionOfSpectralPoints(other._DimensionOfSpectralPoints),
  _DiffeomorphicSpectralDecomposition(other._DiffeomorphicSpectralDecomposition),
  _UpdateSpectralPoints (other._UpdateSpectralPoints),
  _FromTargetToSource(other._FromTargetToSource),
  _FromSourceToTarget(other._FromSourceToTarget),
  _DefaultDirection(other._DefaultDirection)
{
}

// -----------------------------------------------------------------------------
PointCorrespondence::~PointCorrespondence()
{
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
int PointCorrespondence::GetPointDataIndexByCaseInsensitiveName(vtkPointData *pd, const string &name)
{
  string lname = ToLower(name);
  for (int i = 0; i < pd->GetNumberOfArrays(); ++i) {
    const char *array_name = pd->GetArrayName(i);
    if (array_name == NULL) continue;
    string lower_name = ToLower(array_name);
    if (lower_name == lname) return i;
  }
  return -1;
}

// -----------------------------------------------------------------------------
bool PointCorrespondence::Set(const char *param, const char *value)
{
  if (strcmp(param, "No. of spectral coordinates")    == 0 ||
      strcmp(param, "Number of spectral coordinates") == 0 ||
      strcmp(param, "Dimension of spectral points")   == 0 ||
      strcmp(param, "Spectral coordinates number")    == 0) {
    return FromString(value, _DimensionOfSpectralPoints);
  }
  if (strcmp(param, "Diffeomorphic spectral matching")      == 0 ||
      strcmp(param, "Diffeomorphic spectral coordinates")   == 0 ||
      strcmp(param, "Diffeomorphic spectral decomposition") == 0) {
    return FromString(value, _DiffeomorphicSpectralDecomposition);
  }
  if (strcmp(param, "Spectral coordinates update") == 0) {
    return FromString(value, _UpdateSpectralPoints);
  }

  size_t len = strlen(param);
  if (len > 7 && strcmp(param + len - 7, " weight") == 0) {
    double weight;
    if (!FromString(value, weight) || weight < .0) return false;
    string feature_name(param, param + len - 7);
    AddFeature(feature_name.c_str(), weight);
    return true;
  }

  return false;
}

// -----------------------------------------------------------------------------
ParameterList PointCorrespondence::Parameter() const
{
  ParameterList params = Observable::Parameter();
  Insert(params, "Diffeomorphic spectral decompositon", ToString(_DiffeomorphicSpectralDecomposition));
  Insert(params, "No. of spectral coordinates", ToString(_DimensionOfSpectralPoints));
  Insert(params, "Spectral coordinates update", ToString(_UpdateSpectralPoints));
  string name;
  for (size_t i = 0; i < _TargetFeatures.size(); ++i) {
    name = _TargetFeatures[i]._Name;
    if (name.empty()) continue;
    if (name == "eigenmodes") name = "spectral coordinates";
    name += " weight";
    Insert(params, name, ToString(_TargetFeatures[i]._Slope));
  }
  return params;
}

// -----------------------------------------------------------------------------
static string TransformFeatureName(const char *name)
{
  string feature_name = ToLower(name);
  if (feature_name == "spatial normals") {
    feature_name = "normals";
  } else if (feature_name == "spectral coordinates" ||
             feature_name == "spectral points") {
    feature_name = "eigenmodes";
  } else if (feature_name == "spatial coordinates" ||
             feature_name == "spatial point"       ||
             feature_name == "spatial points"      ||
             feature_name == "points") {
    feature_name = "spatial coordinates";
  }
  return feature_name;
}

// -----------------------------------------------------------------------------
bool PointCorrespondence::AddFeature(const char *name, double weight, double slope, double intercept)
{
  const string feature_name = TransformFeatureName(name);
  int index = -2; // < -1: determined after data is set by CompleteFeatureInfo
  if (feature_name == "spatial coordinates") index = -1;
  _TargetFeatures.push_back(FeatureInfo(feature_name, index, weight, slope, intercept));
  _SourceFeatures.push_back(FeatureInfo(feature_name, index, weight, slope, intercept));
  return true;
}

// -----------------------------------------------------------------------------
void PointCorrespondence::RemoveFeature(const char *name)
{
  const string feature_name = TransformFeatureName(name);
  for (FeatureList::iterator i = _TargetFeatures.begin(); i != _TargetFeatures.end(); ++i) {
    if (i->_Name == feature_name) {
      FeatureList::iterator pos = i; --i;
      _TargetFeatures.erase(pos);
    }
  }
  for (FeatureList::iterator i = _SourceFeatures.begin(); i != _SourceFeatures.end(); ++i) {
    if (i->_Name == feature_name) {
      FeatureList::iterator pos = i; --i;
      _SourceFeatures.erase(pos);
    }
  }
}

// -----------------------------------------------------------------------------
void PointCorrespondence::CompleteFeatureInfo(const RegisteredPointSet *input, FeatureList &feature)
{
  vtkPointData * const inputPD = input->PointSet()->GetPointData();
  for (FeatureList::iterator i = feature.begin(); i != feature.end(); ++i) {
    if (i->_Index < -1) {
      if (i->_Name.empty()) {
        cerr << "PointCorrespondence::CompleteFeatureInfo: Encountered feature without name" << endl;
        exit(1);
      }
      i->_Index = GetPointDataIndexByCaseInsensitiveName(inputPD, i->_Name.c_str());
      // Spectral features are computed and added by UpdateEigenmodesIfUsed when missing
      if (i->_Index < 0 && i->_Name != "eigenmodes" && i->_Name != "spectral normals") {
        cerr << "PointCorrespondence::CompleteFeatureInfo: Missing point data " << i->_Name << endl;
        exit(1);
      }
    } else if (i->_Index >= inputPD->GetNumberOfArrays()) {
      cerr << "PointCorrespondence::CompleteFeatureInfo: Feature index is out of bounds" << endl;
      exit(1);
    }
  }
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void PointCorrespondence::Initialize()
{
  // Check that all inputs are set
  if (_Target == NULL) {
    cerr << "PointCorrespondence::Initialize: Target data set not set" << endl;
    exit(1);
  }
  if (_Source == NULL) {
    cerr << "PointCorrespondence::Initialize: Source data set not set" << endl;
    exit(1);
  }

  // Ensure that correspondences in at least one direction are requested
  if (!_FromTargetToSource && !_FromSourceToTarget) {
    cerr << "PointCorrespondence::Initialize: At least one direction,"
            " either FromTargetToSource or FromSourceToTarget must be enabled" << endl;
    exit(1);
  }

  // Fill missing feature info which requires inspection of the input dataset
  CompleteFeatureInfo(_Target, _TargetFeatures);
  CompleteFeatureInfo(_Source, _SourceFeatures);

  // Determine size of feature vectors
  _NumberOfFeatures = GetPointDimension(_Target, &_TargetFeatures);
  if (GetPointDimension(_Source, &_SourceFeatures) != _NumberOfFeatures) {
    cerr << "PointCorrespondence::Initialize: Mismatching feature vector size" << endl;
    exit(1);
  }

  // Initialize this class
  PointCorrespondence::Init();
}

// -----------------------------------------------------------------------------
void PointCorrespondence::Reinitialize()
{
  // Reinitialize this class
  PointCorrespondence::Init();
}

// -----------------------------------------------------------------------------
void PointCorrespondence::Init()
{
  // Determine number of points (samples)
  _M = GetNumberOfPoints(_Target, _TargetSample);
  _N = GetNumberOfPoints(_Source, _SourceSample);
  if (_M == 0) {
    cerr << "PointCorrespondence::Initialize: Target data set has no points!" << endl;
    exit(1);
  }
  if (_N == 0) {
    cerr << "PointCorrespondence::Initialize: Source data set has no points!" << endl;
    exit(1);
  }

  // Compute spectral coordinates if needed
  UpdateEigenmodesIfUsed(_Target, _TargetEigenvalues, _TargetFeatures,
                         _Source, _SourceEigenvalues, _SourceFeatures,
                         _DimensionOfSpectralPoints, _DiffeomorphicSpectralDecomposition);
}

// -----------------------------------------------------------------------------
void PointCorrespondence::Update()
{
  // Recompute spectral coordinates (optional/experimental)
  if (_UpdateSpectralPoints) {
    UpdateEigenmodesIfUsed(_Target, _TargetEigenvalues, _TargetFeatures,
                           _Source, _SourceEigenvalues, _SourceFeatures,
                           _DimensionOfSpectralPoints, _DiffeomorphicSpectralDecomposition);
  }
}

// -----------------------------------------------------------------------------
bool PointCorrespondence::Upgrade()
{
  return false;
}

// -----------------------------------------------------------------------------
int PointCorrespondence::GetTargetIndex(int) const
{
  return -1;
}

// -----------------------------------------------------------------------------
int PointCorrespondence::GetSourceIndex(int) const
{
  return -1;
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void PointCorrespondence::WriteDataSets(const char *prefix, const char *suffix, bool all) const
{
  if (all || debug >= 4) {
    const int sz = 1024;
    char      fname[sz];
    snprintf(fname, sz, "%sspectral_target_points%s.vtp", prefix, suffix);
    this->WriteSpectralPoints(fname, _Target->PointSet());
    snprintf(fname, sz, "%sspectral_source_points%s.vtp", prefix, suffix);
    this->WriteSpectralPoints(fname, _Source->PointSet());
  }
}

// -----------------------------------------------------------------------------
void PointCorrespondence::WriteSpectralPoints(const char *fname, vtkPointSet *d) const
{
  vtkDataArray *modes = d->GetPointData()->GetArray("eigenmodes");
  if (!modes) return;

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  points->SetNumberOfPoints(d->GetNumberOfPoints());

  double *p = new double[max(3, int(modes->GetNumberOfComponents()))];
  for (int c = 0; c < 3; ++c) p[c] = .0;
  for (vtkIdType i = 0; i < points->GetNumberOfPoints(); ++i) {
    modes ->GetTuple(i, p);
    points->SetPoint(i, p);
  }
  delete[] p;

  vtkSmartPointer<vtkPointSet> output;
  output.TakeReference(d->NewInstance());
  output->ShallowCopy(d);
  output->SetPoints(points);
  WritePointSet(fname, output);
}


} // namespace mirtk
