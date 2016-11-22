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

#include "mirtk/RobustClosestPoint.h"

#include "mirtk/Math.h"
#include "mirtk/Array.h"
#include "mirtk/Pair.h"
#include "mirtk/RegisteredPointSet.h"
#include "mirtk/PointLocator.h"
#include "mirtk/SparseMatrix.h"
#include "mirtk/Profiling.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
RobustClosestPoint::RobustClosestPoint()
:
  _Sigma(3.0)
{
}

// -----------------------------------------------------------------------------
RobustClosestPoint::RobustClosestPoint(const RegisteredPointSet *target,
                                       const RegisteredPointSet *source)
:
  _Sigma(3.0)
{
  Target(target);
  Source(source);
  Initialize();
}

// -----------------------------------------------------------------------------
RobustClosestPoint::RobustClosestPoint(const RobustClosestPoint &other)
:
  FuzzyCorrespondence(other),
  _Sigma(other._Sigma)
{
}

// -----------------------------------------------------------------------------
PointCorrespondence *RobustClosestPoint::NewInstance() const
{
  return new RobustClosestPoint(*this);
}

// -----------------------------------------------------------------------------
RobustClosestPoint::~RobustClosestPoint()
{
}

// -----------------------------------------------------------------------------
RobustClosestPoint::TypeId RobustClosestPoint::Type() const
{
  return TypeId::RobustClosestPoint;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool RobustClosestPoint::Set(const char *name, const char *value)
{
  // Standard deviation of outlier rejection
  if (strcmp(name, "Sigma") == 0) {
    return FromString(value, _Sigma);
  }
  return FuzzyCorrespondence::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList RobustClosestPoint::Parameter() const
{
  ParameterList params = FuzzyCorrespondence::Parameter();
  Insert(params, "Sigma", _Sigma);
  return params;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void RobustClosestPoint::CalculateWeights()
{
  MIRTK_START_TIMING();

  // Initialize weight matrix
  _Weight.Initialize(_M, _N);

  // Find closest points
  Array<int   > corr12, corr21;
  Array<double> dist12, dist21;

  UniquePtr<PointLocator> source_locator;
  source_locator.reset(PointLocator::New(_Source->PointSet(), _SourceSample, &_SourceFeatures));
  source_locator->GlobalIndices(true);
  corr12 = source_locator->FindClosestPoint(_Target->PointSet(), _TargetSample, &_TargetFeatures, &dist12);
  source_locator.reset();

  UniquePtr<PointLocator> target_locator;
  target_locator.reset(PointLocator::New(_Target->PointSet(), _TargetSample, &_TargetFeatures));
  target_locator->GlobalIndices(true);
  corr21 = target_locator->FindClosestPoint(_Source->PointSet(), _SourceSample, &_SourceFeatures, &dist21);
  target_locator.reset();

  // Allocate lists for non-zero weight entries
  const int nentries = (_Weight.Layout() == WeightMatrix::CRS ? _M : _N);
  WeightMatrix::Entries *entries = new WeightMatrix::Entries[nentries];
  for (int i = 0; i < nentries; ++i) entries[i].reserve(2);

  // Symmetric closest point matching with outlier rejection
  const WeightMatrix::EntryType one(1.0);
  const int ncorr12 = static_cast<int>(corr12.size());
  const int ncorr21 = static_cast<int>(corr21.size());

  if (_Sigma > .0) {

    double max_dist, dist_mean = .0, dist_mean2 = .0;
    for (Array<double>::iterator dist = dist12.begin(); dist != dist12.end(); ++dist) {
      dist_mean  += (*dist);
      dist_mean2 += (*dist) * (*dist);
    }
    for (Array<double>::iterator dist = dist21.begin(); dist != dist21.end(); ++dist) {
      dist_mean  += (*dist);
      dist_mean2 += (*dist) * (*dist);
    }
    dist_mean  /= (dist12.size() + dist21.size());
    dist_mean2 /= (dist12.size() + dist21.size());

    max_dist = dist_mean + _Sigma * sqrt(dist_mean2 - dist_mean * dist_mean);

    if (_Weight.Layout() == WeightMatrix::CRS) {
      for (int r = 0; r < ncorr12; ++r) {
        if (dist12[r] <= max_dist) entries[r].push_back(MakePair(corr12[r], one));
      }
      for (int c = 0; c < ncorr21; ++c) {
        if (dist21[c] <= max_dist) entries[corr21[c]].push_back(MakePair(c, one));
      }
    } else {
      for (int r = 0; r < ncorr12; ++r) {
        if (dist12[r] <= max_dist) entries[corr12[r]].push_back(MakePair(r, one));
      }
      for (int c = 0; c < ncorr21; ++c) {
        if (dist21[c] <= max_dist) entries[c].push_back(MakePair(corr21[c], one));
      }
    }

  // Symmetric closest point matching without outlier rejection
  } else {

    if (_Weight.Layout() == WeightMatrix::CRS) {
      for (int r = 0; r < ncorr12; ++r) {
        entries[r].push_back(MakePair(corr12[r], one));
      }
      for (int c = 0; c < ncorr21; ++c) {
        entries[corr21[c]].push_back(MakePair(c, one));
      }
    } else {
      for (int r = 0; r < ncorr12; ++r) {
        entries[corr12[r]].push_back(MakePair(r, one));
      }
      for (int c = 0; c < ncorr21; ++c) {
        entries[c].push_back(MakePair(corr21[c], one));
      }
    }

  }

  // Initialize weight matrix
  _Weight.Initialize(_M, _N, entries);

  // Normalize weights to sum up to 1
  _Weight *= 0.5;

  MIRTK_DEBUG_TIMING(6, "calculating correspondence weights");
}


} // namespace mirtk
