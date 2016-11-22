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

#include "mirtk/ClosestPoint.h"

#include "mirtk/Math.h"
#include "mirtk/Array.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
ClosestPoint::ClosestPoint()
:
  _Sigma(-1.0),
  _MaxDistance(numeric_limits<double>::infinity()),
  _MaxSquaredDistance(numeric_limits<double>::infinity())
{
}

// -----------------------------------------------------------------------------
ClosestPoint::ClosestPoint(const ClosestPoint &other)
:
  PointCorrespondence(other),
  _Sigma(other._Sigma),
  _MaxDistance(other._MaxDistance),
  _MaxSquaredDistance(other._MaxSquaredDistance),
  _TargetIndex(other._TargetIndex),
  _SourceIndex(other._SourceIndex)
{
}

// -----------------------------------------------------------------------------
PointCorrespondence *ClosestPoint::NewInstance() const
{
  return new ClosestPoint(*this);
}

// -----------------------------------------------------------------------------
ClosestPoint::~ClosestPoint()
{
}

// -----------------------------------------------------------------------------
ClosestPoint::TypeId ClosestPoint::Type() const
{
  return TypeId::ClosestPoint;
}

// =============================================================================
// Parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool ClosestPoint::Set(const char *name, const char *value)
{
  if (strcmp(name, "Sigma") == 0) {
    return FromString(value, _Sigma);
  }
  if (strcmp(name, "Maximum distance") == 0) {
    return FromString(value, _MaxDistance);
  }
  return PointCorrespondence::Set(name, value);
}

// -----------------------------------------------------------------------------
ParameterList ClosestPoint::Parameter() const
{
  ParameterList params = PointCorrespondence::Parameter();
  if (_Sigma >= .0) {
    Insert(params, "Sigma", _Sigma);
  } else if (_MaxDistance > .0 && !IsInf(_MaxDistance)) {
    Insert(params, "Maximum distance", _MaxDistance);
  }
  return params;
}

// =============================================================================
// Correspondences
// =============================================================================

// -----------------------------------------------------------------------------
void ClosestPoint::Initialize()
{
  // Initialize base class
  PointCorrespondence::Initialize();

  // Set maximum squared distance threshold
  if (_Sigma < .0) {
    if (_MaxDistance > .0 && !IsInf(_MaxDistance)) {
      _MaxSquaredDistance = _MaxDistance * _MaxDistance;
    } else {
      _MaxSquaredDistance = numeric_limits<double>::infinity();
    }
  }
}

// -----------------------------------------------------------------------------
void ClosestPoint::Update()
{
  // Update point features
  PointCorrespondence::Update();

  // Find closest points
  if (_FromTargetToSource) {
    UniquePtr<PointLocator> locator(PointLocator::New(_Source->PointSet(), _SourceSample, &_SourceFeatures));
    locator->GlobalIndices(true);
    _SourceIndex = locator->FindClosestPoint(_Target->PointSet(), _TargetSample, &_TargetFeatures, &_SourceDistance);
  } else {
    _SourceIndex   .clear();
    _SourceDistance.clear();
  }

  if (_FromSourceToTarget) {
    UniquePtr<PointLocator> locator(PointLocator::New(_Target->PointSet(), _TargetSample, &_TargetFeatures));
    locator->GlobalIndices(true);
    _TargetIndex = locator->FindClosestPoint(_Source->PointSet(), _SourceSample, &_SourceFeatures, &_TargetDistance);
  } else {
    _TargetIndex   .clear();
    _TargetDistance.clear();
  }

  // Update maximum distance threshold
  if (_Sigma >= .0) {
    Array<double>::const_iterator d;
    // Robust evaluation of variance of corresponding squared point distances
    int    m = 0;
    double mean = .0, var = .0, delta;
    for (d = _TargetDistance.begin(); d != _TargetDistance.end(); ++d) {
      ++m;
      delta  = (*d) - mean;
      mean  += delta / m;
      var   += delta * ((*d) - mean);
    }
    for (d = _SourceDistance.begin(); d != _SourceDistance.end(); ++d) {
      ++m;
      delta  = (*d) - mean;
      mean  += delta / m;
      var   += delta * ((*d) - mean);
    }
    if (m > 1) var /= m - 1;
    else       var  = .0;
    // Set maximum distance to mean plus _Sigma times standard deviation
    _MaxSquaredDistance = mean + _Sigma * sqrt(var);
  }
}

// -----------------------------------------------------------------------------
bool ClosestPoint::Upgrade()
{
  Array<int> target_index(_TargetIndex), source_index(_SourceIndex);
  this->Update();
  return _TargetIndex != target_index || _SourceIndex != source_index;
}

// -----------------------------------------------------------------------------
bool ClosestPoint::GetInputTargetPoint(int i, Point &p) const
{
  if (_TargetIndex[i] < 0) return false;
  _Target->GetInputPoint(_TargetIndex[i], p);
  return _TargetDistance[i] <= _MaxSquaredDistance;
}

// -----------------------------------------------------------------------------
bool ClosestPoint::GetInputSourcePoint(int i, Point &p) const
{
  if (_SourceIndex[i] < 0) return false;
  _Source->GetInputPoint(_SourceIndex[i], p);
  return _SourceDistance[i] <= _MaxSquaredDistance;
}

// -----------------------------------------------------------------------------
bool ClosestPoint::GetTargetPoint(int i, Point &p) const
{
  if (_TargetIndex[i] < 0) return false;
  _Target->GetPoint(_TargetIndex[i], p);
  return _TargetDistance[i] <= _MaxSquaredDistance;
}

// -----------------------------------------------------------------------------
bool ClosestPoint::GetSourcePoint(int i, Point &p) const
{
  if (_SourceIndex[i] < 0) return false;
  _Source->GetPoint(_SourceIndex[i], p);
  return _SourceDistance[i] <= _MaxSquaredDistance;
}

// -----------------------------------------------------------------------------
int ClosestPoint::GetTargetIndex(int i) const
{
  return _TargetIndex[i];
}

// -----------------------------------------------------------------------------
int ClosestPoint::GetSourceIndex(int i) const
{
  return _SourceIndex[i];
}


} // namespace mirtk
