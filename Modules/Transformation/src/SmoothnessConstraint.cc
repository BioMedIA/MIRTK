/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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

#include "mirtk/SmoothnessConstraint.h"

#include "mirtk/String.h"
#include "mirtk/Memory.h"
#include "mirtk/FreeFormTransformation.h"
#include "mirtk/MultiLevelTransformation.h"
#include "mirtk/ObjectFactory.h"


namespace mirtk {


// Register energy term with object factory during static initialization
mirtkAutoRegisterEnergyTermMacro(SmoothnessConstraint);


// -----------------------------------------------------------------------------
SmoothnessConstraint::SmoothnessConstraint(const char *name, double weight)
:
  TransformationConstraint(name, weight),
  _WithRespectToWorld(true), _UseLatticeSpacing(false),
  _AnnealingRate(.0), _AnnealingWeight(1.0)
{
  _ParameterPrefix.push_back("Smoothness ");
  _ParameterPrefix.push_back("Bending energy ");
  _ParameterPrefix.push_back("Bending ");
}

// -----------------------------------------------------------------------------
bool SmoothnessConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strstr(param, "W.r.t world"          ) == param ||
      strstr(param, "W.r.t. world"         ) == param ||
      strstr(param, "Wrt world"            ) == param ||
      strstr(param, "With respect to world") == param) {
    return FromString(value, _WithRespectToWorld);
  }
  if (strcmp(param, "Use lattice spacing") == 0 ||
      strcmp(param, "Use grid spacing")    == 0 ||
      strcmp(param, "Use spacing")         == 0) {
    return FromString(value, _UseLatticeSpacing);
  }
  if (strstr(param, "Annealing rate") == param) {
    return FromString(value, _AnnealingRate);
  }
  return TransformationConstraint::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList SmoothnessConstraint::Parameter() const
{
  ParameterList params = TransformationConstraint::Parameter();
  InsertWithPrefix(params, "W.r.t. world", _WithRespectToWorld);
  InsertWithPrefix(params, "Use lattice spacing", _UseLatticeSpacing);
  InsertWithPrefix(params, "Annealing rate (experimental)", _AnnealingRate);
  return params;
}

// -----------------------------------------------------------------------------
void SmoothnessConstraint::Initialize()
{
  TransformationConstraint::Initialize();
  _AnnealingWeight = 1.0;
}

// -----------------------------------------------------------------------------
bool SmoothnessConstraint::Upgrade()
{
  TransformationConstraint::Upgrade();
  if (_AnnealingRate != .0) _AnnealingWeight *= _AnnealingRate;
  return false; // whether to continue depends on similarity term, not penalty
}

// -----------------------------------------------------------------------------
double SmoothnessConstraint::Evaluate()
{
  const MultiLevelTransformation *mffd = NULL;
  const FreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  double bending = .0;
  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      bending += ffd->BendingEnergy(_ConstrainPassiveDoFs, _WithRespectToWorld && _UseLatticeSpacing);
    }
  } else if (ffd) {
    bending = ffd->BendingEnergy(_ConstrainPassiveDoFs, _WithRespectToWorld && _UseLatticeSpacing);
  }
  return _AnnealingWeight * bending;
}

// -----------------------------------------------------------------------------
void SmoothnessConstraint::EvaluateGradient(double *gradient, double, double weight)
{
  weight *= _AnnealingWeight;
  const MultiLevelTransformation *mffd = NULL;
  const FreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      ffd->BendingEnergyGradient(gradient, weight, _ConstrainPassiveDoFs, _WithRespectToWorld, _UseLatticeSpacing);
      gradient += ffd->NumberOfDOFs();
    }
  } else if (ffd) {
    ffd->BendingEnergyGradient(gradient, weight, _ConstrainPassiveDoFs, _WithRespectToWorld, _UseLatticeSpacing);
  }
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void SmoothnessConstraint::WriteGradient(const char *p, const char *suffix) const
{
  const int   sz = 1024;
  char        fname[sz];
  string _prefix = Prefix(p);
  const char  *prefix = _prefix.c_str();

  const MultiLevelTransformation *mffd = NULL;
  const FreeFormTransformation   *ffd  = NULL;

  (mffd = MFFD()) || (ffd = FFD());

  if (mffd) {
    for (int l = 0; l < mffd->NumberOfLevels(); ++l) {
      if (!mffd->LocalTransformationIsActive(l)) continue;
      ffd = mffd->GetLocalTransformation(l);
      double *gradient = CAllocate<double>(ffd->NumberOfDOFs());
      ffd->BendingEnergyGradient(gradient, 1.0, _ConstrainPassiveDoFs, _WithRespectToWorld, _UseLatticeSpacing);
      if (mffd->NumberOfActiveLevels() == 1) {
        snprintf(fname, sz, "%sgradient%s", prefix, suffix);
      } else {
        snprintf(fname, sz, "%sgradient_of_ffd_at_level_%d%s", prefix, l+1, suffix);
      }
      WriteFFDGradient(fname, ffd, gradient);
      Deallocate(gradient);
    }
  } else if (ffd) {
    snprintf(fname, sz, "%sgradient%s", prefix, suffix);
    double *gradient = CAllocate<double>(ffd->NumberOfDOFs());
    ffd->BendingEnergyGradient(gradient, 1.0, _ConstrainPassiveDoFs, _WithRespectToWorld, _UseLatticeSpacing);
    WriteFFDGradient(fname, ffd, gradient);
    Deallocate(gradient);
  }
}


} // namespace mirtk
