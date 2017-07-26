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

#include "mirtk/TransformationConstraint.h"


namespace mirtk {


// =============================================================================
// Factory method
// =============================================================================

// -----------------------------------------------------------------------------
TransformationConstraint *TransformationConstraint::New(ConstraintMeasure cm, const char *name, double w)
{
  enum EnergyMeasure em = static_cast<enum EnergyMeasure>(cm);
  if (CM_Begin < em && em < CM_End) {
    EnergyTerm *term = EnergyTerm::TryNew(em, name, w);
    if (term) return dynamic_cast<TransformationConstraint *>(term);
    cerr << NameOfType() << "::New: Transformation constraint not available: ";
  } else {
    cerr << NameOfType() << "::New: Energy term is no transformation constraint: ";
  }
  cerr << ToString(em) << " (" << em << ")" << endl;
  exit(1);
  return NULL;
}

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
TransformationConstraint::TransformationConstraint(const char *name, double weight)
:
  EnergyTerm(name, weight),
  _ConstrainPassiveDoFs(false)
{
}

// -----------------------------------------------------------------------------
TransformationConstraint::TransformationConstraint(const TransformationConstraint &other)
:
  EnergyTerm(other),
  _ConstrainPassiveDoFs(other._ConstrainPassiveDoFs)
{
}

// -----------------------------------------------------------------------------
TransformationConstraint &TransformationConstraint::operator =(const TransformationConstraint &other)
{
  if (this != &other) {
    EnergyTerm::operator =(other);
    _ConstrainPassiveDoFs = other._ConstrainPassiveDoFs;
  }
  return *this;
}

// -----------------------------------------------------------------------------
TransformationConstraint::~TransformationConstraint()
{
}

// =============================================================================
// Configuration
// =============================================================================

// -----------------------------------------------------------------------------
bool TransformationConstraint::SetWithPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Constrain passive DoFs")           == 0 ||
      strcmp(param, "Constrain passive CPs")            == 0 ||
      strcmp(param, "Constrain passive control points") == 0 ||
      strcmp(param, "Constrain passive parameters")     == 0) {
    return FromString(value, _ConstrainPassiveDoFs);
  }
  return EnergyTerm::SetWithPrefix(param, value);
}

// -----------------------------------------------------------------------------
bool TransformationConstraint::SetWithoutPrefix(const char *param, const char *value)
{
  if (strcmp(param, "Of passive DoFs")           == 0 ||
      strcmp(param, "Of passive parameters")     == 0 ||
      strcmp(param, "At passive control points") == 0) {
    return FromString(value, _ConstrainPassiveDoFs);
  }
  return EnergyTerm::SetWithoutPrefix(param, value);
}

// -----------------------------------------------------------------------------
ParameterList TransformationConstraint::Parameter() const
{
  ParameterList params = EnergyTerm::Parameter();
  InsertWithPrefix(params, "At passive control points", _ConstrainPassiveDoFs);
  return params;
}

// =============================================================================
// Debugging
// =============================================================================

// -----------------------------------------------------------------------------
void TransformationConstraint
::WriteFFDGradient(const char *fname, const FreeFormTransformation *ffd, const double *g) const
{
  typedef FreeFormTransformation::CPValue CPValue;
  typedef GenericImage<CPValue>           CPImage;
  CPValue *data = reinterpret_cast<CPValue *>(const_cast<double *>(g));
  CPImage gradient(ffd->Attributes(), data);
  gradient.Write(fname);
}


} // namespace mirtk
