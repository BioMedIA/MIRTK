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

#include "mirtk/PartialAffineTransformation.h"

#include "mirtk/Assert.h"
#include "mirtk/Math.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
PartialAffineTransformation::PartialAffineTransformation(HomogeneousTransformation *t, double f)
:
  _Transformation(NULL),
  _Fraction(f)
{
  _TransformationObserver.Bind(
    ModifiedEvent,
    MakeDelegate(this, &PartialAffineTransformation::OnTransformationChanged)
  );
  Transformation(t);
}

// -----------------------------------------------------------------------------
PartialAffineTransformation::~PartialAffineTransformation()
{
  Transformation(NULL);
}

// -----------------------------------------------------------------------------
void PartialAffineTransformation::Transformation(HomogeneousTransformation *t)
{
  mirtkAssert(!t || t->NumberOfDOFs() <= 12, "Transformation must be rigid, similarity, or affine");
  if (_Transformation) _Transformation->DeleteObserver(_TransformationObserver);
  _Transformation = t;
  if (_Transformation) _Transformation->AddObserver(_TransformationObserver);
  this->OnTransformationChanged();
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool PartialAffineTransformation::HasSameDOFsAs(const class Transformation *t) const
{
  return HaveSameDOFs(_Transformation, t);
}

// -----------------------------------------------------------------------------
void PartialAffineTransformation::OnTransformationChanged()
{
  if (!_Transformation) return;
  // Copy status of transformation parameters
  const int ndofs = _Transformation->NumberOfDOFs();
  if (ndofs <= 12) {
    _Status[TX] = _Transformation->GetStatus(TX);
    _Status[TY] = _Transformation->GetStatus(TY);
    _Status[TZ] = _Transformation->GetStatus(TZ);
    _Status[RX] = _Transformation->GetStatus(RX);
    _Status[RY] = _Transformation->GetStatus(RY);
    _Status[RZ] = _Transformation->GetStatus(RZ);
  }
  if (ndofs == 7) {
    _Status[SX] = _Status[SY] = _Status[SZ] = _Transformation->GetStatus(SG);
  } else if (ndofs == 9 || ndofs == 12) {
    _Status[SX] = _Transformation->GetStatus(SX);
    _Status[SY] = _Transformation->GetStatus(SY);
    _Status[SZ] = _Transformation->GetStatus(SZ);
  }
  if (ndofs == 12) {
    _Status[SXY] = _Transformation->GetStatus(SXY);
    _Status[SXZ] = _Transformation->GetStatus(SXZ);
    _Status[SYZ] = _Transformation->GetStatus(SYZ);
  }
  // Compute partial transformation matrix
  this->PutMatrix((_Transformation->GetMatrix().Log() *= _Fraction).Exp());
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void PartialAffineTransformation
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  double param[12];

  // Compute gradient w.r.t parameters of partial transformation
  AffineTransformation::ParametricGradient(in, param, i2w, wc, t0);

  // Compute normalization factor -- needed for numerical stability
  double max = .0;
  for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
    if (abs(param[dof]) > max) max = abs(param[dof]);
  }
  if (max == .0) return;

  // Convert to gradient w.r.t parameters of decorated transformation
  for (int dof = 0; dof < _NumberOfDOFs; ++dof) {
    param[dof] = _Param[dof] + param[dof] / max;
  }
  Matrix2DOFs(expm(logm(DOFs2Matrix(param)) /= _Fraction), param);
  for (int dof = 0; dof < _Transformation->NumberOfDOFs(); ++dof) {
    if (_Transformation->GetStatus(dof) == Active) {
      out[dof] += w * (param[dof] - _Transformation->Get(dof));
    }
  }
}


} // namespace mirtk
