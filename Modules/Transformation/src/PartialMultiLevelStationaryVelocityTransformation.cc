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

#include "mirtk/PartialMultiLevelStationaryVelocityTransformation.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
PartialMultiLevelStationaryVelocityTransformation
::PartialMultiLevelStationaryVelocityTransformation(MultiLevelStationaryVelocityTransformation *t, double f)
:
  _Transformation(t),
  _Fraction      (f)
{
}

// -----------------------------------------------------------------------------
PartialMultiLevelStationaryVelocityTransformation
::~PartialMultiLevelStationaryVelocityTransformation()
{
}

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
int PartialMultiLevelStationaryVelocityTransformation::NumberOfLevels() const
{
  return _Transformation->NumberOfLevels();
}

// -----------------------------------------------------------------------------
AffineTransformation *PartialMultiLevelStationaryVelocityTransformation
::GetGlobalTransformation()
{
  return _Transformation->GetGlobalTransformation();
}

// -----------------------------------------------------------------------------
const AffineTransformation *PartialMultiLevelStationaryVelocityTransformation
::GetGlobalTransformation() const
{
  return _Transformation->GetGlobalTransformation();
}

// -----------------------------------------------------------------------------
FreeFormTransformation *PartialMultiLevelStationaryVelocityTransformation
::GetLocalTransformation(int pos)
{
  return _Transformation->GetLocalTransformation(pos);
}

// -----------------------------------------------------------------------------
const FreeFormTransformation *PartialMultiLevelStationaryVelocityTransformation
::GetLocalTransformation(int pos) const
{
  return _Transformation->GetLocalTransformation(pos);
}

// ---------------------------------------------------------------------------
FreeFormTransformation *PartialMultiLevelStationaryVelocityTransformation
::PutLocalTransformation(FreeFormTransformation *ffd, int pos, bool transfer_ownership)
{
  return _Transformation->PutLocalTransformation(ffd, pos, transfer_ownership);
}

// ---------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation
::PushLocalTransformation(FreeFormTransformation *ffd, bool transfer_ownership)
{
  _Transformation->PushLocalTransformation(ffd, transfer_ownership);
}

// ---------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation
::InsertLocalTransformation(FreeFormTransformation *ffd, int pos, bool transfer_ownership)
{
  _Transformation->InsertLocalTransformation(ffd, pos, transfer_ownership);
}

// ---------------------------------------------------------------------------
FreeFormTransformation *PartialMultiLevelStationaryVelocityTransformation
::PopLocalTransformation()
{
  return _Transformation->PopLocalTransformation();
}

// ---------------------------------------------------------------------------
FreeFormTransformation *PartialMultiLevelStationaryVelocityTransformation
::RemoveLocalTransformation(int pos)
{
  return _Transformation->RemoveLocalTransformation(pos);
}

// ---------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation
::CombineLocalTransformation()
{
  _Transformation->CombineLocalTransformation();
}

// ---------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation
::MergeGlobalIntoLocalDisplacement()
{
  _Transformation->MergeGlobalIntoLocalDisplacement();
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool PartialMultiLevelStationaryVelocityTransformation::CopyFrom(const class Transformation *t)
{
  const PartialMultiLevelStationaryVelocityTransformation *part = NULL;
  if ((part = dynamic_cast<const PartialMultiLevelStationaryVelocityTransformation *>(t))) {
    return _Transformation->CopyFrom(part->Transformation());
  } else {
    return _Transformation->CopyFrom(t);
  }
}

// -----------------------------------------------------------------------------
int PartialMultiLevelStationaryVelocityTransformation::NumberOfDOFs() const
{
  return _Transformation->NumberOfDOFs();
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::Put(int dof, double x)
{
  _Transformation->Put(dof, x);
}

// -----------------------------------------------------------------------------
double PartialMultiLevelStationaryVelocityTransformation::Get(int dof) const
{
  return _Transformation->Get(dof);
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::Put(const DOFValue *x)
{
  _Transformation->Put(x);
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::Add(const DOFValue *dx)
{
  _Transformation->Add(dx);
}

// -----------------------------------------------------------------------------
double PartialMultiLevelStationaryVelocityTransformation::Update(const DOFValue *dx)
{
  return _Transformation->Update(dx);
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::Get(DOFValue *x) const
{
  _Transformation->Get(x);
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::PutStatus(int dof, DOFStatus s)
{
  _Transformation->PutStatus(dof, s);
}

// -----------------------------------------------------------------------------
PartialMultiLevelStationaryVelocityTransformation::DOFStatus
PartialMultiLevelStationaryVelocityTransformation::GetStatus(int dof) const
{
  return _Transformation->GetStatus(dof);
}

// -----------------------------------------------------------------------------
bool PartialMultiLevelStationaryVelocityTransformation::HasSameDOFsAs(const class Transformation *t) const
{
  return HaveSameDOFs(_Transformation, t);
}

// -----------------------------------------------------------------------------
bool PartialMultiLevelStationaryVelocityTransformation::IsIdentity() const
{
  return (_Fraction == .0) || _Transformation->IsIdentity();
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool PartialMultiLevelStationaryVelocityTransformation::Set(const char *param, const char *value)
{
  return _Transformation->Set(param, value);
}

// -----------------------------------------------------------------------------
ParameterList PartialMultiLevelStationaryVelocityTransformation::Parameter() const
{
  return _Transformation->Parameter();
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
bool PartialMultiLevelStationaryVelocityTransformation::RequiresCachingOfDisplacements() const
{
  return _Transformation->RequiresCachingOfDisplacements();
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::GlobalTransform(double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->GlobalTransform(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::LocalTransform(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->RKE1((m < 0 ? 0 : m), n, x, y, z, + _Fraction * _Transformation->UpperIntegrationLimit(t, t0));
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::Transform(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->RKE1(m, n, x, y, z, + _Fraction * _Transformation->UpperIntegrationLimit(t, t0));
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::GlobalInverse(double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->GlobalInverse(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
bool PartialMultiLevelStationaryVelocityTransformation::LocalInverse(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->RKE1((m < 0 ? 0 : m), n, x, y, z, - _Fraction * _Transformation->UpperIntegrationLimit(t, t0));
  return true;
}

// -----------------------------------------------------------------------------
bool PartialMultiLevelStationaryVelocityTransformation::Inverse(int m, int n, double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->RKE1(m, n, x, y, z, - _Fraction * _Transformation->UpperIntegrationLimit(t, t0));
  return true;
}

// ---------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::Displacement(int m, int n, GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Displacement(m, n, disp, T, .0, wc);
}

// ---------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::Displacement(int m, int n, GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Displacement(m, n, disp, T, .0, wc);
}

// ---------------------------------------------------------------------------
int PartialMultiLevelStationaryVelocityTransformation::InverseDisplacement(int m, int n, GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->InverseDisplacement(m, n, disp, T, .0, wc);
  return 0;
}

// ---------------------------------------------------------------------------
int PartialMultiLevelStationaryVelocityTransformation::InverseDisplacement(int m, int n, GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->InverseDisplacement(m, n, disp, T, .0, wc);
  return 0;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  double t = in->GetTOrigin();
  double T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
  if (T) _Transformation->ParametricGradient(in, out, i2w, wc, T, .0, w);
}

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation
::ParametricGradient(const GenericImage<double> **in, int n, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     const double *t0, double w) const
{
  double t, T;
  for (int i = 0; i < n; ++i) {
    t = in[i]->GetTOrigin();
    T = _Fraction * _Transformation->UpperIntegrationLimit(t, t0 ? t0[i] : t);
    if (T) _Transformation->ParametricGradient(in[i], out, i2w, wc, T, .0, w);
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void PartialMultiLevelStationaryVelocityTransformation::Print(ostream &os, Indent indent) const
{
  os << indent << "Partial Multi-level B-spline SV FFD" << endl;
  _Transformation->Print(os, indent + 1);
  os << indent << "Fraction of transformation applied: " << _Fraction << endl;
}

// -----------------------------------------------------------------------------
Cifstream &PartialMultiLevelStationaryVelocityTransformation::Read(Cifstream &from)
{
  return _Transformation->Read(from);
}

// -----------------------------------------------------------------------------
Cofstream &PartialMultiLevelStationaryVelocityTransformation::Write(Cofstream &to) const
{
  return _Transformation->Write(to);
}


} // namespace mirtk
