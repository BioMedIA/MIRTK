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

#include "mirtk/PartialBSplineFreeFormTransformationSV.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
PartialBSplineFreeFormTransformationSV
::PartialBSplineFreeFormTransformationSV(BSplineFreeFormTransformationSV *t, double f)
:
  _Transformation(t),
  _Fraction      (f)
{
}

// -----------------------------------------------------------------------------
PartialBSplineFreeFormTransformationSV
::~PartialBSplineFreeFormTransformationSV()
{
}

// =============================================================================
// Transformation parameters
// =============================================================================

// -----------------------------------------------------------------------------
bool PartialBSplineFreeFormTransformationSV::CopyFrom(const class Transformation *t)
{
  const PartialBSplineFreeFormTransformationSV *part = NULL;
  if ((part = dynamic_cast<const PartialBSplineFreeFormTransformationSV *>(t))) {
    return _Transformation->CopyFrom(part->Transformation());
  } else {
    return _Transformation->CopyFrom(t);
  }
}

// -----------------------------------------------------------------------------
int PartialBSplineFreeFormTransformationSV::NumberOfDOFs() const
{
  return _Transformation->NumberOfDOFs();
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Put(int dof, double x)
{
  _Transformation->Put(dof, x);
}

// -----------------------------------------------------------------------------
double PartialBSplineFreeFormTransformationSV::Get(int dof) const
{
  return _Transformation->Get(dof);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Put(const DOFValue *x)
{
  _Transformation->Put(x);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Add(const DOFValue *dx)
{
  _Transformation->Add(dx);
}

// -----------------------------------------------------------------------------
double PartialBSplineFreeFormTransformationSV::Update(const DOFValue *dx)
{
  return _Transformation->Update(dx);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Get(DOFValue *x) const
{
  _Transformation->Get(x);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::PutStatus(int dof, DOFStatus s)
{
  _Transformation->PutStatus(dof, s);
}

// -----------------------------------------------------------------------------
PartialBSplineFreeFormTransformationSV::DOFStatus
PartialBSplineFreeFormTransformationSV::GetStatus(int dof) const
{
  return _Transformation->GetStatus(dof);
}

// -----------------------------------------------------------------------------
bool PartialBSplineFreeFormTransformationSV::HasSameDOFsAs(const class Transformation *t) const
{
  return HaveSameDOFs(_Transformation, t);
}

// -----------------------------------------------------------------------------
bool PartialBSplineFreeFormTransformationSV::IsIdentity() const
{
  return (_Fraction == .0) || _Transformation->IsIdentity();
}

// =============================================================================
// Parameters (non-DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
bool PartialBSplineFreeFormTransformationSV::Set(const char *param, const char *value)
{
  return _Transformation->Set(param, value);
}

// -----------------------------------------------------------------------------
ParameterList PartialBSplineFreeFormTransformationSV::Parameter() const
{
  return _Transformation->Parameter();
}

// =============================================================================
// Point transformation
// =============================================================================

// -----------------------------------------------------------------------------
bool PartialBSplineFreeFormTransformationSV::RequiresCachingOfDisplacements() const
{
  return _Transformation->RequiresCachingOfDisplacements();
}

// -----------------------------------------------------------------------------
inline double PartialBSplineFreeFormTransformationSV::UpperIntegrationLimit(double t, double t0) const
{
  return _Fraction * _Transformation->UpperIntegrationLimit(t, t0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::GlobalTransform(double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->GlobalTransform(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::LocalTransform(double &x, double &y, double &z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->LocalTransform(x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Transform(double &x, double &y, double &z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Transform(x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::GlobalInverse(double &x, double &y, double &z, double t, double t0) const
{
  _Transformation->GlobalInverse(x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
bool PartialBSplineFreeFormTransformationSV::LocalInverse(double &x, double &y, double &z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->LocalInverse(x, y, z, T, .0);
  return true;
}

// -----------------------------------------------------------------------------
bool PartialBSplineFreeFormTransformationSV::Inverse(double &x, double &y, double &z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Inverse(x, y, z, T, .0);
  return true;
}

// ---------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Displacement(GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Displacement(disp, T, .0, wc);
}

// ---------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Displacement(GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Displacement(disp, T, .0, wc);
}

// ---------------------------------------------------------------------------
int PartialBSplineFreeFormTransformationSV::InverseDisplacement(GenericImage<double> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->InverseDisplacement(disp, T, .0, wc);
  return 0;
}

// ---------------------------------------------------------------------------
int PartialBSplineFreeFormTransformationSV::InverseDisplacement(GenericImage<float> &disp, double t, double t0, const WorldCoordsImage *wc) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->InverseDisplacement(disp, T, .0, wc);
  return 0;
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::GlobalJacobian(Matrix &jac, double x, double y, double z, double t, double t0) const
{
  _Transformation->GlobalJacobian(jac, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::LocalJacobian(Matrix &jac, double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->LocalJacobian(jac, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Jacobian(Matrix &jac, double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Jacobian(jac, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::GlobalHessian(Matrix h[3], double x, double y, double z, double t, double t0) const
{
  _Transformation->GlobalHessian(h, x, y, z, t, t0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::LocalHessian(Matrix h[3], double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->LocalHessian(h, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Hessian(Matrix h[3], double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->Hessian(h, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::JacobianDOFs(double jac[3], int dof, double x, double y, double z, double t, double t0) const
{
  const double T = UpperIntegrationLimit(t, t0);
  if (T) _Transformation->JacobianDOFs(jac, dof, x, y, z, T, .0);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  const double T = UpperIntegrationLimit(in->GetTOrigin(), t0);
  if (T) _Transformation->ParametricGradient(in, out, i2w, wc, T, .0, w);
}

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV
::ParametricGradient(const GenericImage<double> **in, int n, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     const double *t0, double w) const
{
  double t, T;
  for (int i = 0; i < n; ++i) {
    t = in[i]->GetTOrigin();
    T = UpperIntegrationLimit(t, t0 ? t0[i] : t);
    if (T) _Transformation->ParametricGradient(in[i], out, i2w, wc, T, .0, w);
  }
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void PartialBSplineFreeFormTransformationSV::Print(ostream &os, Indent indent) const
{
  os << indent << "Partial B-spline SV FFD" << endl;
  _Transformation->Print(os, indent + 1);
  os << indent << "Fraction of transformation applied: " << _Fraction << endl;
}

// -----------------------------------------------------------------------------
Cifstream &PartialBSplineFreeFormTransformationSV::Read(Cifstream &from)
{
  return _Transformation->Read(from);
}

// -----------------------------------------------------------------------------
Cofstream &PartialBSplineFreeFormTransformationSV::Write(Cofstream &to) const
{
  return _Transformation->Write(to);
}


} // namespace mirtk
