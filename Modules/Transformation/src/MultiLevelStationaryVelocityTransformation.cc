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

#include "mirtk/MultiLevelStationaryVelocityTransformation.h"

#include "mirtk/Event.h"
#include "mirtk/Array.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/ImageToInterpolationCoefficients.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
MultiLevelStationaryVelocityTransformation
::MultiLevelStationaryVelocityTransformation()
{
  _GlobalTransformationObserver.Bind(
    ModifiedEvent,
    MakeDelegate(this, &MultiLevelStationaryVelocityTransformation::UpdateLogMatrix)
  );
  _GlobalTransformation.AddObserver(_GlobalTransformationObserver);
  UpdateLogMatrix();
}

// -----------------------------------------------------------------------------
MultiLevelStationaryVelocityTransformation
::MultiLevelStationaryVelocityTransformation(const RigidTransformation &t)
:
  MultiLevelTransformation(t)
{
  _GlobalTransformationObserver.Bind(
    ModifiedEvent,
    MakeDelegate(this, &MultiLevelStationaryVelocityTransformation::UpdateLogMatrix)
  );
  _GlobalTransformation.AddObserver(_GlobalTransformationObserver);
  UpdateLogMatrix();
}

// -----------------------------------------------------------------------------
MultiLevelStationaryVelocityTransformation
::MultiLevelStationaryVelocityTransformation(const AffineTransformation &t)
:
  MultiLevelTransformation(t)
{
  _GlobalTransformationObserver.Bind(
    ModifiedEvent,
    MakeDelegate(this, &MultiLevelStationaryVelocityTransformation::UpdateLogMatrix)
  );
  _GlobalTransformation.AddObserver(_GlobalTransformationObserver);
  UpdateLogMatrix();
}

// -----------------------------------------------------------------------------
MultiLevelStationaryVelocityTransformation
::MultiLevelStationaryVelocityTransformation(const MultiLevelStationaryVelocityTransformation &t)
:
  MultiLevelTransformation(t),
  _LogA(t._LogA)
{
  _GlobalTransformationObserver.Bind(
    ModifiedEvent,
    MakeDelegate(this, &MultiLevelStationaryVelocityTransformation::UpdateLogMatrix)
  );
  _GlobalTransformation.AddObserver(_GlobalTransformationObserver);
}

// ------------------------------------------------------------------------------
MultiLevelStationaryVelocityTransformation
::~MultiLevelStationaryVelocityTransformation()
{
  _GlobalTransformation.ClearObservers();
}

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation
::UpdateLogMatrix()
{
  _LogA = _GlobalTransformation.GetMatrix().Log();
}

// =============================================================================
// Transformation parameters (DoFs)
// =============================================================================

// -----------------------------------------------------------------------------
double MultiLevelStationaryVelocityTransformation::DOFGradientNorm(const double *gradient) const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return svffd ? svffd->DOFGradientNorm(gradient) : .0;
}

// -----------------------------------------------------------------------------
int MultiLevelStationaryVelocityTransformation::NumberOfDOFs() const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return svffd ? svffd->NumberOfDOFs() : 0;
}

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation::Put(int idx, double value)
{
  BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  if (svffd) {
    svffd->Put(idx, value);
    this->Changed(svffd->Changed());
  }
}

// -----------------------------------------------------------------------------
double MultiLevelStationaryVelocityTransformation::Get(int idx) const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return svffd ? svffd->Get(idx) : numeric_limits<double>::quiet_NaN();
}

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation::Put(const DOFValue *x)
{
  BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  if (svffd) {
    svffd->Put(x);
    this->Changed(svffd->Changed());
  }
}

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation::Add(const DOFValue *dx)
{
  BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  if (svffd) {
    svffd->Add(dx);
    this->Changed(svffd->Changed());
  }
}

// -----------------------------------------------------------------------------
double MultiLevelStationaryVelocityTransformation::Update(const DOFValue *dx)
{
  BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  double max_delta = .0;
  if (svffd) {
    max_delta = svffd->Update(dx);
    this->Changed(svffd->Changed());
  }
  return max_delta;
}

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation::Get(DOFValue *x) const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  if (svffd) svffd->Get(x);
}

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation::PutStatus(int idx, DOFStatus status)
{
  BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  if (svffd) {
    svffd->PutStatus(idx, status);
    this->Changed(svffd->Changed());
  }
}

// -----------------------------------------------------------------------------
MultiLevelStationaryVelocityTransformation::DOFStatus
MultiLevelStationaryVelocityTransformation::GetStatus(int idx) const
{
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  return svffd ? svffd->GetStatus(idx) : Passive;
}

// =============================================================================
// Levels
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation::CombineLocalTransformation()
{
  FreeFormTransformation *first = NULL, *second = NULL;
  while (_NumberOfLevels > 1) {
    first  = this->PopLocalTransformation();
    second = this->PopLocalTransformation();
    if (first->Attributes() == second->Attributes()) {
      for(int i = 0; i < first->NumberOfDOFs(); ++i) {
        second->Put(i, first->Get(i) + second->Get(i));
      }
    } else {
      cerr << this->NameOfClass() << "::CombineLocalTransformation: Only implemented for local transformations defined on same lattice" << endl;
      exit(1);
    }
    this->PushLocalTransformation(second);
    delete first;
  }
}

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation::MergeGlobalIntoLocalDisplacement()
{
  // Do nothing if global transformation is the identity
  if (_GlobalTransformation.IsIdentity()) return;

  BSplineFreeFormTransformationSV *ffd = SVFFD(0);

  // Approximate global transformation by SV FFD
  BSplineFreeFormTransformationSV global(ffd->Attributes());
  global.ApproximateAsNew(this->GetGlobalTransformation());

  // Add the calculated coefficients to the coefficients of the first FFD
  for (int i = 0; i < ffd->NumberOfDOFs(); ++i) {
    ffd->Put(i, ffd->Get(i) + global.Get(i));
  }

  // Reset matrix previously used for global transformation to identity
  _GlobalTransformation.Reset();
}

// =============================================================================
// Point transformations
// =============================================================================

// -----------------------------------------------------------------------------
template <class ScalarType>
void MultiLevelStationaryVelocityTransformation
::VelocityComponents(int m, int n, GenericImage<ScalarType> &v, bool coeff) const
{
  MIRTK_START_TIMING();
  if (n < 0 || n > _NumberOfLevels) n = _NumberOfLevels;
  const int l1 = (m < 0 ? 0 : m);
  if (v.T() != 3) {
    cerr << "MultiLevelStationaryVelocityTransformation::VelocityComponents: Vector field must have 3 components (_t)" << endl;
    exit(1);
  }
  const ImageAttributes &grid = v.Attributes();
  // Determine local SV FFDs whose coefficients can be added without prior deconvolution
  Array<bool> add_coeff(n, false);
  bool convert_to_coeff = false;
  if (coeff) {
    for (int l = l1; l < n; ++l) {
      add_coeff[l] = grid.EqualInSpace(SVFFD(l)->Attributes());
    }
  }
  // Add global velocities
  if (m < 0 && !_LogA.IsIdentity()) {
    ParallelForEachVoxel(EvaluateGlobalSVFFD3D(_LogA, &v), grid, v);
    convert_to_coeff = coeff;
  } else {
    v = .0;
  }
  // Add local velocities (which require deconvolution)
  for (int l = l1; l < n; ++l) {
    if (!add_coeff[l]) {
      ParallelForEachVoxel(AddBSplineSVFFD3D(SVFFD(l), &v), grid, v);
      convert_to_coeff = coeff;
    }
  }
  // Convert to cubic B-spline coefficients if needed and requested
  if (convert_to_coeff) ConvertToCubicBSplineCoefficients(v);
  // Add local velocities (which do not require deconvolution)
  for (int l = l1; l < n; ++l) {
    if (add_coeff[l]) {
      const BSplineFreeFormTransformationSV *svffd = SVFFD(l);
      ScalarType *vx = v.Data(0, 0, 0, 0);
      ScalarType *vy = v.Data(0, 0, 0, 1);
      ScalarType *vz = v.Data(0, 0, 0, 2);
      for (int k = 0; k < grid._z; ++k)
      for (int j = 0; j < grid._y; ++j)
      for (int i = 0; i < grid._x; ++i, ++vx, ++vy, ++vz) {
        const BSplineFreeFormTransformationSV::Vector &v = svffd->_CPImage(i, j, k);
        *vx += static_cast<ScalarType>(v._x);
        *vy += static_cast<ScalarType>(v._y);
        *vz += static_cast<ScalarType>(v._z);
      }
    }
  }
  MIRTK_DEBUG_TIMING(3, "MultiLevelStationaryVelocityTransformation::Velocity");
}

// -----------------------------------------------------------------------------
template <> void MultiLevelStationaryVelocityTransformation
::Velocity(int m, int n, GenericImage<float> &v, bool coeff) const
{
  VelocityComponents(m, n, v, coeff);
}

// -----------------------------------------------------------------------------
template <> void MultiLevelStationaryVelocityTransformation
::Velocity(int m, int n, GenericImage<double> &v, bool coeff) const
{
  VelocityComponents(m, n, v, coeff);
}

// =============================================================================
// Derivatives
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const GenericImage<double> *wc,
                     double t, double t0, double w) const
{
  MIRTK_START_TIMING();
  typedef BSplineFreeFormTransformationSV::CPImage CPImage;
  // Gradient w.r.t parameters of last active transformation level
  const BSplineFreeFormTransformationSV *svffd = ActiveSVFFD();
  if (!svffd) return;
  // Upper integration limit for given interval
  const double tau = UpperIntegrationLimit(t, t0);
  // Get current SV MFFD velocity
  CPImage v(svffd->Attributes());
  Velocity(-1, this->NumberOfLevels(), v, true);
  // Apply chain rule to compute gradient w.r.t control point displacements
  CPImage d(svffd->Attributes());
  DOFValue * const grd = reinterpret_cast<DOFValue *>(d.Data());
  svffd->BSplineFreeFormTransformation3D::ParametricGradient(in, grd, i2w, wc, .0, 1.0);
  // Convert displacement gradient to velocity update
  svffd->EvaluateBCHFormula(svffd->NumberOfBCHTerms(), d, tau, v, 1., d, true);
  // Adjust weight as update field is computed for tau * v, i.e.,
  //   exp(tau * v_{i+1}) = exp(tau v_i) o exp(\delta u)
  //   ==> v_{i+1} = log(exp(tau * v_{i+1})) / tau
  w /= tau;
  // Add weighted gradient to total energy gradient
  for (int dof = 0; dof < svffd->NumberOfDOFs(); ++dof) out[dof] += w * grd[dof];
  MIRTK_DEBUG_TIMING(2, "parametric gradient computation (SV MFFD)");
}

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation
::ParametricGradient(const GenericImage<double> *in, double *out,
                     const WorldCoordsImage *i2w, const WorldCoordsImage *wc,
                     double t0, double w) const
{
  this->ParametricGradient(in, out, i2w, wc, in->GetTOrigin(), t0, w);
}

// =============================================================================
// Others
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation::Invert()
{
  _GlobalTransformation.Invert();
  for (int n = 0; n < this->NumberOfLevels(); ++n) SVFFD(n)->Invert();
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void MultiLevelStationaryVelocityTransformation::Print(ostream &os, Indent indent) const
{
  os << indent << "Multi-level SV FFD:" << endl;
  MultiLevelTransformation::Print(os, indent + 1);
}


} // namespace mirtk
