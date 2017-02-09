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

#include "mirtk/LinearFreeFormTransformationTD.h"

#include "mirtk/Math.h"
#include "mirtk/BSplineFreeFormTransformationTD.h"
#include "mirtk/DisplacementToVelocityFieldBCH.h"


namespace mirtk {


// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
LinearFreeFormTransformationTD::LinearFreeFormTransformationTD()
:
  _MinTimeStep(0.01),
  _MaxTimeStep(0.5)
{
  _ExtrapolationMode = Extrapolation_NN;
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformationTD
::LinearFreeFormTransformationTD(const ImageAttributes &attr,
                                 double dx, double dy, double dz, double dt)
:
  _MinTimeStep(0.01),
  _MaxTimeStep(0.5)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(attr, dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformationTD
::LinearFreeFormTransformationTD(const BaseImage &target,
                                 double dx, double dy, double dz, double dt)
:
  _MinTimeStep(0.01),
  _MaxTimeStep(0.5)
{
  _ExtrapolationMode = Extrapolation_NN;
  Initialize(target.Attributes(), dx, dy, dz, dt);
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformationTD
::LinearFreeFormTransformationTD(const BSplineFreeFormTransformationTD &t)
:
  LinearFreeFormTransformation4D(t),
  _MinTimeStep(t.MinTimeStep()),
  _MaxTimeStep(t.MaxTimeStep())
{
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformationTD
::LinearFreeFormTransformationTD(const LinearFreeFormTransformationTD &t)
:
  LinearFreeFormTransformation4D(t),
  _MinTimeStep(t._MinTimeStep),
  _MaxTimeStep(t._MaxTimeStep)
{
}

// -----------------------------------------------------------------------------
LinearFreeFormTransformationTD
::~LinearFreeFormTransformationTD()
{
}

// =============================================================================
// Approximation/Interpolation
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformationTD
::ApproximateDOFs(const GenericImage<double> * const *disp,
                  const double *t1, const double *t2, int no,
                  bool smooth, int nterms, int niter)
{
  if (no == 0) return;
 
  // Get last displacement field
  int last = 0;
  for (int n = 1; n < no; n++) {
    if (t2[n] > t2[last]) last = n;
  }

  // Determine number of 3D+t samples
  int npoints = 0;
  for (int n = 0; n < no; n++) {
    npoints += disp[n]->X() * disp[n]->Y() * disp[n]->Z();
  }

  // Allocate memory for computed stationary velocity fields
  double *px = Allocate<double>(npoints);
  double *py = Allocate<double>(npoints);
  double *pz = Allocate<double>(npoints);
  double *pt = Allocate<double>(npoints);
  double *vx = Allocate<double>(npoints);
  double *vy = Allocate<double>(npoints);
  double *vz = Allocate<double>(npoints);

  // Compute stationary velocity fields from displacement fields
  // and store them in 3D+t vector field with associated time t
  DisplacementToVelocityFieldBCH<double> dtov;
  dtov.SmoothVelocities  (smooth);
  dtov.NumberOfTerms     (nterms);
  dtov.NumberOfIterations(niter);

  int p = 0;
  for (int n = 0; n < no; n++) {
    GenericImage<double> v;
    dtov.UpperIntegrationLimit(t2[n] - t1[n]);
    dtov.NumberOfSteps(iround(dtov.UpperIntegrationLimit() / _MinTimeStep));
    dtov.Input(const_cast<GenericImage<double> *>(disp[n]));
    dtov.Output(&v);
    dtov.Run();
    for (int k = 0; k < v.Z(); k++) {
      for (int j = 0; j < v.Y(); j++) {
        for (int i = 0; i < v.X(); i++) {
          px[p] = i;
          py[p] = j;
          pz[p] = k;
          v.ImageToWorld(px[p], py[p], pz[p]);
          pt[p] = t1[n];
          vx[p] = v.Get(i, j, k, 0);
          vy[p] = v.Get(i, j, k, 1);
          vz[p] = v.Get(i, j, k, 2);
          p++;
        }
        if (n == last) {
          px[p] = px[p - 1];
          py[p] = py[p - 1];
          pz[p] = pz[p - 1];
          pt[p] = t2[last ];
          vx[p] = vx[p - 1];
          vy[p] = vy[p - 1];
          vz[p] = vz[p - 1];
          p++;
        }
      }
    }
  }

  // Approximate 3D+t velocity field by linear FFD
  LinearFreeFormTransformation4D::ApproximateDOFs(px, py, pz, pt, vx, vy, vz, npoints);

  // Free memory
  Deallocate(px);
  Deallocate(py);
  Deallocate(pz);
  Deallocate(pt);
  Deallocate(vx);
  Deallocate(vy);
  Deallocate(vz);
}

// -----------------------------------------------------------------------------
void LinearFreeFormTransformationTD
::ApproximateDOFs(const double *, const double *, const double *, const double *,
                  const double *, const double *, const double *, int)
{
  cerr << this->NameOfClass() << "::ApproximateDOFs: Not implemented" << endl;
  cerr << "  --> Try the other overloaded approximation function" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
void LinearFreeFormTransformationTD
::ApproximateDOFsGradient(const double *, const double *, const double *, const double *,
                          const double *, const double *, const double *, int,
                          double *, double) const
{
  cerr << this->NameOfClass() << "::ApproximateDOFsGradient: Not implemented" << endl;
  exit(1);
}

// -----------------------------------------------------------------------------
double LinearFreeFormTransformationTD
::ApproximateAsNew(GenericImage<double> **disp,
                   const double *t1, const double *t2, int no,
                   bool smooth, int nterms, int niter)
{
  // Approximate 3D+t velocity field
  this->ApproximateDOFs(disp, t1, t2, no, smooth, nterms, niter);

  // Evaluate RMS of approximation error
  double error  = 0;
  int    ntotal = 0;
  double dx, dy, dz;

  for (int n = 0; n < no; n++) {
    GenericImage<double> &d = *disp[n];
    for (int k = 0; k < d.Z(); ++k)
    for (int j = 0; j < d.Y(); ++j)
    for (int i = 0; i < d.X(); ++i) {
      dx = i, dy = j, dz = k;
      d.ImageToWorld(dx, dy, dz);
      this->Displacement(dx, dy, dz, t1[n], t2[n]);
      d(i, j, k, 0) -= dx;
      d(i, j, k, 1) -= dy;
      d(i, j, k, 2) -= dz;
      error += sqrt(d(i, j, k, 0) * d(i, j, k, 0) +
                         d(i, j, k, 1) * d(i, j, k, 1) +
                         d(i, j, k, 2) * d(i, j, k, 2));
      ++ntotal;
    }
  }
  if (ntotal > 0) error /= ntotal;
  return error;
}

// -----------------------------------------------------------------------------
void LinearFreeFormTransformationTD::Interpolate(const double *, const double *, const double *)
{
  cerr << this->NameOfClass() << "::Interpolate: Not implemented" << endl;
  exit(1);
}

// =============================================================================
// I/O
// =============================================================================

// -----------------------------------------------------------------------------
void LinearFreeFormTransformationTD::Print(ostream &os, Indent indent) const
{
  os << indent << "Linear TD FFD:" << endl;
  ++indent;
  FreeFormTransformation4D::Print(os, indent);
  os << indent << "Minimum length of steps: " << _MinTimeStep << endl;
  os << indent << "Maximum length of steps: " << _MaxTimeStep << endl;
}

// -----------------------------------------------------------------------------
bool LinearFreeFormTransformationTD::CanRead(TransformationType format) const
{
  switch (format) {
    case TRANSFORMATION_LINEAR_FFD_TD_v1:
    case TRANSFORMATION_LINEAR_FFD_TD_v2:
    case TRANSFORMATION_LINEAR_FFD_TD_v3:
      return true;
    default:
      return false;
  }
}

// -----------------------------------------------------------------------------
Cifstream &LinearFreeFormTransformationTD::ReadDOFs(Cifstream &from, TransformationType format)
{
  // Read FFD data
  if (format < TRANSFORMATION_LINEAR_FFD_TD_v2) {
    LinearFreeFormTransformation4D::ReadDOFs(from, TRANSFORMATION_LINEAR_FFD_4D_v1);
  } else {
    LinearFreeFormTransformation4D::ReadDOFs(from, TRANSFORMATION_LINEAR_FFD_4D_v2);
  }

  // Read minimum/maximum time step length
  from.ReadAsDouble(&_MinTimeStep, 1);
  from.ReadAsDouble(&_MaxTimeStep, 1);

  return from;
}

// -----------------------------------------------------------------------------
Cofstream &LinearFreeFormTransformationTD::WriteDOFs(Cofstream &to) const
{
  // Write FFD data
  LinearFreeFormTransformation4D::WriteDOFs(to);

  // Write minimum/maximum time step length
  to.WriteAsDouble(&_MinTimeStep, 1);
  to.WriteAsDouble(&_MaxTimeStep, 1);

  return to;
}


} // namespace mirtk
