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

#include <Transformation.h>

#include <DisplacementToVelocityField.h>


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
    npoints += disp[n]->GetX() * disp[n]->GetY() * disp[n]->GetZ();
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
  dtov.SetSmoothVelocities  (smooth);
  dtov.SetNumberO