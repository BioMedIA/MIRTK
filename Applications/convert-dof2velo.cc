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

#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkTransformation.h>
#include <mirtkFFDIntegrationMethod.h>
#include <mirtkBSplineFreeFormTransformationSV.h>
#include <mirtkBSplineFreeFormTransformationTD.h>
#include <mirtkLinearFreeFormTransformationTD.h>
#include <mirtkHomogeneousTransformation.h>
#include <mirtkFreeFormTransformation.h>
#include <mirtkMultiLevelTransformation.h>
#include <mirtkMultiLevelFreeFormTransformation.h>

using namespace mirtk;


// =============================================================================
// Global variables
// =============================================================================

FFDIntegrationMethod IntegrationMethod           = FFDIM_RKE2;
int                  MinNumberOfIntegrationSteps = 10;
int                  MaxNumberOfIntegrationSteps = 100;
double               Tolerance                   = 1.0e-3;
int                  NumberOfBCHSteps            = 8;
int                  NumberOfBCHTerms            = 3;
bool                 SmoothBCHApproximation      = false;

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "Usage: " << name << " <dofin> <dofout> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Reads a 3D transformation and writes a 3D(+t) free-form transformation" << endl;
  cout << "  which is parameterized by a (non-)stationary velocity field and approximates" << endl;
  cout << "  the displacements of the input transformation." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  dofin             Input transformation." << endl;
  cout << "  dofout            Output transformation parameterized by velocities." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -target <image>   Target image. Required in case of non-FFD input transformation. (default: none)" << endl;
  cout << "  -t1 <double>      Time point of first frame. (default: temporal origin of FFD/target)" << endl;
  cout << "  -t2 <double>      Time point of second frame. (default: t1 + 1)" << endl;
  cout << "  -stationary       Use stationary velocity field. (default: off)" << endl;
  cout << "  -int <string>     Integration method: RKE1, RKE2, RKBS2, RK4,... (default: " << ToString(IntegrationMethod) << ")" << endl;
  cout << "  -minsteps <int>   Minimum number of integration steps. (default: " << MinNumberOfIntegrationSteps << ")" << endl;
  cout << "  -maxsteps <int>   Maximum number of integration steps. (default: " << MaxNumberOfIntegrationSteps << ")" << endl;
  cout << "  -steps <int>      Number of integration steps. Sets both :option:`-minsteps` and :option:`-maxsteps`." << endl;
  cout << "  -tol <double>     Error tolerance of integration method with adaptive step size. (default: " << Tolerance << ")" << endl;
  cout << "  -terms <int>      Number of BCH approximation terms (either 2 or 3). (default: " << NumberOfBCHTerms << ")" << endl;
  cout << "  -iters <int>      Number of BCH update steps. (default: " << NumberOfBCHSteps << ")" << endl;
  cout << "  -smooth           Smooth velocities before each update. (default: off)" << endl;
  cout << "  -merge            Merge global into local displacements. (default: off)" << endl;
  cout << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
void read_dofin(const char *name_in, const char *target_name, Matrix &global, GenericImage<double> &local, bool merge)
{
  global.Initialize(4, 4);
  global.Ident();

  unique_ptr<Transformation> T(Transformation::New(name_in));
  MultiLevelTransformation  *mffd = dynamic_cast<MultiLevelTransformation *> (T.get());
  FreeFormTransformation    *ffd  = dynamic_cast<FreeFormTransformation *>   (T.get());
  HomogeneousTransformation *aff  = dynamic_cast<HomogeneousTransformation *>(T.get());
  double                     dx, dy, dz;
  ImageAttributes            attr;

  if (mffd != nullptr) {
    ffd = mffd->GetLocalTransformation(mffd->NumberOfLevels() - 1);
    if (merge) {
      mffd->MergeGlobalIntoLocalDisplacement();
    } else {
      global = mffd->GetGlobalTransformation()->GetMatrix();
    }
  }
  if (ffd != nullptr) {
    attr = ffd->Attributes();
    if (attr._t > 1) {
      FatalError("Input FFD " << name_in << " is 4D, but must be 3D.");
    }
    attr._t = 3;
    local.Initialize(attr);
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      dx = i;
      dy = j;
      dz = k;
      local.ImageToWorld(dx, dy, dz);
      ffd->LocalDisplacement(dx, dy, dz);
      local(i, j, k, 0) = dx;
      local(i, j, k, 1) = dy;
      local(i, j, k, 2) = dz;
    }
  } else if (aff != NULL && target_name != NULL) {
    GreyImage target(target_name);
    attr    = target.Attributes();
    attr._t = 3;
    local.Initialize(attr);
    if (merge) {
      for (int k = 0; k < attr._z; ++k)
      for (int j = 0; j < attr._y; ++j)
      for (int i = 0; i < attr._x; ++i) {
        dx = i;
        dy = j;
        dz = k;
        local.ImageToWorld(dx, dy, dz);
        aff->GlobalDisplacement(dx, dy, dz);
        local(i, j, k, 0) = dx;
        local(i, j, k, 1) = dy;
        local(i, j, k, 2) = dz;
      }
    } else {
      global = aff->GetMatrix();
    }
  } else {
    FatalError("Input transformation is no FFD. Specify target image using -target option to define the domain of the output FFD.");
  }
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(2);

  // parse arguments
  const char *dofin      = POSARG(1);
  const char *dofout     = POSARG(2);
  const char *target     = nullptr;
  double      t1         = numeric_limits<double>::quiet_NaN();
  double      t2         = numeric_limits<double>::quiet_NaN();
  bool        stationary = false;
  bool        linear     = false;
  bool        merge      = false;

  for (ALL_OPTIONS) {
    if      (OPTION("-target")) target = ARGUMENT;
    else if (OPTION("-t1")) PARSE_ARGUMENT(t1);
    else if (OPTION("-t2")) PARSE_ARGUMENT(t2);
    else if (OPTION("-int")) PARSE_ARGUMENT(IntegrationMethod);
    else if (OPTION("-minsteps")) PARSE_ARGUMENT(MinNumberOfIntegrationSteps);
    else if (OPTION("-maxsteps")) PARSE_ARGUMENT(MaxNumberOfIntegrationSteps);
    else if (OPTION("-steps")) {
      PARSE_ARGUMENT(MinNumberOfIntegrationSteps);
      MaxNumberOfIntegrationSteps = MinNumberOfIntegrationSteps;
    }
    else if (OPTION("-tol")) PARSE_ARGUMENT(Tolerance);
    else if (OPTION("-iters")) PARSE_ARGUMENT(NumberOfBCHSteps);
    else if (OPTION("-terms")) PARSE_ARGUMENT(NumberOfBCHTerms);
    else if (OPTION("-smooth")) SmoothBCHApproximation = true;
    else if (OPTION("-merge")) merge = true;
    else if (OPTION("-stationary")) stationary = true;
    else if (OPTION("-linear")) linear = true;
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  // Read input FFD
  Matrix               global;
  GenericImage<double> local;

  read_dofin(dofin, target, global, local, merge);

  // Attributes of output FFD
  ImageAttributes attr = local.GetImageAttributes();

  if (IsNaN(t1)) t1 = local.ImageToTime(.0);
  if (IsNaN(t2)) t2 = t1 + 1.0;

  attr._torigin = t1;
  attr._dt      = t2 - t1;

  // Approximate (local) displacements by velocities
  MultiLevelFreeFormTransformation mffd;

  mffd.GetGlobalTransformation()->PutMatrix(global);
  GenericImage<double> *disp = new GenericImage<double>(local);

  if (stationary) {
    if (linear) {
      delete disp;
      FatalError("LinearFreeFormTransformationSV not implemented!");
    }
    attr._t = 1;
    BSplineFreeFormTransformationSV *svffd = new BSplineFreeFormTransformationSV(attr, attr._dx, attr._dy, attr._dz);
    svffd->NumberOfSteps(MinNumberOfIntegrationSteps);
    svffd->ApproximateAsNew(*disp, SmoothBCHApproximation, NumberOfBCHTerms, NumberOfBCHSteps);
    mffd.PushLocalTransformation(svffd);
  } else if (linear) {
    attr._t = 2;
    LinearFreeFormTransformationTD *tdffd = new LinearFreeFormTransformationTD(attr, attr._dx, attr._dy, attr._dz, attr._dt);
    tdffd->MinTimeStep((t2 - t1) / MaxNumberOfIntegrationSteps);
    tdffd->MaxTimeStep((t2 - t1) / MinNumberOfIntegrationSteps);
    tdffd->ApproximateAsNew(&disp, &t1, &t2, 1, SmoothBCHApproximation, NumberOfBCHTerms, NumberOfBCHSteps);
    mffd.PushLocalTransformation(tdffd);
  } else {
    attr._t = 2;
    BSplineFreeFormTransformationTD *tdffd = new BSplineFreeFormTransformationTD(attr, attr._dx, attr._dy, attr._dz, attr._dt);
    tdffd->IntegrationMethod(IntegrationMethod);
    tdffd->MinTimeStep((t2 - t1) / MaxNumberOfIntegrationSteps);
    tdffd->MaxTimeStep((t2 - t1) / MinNumberOfIntegrationSteps);
    tdffd->Tolerance(Tolerance);
    tdffd->ApproximateAsNew(&disp, &t1, &t2, 1, SmoothBCHApproximation, NumberOfBCHTerms, NumberOfBCHSteps);
    mffd.PushLocalTransformation(tdffd);
  }

  delete disp;

  // Write output FFD
  mffd.Write(dofout);

  // Evaluate error over whole image domain
  double dx1, dy1, dz1;
  double dx2, dy2, dz2;
  double error, mag;
  double avg_dispin  = 0;
  double max_dispin  = 0;
  double avg_dispout = 0;
  double max_dispout = 0;
  double avg_error   = 0;
  double max_error   = 0;

  for (int k = 0; k < local.GetZ(); k++) {
    for (int j = 0; j < local.GetY(); j++) {
      for (int i = 0; i < local.GetX(); i++) {
        // Get local input FFD
        dx1 = local(i, j, k, 0);
        dy1 = local(i, j, k, 1);
        dz1 = local(i, j, k, 2);
        mag = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
        avg_dispin += mag;
        if (mag > max_dispin) max_dispin = mag;
        // World coordinates of control point
        dx2 = i;
        dy2 = j;
        dz2 = k;
        local.ImageToWorld(dx2, dy2, dz2);
        // Local displacement
        mffd.GetLocalTransformation(0)->Displacement(dx2, dy2, dz2, t1, t2);
        mag = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
        avg_dispout += mag;
        if (mag > max_dispout) max_dispout = mag;
        // Compute error
        error = sqrt((dx2 - dx1) * (dx2 - dx1) +
                     (dy2 - dy1) * (dy2 - dy1) +
                     (dz2 - dz1) * (dz2 - dz1));
        avg_error += error;
        if (error > max_error) max_error = error;
      }
    }
  }
  avg_dispout /= static_cast<double>(local.GetX() * local.GetY() * local.GetZ());
  avg_dispin  /= static_cast<double>(local.GetX() * local.GetY() * local.GetZ());
  avg_error   /= static_cast<double>(local.GetX() * local.GetY() * local.GetZ());

  cout << "Approximation error (incl. boundary):" << endl;
  cout << "  Average input displacement:  " << avg_dispin  << endl;
  cout << "  Maximum input displacement:  " << max_dispin  << endl;
  cout << "  Average output displacement: " << avg_dispout << endl;
  cout << "  Maximum output displacement: " << max_dispout << endl;
  cout << "  Average RMS error:           " << avg_error  << endl;
  cout << "  Maximum RMS error:           " << max_error  << endl;

  // Evaluate error ignoring boundary
  avg_dispin  = 0;
  max_dispin  = 0;
  avg_dispout = 0;
  max_dispout = 0;
  avg_error   = 0;
  max_error   = 0;

  int mini = 2;
  int maxi = local.GetX() - 3;
  int minj = 2;
  int maxj = local.GetY() - 3;
  int mink = 2;
  int maxk = local.GetZ() - 3;

  if (mini <= maxi || minj <= maxj || mink <= maxk) {
    for (int k = mink; k <= maxk; k++) {
      for (int j = minj; j <= maxj; j++) {
        for (int i = mini; i <= maxi; i++) {
          // Get local input FFD
          dx1 = local(i, j, k, 0);
          dy1 = local(i, j, k, 1);
          dz1 = local(i, j, k, 2);
          mag = sqrt(dx1 * dx1 + dy1 * dy1 + dz1 * dz1);
          avg_dispin += mag;
          if (mag > max_dispin) max_dispin = mag;
          // World coordinates of control point
          dx2 = i;
          dy2 = j;
          dz2 = k;
          local.ImageToWorld(dx2, dy2, dz2);
          // Local displacement
          mffd.GetLocalTransformation(0)->Displacement(dx2, dy2, dz2, t1, t2);
          mag = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
          avg_dispout += mag;
          if (mag > max_dispout) max_dispout = mag;
          // Compute error
          error = sqrt((dx2 - dx1) * (dx2 - dx1) +
                       (dy2 - dy1) * (dy2 - dy1) +
                       (dz2 - dz1) * (dz2 - dz1));
          avg_error += error;
          if (error > max_error) max_error = error;
        }
      }
    }
    avg_dispout /= static_cast<double>((maxi - mini) * (maxj - minj) * (maxk - mink));
    avg_dispin  /= static_cast<double>((maxi - mini) * (maxj - minj) * (maxk - mink));
    avg_error   /= static_cast<double>((maxi - mini) * (maxj - minj) * (maxk - mink));

    cout << endl;
    cout << "Approximation error (excl. boundary):" << endl;
    cout << "  Average input displacement:  " << avg_dispin  << endl;
    cout << "  Maximum input displacement:  " << max_dispin  << endl;
    cout << "  Average output displacement: " << avg_dispout << endl;
    cout << "  Maximum output displacement: " << max_dispout << endl;
    cout << "  Average RMS error:           " << avg_error  << endl;
    cout << "  Maximum RMS error:           " << max_error  << endl;
  }

  return 0;
}
