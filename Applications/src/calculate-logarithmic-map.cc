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

#include "mirtk/Common.h"
#include "mirtk/Options.h"

#include "mirtk/IOConfig.h"

#include "mirtk/VelocityToDisplacementField.h"
#include "mirtk/VelocityToDisplacementFieldEuler.h"
#include "mirtk/VelocityToDisplacementFieldSS.h"
#include "mirtk/DisplacementToVelocityFieldBCH.h"

using namespace mirtk;


// =============================================================================
// Global variables with default arguments
// =============================================================================

char IntegrationMethod[32]  = {"SS"};
int  NumberOfSteps          = 64;
int  NumberOfBCHSteps       = 8;
int  NumberOfBCHTerms       = 6;
bool SmoothBCHApproximation = false;
bool UseJacobian            = false;

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <disp> <velo> [options]" << endl;
  cout << "       " << name << " <dx> <dy> [<dz>] <vx> <vy> [<vz>] [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Reads a dense 3D (2D) displacement field, computes the corresponding" << endl;
  cout << "  stationary 3D (2D) velocity field, and writes the resulting vector field." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -ss             Use scaling and squaring (SS) method for exponentiation. (default)" << endl;
  cout << "  -euler          Use forward Euler method for exponentiation." << endl;
  cout << "  -steps <int>    Number of integration steps (2^n in case of SS). "
                             << "(default: " << NumberOfSteps << ")" << endl;
  cout << "  -terms <int>    Number of BCH approximation terms (either 2 or 3). (default: " << NumberOfBCHTerms << ")" << endl;
  cout << "  -iters <int>    Number of BCH update steps. (default: " << NumberOfBCHSteps << ")" << endl;
  cout << "  -smooth         Smooth velocities before each update. (default: off)" << endl;
  cout << "  -jac            Use Jacobian for Lie bracket computation. (default: off)" << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Print help or version and exit if requested
  HANDLE_HELP_OR_VERSION();

  InitializeIOLibrary();

  // Parse positional arguments
  const char *disp_fname = NULL;
  const char *velo_fname = NULL;
  const char *dx_fname   = NULL;
  const char *dy_fname   = NULL;
  const char *dz_fname   = NULL;
  const char *vx_fname   = NULL;
  const char *vy_fname   = NULL;
  const char *vz_fname   = NULL;

  int posargc = 0;
  while (posargc+1 < argc && argv[posargc+1][0] != '-') posargc++;

  switch (posargc) {
    case 2:
      disp_fname = argv[1];
      velo_fname = argv[2];
      break;
    case 4:
      dx_fname = argv[1];
      dy_fname = argv[2];
      vx_fname = argv[3];
      vy_fname = argv[4];
      break;
    case 6:
      dx_fname = argv[1];
      dy_fname = argv[2];
      dz_fname = argv[3];
      vx_fname = argv[4];
      vy_fname = argv[5];
      vz_fname = argv[6];
      break;
    default:
      PrintHelp(EXECNAME);
      exit(1);
  }

  // Parse optional arguments
  for (ARGUMENTS_AFTER(posargc)) {
    if      (OPTION("-steps" )) NumberOfSteps          = atoi(ARGUMENT);
    else if (OPTION("-iters" )) NumberOfBCHSteps       = atoi(ARGUMENT);
    else if (OPTION("-terms" )) NumberOfBCHTerms       = atoi(ARGUMENT);
    else if (OPTION("-smooth")) SmoothBCHApproximation = true;
    else if (OPTION("-jac"))    UseJacobian            = true;
    else if (OPTION("-int")) {
      const char *arg = ARGUMENT;
      size_t len = strlen(arg);
      if (len > 99u) len = 99u;
      for (size_t i = 0u; i < len; ++i) {
        if ('a' <= arg[i] && arg[i] <= 'z') IntegrationMethod[i] = 'A' + (arg[i] - 'a');
        else                                IntegrationMethod[i] = arg[i];
      }
      IntegrationMethod[len] = '\0';
    } else if (OPTION("-ss")) {
      IntegrationMethod[0] = 'S';
      IntegrationMethod[1] = 'S';
      IntegrationMethod[2] = '\0';
    } else if (OPTION("-euler")) {
      IntegrationMethod[0] = 'R';
      IntegrationMethod[1] = 'K';
      IntegrationMethod[2] = 'E';
      IntegrationMethod[3] = '1';
      IntegrationMethod[4] = '\0';
    } else {
      HANDLE_COMMON_OR_UNKNOWN_OPTION();
    }
  }

  // Read input displacements
  UniquePtr<RealImage> d;

  if (disp_fname) {
    d.reset(new RealImage(disp_fname));
  } else {
    RealImage dx(dx_fname);
    RealImage dy(dy_fname);
    UniquePtr<RealImage> dz(dz_fname ? new RealImage(dz_fname) : nullptr);

    ImageAttributes attr = dx.Attributes();
    if (attr._t > 1) {
      FatalError("Separate input component images must be scalar displacement fields.\n"
                 "Try calling this program with single vector field input and output images instead.");
    }
    if (dy.Attributes() != attr || (dz.get() && dz->Attributes() != attr)) {
      FatalError("Image attributes of input displacement component images do not match!");
    }

    d.reset(new RealImage(attr, dz.get() ? 3 : 2));
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      d->Put(i, j, k, 0, dx.Get(i, j, k));
      d->Put(i, j, k, 1, dy.Get(i, j, k));
      if (dz.get()) d->Put(i, j, k, 2, dz->Get(i, j, k));
    }
  }

  // Instantiate exponential filter
  UniquePtr<VelocityToDisplacementField<RealPixel> > exp;
  if (strcmp(IntegrationMethod, "SS") == 0) {
    exp.reset(new VelocityToDisplacementFieldSS<RealPixel>());
  } else if (strcmp(IntegrationMethod, "RKE1")         == 0 ||
             strcmp(IntegrationMethod, "FORWARDEULER") == 0 ||
             strcmp(IntegrationMethod, "EULER")        == 0) {
    exp.reset(new VelocityToDisplacementFieldEuler<RealPixel>());
  } else {
    FatalError("Unknown integration method: " << IntegrationMethod);
  }

  // Compute stationary velocity field
  DisplacementToVelocityFieldBCH<RealPixel> dtov;
  RealImage                                 v;

  dtov.ExponentialFilter (exp.get());
  dtov.NumberOfSteps     (NumberOfSteps);
  dtov.NumberOfIterations(NumberOfBCHSteps);
  dtov.NumberOfTerms     (NumberOfBCHTerms);
  dtov.SmoothVelocities  (SmoothBCHApproximation);
  dtov.UseJacobian       (UseJacobian);

  dtov.Input(d.get());
  dtov.Output(&v);

  dtov.Run();

  // Save output velocities
  if (velo_fname) {
    v.Write(velo_fname);
  } else {
    v.GetFrame(0).Write(vx_fname);
    v.GetFrame(1).Write(vy_fname);
    if (vz_fname) v.GetFrame(2).Write(vz_fname);
  }

  return 0;
}
