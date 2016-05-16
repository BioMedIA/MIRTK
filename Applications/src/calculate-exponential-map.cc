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

using namespace mirtk;


// =============================================================================
// Global variables with default arguments
// =============================================================================

char IntegrationMethod[32] = {"SS"};
int  NumberOfSteps         = 128;

// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <velo> <disp> [options]" << endl;
  cout << "       " << name << " <vx> <vy> [<vz>] <dx> <dy> [<dz>] [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Reads a dense 3D (2D) stationary velocity field, computes the corresponding" << endl;
  cout << "  dense 3D (2D) displacement field, and writes the resulting vector field." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -ss             Use scaling and squaring (SS) method. (default)" << endl;
  cout << "  -euler          Use forward Euler method." << endl;
  cout << "  -steps <int>    Number of integration steps (2^n in case of SS). "
                          << "(default: " << NumberOfSteps << ")" << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Print help or version and exit if requested
  HANDLE_HELP_OR_VERSION();

  InitializeIOLibrary();

  // Parse positional arguments
  const char *velo_fname = NULL;
  const char *disp_fname = NULL;
  const char *vx_fname   = NULL;
  const char *vy_fname   = NULL;
  const char *vz_fname   = NULL;
  const char *dx_fname   = NULL;
  const char *dy_fname   = NULL;
  const char *dz_fname   = NULL;

  int posargc = 0;
  while (posargc+1 < argc && argv[posargc+1][0] != '-') posargc++;

  switch (posargc) {
    case 2:
      velo_fname = argv[1];
      disp_fname = argv[2];
      break;
    case 4:
      vx_fname = argv[1];
      vy_fname = argv[2];
      dx_fname = argv[3];
      dy_fname = argv[4];
      break;
    case 6:
      vx_fname = argv[1];
      vy_fname = argv[2];
      vz_fname = argv[3];
      dx_fname = argv[4];
      dy_fname = argv[5];
      dz_fname = argv[6];
      break;
    default:
      PrintHelp(EXECNAME);
      exit(1);
  }

  // Parse optional arguments
  for (ARGUMENTS_AFTER(posargc)) {
    if (OPTION("-int")) {
      const char *arg = ARGUMENT;
      size_t len = strlen(arg);
      if (len > 99u) len = 99u;
      for (size_t i = 0u; i < len; i++) {
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
    } else if (OPTION("-steps")) {
      NumberOfSteps = atoi(ARGUMENT);
    } else {
      HANDLE_COMMON_OR_UNKNOWN_OPTION();
    }
  }

  // Read input velocities
  UniquePtr<RealImage> v;

  if (velo_fname) {
    v.reset(new RealImage(velo_fname));
  } else {
    RealImage vx(vx_fname);
    RealImage vy(vy_fname);
    UniquePtr<RealImage> vz(vz_fname ? new RealImage(vz_fname) : nullptr);

    ImageAttributes attr = vx.Attributes();
    if (attr._t > 1) {
      FatalError("Separate input component images must be scalar velocity fields.\n"
                 "Try calling this program with single vector field input and output images instead.");
    }
    if (vy.Attributes() != attr || (vz.get() && vz->Attributes() != attr)) {
      FatalError("Image attributes of input velocity component images do not match!");
    }
    attr._t = (vz.get() ? 3 : 2);

    v.reset(new RealImage(attr));
    for (int k = 0; k < attr._z; ++k)
    for (int j = 0; j < attr._y; ++j)
    for (int i = 0; i < attr._x; ++i) {
      v->Put(i, j, k, 0, vx.Get(i, j, k));
      v->Put(i, j, k, 1, vy.Get(i, j, k));
      if (vz.get()) v->Put(i, j, k, 2, vz->Get(i, j, k));
    }
  }

  // Instantiate filter which implements the desired integration method
  UniquePtr<VelocityToDisplacementField<RealPixel> > vtod;

  if (strcmp(IntegrationMethod, "SS") == 0) {
    vtod.reset(new VelocityToDisplacementFieldSS<RealPixel>());
  } else if (strcmp(IntegrationMethod, "RKE1") == 0 || strcmp(IntegrationMethod, "FORWARDEULER") == 0 || strcmp(IntegrationMethod, "EULER") == 0) {
    vtod.reset(new VelocityToDisplacementFieldEuler<RealPixel>());
  } else {
    FatalError("Unknown integration method: " << IntegrationMethod);
  }

  // Compute displacement field
  RealImage d;

  vtod->Input(v.get());
  vtod->Output(&d);
  vtod->NumberOfSteps(NumberOfSteps);
  vtod->Run();

  // Save output velocities
  if (disp_fname) {
    d.Write(disp_fname);
  } else {
    d.GetFrame(0).Write(dx_fname);
    d.GetFrame(1).Write(dy_fname);
    if (dz_fname) d.GetFrame(2).Write(dz_fname);
  }

  return 0;
}
