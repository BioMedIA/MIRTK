/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
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
#include "mirtk/Transformations.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <dofin> <dofout> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Inverts any transformation. In case of a non-rigid transformation\n";
  cout << "  the output transformation only approximates the true inverse.\n";
  cout << "  When the inverse mapping is not defined at a given point, the\n";
  cout << "  output transformation at this point depends on the interpolation\n";
  cout << "  numerical approximate solution found.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  dofin    Rigid, affine, or non-rigid input transformation.\n";
  cout << "  dofout   Inverse of input transformation. (approximation)\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Parse arguments
  REQUIRES_POSARGS(2);

  const char *dofin_name  = POSARG(1);
  const char *dofout_name = POSARG(2);

  for (ALL_OPTIONS) {
    HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read transformation
  UniquePtr<Transformation> dof(Transformation::New(dofin_name));

  HomogeneousTransformation                  *lin    = NULL;
  BSplineFreeFormTransformationSV            *svffd  = NULL;
  BSplineFreeFormTransformationTD            *tdffd  = NULL;
  MultiLevelFreeFormTransformation           *mffd   = NULL;
  MultiLevelStationaryVelocityTransformation *msvffd = NULL;

  (lin    = dynamic_cast<HomogeneousTransformation                  *>(dof.get())) ||
  (tdffd  = dynamic_cast<BSplineFreeFormTransformationTD            *>(dof.get())) ||
  (svffd  = dynamic_cast<BSplineFreeFormTransformationSV            *>(dof.get())) ||
  (mffd   = dynamic_cast<MultiLevelFreeFormTransformation           *>(dof.get())) ||
  (msvffd = dynamic_cast<MultiLevelStationaryVelocityTransformation *>(dof.get()));

  FreeFormTransformation *ffd1 = NULL;
  FreeFormTransformation *ffd2 = NULL;

  if      (lin)    lin   ->Invert();
  else if (tdffd)  tdffd ->Invert();
  else if (svffd)  svffd ->Invert();
  else if (msvffd) msvffd->Invert();
  else if (mffd) {

    // Check number of levels
    if (mffd->NumberOfLevels() > 1) {
      cerr << EXECNAME << ": Inverting of FFDs with more than one level currently not implemented" << endl;
      exit(1);
    }

    // Invert global transformation
    mffd->GetGlobalTransformation()->Invert();

    // Invert local FFD
    if ((ffd1 = mffd->PopLocalTransformation())) {
      if (verbose > 1) {
        cout << "Invert local transformation of type " << ffd1->NameOfClass() << endl;
      }
      if (strcmp(ffd1->NameOfClass(), "BSplineFreeFormTransformationTD") == 0) {

        BSplineFreeFormTransformationTD *affd1 = dynamic_cast<BSplineFreeFormTransformationTD *>(ffd1);
        BSplineFreeFormTransformationTD *affd2 = new BSplineFreeFormTransformationTD(*affd1);

        affd2->Invert();

        ffd2 = affd2;

      } else if (strcmp(ffd1->NameOfClass(), "BSplineFreeFormTransformationSV") == 0) {

        BSplineFreeFormTransformationSV *affd1 = dynamic_cast<BSplineFreeFormTransformationSV *>(ffd1);
        BSplineFreeFormTransformationSV *affd2 = new BSplineFreeFormTransformationSV(*affd1);

        affd2->Invert();

        ffd2 = affd2;

      } else if (strcmp(ffd1->NameOfClass(), "BSplineFreeFormTransformation3D") == 0) {

        BSplineFreeFormTransformation3D *affd1 = dynamic_cast<BSplineFreeFormTransformation3D *>(ffd1);
        BSplineFreeFormTransformation3D *affd2 = new BSplineFreeFormTransformation3D(*affd1);

        // Evaluate inverse displacements
        double x, y, z;

        double *dx = new double[affd1->NumberOfDOFs() / 3];
        double *dy = new double[affd1->NumberOfDOFs() / 3];
        double *dz = new double[affd1->NumberOfDOFs() / 3];

        int index = 0;
        for (int k = 0; k < affd1->Z(); ++k)
        for (int j = 0; j < affd1->Y(); ++j)
        for (int i = 0; i < affd1->X(); ++i, ++index) {
          x = i, y = j, z = k;
          affd1->LatticeToWorld(x, y, z);
          affd1->InverseDisplacement(x, y, z);
          dx[index] = x;
          dy[index] = y;
          dz[index] = z;
        }

        // Interpolate inverse displacements
        affd2->Interpolate(dx, dy, dz);

        delete[] dx;
        delete[] dy;
        delete[] dz;

        ffd2 = affd2;

      } else {
        cerr << EXECNAME << ": Local transformation is of unsupported type: " << ffd1->NameOfClass() << endl;
        exit(1);
      }

      // Push inverted FFD
      mffd->PushLocalTransformation(ffd2);
    }

  } else {
    cerr << EXECNAME << ": Cannot invert transformation of type " << dof->NameOfClass() << endl;
    exit(1);
  }

  // Write inverted transformation
  dof->Write(dofout_name);

  // Evaluate error
  if (verbose > 0 && ffd1 && ffd2) {
    double x1, y1, z1, x2, y2, z2;

    double error     = 0.0;
    double max_error = 0.0;
    double rms_error = 0.0;

    double t1 = ffd1->LatticeToTime(0);
    double t2 = ffd1->LatticeToTime(ffd1->T() - 1);

    for (int k = 0; k < ffd1->Z(); ++k)
    for (int j = 0; j < ffd1->Y(); ++j)
    for (int i = 0; i < ffd1->X(); ++i) {
      x1 = i, y1 = j, z1 = k;
      ffd1->LatticeToWorld(x1, y1, z1);
      x2 = x1, y2 = y1, z2 = z1;
      ffd1->Transform(x2, y2, z2, t1, t2);
      ffd2->Transform(x2, y2, z2, t1, t2);
      error = sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
      rms_error += error;
      if (error > max_error) max_error = error;
    }

    rms_error /= static_cast<double>(ffd1->NumberOfDOFs() / 3);
    cerr << "RMS error for inverse FFD is " << rms_error << endl;
    cerr << "Max error for inverse FFD is " << max_error << endl;
  }

  // Clean up
  if (ffd1 != ffd2) delete ffd1;

  return 0;
}
