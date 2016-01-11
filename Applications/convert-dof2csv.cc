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

#include <mirtkGenericImage.h>
#include <mirtkTransformation.h>
#include <mirtkHomogeneousTransformation.h>
#include <mirtkBSplineFreeFormTransformation3D.h>
#include <mirtkMultiLevelFreeFormTransformation.h>

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "Usage: " << name << " <dof> <csv> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Converts a transformation to a CSV file which can be loaded into STAR-CCM+ GUI." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -delimiter <str>   Delimiting string. (default: ' ')" << endl;
  cout << "  -target <image>    Lattice points for which to store the displacements. (default: FFD lattice)" << endl;
  cout << "  -points <file>     Point set for which to store the displacements." << endl;
  cout << "  -t1 <t>            Lower time interval limit. (default: -inf)" << endl;
  cout << "  -t2 <t>            Upper time interval limit. (default:  inf)" << endl;
  cout << "  -dt <t>            Temporal resolution. (default: lattice resolution)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // ---------------------------------------------------------------------------
  // Parse arguments
  EXPECTS_POSARGS(2);

  const char *input_name    = POSARG(1);
  const char *output_name   = POSARG(2);
  const char *target_name   = NULL;
  const char *points_name   = NULL;
  const char *delimiter     = " ";
  double      t, t0         = .0;
  double      tmin          = -numeric_limits<double>::infinity();
  double      tmax          = +numeric_limits<double>::infinity();
  double      dt            = .0;
  bool        displacements = true;

  for (ALL_OPTIONS) {
    if (OPTION("-delimiter") || OPTION("-delim") || OPTION("-d")) delimiter = ARGUMENT;
    else if (OPTION("-target")) target_name = ARGUMENT;
    else if (OPTION("-points")) points_name = ARGUMENT;
    else if (OPTION("-Tt")) t0   = atof(ARGUMENT);
    else if (OPTION("-t1")) tmin = atof(ARGUMENT);
    else if (OPTION("-t2")) tmax = atof(ARGUMENT);
    else if (OPTION("-dt")) dt   = atof(ARGUMENT);
    else if (OPTION("-pos")) displacements = false;
    else if (OPTION("-disp") || OPTION("-dx")) displacements = true;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // ---------------------------------------------------------------------------
  // Read input transformation
  unique_ptr<Transformation> dof(Transformation::New(input_name));

  HomogeneousTransformation        *aff  = NULL;
  FreeFormTransformation           *ffd  = NULL;
  MultiLevelFreeFormTransformation *mffd = NULL;

  (aff  = dynamic_cast<HomogeneousTransformation        *>(dof.get())) ||
  (ffd  = dynamic_cast<BSplineFreeFormTransformation3D  *>(dof.get())) ||
  (mffd = dynamic_cast<MultiLevelFreeFormTransformation *>(dof.get()));
  if (mffd) ffd = mffd->GetLocalTransformation(-1);

  ImageAttributes domain;
  if (target_name) {
    BinaryImage target(target_name);
    domain = target.Attributes();
  } else if (ffd) {
    domain = ffd->Attributes();
  }
  if (dt > .0) {
    double T = domain.LatticeToTime(domain._t - 1) - domain.LatticeToTime(0);
    domain._t  = floor(T / dt);
    domain._dt = dt;
  }

  // ---------------------------------------------------------------------------
  // Open CSV file
  ofstream csv(output_name);

  // ---------------------------------------------------------------------------
  // Write header
  if (displacements) {
    csv << "X" << delimiter << "Y" << delimiter << "Z";
    if (domain._t == 1) {
      csv << delimiter << "DX" << delimiter << "DY" << delimiter << "DZ";
    } else {
      for (int l = 0; l < domain._t; ++l) {
        t = domain.LatticeToTime(l);
        if (tmin <= t && t <= tmax) {
          csv << delimiter << "DX[t=" << domain.LatticeToTime(l) << "ms]";
          csv << delimiter << "DY[t=" << domain.LatticeToTime(l) << "ms]";
          csv << delimiter << "DZ[t=" << domain.LatticeToTime(l) << "ms]";
        }
      }
    }
  } else {
    if (domain._t == 1) {
      csv << "X[t=0ms]" << delimiter << "Y[t=0ms]" << delimiter << "Z[t=0ms]";
      csv << delimiter;
      csv << "X[t=1ms]" << delimiter << "Y[t=1ms]" << delimiter << "Z[t=1ms]";
    } else {
      for (int l = 0, c = 0; l < domain._t; ++l) {
        t = domain.LatticeToTime(l);
        if (tmin <= t && t <= tmax) {
          if (++c > 1) csv << delimiter;
          csv << "X[t=" << domain.LatticeToTime(l) << "ms]";
          csv << delimiter << "Y[t=" << domain.LatticeToTime(l) << "ms]";
          csv << delimiter << "Z[t=" << domain.LatticeToTime(l) << "ms]";
        }
      }
    }
  }
  csv << "\n";

  // ---------------------------------------------------------------------------
  // Either write displacement/new position for each input point
  if (points_name) {

    PointSet points;
    points.Read(points_name);

    double dx, dy, dz;
    for (int r = 0; r < points.Size(); ++r) {
      if (displacements) {
        csv << points(r)._x << delimiter << points(r)._y << delimiter << points(r)._z;
        csv << delimiter;
      }
      for (int l = 0, c = 0; l < domain._t; ++l) {
        t = domain.LatticeToTime(l);
        if (domain._t == 1 || (tmin <= t && t <= tmax)) {
          if (++c > 1) csv << delimiter;
          dx = points(r)._x, dy = points(r)._y, dz = points(r)._z;
          if (displacements) {
            dof->Displacement(dx, dy, dz, t, t0);
          } else {
            dof->Transform(dx, dy, dz, t, t0);
          }
          csv << dx << delimiter << dy << delimiter << dz;
        }
      }
      csv << "\n";
    }

  // ---------------------------------------------------------------------------
  // or write displacement/new position for each lattice point
  } else {

    // Evaluate displacements
    if (!domain) {
      FatalError("-target image required for converting a linear transformation!");
    }
    WorldCoordsImage wc;
    RealImage *disp = new RealImage[domain._t];
    for (int l = 0; l < domain._t; ++l) {
      t = domain.LatticeToTime(l);
      if (domain._t == 1 || (tmin <= t && t <= tmax)) {
        disp[l].Initialize(ffd->Attributes(), 3);
        disp[l].PutTOrigin(domain.LatticeToTime(l));
        if (wc.IsEmpty()) disp[l].ImageToWorld(wc);
        dof->Displacement(disp[l], t0, &wc);
      }
    }
    if (wc.IsEmpty()) {
      FatalError("Invalid time interval [" << tmin << " " << tmax << "]");
    }

    // Write rows
    double dx, dy, dz;
    const int nvox = domain.NumberOfSpatialPoints();
    const WorldCoordsImage::VoxelType *x = wc.Data();
    const WorldCoordsImage::VoxelType *y = x + nvox;
    const WorldCoordsImage::VoxelType *z = y + nvox;

    for (int r = 0; r < nvox; ++r, ++x, ++y, ++z) {
      if (displacements) {
        csv << *x << delimiter << *y << delimiter << *z;
        csv << delimiter;
      }
      for (int l = 0, c = 0; l < domain._t; ++l) {
        t = domain.LatticeToTime(l);
        if (domain._t == 1 || (tmin <= t && t <= tmax)) {
          if (++c > 1) csv << delimiter;
          dx = *(disp[l].Data() + r           );
          dy = *(disp[l].Data() + r +     nvox);
          dz = *(disp[l].Data() + r + 2 * nvox);
          if (!displacements) {
            dx += *x, dy += *y, dz += *z;
          }
          csv << dx << delimiter << dy << delimiter << dz;
        }
      }
      csv << "\n";
    }

    // Free displacement fields
    delete[] disp;

  }

  // ---------------------------------------------------------------------------
  // Close CSV file
  csv.close();
  return 0;
}
