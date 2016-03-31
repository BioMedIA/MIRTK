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

#include "mirtk/GenericImage.h"
#include "mirtk/Transformation.h"
#include "mirtk/Noise.h"
#include "mirtk/GaussianBlurring.h"
#include "mirtk/GaussianBlurring2D.h"
#include "mirtk/ConvolutionWithGaussianDerivative.h"
#include "mirtk/ConvolutionWithGaussianDerivative2.h"
#include "mirtk/DisplacementToVelocityField.h"
#include "mirtk/VelocityToDisplacementFieldEuler.h"

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace mirtk;


// ===========================================================================
// Initialization
// ===========================================================================

// ---------------------------------------------------------------------------
void init(char *name_in, ImageAttributes &grid,
          double sigma, double max,
          double *&x,  double *&y,  double *&z,
          double *&dx, double *&dy, double *&dz)
{
  GenericImage<RealPixel> disp;

  if (name_in == NULL) {
    // create zero displacement field
    ImageAttributes attr = grid;
    attr._t = 3;
    disp.Initialize(attr);

    // add random noise to displacement field
    GaussianNoise<RealPixel> noise(0.0, sigma, -fabs(max), +fabs(max));
    noise.SetInput (&disp);
    noise.SetOutput(&disp);
    noise.Run();

    // smooth displacements
    if (grid._z > 1) {
      GaussianBlurring<RealPixel> blur(1.25);
      blur.SetInput (&disp);
      blur.SetOutput(&disp);
      blur.Run();
    } else {
      GaussianBlurring2D<RealPixel> blur(1.25);
      blur.SetInput (&disp);
      blur.SetOutput(&disp);
      blur.Run();
    }
  } else if (strstr(name_in, ".dof") != NULL) {
    // read transformation from file
    Transformation           *t    = Transformation::New(name_in);
    MultiLevelTransformation *mffd = dynamic_cast<MultiLevelTransformation *>(t);
    FreeFormTransformation   *ffd  = dynamic_cast<FreeFormTransformation *>(t);
    if (mffd != NULL && mffd->NumberOfLevels() > 0) {
      ffd = mffd->GetLocalTransformation(mffd->NumberOfLevels() - 1);
    }
    if (ffd != NULL) grid = ffd->Attributes();
    ImageAttributes attr = grid;
    attr._t = 3;
    disp.Initialize(attr);
    t->Displacement(disp);
  } else {
    // read displacements from file
    disp.Read(name_in);
    if (disp.GetT() != 3) {
      cerr << "Invalid input displacement field: " << name_in << endl;
      exit(1);
    }
    grid = disp.GetImageAttributes();
    grid._t = 1;
  }

  // copy to output arrays
  const int n = grid._x * grid._y * grid._z;

  x  = new double[n];
  y  = new double[n];
  z  = new double[n];
  dx = new double[n];
  dy = new double[n];
  dz = new double[n];

  int p = 0;
  for (int k = 0; k < grid._z; k++) {
    for (int j = 0; j < grid._y; j++) {
      for (int i = 0; i < grid._x; i++) {
        x[p] = i;
        y[p] = j;
        z[p] = k;
        disp.ImageToWorld(x[p], y[p], z[p]);
        dx[p] = disp(i, j, k, 0);
        dy[p] = disp(i, j, k, 1);
        dz[p] = disp(i, j, k, 2);
        p++;
      }
    }
  }
}

// ===========================================================================
// Auxiliary functions
// ===========================================================================

// ---------------------------------------------------------------------------
inline void jac(Matrix &J, GenericImage<double> &img, int i, int j, int k)
{
  J.Initialize(3, 3);

  double sx = 1.0 / img.GetXSize();
  double sy = 1.0 / img.GetYSize();
  double sz = 1.0 / img.GetZSize();

  if (i <= 0 || i >= img.GetX() - 1) {
    J(0, 0) = 0;
    J(0, 1) = 0;
    J(0, 2) = 0;
  } else {
    // central difference
    J(0, 0) = 0.5 * sx * (img(i + 1, j, k, 0) - img(i - 1, j, k, 0));
    J(0, 1) = 0.5 * sx * (img(i + 1, j, k, 1) - img(i - 1, j, k, 1));
    J(0, 2) = 0.5 * sx * (img(i + 1, j, k, 2) - img(i - 1, j, k, 2));
  }

  if (j <= 0 || j >= img.GetY() - 1) {
    J(1, 0) = 0;
    J(1, 1) = 0;
    J(1, 2) = 0;
  } else {
    // central difference
    J(1, 0) = 0.5 * sy * (img(i, j + 1, k, 0) - img(i, j - 1, k, 0));
    J(1, 1) = 0.5 * sy * (img(i, j + 1, k, 1) - img(i, j - 1, k, 1));
    J(1, 2) = 0.5 * sy * (img(i, j + 1, k, 2) - img(i, j - 1, k, 2));
  }

  if (k <= 0 || k >= img.GetZ() - 1) {
    J(2, 0) = 0;
    J(2, 1) = 0;
    J(2, 2) = 0;
  } else {
    // central difference
    J(2, 0) = 0.5 * sz * (img(i, j, k + 1, 0) - img(i, j, k - 1, 0));
    J(2, 1) = 0.5 * sz * (img(i, j, k + 1, 1) - img(i, j, k - 1, 1));
    J(2, 2) = 0.5 * sz * (img(i, j, k + 1, 2) - img(i, j, k - 1, 2));
  }
}

// ===========================================================================
// Exponential map
// ===========================================================================

// ---------------------------------------------------------------------------
inline void expEuler(GenericImage<double> &v, double &x, double &y, double &z, bool inv)
{
  const int    NumberOfSteps = 100; // Vercauteren et al. use 500 Euler steps for log()
  const double dt            = (inv ? -1.0 : 1.0) / static_cast<double>(NumberOfSteps);

  LinearInterpolateImageFunction vf;
  vf.SetInput(&v);
  vf.Initialize();

  for (int s = 0; s < NumberOfSteps; s++) {
    double i = x;
    double j = y;
    double k = z;
    v.WorldToImage(i, j, k);
    x += vf.Evaluate(i, j, k, 0) * dt;
    y += vf.Evaluate(i, j, k, 1) * dt;
    z += vf.Evaluate(i, j, k, 2) * dt;
  }
}

// ===========================================================================
// Logarithmic map
// ===========================================================================

// ---------------------------------------------------------------------------
void testLogWithExpEuler(const ImageAttributes &grid,
                         const double *x,  const double *y,  const double *z,
                         const double *dx, const double *dy, const double *dz)
{
  const int nsteps = 5;
  const int n      = grid._x * grid._y * grid._z;

  double x1, y1, z1;
  double x2, y2, z2;
  double lx, ly, lz;
  double vx, vy, vz;
  double dvx, dvy, dvz;
  Matrix Jv, Jdv;
  double error, avg_error, max_error;
  int p;

  ImageAttributes attr = grid;
  attr._t = 3;

  GenericImage<double> phi(attr); // displacement field
  GenericImage<double> v  (attr); // current velocity field
  GenericImage<double> vn (attr); // next velocity field
  GenericImage<double> dv (attr); // dv = exp(-vp) ° Phi(x) - x

  // Initialize:
  // v_0 = 0
  // v_1 = exp(-v_0) ° Phi(x) - x = Phi(x) - x
  p = 0;
  for (int k = 0; k < grid._z; k++) {
    for (int j = 0; j < grid._y; j++) {
      for (int i = 0; i < grid._x; i++) {
        phi(i, j, k, 0) = dx[p];
        phi(i, j, k, 1) = dy[p];
        phi(i, j, k, 2) = dz[p];
        v  (i, j, k, 0) = dx[p];
        v  (i, j, k, 1) = dy[p];
        v  (i, j, k, 2) = dz[p];
        p++;
      }
    }
  }

  for (int s = 0; s < nsteps; s++) {
    // Compute dv = exp(-v) ° Phi(x) - x
    for (int k = 0; k < grid._z; k++) {
      for (int j = 0; j < grid._y; j++) {
        for (int i = 0; i < grid._x; i++) {
          // Convert to world coordinates
          x1 = i;
          y1 = j;
          z1 = k;
          phi.ImageToWorld(x1, y1, z1);
          // Transform by Phi
          x2 = x1 + phi(i, j, k, 0);
          y2 = y1 + phi(i, j, k, 1);
          z2 = z1 + phi(i, j, k, 2);
          // Transform by inverse of current transformation [= exp(-v)]
          expEuler(v, x2, y2, z2, true);
          // Compute delta values
          dv(i, j, k, 0) = x2 - x1;
          dv(i, j, k, 1) = y2 - y1;
          dv(i, j, k, 2) = z2 - z1;
        }
      }
    }
    // Smooth dv to stabilize computation (see Vercauteren et al. ITK filter)
    GaussianBlurring2D<double> blur(2.0 * (grid._dx > grid._dy ? grid._dx : grid._dy));
    blur.SetInput (&dv);
    blur.SetOutput(&dv);
    blur.Run();
    // Update velocities using the Baker-Campbell-Hausdorff (BCH) formula
    for (int k = 0; k < grid._z; k++) {
      for (int j = 0; j < grid._y; j++) {
        for (int i = 0; i < grid._x; i++) {
          vx  = v (i, j, k, 0);
          vy  = v (i, j, k, 1);
          vz  = v (i, j, k, 2);
          dvx = dv(i, j, k, 0);
          dvy = dv(i, j, k, 1);
          dvz = dv(i, j, k, 2);
          // Lie bracket [v, dv]
#if 1
          jac(Jv,  v,  i, j, k);
          jac(Jdv, dv, i, j, k);
          lx = dvx * Jv(0, 0) - vx * Jdv(0, 0) +
               dvy * Jv(0, 1) - vy * Jdv(0, 1) +
               dvz * Jv(0, 2) - vz * Jdv(0, 2);
          ly = dvx * Jv(1, 0) - vx * Jdv(1, 0) +
               dvy * Jv(1, 1) - vy * Jdv(1, 1) +
               dvz * Jv(1, 2) - vz * Jdv(1, 2);
          lz = dvx * Jv(2, 0) - vx * Jdv(2, 0) +
               dvy * Jv(2, 1) - vy * Jdv(2, 1) +
               dvz * Jv(2, 2) - vz * Jdv(2, 2);
          dvx += 0.5 * lx;
          dvy += 0.5 * ly;
          dvz += 0.5 * lz;
#endif
          // BCH update step
          vn(i, j, k, 0) = vx + dvx;
          vn(i, j, k, 1) = vy + dvy;
          vn(i, j, k, 2) = vz + dvz;
        }
      }
    }
    // Set velocities to updated ones
    v = vn;
  }

  // The stationary velocity field v = log(Phi) is now computed.
  // In the following, the dense vector field is approximated/interpolated
  // by a 3D B-spline FFD with a control point at each grid point (voxel).

  // Store integrated displacements in [rx, ry, rz] while also calculating
  // the RMS error of the resulting (dense) displacement field.

  double *rx = new double[n];
  double *ry = new double[n];
  double *rz = new double[n];

  avg_error = 0;
  max_error = 0;
  p = 0;
  for (int k = 0; k < grid._z; k++) {
    for (int j = 0; j < grid._y; j++) {
      for (int i = 0; i < grid._x; i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        v.ImageToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        expEuler(v, x2, y2, z2, false);
        rx[p] = x2 - x1;
        ry[p] = y2 - y1;
        rz[p] = z2 - z1;
        error = sqrt(pow(rx[p] - phi(i, j, k, 0), 2) +
                     pow(ry[p] - phi(i, j, k, 1), 2) +
                     pow(rz[p] - phi(i, j, k, 2), 2));
        if (max_error < error) max_error = error;
        avg_error += error;
        p++;
      }
    }
  }
  avg_error /= static_cast<double>(n);

  cout << "Average RMS error SV (grid, expEuler): " << avg_error << endl;
  cout << "Maximum RMS error SV (grid, expEuler): " << max_error << endl;
  cout.flush();

  // Now approximate a B-spline FFD from the given displacements
  BSplineFreeFormTransformation3D ffd(const_cast<ImageAttributes&>(grid), grid._dx, grid._dy, grid._dz);

  double *tmpx = new double[n];
  double *tmpy = new double[n];
  double *tmpz = new double[n];

  memcpy(tmpx, rx, sizeof(double) * n);
  memcpy(tmpy, ry, sizeof(double) * n);
  memcpy(tmpz, rz, sizeof(double) * n);

  ffd.ApproximateAsNew(x, y, z, tmpx, tmpy, tmpz, n);
  ffd.Write("test1-vffd1.dof.gz");

  // Error of the approximation, thus compare FFD displacement to
  // displacement [rx, ry, rz] which was used as input to ApproximateAsNew()
  avg_error = 0;
  max_error = 0;
  for (int k = 2; k < grid._z - 2; k++) {
    for (int j = 2; j < grid._y - 2; j++) {
      for (int i = 2; i < grid._x - 2; i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        v.ImageToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        ffd.Transform(x2, y2, z2);
        p = k * grid._x * grid._y + j * grid._x + i;
        error = sqrt(pow((x2 - x1) - rx[p], 2) +
                     pow((y2 - y1) - ry[p], 2) +
                     pow((z2 - z1) - rz[p], 2));
        if (max_error < error) max_error = error;
        avg_error += error;
        p++;
      }
    }
  }
  avg_error /= static_cast<double>(n);

  cout << "Average RMS error SV (grid, approx. FFD): " << avg_error << endl;
  cout << "Maximum RMS error SV (grid, approx. FFD): " << max_error << endl;
  cout.flush();

  // Error of the approximated FFD relative to the initial displacement Phi
  avg_error = 0;
  max_error = 0;
  for (int k = 2; k < grid._z - 2; k++) {
    for (int j = 2; j < grid._y - 2; j++) {
      for (int i = 2; i < grid._x - 2; i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        v.ImageToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        ffd.Transform(x2, y2, z2);
        error = sqrt(pow((x2 - x1) - phi(i, j, k, 0), 2) +
                     pow((y2 - y1) - phi(i, j, k, 1), 2) +
                     pow((z2 - z1) - phi(i, j, k, 2), 2));
        if (max_error < error) max_error = error;
        avg_error += error;
        p++;
      }
    }
  }
  avg_error /= static_cast<double>(n);

  cout << "Average RMS error SV (grid, expEuler, approx. FFD): " << avg_error << endl;
  cout << "Maximum RMS error SV (grid, expEuler, approx. FFD): " << max_error << endl;
  cout.flush();

  memcpy(tmpx, rx, sizeof(double) * n);
  memcpy(tmpy, ry, sizeof(double) * n);
  memcpy(tmpz, rz, sizeof(double) * n);

  // Now do the same using the B-spline interpolation instead
  ffd.Interpolate(tmpx, tmpy, tmpz);
  ffd.Write("test1-vffd2.dof.gz");

  avg_error = 0;
  max_error = 0;
  for (int k = 2; k < grid._z - 2; k++) {
    for (int j = 2; j < grid._y - 2; j++) {
      for (int i = 2; i < grid._x - 2; i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        v.ImageToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        ffd.Transform(x2, y2, z2);
        p = k * grid._x * grid._y + j * grid._x + i;
        error = sqrt(pow((x2 - x1) - rx[p], 2) +
                     pow((y2 - y1) - ry[p], 2) +
                     pow((z2 - z1) - rz[p], 2));
        if (max_error < error) max_error = error;
        avg_error += error;
        p++;
      }
    }
  }
  avg_error /= static_cast<double>(n);

  cout << "Average RMS error SV (grid, interp. FFD): " << avg_error << endl;
  cout << "Maximum RMS error SV (grid, interp. FFD): " << max_error << endl;
  cout.flush();

  avg_error = 0;
  max_error = 0;
  for (int k = 2; k < grid._z - 2; k++) {
    for (int j = 2; j < grid._y - 2; j++) {
      for (int i = 2; i < grid._x - 2; i++) {
        x1 = i;
        y1 = j;
        z1 = k;
        v.ImageToWorld(x1, y1, z1);
        x2 = x1;
        y2 = y1;
        z2 = z1;
        ffd.Transform(x2, y2, z2);
        error = sqrt(pow((x2 - x1) - phi(i, j, k, 0), 2) +
                     pow((y2 - y1) - phi(i, j, k, 1), 2) +
                     pow((z2 - z1) - phi(i, j, k, 2), 2));
        if (max_error < error) max_error = error;
        avg_error += error;
        p++;
      }
    }
  }
  avg_error /= static_cast<double>(n);

  cout << "Average RMS error SV (grid, expEuler, interp. FFD): " << avg_error << endl;
  cout << "Maximum RMS error SV (grid, expEuler, interp. FFD): " << max_error << endl;
  cout.flush();

  delete[] tmpx;
  delete[] tmpy;
  delete[] tmpz;

  delete[] rx;
  delete[] ry;
  delete[] rz;
}

// ===========================================================================
// Main
// ===========================================================================

int DisplacementToVelocityFieldTest(int ac, char *av[])
{
  // extract/discard program name
  ac--;
  av++;

  // parse arguments
  char  *dofin  = NULL;
  char  *dofout = NULL;
  double std    = 2.0;
  double max    = 5.0;

  ImageAttributes grid;
  grid._x  = 64;
  grid._y  = 64;
  grid._z  = 1;
  grid._dx = 1.0;
  grid._dy = 1.0;
  grid._dz = 1.0;

  for (int i = 0; i < ac; i++) {
    char *opt = av[i];   // option/flag name
    char *arg = av[i+1]; // argument for option
    bool flag = false;     // whether option has no argument

    if (strcmp(opt, "-dofin") == 0) {
      dofin = arg;
    } else if (strcmp(opt, "-dofout") == 0) {
      dofout = arg;
    } else if (strcmp(opt, "-std") == 0) {
      std = atof(arg);
    } else if (strcmp(opt, "-max") == 0) {
      max = atof(arg);
    } else if (strcmp(opt, "-x") == 0) {
      grid._x = atoi(arg);
    } else if (strcmp(opt, "-y") == 0) {
      grid._y = atoi(arg);
    } else if (strcmp(opt, "-z") == 0) {
      grid._z = atoi(arg);
    } else if (strcmp(opt, "-dx") == 0) {
      grid._dx = atoi(arg);
    } else if (strcmp(opt, "-dy") == 0) {
      grid._dy = atoi(arg);
    } else if (strcmp(opt, "-dz") == 0) {
      grid._dz = atoi(arg);
    } else {
      cerr << "Unknown option: " << opt << endl;
      exit(1);
    }

    if (!flag) i++; // skip option argument
  }

  // initial displacement field
  double *x, *y, *z, *dx, *dy, *dz;
  init(dofin, grid, std, max, x, y, z, dx, dy, dz);

  ImageAttributes attr = grid; attr._t = 3;
  GenericImage<double> din(attr);

  int p = 0;
  for (int k = 0; k < grid._z; k++) {
    for (int j = 0; j < grid._y; j++) {
      for (int i = 0; i < grid._x; i++) {
        din(i, j, k, 0) = dx[p];
        din(i, j, k, 1) = dy[p];
        din(i, j, k, 2) = dz[p];
        p++;
      }
    }
  }

  DisplacementToVelocityFieldBCH<double> dtov;
  GenericImage<double>                   v;

  dtov.SetInput (&din);
  dtov.SetOutput(&v);

  dtov.Run();

  VelocityToDisplacementFieldEuler<double> vtod;
  GenericImage<double>                     dout;

  vtod.SetInput (&v);
  vtod.SetOutput(&dout);

  vtod.Run();

  if (dofout != NULL) {
    if (strstr(dofout, ".dof") != NULL) {
      MultiLevelFreeFormTransformation mffd;
      FreeFormTransformation3D *ffd = new LinearFreeFormTransformation(dout, grid._dx, grid._dy, grid._dz);
      ffd->Interpolate(dout.GetPointerToVoxels(0, 0, 0, 0),
                       dout.GetPointerToVoxels(0, 0, 0, 1),
                       dout.GetPointerToVoxels(0, 0, 0, 2));
      mffd.PushLocalTransformation(ffd);
      mffd.Write(dofout);
    } else {
      dout.Write(dofout);
    }
  }

  // perform tests
  //testLogWithExpEuler(grid, x, y, z, dx, dy, dz);

  // clean up
  delete[] x;
  delete[] y;
  delete[] z;
  delete[] dx;
  delete[] dy;
  delete[] dz;

  exit(0);
}
