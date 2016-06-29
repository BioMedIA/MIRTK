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

#include "mirtk/Vector3.h"
#include "mirtk/Matrix3x3.h"
#include "mirtk/PointSet.h"
#include "mirtk/GenericImage.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
/// Print help screen
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <image> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Cut input brain volume/mask into left and/or right hemisphere(s)." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  image   Input brain image, mask, or (tissue) probability map." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -plane <nx> <ny> <nz> <b>   Cutting plane parameters. (default: 1 0 0 0)" << endl;
  cout << "  -hemispheres <file>         Hemispheres mask (0: outside, 1: right, 2: left) from" << endl;
  cout << "                              which cutting plane is computed. (default: none)" << endl;
  cout << "  -subcortical <file>         Mask of subcortical structures which may be cut in half" << endl;
  cout << "                              no matter whether they are labeled left or right. (default: none)" << endl;
  cout << "  -left  <file>               Output name of left hemisphere image. (default: none)" << endl;
  cout << "  -right <file>               Output name of right hemisphere image. (default: none)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Labels of brain hemispheres in input mask image
///
/// @sa NeoSeg's segmentations-data/cortmrf
enum Hemisphere
{
  UH = 0, ///< Unspecified hemishpere
  RH = 1, ///< Right brain hemisphere
  LH = 2  ///< Left  brain hemisphere
};

// -----------------------------------------------------------------------------
/// Naive cutting plane calculation assuming symmetric brain hemispheres mask
double SymmetricCuttingPlane(const ByteImage &mask, double n[3])
{
  BytePixel h;
  Point center[2], p;
  int   numvox[2] = {0};
  for (int k = 0; k < mask.Z(); ++k)
  for (int j = 0; j < mask.Y(); ++j)
  for (int i = 0; i < mask.X(); ++i) {
    h = mask(i, j, k) - 1;
    if (h < 0 || h > 1) continue;
    p._x = i, p._y = j, p._z = k;
    mask.ImageToWorld(p);
    center[h] += p;
    numvox[h] += 1;
  }
  if (numvox[0] == 0 || numvox[1] == 0) {
    FatalError("Expected mask to contain label " << RH
                 << " for right hemisphere and " << LH << " for left hemisphere");
  }
  center[0] /= numvox[0];
  center[1] /= numvox[1];
  p = (center[0] + center[1]) / 2.0;
  n[0] = center[0]._x - center[1]._x;
  n[1] = center[0]._y - center[1]._y;
  n[2] = center[0]._z - center[1]._z;
  double m = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] /= m, n[1] /= m, n[2] /= m;
  return - (p[0]*n[0] + p[1]*n[1] + p[2]*n[2]);
}

// -----------------------------------------------------------------------------
/// Cutting plane separating left/right brain hemisphere clusters
double CuttingPlane(const ByteImage &mask, double n[3])
{
  double p[3], c[3] = {.0};
  bool   include;

  PointSet points;
  for (int k = 0; k < mask.Z(); ++k)
  for (int j = 0; j < mask.Y(); ++j)
  for (int i = 0; i < mask.X(); ++i) {
    if (mask(i, j, k) == RH || mask(i, j, k) == LH) {
      include = false;
      for (int nk = k-1; nk <= k+1; ++nk) {
        if (nk < 0 || nk >= mask.Z()) continue;
        for (int nj = j-1; nj <= j+1; ++nj) {
          if (nj < 0 || nj >= mask.Y()) continue;
          for (int ni = i-1; ni <= i+1; ++ni) {
            if (ni < 0 || ni >= mask.X()) continue;
            if ((mask(ni, nj, nk) == RH || mask(ni, nj, nk) == LH) && mask(ni, nj, nk) != mask(i, j, k)) {
              include = true;
              break;
            }
          }
          if (include) break;
        }
        if (include) break;
      }
      if (include) {
        p[0] = i, p[1] = j, p[2] = k;
        mask.ImageToWorld(p[0], p[1], p[2]);
        points.Add(p);
        c[0] += p[0], c[1] += p[1], c[2] += p[2];
      }
    }
  }

  if (points.Size() == 0) return SymmetricCuttingPlane(mask, n);

  c[0] /= points.Size();
  c[1] /= points.Size();
  c[2] /= points.Size();

  Matrix3x3 covar(.0);
  for (int i = 0; i < points.Size(); ++i) {
    Point &point = points(i);
    p[0] = point._x - c[0];
    p[1] = point._y - c[1];
    p[2] = point._z - c[2];
    for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c) {
      covar[r][c] += p[r] * p[c];
    }
  }

  double  eigen[3];
  Vector3 axis [3];
  covar.EigenSolveSymmetric(eigen, axis);

  int i = abs(eigen[0]) < abs(eigen[1]) ? 0 : 1;
  if (abs(eigen[2]) < abs(eigen[i])) i = 2;
  n[0] = axis[i][0];
  n[1] = axis[i][1];
  n[2] = axis[i][2];

  double m = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
  n[0] /= m, n[1] /= m, n[2] /= m;
  return - (c[0]*n[0] + c[1]*n[1] + c[2]*n[2]);
}

// -----------------------------------------------------------------------------
/// Pad voxels of specified brain hemisphere
BaseImage *Pad(const BaseImage *image,
               const ByteImage &hemi, const ByteImage &mask,
               double b, double n[3], Hemisphere hemisphere, double padding = .0)
{
  double wx, wy, wz, d;
  BaseImage *output = image->Copy();

  for (int l = 0; l < image->T(); ++l)
  for (int k = 0; k < image->Z(); ++k)
  for (int j = 0; j < image->Y(); ++j)
  for (int i = 0; i < image->X(); ++i) {
    if (!hemi.IsEmpty() && hemi(i, j, k) && (mask.IsEmpty() || !mask(i, j, k))) {
      // trust input mask where defined as cutting plane may be inaccurate
      if (hemi(i, j, k) == hemisphere) {
        output->PutAsDouble(i, j, k, l, padding);
      }
    } else {
      wx = i, wy = j, wz = k;
      image->ImageToWorld(wx, wy, wz);
      d = wx * n[0] + wy * n[1] + wz * n[2] + b;
      if ((hemisphere == LH && d <  .0) ||
          (hemisphere == RH && d >= .0)) {
        output->PutAsDouble(i, j, k, l, padding);
      }
    }
  }

  return output;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(1);

  const char *image_name = POSARG(1);
  const char *hemi_name  = NULL;
  const char *mask_name  = NULL;
  const char *left_name  = NULL;
  const char *right_name = NULL;

  double n[3] = {1.0, .0, .0}, b = .0;
  double padding = .0;

  for (ALL_OPTIONS) {
    if (OPTION("-plane")) {
      n[0] = atof(ARGUMENT);
      n[1] = atof(ARGUMENT);
      n[2] = atof(ARGUMENT);
      b    = atof(ARGUMENT);
    }
    else if (OPTION("-hemispheres")) hemi_name = ARGUMENT;
    else if (OPTION("-subcortical")) mask_name = ARGUMENT;
    else if (OPTION("-left"))    left_name   = ARGUMENT;
    else if (OPTION("-right"))   right_name  = ARGUMENT;
    else if (OPTION("-padding")) padding     = atof(ARGUMENT);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input image
  InitializeIOLibrary();
  UniquePtr<BaseImage> image(BaseImage::New(image_name));

  // Compute parameters of cutting plane
  BinaryImage hemi;
  if (hemi_name) {
    hemi.Read(hemi_name);
    if (!hemi.Attributes().EqualInSpace(image->Attributes())) {
      FatalError("Input hemispheres mask has differing spatial attributes");
    }
    b = CuttingPlane(hemi, n);
  }
  if (verbose) {
    cout << "Cutting plane equation: x^T [" << n[0] << " " << n[1] << " " << n[2] << "] + " << b << " = 0" << endl;
  }

  // Optional structures mask to be cut
  BinaryImage mask;
  if (mask_name) {
    mask.Read(mask_name);
    if (!mask.Attributes().EqualInSpace(image->Attributes())) {
      FatalError("Error: Input subcortical structures mask has differing spatial attributes");
    }
  }

  // Write left/right hemisphere
  if (left_name) {
    UniquePtr<BaseImage> left(Pad(image.get(), hemi, mask, b, n, RH, padding));
    left->Write(left_name);
  }
  if (right_name) {
    UniquePtr<BaseImage> right(Pad(image.get(), hemi, mask, b, n, LH, padding));
    right->Write(right_name);
  }

  return 0;
}
