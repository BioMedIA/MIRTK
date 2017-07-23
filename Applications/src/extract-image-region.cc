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

#include "mirtk/PointSet.h"
#include "mirtk/BaseImage.h"
#include "mirtk/IOConfig.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <input> <output> [options]\n";
  cout << "\n";
  cout << "Description:" << endl;
  cout << "  Crop/pad image by extracting a region of interest and optionally split";
  cout << "  the extracted region into separate image files, e.g., individual slices";
  cout << "  of a volume saved as individual image files. The output image region\n";
  cout << "  is chosen such that it contains the union of all specified axis-aligned\n";
  cout << "  rectangular input image regions. In case of :option:`-pad`, the output\n";
  cout << "  region does not have to be fully contained within the input image.\n";
  cout << "  Values outside are then set to the specified padding value.\n";
  cout << "\n";
  cout << "Positional arguments:\n";
  cout << "  input    Input image file path.\n";
  cout << "  output   Output file path.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -Rx1 <int>\n";
  cout << "      Leftmost input voxel index along x axis.\n";
  cout << "  -Rx2 <int>\n";
  cout << "      Rightmost input voxel index along x axis.\n";
  cout << "  -Ry1 <int>\n";
  cout << "      Leftmost input voxel index along y axis.\n";
  cout << "  -Ry2 <int>\n";
  cout << "      Rightmost input voxel index along y axis.\n";
  cout << "  -Rz1 <int>\n";
  cout << "      Leftmost input voxel index along z axis.\n";
  cout << "  -Rz2 <int>\n";
  cout << "      Rightmost input voxel index along z axis.\n";
  cout << "  -Rt1 <int>\n";
  cout << "      Leftmost input voxel index along t axis.\n";
  cout << "  -Rt2 <int>\n";
  cout << "      Rightmost input voxel index along t axis.\n";
  cout << "  -patch <i> <j> <k> <nx> [<ny> [<nz>]]\n";
  cout << "      Extract image patch of given size centered at the specified voxel.\n";
  cout << "  -closest-patch <x> <y> <z> <nx> [<ny> [<nz>]]\n";
  cout << "      Extract image patch of given size centered at nearest voxel to specified point [mm].\n";
  cout << "  -landmarks <file>...\n";
  cout << "      Extract minimum bounding box containing the landmark points.\n";
  cout << "  -ref <file>\n";
  cout << "      Extract region specified by discrete reference image domain.\n";
  cout << "  -margin <int>\n";
  cout << "      Add fixed-width margin to union of image regions. (default: 0)\n";
  cout << "  -scale <float>\n";
  cout << "      Scale resulting region by specified factor. (default: 1)\n";
  cout << "  -crop [value]\n";
  cout << "      Crop background with intensity below or equal specified value. (default: 0)\n";
  cout << "  -pad [value]\n";
  cout << "      Pad output image by the specified value. (default: 0)\n";
  cout << "  -split <dim>|x|y|z|t\n";
  cout << "      Split extracted region along specified dimension into separate images.\n";
  cout << "      For example, use '-split z' to save individual slices of a 3D volume.\n";
  cout << "      The <output> file path is appended with a format string '_%03d' before\n";
  cout << "      the file name extension to create unique output file paths, unless such\n";
  cout << "      format specification is part of the given file path already.\n";
  cout << "  -swap [on|off]\n";
  cout << "      When :option:`-split` is used, first swap the dimension along which the\n";
  cout << "      extracted region is split with the last dimension such that when splitting\n";
  cout << "      a volume along the first dimension, the output image files have z=1. (default: off)\n";
  cout << "  -noswap\n";
  cout << "      Set :option:`-swap` to off.\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

const int MIN_INDEX = numeric_limits<int>::min();
const int MAX_INDEX = numeric_limits<int>::max();

// -----------------------------------------------------------------------------
int ToVoxelIndex(const char *arg)
{
  int i;
  if (!FromString(arg, i)) {
    cerr << "Invalid voxel index argument: " << arg << endl;
    exit(1);
  }
  return i;
}

// -----------------------------------------------------------------------------
int ToRegionMargin(const char *arg)
{
  int i;
  if (!FromString(arg, i)) {
    cerr << "Invalid region margin argument: " << arg << endl;
    exit(1);
  }
  return i;
}

// -----------------------------------------------------------------------------
double ToScaleFactor(const char *arg)
{
  double s;
  if (!FromString(arg, s)) {
    cerr << "Invalid scale factor argument: " << arg << endl;
    exit(1);
  }
  return s;
}

// -----------------------------------------------------------------------------
int MinIndex(int i, int j)
{
  return (i <= MIN_INDEX ? j : min(i, j));
}

// -----------------------------------------------------------------------------
int MaxIndex(int i, int j)
{
  return (i >= MAX_INDEX ? j : max(i, j));
}

// -----------------------------------------------------------------------------
void ScaleInterval(int &i, int &j, double s)
{
  if (s == 1.0) return;
  double m = (i + j) / 2.0;
  double r = s * abs(j - i) / 2.0;
  i = ifloor(m - r);
  j = iceil (m + r);
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  // Read input image
  InitializeIOLibrary();
  UniquePtr<BaseImage> in(BaseImage::New(input_name));

  // Parse optional arguments and adjust image region accordingly
  int    xmargin = 0, ymargin = 0, zmargin = 0, tmargin = 0;
  double xscale  = 1, yscale  = 1, zscale  = 1, tscale  = 1;
  bool   pad           = false;
  double padding_value = .0;
  int    split         = -1;
  bool   swap_dims     = false;

  int x1 = MIN_INDEX, x2 = MAX_INDEX;
  int y1 = MIN_INDEX, y2 = MAX_INDEX;
  int z1 = MIN_INDEX, z2 = MAX_INDEX;
  int t1 = MIN_INDEX, t2 = MAX_INDEX;

  for (ALL_OPTIONS) {
    if      (OPTION("-Rx1")) x1 = MinIndex(x1, ToVoxelIndex(ARGUMENT));
    else if (OPTION("-Rx2")) x2 = MaxIndex(x2, ToVoxelIndex(ARGUMENT));
    else if (OPTION("-Ry1")) y1 = MinIndex(y1, ToVoxelIndex(ARGUMENT));
    else if (OPTION("-Ry2")) y2 = MaxIndex(y2, ToVoxelIndex(ARGUMENT));
    else if (OPTION("-Rz1")) z1 = MinIndex(z1, ToVoxelIndex(ARGUMENT));
    else if (OPTION("-Rz2")) z2 = MaxIndex(z2, ToVoxelIndex(ARGUMENT));
    else if (OPTION("-Rt1")) t1 = MinIndex(t1, ToVoxelIndex(ARGUMENT));
    else if (OPTION("-Rt2")) t2 = MaxIndex(t2, ToVoxelIndex(ARGUMENT));
    else if (OPTION("-patch")) {
      int i, j, k, nx, ny, nz;
      PARSE_ARGUMENT(i);
      PARSE_ARGUMENT(j);
      PARSE_ARGUMENT(k);
      PARSE_ARGUMENT(nx);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(ny);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(nz);
        else              nz = 1;
      } else ny = nz = nx;
      int rx = (nx - 1) / 2;
      int ry = (ny - 1) / 2;
      int rz = (nz - 1) / 2;
      x1 = i - rx, x2 = i + rx;
      y1 = j - ry, y2 = j + ry;
      z1 = k - rz, z2 = k + rz;
    }
    else if (OPTION("-closest-patch")) {
      double x, y, z;
      PARSE_ARGUMENT(x);
      PARSE_ARGUMENT(y);
      PARSE_ARGUMENT(z);
      int nx, ny, nz;
      PARSE_ARGUMENT(nx);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(ny);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(nz);
        else              nz = 1;
      } else ny = nz = nx;
      in->WorldToImage(x, y, z);
      int i  = iround(x);
      int j  = iround(y);
      int k  = iround(z);
      int rx = (nx - 1) / 2;
      int ry = (ny - 1) / 2;
      int rz = (nz - 1) / 2;
      x1 = i - rx, x2 = i + rx;
      y1 = j - ry, y2 = j + ry;
      z1 = k - rz, z2 = k + rz;
    }
    else if (OPTION("-landmarks") || OPTION("-points") || OPTION("-poly")) {
      PointSet landmarks;
      do {
        landmarks.Read(ARGUMENT);
        for (int i = 0; i < landmarks.Size(); ++i) {
          Point &p = landmarks(i);
          in->WorldToImage(p);
          x1 = MinIndex(x1, ifloor(p._x));
          x2 = MaxIndex(x2, iceil (p._x));
          y1 = MinIndex(y1, ifloor(p._y));
          y2 = MaxIndex(y2, iceil (p._y));
          z1 = MinIndex(z1, ifloor(p._z));
          z2 = MaxIndex(z2, iceil (p._z));
        }
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-ref")) {
      do {
        BinaryImage ref(ARGUMENT);
        double y, z;
        Point p;
        for (int k = 0; k < 2; ++k) {
          z = k * ref.Z();
          for (int j = 0; j < 2; ++j) {
            y = j * ref.Y();
            for (int i = 0; i < 2; ++i) {
              p._x = i * ref.X(), p._y = y, p._z = z;
              ref.ImageToWorld(p);
              in->WorldToImage(p);
              x1 = MinIndex(x1, ifloor(p._x));
              x2 = MaxIndex(x2, iceil (p._x));
              y1 = MinIndex(y1, ifloor(p._y));
              y2 = MaxIndex(y2, iceil (p._y));
              z1 = MinIndex(z1, ifloor(p._z));
              z2 = MaxIndex(z2, iceil (p._z));
            }
          }
        }
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-margin")) {
      xmargin = ymargin = ToRegionMargin(ARGUMENT);
      zmargin = tmargin = 0;
      if (HAS_ARGUMENT) ymargin = ToRegionMargin(ARGUMENT);
      if (HAS_ARGUMENT) zmargin = ToRegionMargin(ARGUMENT);
      if (HAS_ARGUMENT) tmargin = ToRegionMargin(ARGUMENT);
    }
    else if (OPTION("-scale")) {
      xscale = yscale = ToScaleFactor(ARGUMENT);
      zscale = tscale = 1.0;
      if (HAS_ARGUMENT) yscale = ToScaleFactor(ARGUMENT);
      if (HAS_ARGUMENT) zscale = ToScaleFactor(ARGUMENT);
      if (HAS_ARGUMENT) tscale = ToScaleFactor(ARGUMENT);
    }
    else if (OPTION("-crop")) {
      double bg = padding_value;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(bg);
      auto bbox = in->ForegroundDomain(bg, false);
      Point p(bbox._xorigin, bbox._yorigin, bbox._zorigin);
      in->WorldToImage(p);
      Point p1 = p - Point(bbox._x / 2, bbox._y / 2, bbox._z / 2);
      Point p2 = p + Point(bbox._x / 2, bbox._y / 2, bbox._z / 2);
      x1 = MinIndex(x1, ifloor(p1._x));
      x2 = MaxIndex(x2, iceil (p2._x));
      y1 = MinIndex(y1, ifloor(p1._y));
      y2 = MaxIndex(y2, iceil (p2._y));
      z1 = MinIndex(z1, ifloor(p1._z));
      z2 = MaxIndex(z2, iceil (p2._z));
    }
    else if (OPTION("-pad")) {
      pad = true;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(padding_value);
    }
    else if (OPTION("-nopad")) {
      pad = false;
    }
    else if (OPTION("-split")) {
      const string arg = ToLower(ARGUMENT);
      if      (arg == "x") split = 0;
      else if (arg == "y") split = 1;
      else if (arg == "z") split = 2;
      else if (arg == "t") split = 3;
      else if (!FromString(arg, split) || split < 0 || split > 3) {
        FatalError("Invalid -split option argument: " << arg);
      }
    }
    else HANDLE_BOOLEAN_OPTION("swap", swap_dims);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Full input image region by default
  if (x1 <= MIN_INDEX) x1 = 0;
  if (x2 >= MAX_INDEX) x2 = in->X() - 1;
  if (y1 <= MIN_INDEX) y1 = 0;
  if (y2 >= MAX_INDEX) y2 = in->Y() - 1;
  if (z1 <= MIN_INDEX) z1 = 0;
  if (z2 >= MAX_INDEX) z2 = in->Z() - 1;
  if (t1 <= MIN_INDEX) t1 = 0;
  if (t2 >= MAX_INDEX) t2 = in->T() - 1;

  // Add fixed-width margin
  x1 -= xmargin, x2 += xmargin;
  y1 -= ymargin, y2 += ymargin;
  z1 -= zmargin, z2 += zmargin;
  t1 -= tmargin, t2 += tmargin;

  // Swap indices if necessary
  if (x1 > x2) swap(x1, x2);
  if (y1 > y2) swap(y1, y2);
  if (z1 > z2) swap(z1, z2);
  if (t1 > t2) swap(t1, t2);

  // Scale region
  ScaleInterval(x1, x2, xscale);
  ScaleInterval(y1, y2, yscale);
  ScaleInterval(z1, z2, zscale);
  ScaleInterval(t1, t2, tscale);

  // Limit to extend of input region if padding not requested
  if (!pad) {
    x1 = max(x1, 0);
    y1 = max(y1, 0);
    z1 = max(z1, 0);
    t1 = max(t1, 0);
    x2 = min(x2, in->X() - 1);
    y2 = min(y2, in->Y() - 1);
    z2 = min(z2, in->Z() - 1);
    t2 = min(t2, in->T() - 1);
  }

  // Verbose reporting/logging of resulting region
  if (verbose) {
    cout <<  "-Rx1 " << x1 << " -Rx2 " << x2 << " -Ry1 " << y1 << " -Ry2 " << y2;
    cout << " -Rz1 " << z1 << " -Rz2 " << z2 << " -Rt1 " << t1 << " -Rt2 " << t2;
    cout << endl;
  }

  // Allocate output image
  ImageAttributes attr = in->Attributes();
  attr._x = abs(x2 - x1 + 1);
  attr._y = abs(y2 - y1 + 1);
  attr._z = abs(z2 - z1 + 1);
  attr._t = abs(t2 - t1 + 1);
  attr._xorigin = .0;
  attr._yorigin = .0;
  attr._zorigin = .0;
  attr._torigin = in->ImageToTime(t1);
  UniquePtr<BaseImage> out(BaseImage::New(in->GetDataType()));
  out->Initialize(attr);

  // Adjust spatial origin (i.e., image center)
  double o1[3] = {.0};
  double o2[3] = {static_cast<double>(x1), static_cast<double>(y1), static_cast<double>(z1)};
  out->ImageToWorld(o1[0], o1[1], o1[2]);
  in ->ImageToWorld(o2[0], o2[1], o2[2]);
  out->PutOrigin(o2[0] - o1[0], o2[1] - o1[1], o2[2] - o1[2]);

  // Copy image region
  int idx = 0;
  for (int l = t1; l <= t2; ++l)
  for (int k = z1; k <= z2; ++k)
  for (int j = y1; j <= y2; ++j)
  for (int i = x1; i <= x2; ++i, ++idx) {
    if (in->IsInside(i, j, k, l)) {
      out->PutAsDouble(idx, in->Get(i, j, k, l));
    } else {
      out->PutAsDouble(idx, padding_value);
    }
  }

  // Write output region
  if (split == -1) {
    out->Write(output_name);
  } else {
    // Format string for output file paths
    string output_path(output_name);
    auto pos = output_path.find('%');
    auto end = pos;
    if (pos != string::npos) {
      end = output_path.find('d', pos);
    }
    if (pos == end) {
      const int dims[] = {out->X(), out->Y(), out->Z(), out->T()};
      auto name = FilePrefix(output_path);
      auto ext  = Extension(output_path);
      string fmt = "%03d";
      if      (dims[split] <  10) fmt = "%d";
      else if (dims[split] < 100) fmt = "%02d";
      else if (dims[split] > 999) fmt = "%04d";
      output_path = name + string("_") + fmt + ext;
    }
    // Buffer for output file paths
    size_t max_image_path_length = output_path.length() + 10;
    UniquePtr<char> image_path(new char[max_image_path_length + 1]);
    // Swap image dimensions
    if (swap_dims && split < 2) {
      if (split == 0) out->SwapXY();
      if (split <= 1) out->SwapYZ();
      split = 2;
    }
    // Allocate memory for sub-regions
    double origin[3];
    int dims[] = {out->X(), out->Y(), out->Z(), out->T()};
    auto attr = out->Attributes();
    if      (split == 0) attr._x = 1;
    else if (split == 1) attr._y = 1;
    else if (split == 2) attr._z = 1;
    else if (split == 3) attr._t = 1;
    UniquePtr<BaseImage> sub(BaseImage::New(out->GetDataType()));
    sub->Initialize(attr);
    // Save individual sub-regions
    int outidx[4], subidx[4];
    int perm[] = {0, 1, 2, 3};  // permutation of for loops
    if (split == 0) swap(perm[0], perm[1]);
    if (split <= 1) swap(perm[1], perm[2]);
    if (split <= 2) swap(perm[2], perm[3]);
    for (outidx[perm[3]] = 0; outidx[perm[3]] < dims[perm[3]]; ++outidx[perm[3]]) {
      // Set origin of sub-region image
      if (split < 3) {
        origin[0] = .5 * static_cast<double>(dims[0] - 1);
        origin[1] = .5 * static_cast<double>(dims[1] - 1);
        origin[2] = .5 * static_cast<double>(dims[2] - 1);
        origin[split] = static_cast<double>(outidx[perm[3]]);
        out->ImageToWorld(origin[0], origin[1], origin[2]);
        sub->PutOrigin(origin[0], origin[1], origin[2]);
      } else {
        sub->PutTOrigin(out->ImageToTime(outidx[perm[3]]));
      }
      // Copy sub-region to output image
      double value;
      for (outidx[perm[2]] = 0; outidx[perm[2]] < dims[perm[2]]; ++outidx[perm[2]])
      for (outidx[perm[1]] = 0; outidx[perm[1]] < dims[perm[1]]; ++outidx[perm[1]])
      for (outidx[perm[0]] = 0; outidx[perm[0]] < dims[perm[0]]; ++outidx[perm[0]]) {
        value = out->GetAsDouble(outidx[0], outidx[1], outidx[2], outidx[3]);
        memcpy(subidx, outidx, 4 * sizeof(int));
        subidx[split] = 0;
        sub->PutAsDouble(subidx[0], subidx[1], subidx[2], subidx[3], value);
      }
      // Write sub-region to output image
      snprintf(image_path.get(), max_image_path_length, output_path.c_str(), outidx[perm[3]]);
      sub->Write(image_path.get());

    }
  }

  return 0;
}
