/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2008-2015 Imperial College London
 * Copyright 2008-2013 Daniel Rueckert, Julia Schnabel
 * Copyright 2013-2016 Andreas Schuh
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
#include "mirtk/Resampling.h"
#include "mirtk/LinearInterpolateImageFunction.h"

#include "mirtk/CityBlockDistanceTransform.h"
#include "mirtk/EuclideanDistanceTransform.h"

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
  cout << "Description:\n";
  cout << "  Computes distance transformation of binary object mask.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Binary object mask.\n";
  cout << "  output   Distance transform.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -distance <name>   Name of the distance transform. (default: euclidean)\n";
  cout << "                     - ``cityblock``: Manhatten city block (L1) distance.\n";
  cout << "                     - ``euclidean``: Euclidean distance.\n";
  cout << "                     - ``laplacian``: Laplacian distance.\n";
  cout << "\n";
  cout << "Euclidean distance transform options:\n";
  cout << "  -2D                Compute distance transform in 2D.\n";
  cout << "  -3D                Compute distance transform in 3D. (default)\n";
  cout << "  -radial            Compute radial distance transform. (default: off)\n";
  cout << "  -tradial           Compute radial distance transform in t dimension only.\n";
  cout << "  -isotropic [<s>]   Resample distance transform to isotropic voxel size,\n";
  cout << "                     where the voxel size is s times the minimum voxel size.\n";
  cout << "                     When no argument given, s defaults to 1. (default: 0/off)\n";
  cout << "  -outside <file>    Binary mask with non-zero values at voxels which are outside\n";
  cout << "                     the region for which distance values are to be computed.\n";
  cout << "                     Required by ``laplacian`` :option:`-distance`, where it defines\n";
  cout << "                     the boundary voxels with a maximum distance value of 1.\n";
  PrintCommonOptions(cout);
  cout << "\n";
}

// =============================================================================
// Auxiliaries
// =============================================================================

namespace mirtk {


/// Type of distance transform
enum DistanceTransformType
{
  DT_Unknown,
  DT_Euclidean,
  DT_CityBlock,
  DT_Laplacian,
  DT_Last
};

// -----------------------------------------------------------------------------
/// Convert distance transform type enumeration value to string
template <>
inline string ToString(const DistanceTransformType &value, int w, char c, bool left)
{
  const char *str;
  switch (value) {
    case DT_Euclidean: str = "euclidean"; break;
    case DT_CityBlock: str = "cityblock"; break;
    case DT_Laplacian: str = "laplacian"; break;
    default:           str = "unknown";   break;
  }
  return ToString(str, w, c, left);
}

// -----------------------------------------------------------------------------
/// Convert string to distance transform type enumeration value
template <>
inline bool FromString(const char *str, DistanceTransformType &value)
{
  const string lstr = ToLower(str);
  value = static_cast<DistanceTransformType>(DT_Last - 1);
  while (value != DT_Unknown) {
    if (ToString(value) == lstr) break;
    value = static_cast<DistanceTransformType>(value - 1);
  }
  return value != DT_Unknown;
}

// =============================================================================
// Laplacian distance map -- TODO: Should be a LaplacianDistanceTransform class
// =============================================================================

// -----------------------------------------------------------------------------
/// Voxel function for parallel iterative application of discrete Laplace operator
class LaplacianFunction : public VoxelReduction
{
  RealImage *_Distance;
  int        _Changed;
  double     _SumDelta;
  double     _MaxDelta;

public:

  LaplacianFunction(RealImage *distance)
  :
    _Distance(distance), _Changed(0), _SumDelta(0.), _MaxDelta(.0)
  {}

  void split(const LaplacianFunction &other)
  {
    _Changed  = 0;
    _SumDelta = 0.;
    _MaxDelta = 0.;
  }

  void join(const LaplacianFunction &other)
  {
    _Changed  += other._Changed;
    _SumDelta += other._SumDelta;
    _MaxDelta = max(_MaxDelta, other._MaxDelta);
  }

  double AvgDelta() const
  {
    return (_Changed > 0 ? _SumDelta / _Changed : 0.);
  }

  double MaxDelta() const
  {
    return _MaxDelta;
  }

  void Reset()
  {
    _Changed  = 0;
    _SumDelta = 0.;
    _MaxDelta = 0.;
  }

  void operator ()(int i, int j, int k, int, const BinaryPixel *isobj, const BinaryPixel *isout, RealPixel *distance)
  {
    if (*isobj || *isout) return;

    RealPixel val(0);
    int       num = 0;

    if (!IsNaN(_Distance->Get(i-1, j, k))) val += _Distance->Get(i-1, j, k), ++num;
    if (!IsNaN(_Distance->Get(i+1, j, k))) val += _Distance->Get(i+1, j, k), ++num;
    if (!IsNaN(_Distance->Get(i, j-1, k))) val += _Distance->Get(i, j-1, k), ++num;
    if (!IsNaN(_Distance->Get(i, j+1, k))) val += _Distance->Get(i, j+1, k), ++num;
    if (!IsNaN(_Distance->Get(i, j, k-1))) val += _Distance->Get(i, j, k-1), ++num;
    if (!IsNaN(_Distance->Get(i, j, k+1))) val += _Distance->Get(i, j, k+1), ++num;

    if (num > 0) {
      val /= num;
      RealPixel prev = *distance;
      if (IsNaN(prev)) prev = RealPixel(0);
      *distance = val;
      double delta = abs(prev - val);
      if (delta > _MaxDelta) _MaxDelta = delta;
      _SumDelta += delta;
      ++_Changed;
    }
  }
};

// -----------------------------------------------------------------------------
/// Computes Laplacian distance map from object boundary to outside voxels
RealImage ComputeLaplacianTransform(const BinaryImage &inside, BinaryImage &outside, int maxiter)
{
  MIRTK_START_TIMING();

  const auto prev_flags     = cout.flags(ios::fixed);
  const auto prev_precision = cout.precision(6);

  const ImageAttributes &attr = inside.Attributes();
  if (!outside.HasSpatialAttributesOf(&inside)) {
    FatalError("Outside mask image must have attributes identical to input object mask!");
  }

  // Initialize distance image
  RealImage init(attr);
  for (int k = 0; k < attr._z; ++k)
  for (int j = 0; j < attr._y; ++j)
  for (int i = 0; i < attr._x; ++i) {
    if (i == 0 || j == 0 || k == 0 || i == attr._x-1 || j == attr._y-1 || k == attr._z-1) {
      outside(i, j, k) = BinaryPixel(1);
    }
    if      (inside (i, j, k)) init(i, j, k) = RealPixel(0.);
    else if (outside(i, j, k)) init(i, j, k) = RealPixel(1.);
    else                       init(i, j, k) = RealPixel(NaN);
  }

  // Iteratively apply discrete Laplace operator
  RealImage temp(attr);
  RealImage *input = &init, *output = &temp;
  for (int iter = 0; iter < maxiter; ++iter) {
    LaplacianFunction op(input);
    ParallelForEachVoxel(attr, inside, outside, *output, op);
    swap(input, output); // before breaking the loop!
    if (verbose && (iter+1) % 10 == 0) {
      cout << setw(3)<< right << (iter+1) << left
           << " Delta: max = " << op.MaxDelta()
           << ", avg = " << op.AvgDelta() << endl;
    }
    if (op.MaxDelta() < 1e-12) break;
  }
  swap(input, output);

  cout.flags(prev_flags);
  cout.precision(prev_precision);

  MIRTK_DEBUG_TIMING(1, "computing Laplacian map");
  return *output;
}


} // namespace mirtk

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  typedef EuclideanDistanceTransform<RealPixel> EuclideanDistanceTransformType;
  typedef CityBlockDistanceTransform<RealPixel> CityBlockDistanceTransformType;

  InitializeIOLibrary();

  EXPECTS_POSARGS(2);

  const char *input_name   = POSARG(1);
  const char *output_name  = POSARG(2);
  const char *outside_name = nullptr;
 
  DistanceTransformType type = DT_Euclidean;

  EuclideanDistanceTransformType::Mode euclidean_mode;
  euclidean_mode = EuclideanDistanceTransformType::DT_3D;
  
  int    radial    = 0;
  double isotropic = .0;
  int    maxiter   = 5000;

  for (ALL_OPTIONS) {
    if (OPTION("-distance") || OPTION("-mode")) {
      PARSE_ARGUMENT(type);
    }
    else if (OPTION("-2D")) {
      euclidean_mode = EuclideanDistanceTransformType::DT_2D;
    }
    else if (OPTION("-3D")) {
      euclidean_mode = EuclideanDistanceTransformType::DT_3D;
    }
    else if (OPTION("-radial")) {
      radial = 1;
    }
    else if (OPTION("-tradial")) {
      radial = 2;
    }
    else if (OPTION("-isotropic")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(isotropic);
      else isotropic = 1.0;
    }
    else if (OPTION("-outside")) {
      outside_name = ARGUMENT;
    }
    else if (OPTION("-max-iterations") || OPTION("-iterations") ||
             OPTION("-max-iter")       || OPTION("-iter")) {
      PARSE_ARGUMENT(maxiter);
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // No linear isotropic resampling for city block distance
  if (type == DT_CityBlock) isotropic = .0;

  // Read object image
  if (verbose) cout << "Reading image " << input_name << "...", cout.flush();
  RealImage image(input_name);
  if (verbose) cout << " done" << endl;

  // Compute distance transform
  if (verbose) cout << "Computing " << ToString(type) << " distance transform..." << endl;
  RealImage dmap;

  switch (type) {

    // City Block distance transform
    case DT_CityBlock: {
      CityBlockDistanceTransformType cbdt;
      RealImage &outputA = dmap, outputB;

      if (verbose) cout << "  Computing distance transform of interior...", cout.flush();

      cbdt.Input(&image);
      cbdt.Output(&outputA);
      cbdt.Run();

      if (verbose) cout << " done\n  Computing distance transform of exterior...", cout.flush();

      for (int idx = 0; idx < image.NumberOfVoxels(); ++idx) {
        image(idx) = (image(idx) > .0f ? 0 : 1);
      }

      cbdt.Output(&outputB);
      cbdt.Run();

      if (verbose) cout << " done" << endl;

      for (int idx = 0; idx < image.NumberOfVoxels(); ++idx) {
        dmap(idx) = outputB(idx) - outputA(idx);
      }

    } break;

    // Euclidean distance transform
    case DT_Euclidean: {
      EuclideanDistanceTransformType edt(euclidean_mode);
      RealImage &inputA  = image, inputB(image.Attributes());
      RealImage &outputA = dmap,  outputB;

      for (int idx = 0; idx < image.NumberOfVoxels(); ++idx) {
        if (image(idx) > .5) {
          inputA(idx) = 1.0;
          inputB(idx) = 0.0;
        } else {
          inputA(idx) = 0.0;
          inputB(idx) = 1.0;
        }
      }

      if (verbose) cout << "  Computing distance transform of interior...", cout.flush();

      edt.Input (& inputA);
      edt.Output(&outputA);
      edt.Run();

      if (verbose) cout << " done\n  Computing distance transform of exterior...", cout.flush();

      edt.Input (& inputB);
      edt.Output(&outputB);
      edt.Run();

      if (verbose) cout << " done" << endl;

      for (int idx = 0; idx < image.NumberOfVoxels(); ++idx) {
        dmap(idx) = sqrt(outputA(idx)) - sqrt(outputB(idx));
      }

      if (radial > 0) {
        edt.Input (&dmap);
        edt.Output(&dmap);
        switch (radial) {
          case 1: edt.Radial(); break;
          case 2: edt.TRadial(); break;
        }
      }
    } break;
 
    // Laplacian map
    case DT_Laplacian: {
      BinaryImage inside(image), outside;
      if (outside_name) {
        outside.Read(outside_name);
      } else {
        outside.Initialize(image.Attributes());
      }
      dmap = ComputeLaplacianTransform(inside, outside, maxiter);
    } break;

    default:
      FatalError("Unknown distance transform type: " << type);
  }

  // Resample image to isotropic voxels (smallest voxel dimension)
  if (isotropic > .0) {
    double ds = min(min(dmap.GetXSize(), dmap.GetYSize()), dmap.GetZSize()) * isotropic;
    if (verbose) cout << "Resampling distance map to isotropic voxel size: " << ds << " [mm]...", cout.flush();
    Resampling<RealPixel> resampling(ds, ds, ds);
    GenericLinearInterpolateImageFunction<GenericImage<RealPixel> > interpolator;
    resampling.Interpolator(&interpolator);
    resampling.Input (&dmap);
    resampling.Output(&dmap);
    resampling.Run();
    if (verbose) cout << " done" << endl;
  }

  if (verbose) cout << "Writing distance transform to " << output_name << "...";
  dmap.Write(output_name);
  if (verbose) cout << " done" << endl;

  return 0;
}
