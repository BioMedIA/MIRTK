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

#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkImageIOConfig.h>
#include <mirtkResampling.h>
#include <mirtkLinearInterpolateImageFunction.h>

#include <mirtkCityBlockDistanceTransform.h>
#include <mirtkEuclideanDistanceTransform.h>

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
  cout << "\n";
  cout << "Euclidean distance transform options:\n";
  cout << "  -2D                Compute distance transform in 2D.\n";
  cout << "  -3D                Compute distance transform in 3D. (default)\n";
  cout << "  -radial            Compute radial distance transform. (default: off)\n";
  cout << "  -tradial           Compute radial distance transform in t dimension only.\n";
  cout << "  -isotropic [<s>]   Resample distance transform to isotropic voxel size,\n";
  cout << "                     where the voxel size is s times the minimum voxel size.\n";
  cout << "                     When no argument given, s defaults to 1. (default: 0/off)\n";
  PrintStandardOptions(cout);
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
    default:           str = "unknown"; break;
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


} // namespace mirtk

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  typedef EuclideanDistanceTransform<RealPixel> EuclideanDistanceTransformType;
  typedef CityBlockDistanceTransform<RealPixel> CityBlockDistanceTransformType;

  InitializeImageIOLibrary();

  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);
 
  DistanceTransformType type = DT_Euclidean;

  EuclideanDistanceTransformType::Mode euclidean_mode;
  euclidean_mode = EuclideanDistanceTransformType::DT_3D;
  
  int    radial    = 0;
  double isotropic = .0;

  for (ALL_OPTIONS) {
    if (OPTION("-distance") || OPTION("-mode")) {
      PARSE_ARGUMENT(type);
    }
    else if (OPTION("-2D")) {
      euclidean_mode = EuclideanDistanceTransformType::DT_2D;
    }
    else if (OPTION("-3D")) {
      euclidean_mode = EuclideanDistanceTransformType::DT_2D;
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

      if (verbose) cout << "  Computing distance transform of interior..." << endl;

      cbdt.Input(&image);
      cbdt.Output(&outputA);
      cbdt.Run();

      if (verbose) cout << " done\n  Computing distance transform of exterior..." << endl;

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
      RealImage inputA, inputB, &outputA = dmap, outputB;

      for (int idx = 0; idx < image.NumberOfVoxels(); ++idx) {
        if (image(idx) > .5) {
          inputA(idx) = 1.0;
          inputB(idx) = 0.0;
        } else {
          inputA(idx) = 0.0;
          inputB(idx) = 1.0;
        }
      }

      if (verbose) cout << "  Computing distance transform of interior..." << endl;

      edt.Input (& inputA);
      edt.Output(&outputA);
      edt.Run();

      if (verbose) cout << " done\n  Computing distance transform of exterior..." << endl;

      edt.Input (& inputB);
      edt.Output(&outputB);
      edt.Run();

      if (verbose) cout << " done" << endl;

      for (int idx = 0; idx < image.NumberOfVoxels(); ++idx) {
        dmap(idx) = sqrt(outputA(idx)) - sqrt(outputA(idx));
      }

      if (radial > 0) {
        edt.Input (&dmap);
        edt.Output(&dmap);
        switch (radial) {
          case 1: edt.Radial(); break;
          case 2: edt.TRadial(); break;
        }
      }
    }
 
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
