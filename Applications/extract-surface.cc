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

#include <mirtkIOConfig.h>

#include <mirtkPointSetUtils.h>
#include <mirtkImplicitSurfaceUtils.h>

using namespace mirtk;
using namespace mirtk::ImplicitSurfaceUtils;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Extract the isosurface from an intensity image or segmentation." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input image." << endl;
  cout << "  output   Output surface mesh." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -isovalue <value>     Isovalue of surface intensities. (default: 0)" << endl;
  cout << "  -blur <sigma>         Blur input image with kernel size sigma before running filter. (default: 0)" << endl;
  cout << "  -isotropic            Resample image to an isotropic voxel size (minimum of input voxel size)." << endl;
  cout << "  -normals [on|off]     Choose whether to generate normals (default) or not." << endl;
  cout << "  -gradients [on|off]   Choose whether to generate gradients or not (default)." << endl;
  cout << "  -close [on|off]       Put zeros around the image to generate a closed surface(s)."<<endl;
  cout << "  -[no]compress         Whether to compress output .vtp file. (default: on)" << endl;
  cout << "  -binary               Write binary data when output file name extension is .vtk. (default: on)" << endl;
  cout << "  -ascii                Write ASCII  data when output file name extension is .vtk. (default: off)" << endl;
  PrintCommonOptions(cout);
  cout << endl;
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

  double isovalue  = .0;
  bool   compress  = true;
  bool   ascii     = false;
  bool   isotropic = false;
  bool   close     = false;
  bool   gradients = false;
  bool   normals   = true;
  double blurring  = .0;

  for (ALL_OPTIONS) {
    if (OPTION("-isovalue")) PARSE_ARGUMENT(isovalue);
    else if (OPTION("-close")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(close);
      else close = true;
    }
    else if (OPTION("-isotropic")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(isotropic);
      else isotropic = true;
    }
    else if (OPTION("-gradients")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(gradients);
      else gradients = true;
    }
    else if (OPTION("-nogradients")) {
      gradients = false;
    }
    else if (OPTION("-normals")) {
      if (HAS_ARGUMENT) PARSE_ARGUMENT(normals);
      else normals = true;
    }
    else if (OPTION("-nonormals")) {
      normals = false;
    }
    else if (OPTION("-blur")) PARSE_ARGUMENT(blurring);
    else if (OPTION("-compress"))   compress = true;
    else if (OPTION("-nocompress")) compress = false;
    else if (OPTION("-ascii"))      ascii    = true;
    else if (OPTION("-binary"))     ascii    = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  InitializeIOLibrary();
  DistanceImage image(input_name);
  if (image.Z() == 1 || image.T() > 1) {
    FatalError("Input image must be three-dimensional!");
  }

  vtkSmartPointer<vtkPolyData> isosurface;
  isosurface = Isosurface(image, isovalue, blurring, isotropic, close, normals, gradients);

  if (!WritePolyData(output_name, isosurface, compress, ascii)) {
    FatalError("Failed to write surface to " << output_name);
  }

  return 0;
}
