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

#include "mirtk/IOConfig.h"

#include "mirtk/PointSetIO.h"
#include "mirtk/ImplicitSurfaceUtils.h"

#include "vtkNew.h"
#include "vtkPolyDataNormals.h"

using namespace mirtk;
using namespace mirtk::ImplicitSurfaceUtils;


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
  cout << "  Extract the isosurface from an intensity image or segmentation.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input image.\n";
  cout << "  output   Output surface mesh.\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -isovalue <value>       Isovalue of surface intensities. (default: 0)\n";
  cout << "  -blur, -sigma <sigma>   Blur input image with kernel size sigma before running filter. (default: 0)\n";
  cout << "  -isotropic              Resample image to an isotropic voxel size (minimum of input voxel size). (default: off)\n";
  cout << "  -[no]normals            Whether to calculate surface normals. (default: on)\n";
  cout << "  -[no]gradients          Whether to calculate image gradients. (default: off)\n";
  cout << "  -[no]close              Put zeros around the image to generate a closed surface(s). (default: off)\n";
  cout << "  -[no]compress           Whether to compress output .vtp file. (default: on)\n";
  cout << "  -binary                 Write binary data when output file name extension is .vtk. (default: on)\n";
  cout << "  -ascii                  Write ASCII  data when output file name extension is .vtk. (default: off)\n";
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
  FileOption output_fopt  = FO_Default;

  double isovalue  = .0;
  bool   isotropic = false;
  bool   close     = false;
  bool   gradients = false;
  bool   normals   = true;
  double blurring  = .0;

  for (ALL_OPTIONS) {
    if (OPTION("-isovalue")) {
      PARSE_ARGUMENT(isovalue);
    }
    else if (OPTION("-blur") || OPTION("-sigma")) {
      PARSE_ARGUMENT(blurring);
    }
    else HANDLE_BOOL_OPTION(close);
    else HANDLE_BOOL_OPTION(isotropic);
    else HANDLE_BOOL_OPTION(gradients);
    else HANDLE_BOOL_OPTION(normals);
    else HANDLE_POINTSETIO_OPTION(output_fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  InitializeIOLibrary();
  DistanceImage image(input_name);
  if (image.Z() == 1 || image.T() > 1) {
    FatalError("Input image must be three-dimensional!");
  }

  vtkSmartPointer<vtkPolyData> isosurface;
  isosurface = Isosurface(image, isovalue, blurring, isotropic, close, false, gradients);

  if (normals) {
    vtkNew<vtkPolyDataNormals> filter;
    filter->SplittingOff();
    filter->ConsistencyOn();
    filter->AutoOrientNormalsOn();
    filter->ComputePointNormalsOn();
    filter->ComputeCellNormalsOff();
    filter->SetInputData(isosurface);
    filter->Update();
    isosurface = filter->GetOutput();
  }

  if (gradients) {
    isosurface->GetPointData()->GetVectors()->SetName("Gradient");
  }

  if (!WritePolyData(output_name, isosurface, output_fopt)) {
    FatalError("Failed to write surface to " << output_name);
  }

  return 0;
}
