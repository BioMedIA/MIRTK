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

#include "mirtk/Vtk.h"
#include "mirtk/IOConfig.h"
#include "mirtk/GenericImage.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/Transformation.h"
#include "mirtk/PointSetUtils.h"

#include "vtkSmartPointer.h"
#include "vtkDataSetReader.h"
#include "vtkXMLGenericDataObjectReader.h"
#include "vtkDataSet.h"
#include "vtkImageData.h"
#include "vtkPointSet.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkPolyDataNormals.h"
#include "vtkImageDataToPointSet.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " [options]" << endl;
  cout << "       " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Applies one or more transformations to a set of points. Input point" << endl;
  cout << "  x, y, and z coordinates are either read from STDIN (space separated) or from a" << endl;
  cout << "  point set file. The corresponding transformed coordinates are then written" << endl;
  cout << "  either to STDOUT or an output point set file, respectively. If multiple" << endl;
  cout << "  transformations are specified, these are applied in the order as they appear" << endl;
  cout << "  on the command line." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -dofin <file>...   Transformation file. (default: identity)" << endl;
  cout << "  -invert            Invert preceding transformation. (default: no)" << endl;
  cout << "  -source <image>    Reference source image. (default: none)" << endl;
  cout << "  -target <image>    Reference target image. (default: none)" << endl;
  cout << "  -partial <n>       Transform the first n points and copy the rest. (default: #points)" << endl;
  cout << "  -normals           Whether to recompute normals of input surface mesh. (default: off)" << endl;
  cout << "  -point-normals     Whether to recompute point normals of input surface mesh. (default: off)" << endl;
  cout << "  -cell-normals      Whether to recompute cell  normals of input surface mesh. (default: off)" << endl;
  cout << "  -[no]compress      Enable/disable compression of XML VTK file. (default: on)" << endl;
  cout << "  -ascii             Set file type of output legacy VTK file to ASCII.  (default: input file type)" << endl;
  cout << "  -binary            Set file type of output legacy VTK file to BINARY. (default: input file type)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  REQUIRES_POSARGS(0);

  // Nevertheless, show help if called without arguments
  if (argc == 1) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  // Initialize image I/O library
  InitializeIOLibrary();

  vtkSmartPointer<vtkPointSet> pointset;
  PointSet                     points;

  // Positional arguments
  const char *input_name  = NULL;
  const char *output_name = NULL;

  if (NUM_POSARGS == 2) {
    input_name  = POSARG(1);
    output_name = POSARG(2);
  } else if (NUM_POSARGS != 0) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  // Optional arguments
  Array<const char *> dofin_name;
  Array<bool>         dofin_invert;
  GreyImage           target, source;
  bool compute_point_normals = false;
  bool compute_cell_normals  = false;
  int  pnumber               = 0;
  int  output_file_type      = -1;
  bool compress              = true;
  double tt = .0, ts = 1.0; // Temporal origin, used by velocity based transformations

  for (ALL_OPTIONS) {
    if (OPTION("-dofin")) {
      do {
        dofin_name.push_back(ARGUMENT);
      } while (HAS_ARGUMENT);
      while (dofin_invert.size() < dofin_name.size()) {
        dofin_invert.push_back(false);
      }
    }
    else if (OPTION("-invert")) {
      if (dofin_invert.empty()) dofin_invert.push_back(true);
      else dofin_invert.back() = true;
    }
    else if (OPTION("-Tt")) tt = atof(ARGUMENT);
    else if (OPTION("-St")) ts = atof(ARGUMENT);
    else if (OPTION("-source")) source.Read(ARGUMENT);
    else if (OPTION("-target")) target.Read(ARGUMENT);
    else if (OPTION("-partial")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, pnumber)) {
        cerr << "Failed to convert -partial argument " << arg << " to integer" << endl;
        exit(1);
      }
    }
    else if (OPTION("-normals")) compute_point_normals = compute_cell_normals = true;
    else if (OPTION("-point-normals")) compute_point_normals = true;
    else if (OPTION("-cell-normals")) compute_cell_normals = true;
    else if (OPTION("-compress")) compress = true;
    else if (OPTION("-nocompress")) compress = false;
    else if (OPTION("-ascii"))  output_file_type = VTK_ASCII;
    else if (OPTION("-binary")) output_file_type = VTK_BINARY;
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }
  dofin_invert.resize(dofin_name.size(), false);

  // Read input points
  double p[3];
  if (input_name) {
    pointset = ReadPointSet(input_name, NULL, false);
    if (!pointset) {
      vtkSmartPointer<vtkDataSet> dataset;
      const string ext = Extension(input_name);
      if (ext == ".vtk") {
        vtkNew<vtkDataSetReader> reader;
        reader->SetFileName(input_name);
        reader->Update();
        dataset = reader->GetOutput();
        if (output_file_type == -1) output_file_type = reader->GetFileType();
      } else {
        vtkNew<vtkXMLGenericDataObjectReader> reader;
        reader->SetFileName(input_name);
        reader->Update();
        dataset = vtkDataSet::SafeDownCast(reader->GetOutput());
      }
      vtkImageData *image = vtkImageData::SafeDownCast(dataset);
      if (!image) {
        FatalError("Input VTK dataset must be either of type vtkPointSet or vtkImageData");
      }
      vtkNew<vtkImageDataToPointSet> converter;
      SetVTKInput(converter, image);
      converter->Update();
      pointset = converter->GetOutput();
    }
    if (pointset->GetNumberOfPoints() == 0) {
      FatalError("Input VTK data set has no points");
    }
    points.Resize(pointset->GetNumberOfPoints());
    for (vtkIdType i = 0; i < pointset->GetNumberOfPoints(); ++i) {
      pointset->GetPoint(i, p);
      points(i) = p;
    }
  } else {
    if (output_file_type == -1) output_file_type = VTK_BINARY;
    while (cin) {
      cin >> p[0] >> p[1] >> p[2];
      if (!cin) break;
      points.Add(p);
    }
    points.ShrinkToFit();
  }
  if (pnumber <= 0) pnumber = points.Size();

  // Map all points from target world to source world
  if (!target.IsEmpty() && !source.IsEmpty()) {
    for (int i = 0; i < points.Size(); ++i) {
      Point &pt = points(i);
      target.WorldToImage(pt);
      source.ImageToWorld(pt);
    }
  }

  // Transform first pnumber points
  for (size_t i = 0; i < dofin_name.size(); ++i) {
    unique_ptr<Transformation> dofin(Transformation::New(dofin_name[i]));
    if (dofin_invert[i]) {
      if (verbose) cout << "Apply inverse of " << dofin_name[i] << endl;
      for (int i = 0; i < pnumber; ++i) dofin->Inverse(points(i), ts, tt);
    } else {
      if (verbose) cout << "Apply " << dofin_name[i] << endl;
      for (int i = 0; i < pnumber; ++i) dofin->Transform(points(i), ts, tt);
    }
  }

  if (pointset) {

    // Set output points
    vtkPoints *output_points = pointset->GetPoints();
    for (int i = 0; i < pnumber; ++i) {
      const Point &pt = points(i);
      output_points->SetPoint(i, pt._x, pt._y, pt._z);
    }

    // Generate surface normals
    if (pointset->GetPointData()->HasArray("Normals")) compute_point_normals = true;
    if (pointset->GetCellData ()->HasArray("Normals")) compute_cell_normals  = true;

    if (compute_point_normals || compute_cell_normals) {
      vtkSmartPointer<vtkPolyDataNormals> normals;
      normals = vtkSmartPointer<vtkPolyDataNormals>::New();
      SetVTKInput(normals, pointset);
      normals->SplittingOff();
      normals->ConsistencyOn();
      normals->SetComputePointNormals(compute_point_normals);
      normals->SetComputeCellNormals (compute_cell_normals);
      normals->Update();
      pointset = normals->GetOutput();
    }

    // Write the final dataset
    WritePointSet(output_name, pointset, compress, output_file_type == VTK_ASCII);

  } else {

    // Print transformed points
    for (int i = 0; i < points.Size(); ++i) {
      const Point &pt = points(i);
      cout << pt._x << " " << pt._y << " " << pt._z << endl;
    }

  }

  return 0;
}
