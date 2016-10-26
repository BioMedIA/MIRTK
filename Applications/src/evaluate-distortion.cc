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

#include "mirtk/PointSetIO.h"
#include "mirtk/DataStatistics.h"
#include "mirtk/Numeric.h"

// include of mirtk/vtkMath.h with fix for VTK 6.0 isinf ambiguity
// needed because vtkTriangle.h otherwise includes vtkMath.h
#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkCell.h"
#include "vtkTriangle.h"
#include "vtkFloatArray.h"

using namespace mirtk;
using namespace mirtk::data::statistic;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
/// Print help screen
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <original> <input> [<output>] [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Computes the distortion of a surface mesh under deformation, when\n";
  cout << "  mapped to the surface of another solid with equivalent topology, or\n";
  cout << "  flattened to the plane or sphere, respectively.\n";
  cout << "\n";
  cout << "  The distortion measures are stored as point/cell data of the output\n";
  cout << "  dataset. If no output file name is given, the mean and standard deviation\n";
  cout << "  are reported even when the verbose option -v is not given.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  original   File name of original dataset.\n";
  cout << "  input      File name of deformed/mapped dataset.\n";
  cout << "  output     File name of output dataset.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -angular\n";
  cout << "      Quantify angular distortion by computing the average ratio of normalized\n";
  cout << "      angles at each vertex before and after the mapping. This ratio is 1 when\n";
  cout << "      the angles of the original mesh are preserved in the input mesh.\n";
  cout << "\n";
  cout << "  -areal\n";
  cout << "      Compute areal distortion of cells: log_2(area(deformed)/area(original)).\n";
  cout << "\n";
  cout << "  -abs-areal\n";
  cout << "      Compute absolute value of :option:`-areal` distortion.\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of available distortion measures
enum DistortionMeasure
{
  Distortion_Angular,
  Distortion_Areal,
  Distortion_ArealAbs
};

// -----------------------------------------------------------------------------
double Angle(vtkPolyData *mesh, vtkIdType cellId, vtkIdType ptId)
{
  vtkIdType npts, *pts, i[3] = {0, 1, 2};
  mesh->GetCellPoints(cellId, npts, pts);
  for (vtkIdType j = 0; j < npts; ++j) {
    if (pts[j] == ptId) {
      i[0] = i[j];
      i[j] = 0;
      break;
    }
  }
  double p[3][3], ab, ac, bc;
  mesh->GetPoint(pts[0], p[i[0]]);
  mesh->GetPoint(pts[1], p[i[1]]);
  mesh->GetPoint(pts[2], p[i[2]]);
  ab = pow(p[1][0] - p[0][0], 2) + pow(p[1][1] - p[0][1], 2) + pow(p[1][2] - p[0][2], 2);
  ac = pow(p[2][0] - p[0][0], 2) + pow(p[2][1] - p[0][1], 2) + pow(p[2][2] - p[0][2], 2);
  bc = pow(p[2][0] - p[1][0], 2) + pow(p[2][1] - p[1][1], 2) + pow(p[2][2] - p[1][2], 2);
  return acos((ab + ac - bc) / (2.0 * sqrt(ab) * sqrt(ac)));
}

// -----------------------------------------------------------------------------
/// Evaluate angular distortion, i.e., average ratio of normalized angles at each point
vtkSmartPointer<vtkDataArray>
AngularDistortion(vtkPolyData *orig, vtkPolyData *mesh)
{
  const vtkIdType npoints = orig->GetNumberOfPoints();
  if (mesh->GetNumberOfPoints() != npoints) {
    FatalError("AngularDistortion: The two meshes have differing number of points!");
  }

  vtkSmartPointer<vtkDataArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetName("AngularDistortion");
  array->SetNumberOfComponents(2);
  array->SetNumberOfTuples(npoints);

  unsigned short ncells1, ncells2;
  vtkIdType      *cells1, *cells2;
  double         sum1,    sum2;
  Array<double>  angles1, angles2;
  double         ratio, avg_ratio, max_ratio;

  for (vtkIdType ptId = 0; ptId < npoints; ++ptId) {
    orig->GetPointCells(ptId, ncells1, cells1);
    mesh->GetPointCells(ptId, ncells2, cells2);
    if (ncells1 != ncells2) {
      FatalError("AngularDistortion: Point " << ptId << " belongs to different number of cells in the two meshes!");
    }
    angles1.resize(ncells1);
    angles2.resize(ncells2);
    for (unsigned short i = 0; i < ncells1; ++i) {
      if (cells1[i] != cells2[i]) {
        FatalError("AngularDistortion: Lists of cell IDs to which point " << ptId << " belongs differ!");
      }
      angles1[i] = Angle(orig, cells1[i], ptId);
      angles2[i] = Angle(mesh, cells2[i], ptId);
    }
    sum1 = Accumulate(angles1);
    sum2 = Accumulate(angles2);
    avg_ratio = .0;
    max_ratio = .0;
    for (unsigned short i = 0; i < ncells1; ++i) {
      ratio = (angles1[i] / sum1) / (angles2[i] / sum2);
      max_ratio = max(max_ratio, ratio);
      avg_ratio += ratio;
    }
    avg_ratio /= ncells1;
    array->SetComponent(ptId, 0, avg_ratio);
    array->SetComponent(ptId, 1, max_ratio);
  }

  return array;
}

// -----------------------------------------------------------------------------
/// Evaluate areal distortion for each cell
vtkSmartPointer<vtkDataArray>
ArealDistortion(vtkPolyData *orig, vtkPolyData *mesh)
{
  const double log2 = log(2);

  const vtkIdType ncells = orig->GetNumberOfCells();
  if (mesh->GetNumberOfCells() != ncells) {
    FatalError("ArealDistortion: The two meshes have differing number of cells!");
  }

  vtkSmartPointer<vtkDataArray> array = vtkSmartPointer<vtkFloatArray>::New();
  array->SetName("ArealDistortion");
  array->SetNumberOfComponents(1);
  array->SetNumberOfTuples(ncells);

  double distortion;
  for (vtkIdType cellId = 0; cellId < ncells; ++cellId) {
    vtkCell * const cell1 = orig->GetCell(cellId);
    vtkCell * const cell2 = mesh->GetCell(cellId);
    if (cell1->GetCellType() != cell2->GetCellType()) {
      FatalError("ArealDistortion: Cell " << cellId << " has differing type in the two meshes!");
    }
    if (cell1->GetNumberOfPoints() != cell2->GetNumberOfPoints()) {
      FatalError("ArealDistortion: Cell " << cellId << " has different number of points in the two meshes!");
    }
    switch (cell1->GetCellType()) {
      case VTK_TRIANGLE: {
          vtkTriangle * const triangle1 = dynamic_cast<vtkTriangle *>(cell1);
          vtkTriangle * const triangle2 = dynamic_cast<vtkTriangle *>(cell2);
          distortion = log(triangle2->ComputeArea() / triangle1->ComputeArea()) / log2;
          array->SetComponent(cellId, 0, distortion);
        } break;
      default:
        FatalError("ArealDistortion: Only triangular cells supported");
    }
  }

  return array;
}

// -----------------------------------------------------------------------------
/// Print mean and standard deviation of measured distortions
void PrintMeanAndSigma(vtkDataArray *data, int j = 0, const char *desc = nullptr)
{
  const int n = data->GetNumberOfTuples();
  double mean = .0, var = .0;
  int    m = 0;

  double delta, d;
  for (vtkIdType i = 0; i < n; ++i) {
    ++m;
    d = data->GetComponent(i, j);
    delta = d - mean;
    mean += delta / m;
    var  += delta * (d - mean);
  }

  if (m < 1) {
    mean = var = mirtk::nan;
  } else if (m < 2) {
    var = .0;
  } else {
    var /= m - 1;
  }

  string name(desc ? desc : CamelCaseToPrettyParameterName(data->GetName()));
  cout << name << ": mean = " << mean << ", sigma = " << sqrt(var) << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Positional arguments
  REQUIRES_POSARGS(2);

  const char *orig_name   = POSARG(1);
  const char *input_name  = POSARG(2);
  const char *output_name = (NUM_POSARGS > 2 ? POSARG(3) : nullptr);

  if (NUM_POSARGS > 3) FatalError("Too many positional arguments given");

  // Optional arguments
  Array<DistortionMeasure> measures;

  for (ALL_OPTIONS) {
    if      (OPTION("-angular"))   measures.push_back(Distortion_Angular);
    else if (OPTION("-areal"))     measures.push_back(Distortion_Areal);
    else if (OPTION("-abs-areal")) measures.push_back(Distortion_ArealAbs);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (measures.empty()) {
    measures.push_back(Distortion_Angular);
    measures.push_back(Distortion_Areal);
  }

  if (!output_name && verbose == 0) verbose = 1;

  // Read input files
  if (verbose > 1) cout << "Reading original surface from " << orig_name << endl;
  vtkSmartPointer<vtkPolyData> orig = ReadPolyData(orig_name);
  if (verbose > 1) cout << "Reading input surface from " << input_name << endl;
  vtkSmartPointer<vtkPolyData> mesh = ReadPolyData(input_name);

  // Check input
  vtkIdType npoints = orig->GetNumberOfPoints();
  if (npoints == 0) {
    cerr << "Original surface has no points!" << endl;
    exit(1);
  }
  if (mesh->GetNumberOfPoints() != npoints) {
    cerr << "Surface meshes have differing number of points!" << endl;
    exit(1);
  }
  vtkIdType ncells = orig->GetNumberOfCells();
  if (ncells == 0) {
    cerr << "Original surface has no cells!" << endl;
    exit(1);
  }
  if (mesh->GetNumberOfCells() != ncells) {
    cerr << "Surface meshes have differing number of cells!" << endl;
    exit(1);
  }

  // Build links
  orig->BuildLinks();
  mesh->BuildLinks();

  // Compute distortion measures
  vtkSmartPointer<vtkDataArray> angular_distortion;
  vtkSmartPointer<vtkDataArray> areal_distortion;
  vtkSmartPointer<vtkDataArray> abs_areal_distortion;

  for (auto it = measures.begin(); it != measures.end(); ++it) {
    switch (*it) {
      case Distortion_Angular: {
        if (!angular_distortion) {
          angular_distortion = AngularDistortion(orig, mesh);
        }
      } break;
      case Distortion_Areal: {
        if (!areal_distortion) {
          areal_distortion = ArealDistortion(orig, mesh);
        }
      } break;
      case Distortion_ArealAbs: {
        if (!abs_areal_distortion) {
          if (!areal_distortion) {
            areal_distortion = ArealDistortion(orig, mesh);
          }
          abs_areal_distortion.TakeReference(areal_distortion->NewInstance());
          abs_areal_distortion->SetName("AbsArealDistortion");
          abs_areal_distortion->SetNumberOfComponents(1);
          abs_areal_distortion->SetNumberOfTuples(ncells);
          for (vtkIdType cellId = 0; cellId < ncells; ++cellId) {
            abs_areal_distortion->SetComponent(cellId, 0, abs(areal_distortion->GetComponent(cellId, 0)));
          }
        }
      } break;
    }
  }

  // Add computed data arrays to output mesh and/or report mean/sigma values
  for (auto it = measures.begin(); it != measures.end(); ++it) {
    switch (*it) {
      case Distortion_Angular: {
        if (angular_distortion) {
          if (verbose) PrintMeanAndSigma(angular_distortion);
          mesh->GetPointData()->AddArray(angular_distortion);
          angular_distortion = nullptr;
        }
      } break;
      case Distortion_Areal: {
        if (areal_distortion) {
          if (verbose) PrintMeanAndSigma(areal_distortion);
          mesh->GetPointData()->AddArray(areal_distortion);
          areal_distortion = nullptr;
        }
      } break;
      case Distortion_ArealAbs: {
        if (abs_areal_distortion) {
          if (verbose) PrintMeanAndSigma(abs_areal_distortion);
          mesh->GetPointData()->AddArray(abs_areal_distortion);
          abs_areal_distortion = nullptr;
        }
      } break;
    }
  }

  // Write output dataset
  if (output_name && !WritePolyData(output_name, mesh)) {
    FatalError("Failed to write result to " << output_name);
  }

  return 0;
}
