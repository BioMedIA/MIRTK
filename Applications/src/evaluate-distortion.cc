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

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkCell.h"
#include "vtkTriangle.h"
#include "vtkFloatArray.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
/// Print help screen
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <original> <deformed> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Computes the distortion of a surface mesh under deformation and stores" << endl;
  cout << "  the results in the output dataset. If no output dataset is given, the" << endl;
  cout << "  mean and standard deviation are reported only." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  original   File name of original dataset." << endl;
  cout << "  deformed   File name of deformed dataset." << endl;
  cout << "  output     File name of output dataset." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -areal       Compute areal distortion of cells: log_2(area(deformed)/area(original)). (default)" << endl;
  cout << "  -abs-areal   Compute absolute value of :option:`-areal` distortion." << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Positional arguments
  REQUIRES_POSARGS(2);

  const char *original_name = POSARG(1);
  const char *deformed_name = POSARG(2);
  const char *output_name   = NULL;

  if (NUM_POSARGS == 3) output_name = POSARG(3);
  else if (NUM_POSARGS > 3) {
    cerr << "Too many positional arguments given (see -help for usage)" << endl;
    exit(1);
  }

  // Optional arguments
  bool areal     = false;
  bool areal_abs = false;

  for (ALL_OPTIONS) {
    if      (OPTION("-areal"))     areal     = true;
    else if (OPTION("-abs-areal")) areal_abs = true;
    else HANDLE_STANDARD_OR_UNKNOWN_OPTION();
  }

  if (!areal && !areal_abs) areal = true;

  // Read input files
  if (verbose > 1) cout << "Reading original surface from " << original_name << endl;
  vtkSmartPointer<vtkPolyData> original = ReadPolyData(original_name);
  if (verbose > 1) cout << "Reading deformed surface from " << deformed_name << endl;
  vtkSmartPointer<vtkPolyData> deformed = ReadPolyData(deformed_name);

  // Check input
  vtkIdType npoints = original->GetNumberOfPoints();
  if (npoints == 0) {
    cerr << "Original surface has no points!" << endl;
    exit(1);
  }
  if (deformed->GetNumberOfPoints() != original->GetNumberOfPoints()) {
    cerr << "Surface meshes have differing number of points!" << endl;
    exit(1);
  }
  vtkIdType ncells = original->GetNumberOfCells();
  if (ncells == 0) {
    cerr << "Original surface has no cells!" << endl;
    exit(1);
  }
  if (deformed->GetNumberOfCells() != original->GetNumberOfCells()) {
    cerr << "Surface meshes have differing number of cells!" << endl;
    exit(1);
  }

  // Compute areal distortion
  if (areal || areal_abs) {
    vtkSmartPointer<vtkFloatArray> areal_distortion;
    vtkSmartPointer<vtkFloatArray> abs_areal_distortion;

    if (areal && output_name) {
      areal_distortion = vtkSmartPointer<vtkFloatArray>::New();
      areal_distortion->SetNumberOfComponents(1);
      areal_distortion->SetNumberOfTuples(ncells);
      areal_distortion->SetName("areal_distortion");
      deformed->GetCellData()->AddArray(areal_distortion);
    }
    if (areal_abs && output_name) {
      abs_areal_distortion = vtkSmartPointer<vtkFloatArray>::New();
      abs_areal_distortion->SetNumberOfComponents(1);
      abs_areal_distortion->SetNumberOfTuples(ncells);
      abs_areal_distortion->SetName("abs_areal_distortion");
      deformed->GetCellData()->AddArray(abs_areal_distortion);
    }
  
    const double log2 = log(2);
    double distortion;
    double areal_distortion_sum     = .0, areal_distortion_sum2     = .0;
    double abs_areal_distortion_sum = .0, abs_areal_distortion_sum2 = .0;

    for (vtkIdType i = 0; i < ncells; ++i) {
      vtkCell * const original_cell = original->GetCell(i);
      vtkCell * const deformed_cell = deformed->GetCell(i);
      if (deformed_cell->GetCellType() != original_cell->GetCellType()) {
        cerr << "Cell " << (i+1) << " has differing type in the two meshes!" << endl;
        exit(1);
      }
      if (deformed_cell->GetNumberOfPoints() != original_cell->GetNumberOfPoints()) {
        cerr << "Cell " << (i+1) << " has different number of points in the two meshes!" << endl;
        exit(1);
      }
      switch (original_cell->GetCellType()) {
        case VTK_TRIANGLE: {
            vtkTriangle * const original_triangle = dynamic_cast<vtkTriangle *>(original_cell);
            vtkTriangle * const deformed_triangle = dynamic_cast<vtkTriangle *>(deformed_cell);
            distortion = log(deformed_triangle->ComputeArea() / original_triangle->ComputeArea()) / log2;
            if (areal) {
              areal_distortion_sum  += distortion;
              areal_distortion_sum2 += distortion * distortion;
              if (areal_distortion) areal_distortion    ->SetTuple1(i, distortion);
            }
            if (areal_abs) {
              distortion = abs(distortion);
              abs_areal_distortion_sum  += distortion;
              abs_areal_distortion_sum2 += distortion * distortion;
              if (abs_areal_distortion) abs_areal_distortion->SetTuple1(i, distortion);
            }
          } break;
        default:
          cerr << "Unsupported cell type: " << original_cell->GetCellType() << endl;
          exit(1);
      }
    }

    if (!output_name || verbose) {
      if (areal) {
        double mean = areal_distortion_sum / ncells;
        cout << "Areal distortion: mean = " << mean << ", sigma = "
             << sqrt(areal_distortion_sum2 / ncells - mean * mean) << endl;
      }
      if (areal_abs) {
        double mean = abs_areal_distortion_sum / ncells;
        cout << "Absolute areal distortion: mean = " << mean << ", sigma = "
             << sqrt(abs_areal_distortion_sum2 / ncells - mean * mean) << endl;
      }
    }
  }

  // Write output dataset
  if (output_name) WritePolyData(output_name, deformed);

  return 0;
}
