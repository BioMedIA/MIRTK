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

#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkCellData.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> [<output>] [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Delete point data and/or cell data from input point set." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input point set." << endl;
  cout << "  output   Output point set. (default: input)" << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -all                        Remove all point and cell data arrays." << endl;
  cout << "  -scalars                    Delete SCALARS attribute(s)." << endl;
  cout << "  -vectors                    Delete VECTORS attribute(s)." << endl;
  cout << "  -normals                    Delete NORMALS attribute(s)." << endl;
  cout << "  -tcoords                    Delete TCOORDS attribute(s)." << endl;
  cout << "  -name <name>                Name of point/cell data array to remove." << endl;
  cout << "  -pointdata <name>|<index>   Name of point data array to remove." << endl;
  cout << "  -celldata  <name>|<index>   Name of cell data array to remove." << endl;
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
  REQUIRES_POSARGS(1);

  const char *input_name  = POSARG(1);
  const char *output_name = input_name;
  if (NUM_POSARGS == 2) {
    output_name = POSARG(2);
  } else if (NUM_POSARGS > 2) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  // Handle -help and -version before loading any data
  for (ALL_OPTIONS) HANDLE_HELP_OR_VERSION();

  // Read input point set
  vtkSmartPointer<vtkPointSet> pointset = ReadPointSet(input_name);

  // Delete point/cell data array(s) while processing options
  int index;
  for (ALL_OPTIONS) {
    if (OPTION("-all")) {
      pointset->GetPointData()->Initialize();
      pointset->GetCellData ()->Initialize();
    } else if (OPTION("-name")) {
      const char *name = ARGUMENT;
      pointset->GetPointData()->RemoveArray(name);
      pointset->GetCellData ()->RemoveArray(name);
    } else if (OPTION("-pointdata")) {
      const char *arg = ARGUMENT;
      if (FromString(arg, index)) {
        pointset->GetPointData()->RemoveArray(index);
      } else {
        pointset->GetPointData()->RemoveArray(arg);
      }
    } else if (OPTION("-celldata")) {
      const char *arg = ARGUMENT;
      if (FromString(arg, index)) {
        pointset->GetPointData()->RemoveArray(index);
      } else {
        pointset->GetPointData()->RemoveArray(arg);
      }
    } else if (OPTION("-scalars")) {
      pointset->GetPointData()->SetScalars(NULL);
      pointset->GetCellData ()->SetScalars(NULL);
    } else if (OPTION("-vectors")) {
      pointset->GetPointData()->SetVectors(NULL);
      pointset->GetCellData ()->SetVectors(NULL);
    } else if (OPTION("-normals")) {
      pointset->GetPointData()->SetNormals(NULL);
      pointset->GetCellData ()->SetNormals(NULL);
    } else if (OPTION("-tcoords")) {
      pointset->GetPointData()->SetTCoords(NULL);
      pointset->GetCellData ()->SetTCoords(NULL);
    } else {
      HANDLE_COMMON_OR_UNKNOWN_OPTION();
    }
  }

  // Write output point set
  if (!WritePointSet(output_name, pointset)) {
    FatalError("Failed to write point set to " << output_name);
  }

  return 0;
}
