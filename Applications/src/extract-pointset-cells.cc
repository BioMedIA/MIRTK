/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2016 Imperial College London
 * Copyright 2016 Andreas Schuh
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

#include "mirtk/DataSelection.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"

#include "vtkNew.h"
#include "vtkSmartPointer.h"
#include "vtkExtractCells.h"
#include "vtkPointSet.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkPolyDataNormals.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"

using namespace mirtk;
using namespace mirtk::data;
using namespace mirtk::data::select;


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
  cout << "  Extract point set/surface elements which meet the specified criteria.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input point set.\n";
  cout << "  output   Output point set.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -select <name> [<j>]    Use specified component of cell data array for selection criteria.\n";
  cout << "                          By default, when <j> is not specified, the first component is used. (default: SCALARS)\n";
  cout << "  -and                    Combine current selection criterium to the left of this option\n";
  cout << "                          with the following criteria using a logical AND operator. (default)\n";
  cout << "  -or                     Combine current selection criterium to the left of this option\n";
  cout << "                          with the following criteria using a logical OR operator.\n";
  cout << "  -eq <value>             Select cells whose data value is equal the specified value.\n";
  cout << "  -ne <value>             Select cells whose data value is not equal the specified value.\n";
  cout << "  -lt <value>             Select cells whose data value is less than the specified value.\n";
  cout << "  -le <value>             Select cells whose data value is less than or equal the specified value.\n";
  cout << "  -gt <value>             Select cells whose data value is greater than the specified value.\n";
  cout << "  -ge <value>             Select cells whose data value is greater than or equal the specified value.\n";
  cout << "  -[no]surface [on|off]   Output surface of resulting point set. (default: off)\n";
  cout << "  -[no]normals [on|off]   Recalculate normals of output surface. (default: off)\n";
  cout << "  -origids <name>         Add original cell IDs as output cell data array of given name. (default: off)\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  EXPECTS_POSARGS(2);

  const char * const input_name  = POSARG(1);
  const char * const output_name = POSARG(2);

  vtkSmartPointer<vtkPointSet> input = ReadPointSet(input_name);
  vtkCellData * const cd = input->GetCellData();

  const char *orig_cell_ids = nullptr;
  int  output_surface = -1;
  bool recalc_normals = false;

  double               value;
  Array<double>        values;
  SharedPtr<LogicalOp> selector(new LogicalAnd());

  for (ALL_OPTIONS) {
    if (OPTION("-a") || OPTION("-array") || OPTION("-select") || OPTION("-where")) {
      const char * const name = ARGUMENT;
      vtkDataArray *scalars = cd->GetArray(name);
      if (scalars == nullptr) {
        scalars = GetArrayByCaseInsensitiveName(cd, name);
        if (scalars == nullptr) {
          FatalError("Input point set has no cell data array named: " << name);
        }
      }
      int j = 0;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(j);
      if (j < 0 || j >= scalars->GetNumberOfComponents()) {
        FatalError("Cell data array " << name << " has only " << scalars->GetNumberOfComponents() << " component(s)");
      }
      values.resize(scalars->GetNumberOfTuples());
      for (vtkIdType cellId = 0; cellId < scalars->GetNumberOfTuples(); ++cellId) {
        values[cellId] = scalars->GetComponent(cellId, j);
      }
    }
    else if (OPTION("-and")) {
      if (dynamic_cast<LogicalAnd *>(selector.get()) == nullptr) {
        SharedPtr<LogicalAnd> op(new LogicalAnd());
        if (selector->NumberOfCriteria() > 1) {
          op->Push(selector);
        } else if (selector->NumberOfCriteria() == 1) {
          op->Push(selector->Criterium(0));
        }
        selector = op;
      }
    }
    else if (OPTION("-or")) {
      if (dynamic_cast<LogicalOr *>(selector.get()) == nullptr) {
        SharedPtr<LogicalOr> op(new LogicalOr());
        if (selector->NumberOfCriteria() > 1) {
          op->Push(selector);
        } else if (selector->NumberOfCriteria() == 1) {
          op->Push(selector->Criterium(0));
        }
        selector = op;
      }
    }
    else if (OPTION("-eq")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<Equal>(value));
    }
    else if (OPTION("-ne")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<NotEqual>(value));
    }
    else if (OPTION("-lt")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<LessThan>(value));
    }
    else if (OPTION("-le")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<LessOrEqual>(value));
    }
    else if (OPTION("-gt")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<GreaterThan>(value));
    }
    else if (OPTION("-ge")) {
      PARSE_ARGUMENT(value);
      selector->Push(NewShared<GreaterOrEqual>(value));
    }
    else if (OPTION("-origids")) {
      orig_cell_ids = ARGUMENT;
    }
    else if (OPTION("-surface")) {
      bool bval = true;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(bval);
      output_surface = (bval ? 1 : 0);
    }
    else if (OPTION("-nosurface")) {
      output_surface = 0;
    }
    else HANDLE_BOOLEAN_OPTION("normals", recalc_normals);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }
  if (output_surface == -1) {
    output_surface = (vtkPolyData::SafeDownCast(input) == nullptr ? 0 : 1);
  }

  if (values.empty()) {
    vtkDataArray * const scalars = cd->GetScalars();
    if (scalars) {
      values.resize(scalars->GetNumberOfTuples());
      for (vtkIdType cellId = 0; cellId < scalars->GetNumberOfTuples(); ++cellId) {
        values[cellId] = scalars->GetComponent(cellId, 0);
      }
    } else {
      FatalError("Input point set has no SCALARS cell data array! Use -select option to choose one.");
    }
  }

  auto selection = selector->Evaluate(values);
  vtkNew<vtkIdList> cellIds;
  cellIds->Allocate(static_cast<vtkIdType>(selection.size()));
  for (auto cellId : selection) {
    cellIds->InsertNextId(static_cast<vtkIdType>(cellId));
  }

  vtkNew<vtkExtractCells> extractor;
  SetVTKInput(extractor, input);
  extractor->SetCellList(cellIds.GetPointer());
  extractor->Update();

  vtkSmartPointer<vtkPointSet> output = extractor->GetOutput();
  if (orig_cell_ids && orig_cell_ids[0] != '\0') {
    output->GetCellData()->GetArray("vtkOriginalCellIds")->SetName(orig_cell_ids);
  } else {
    output->GetCellData()->RemoveArray("vtkOriginalCellIds");
  }

  if (output_surface) {
    output = DataSetSurface(output);
    if (recalc_normals) {
      vtkNew<vtkPolyDataNormals> filter;
      SetVTKInput(filter, output);
      filter->ComputeCellNormalsOn();
      filter->ComputePointNormalsOn();
      filter->SplittingOff();
      filter->NonManifoldTraversalOff();
      filter->ConsistencyOn();
      filter->AutoOrientNormalsOff();
      filter->Update();
      output = filter->GetOutput();
    }
  }

  if (!WritePointSet(output_name, output)) {
    FatalError("Failed to write point set to " << output_name);
  }

  return 0;
}
