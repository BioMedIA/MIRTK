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

#include "mirtk/UnorderedSet.h"
#include "mirtk/Algorithm.h"
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
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
class DataSelector
{
public:

  /// Destructor
  virtual ~DataSelector() {};

  /// Run selection query and return IDs of selected data points
  virtual UnorderedSet<int> Evaluate(const Array<double> &values) const = 0;
};

// -----------------------------------------------------------------------------
/// Combine one or more data selection criteria using a logical operator
class LogicalOperator : public DataSelector
{
protected:

  /// Selection criteria
  List<SharedPtr<const DataSelector> > _Criteria;

public:

  /// Number of data criteria
  int NumberOfCriteria() const
  {
    return static_cast<int>(_Criteria.size());
  }

  /// Get n-th criterium
  SharedPtr<const DataSelector> Criterium(int i) const
  {
    auto it = _Criteria.begin();
    for (int pos = 0; pos < i; ++pos) ++it;
    return *it;
  }

  /// Add selection criterium
  void Push(const SharedPtr<const DataSelector> &criterium)
  {
    _Criteria.push_back(criterium);
  }
};

// -----------------------------------------------------------------------------
/// Combine one or more data selection criteria using a logical AND
class LogicalAnd : public LogicalOperator
{
public:

  /// Run selection query and return IDs of selected data points
  virtual UnorderedSet<int> Evaluate(const Array<double> &values) const
  {
    if (_Criteria.empty()) return UnorderedSet<int>();
    auto criterium = _Criteria.begin();
    auto cellIds   = (*criterium)->Evaluate(values);
    while (++criterium != _Criteria.end()) {
      cellIds = Intersection(cellIds, (*criterium)->Evaluate(values));
    }
    return cellIds;
  }
};

// -----------------------------------------------------------------------------
/// Combine one or more data selection criteria using a logical OR
class LogicalOr : public LogicalOperator
{
public:

  /// Run selection query and return IDs of selected data points
  virtual UnorderedSet<int> Evaluate(const Array<double> &values) const
  {
    UnorderedSet<int> cellIds;
    for (auto criterium : _Criteria) {
      auto ids = criterium->Evaluate(values);
      for (auto cellId : ids) {
        cellIds.insert(cellId);
      }
    }
    return cellIds;
  }
};

// -----------------------------------------------------------------------------
class DataCriterium : public DataSelector
{
public:

  /// Run selection query and return IDs of selected data points
  virtual UnorderedSet<int> Evaluate(const Array<double> &values) const
  {
    UnorderedSet<int> cellIds;
    cellIds.reserve(values.size());
    const int n = static_cast<int>(values.size());
    for (int cellId = 0; cellId < n; ++cellId) {
      if (this->Select(values[cellId])) {
        cellIds.insert(cellId);
      }
    }
    return cellIds;
  }

  /// Evaluate criterium for given data value
  virtual bool Select(double) const = 0;
};

// -----------------------------------------------------------------------------
class Equal : public DataCriterium
{
  /// Upper data value threshold
  mirtkPublicAttributeMacro(double, Value);

public:

  /// Constructor
  Equal(double value) : _Value(value) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return fequal(value, _Value);
  }
};

// -----------------------------------------------------------------------------
class NotEqual : public DataCriterium
{
  /// Upper data value threshold
  mirtkPublicAttributeMacro(double, Value);

public:

  /// Constructor
  NotEqual(double value) : _Value(value) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return !fequal(value, _Value);
  }
};

// -----------------------------------------------------------------------------
class LessThan : public DataCriterium
{
  /// Upper data value threshold
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  LessThan(double threshold) : _Threshold(threshold) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return value < _Threshold;
  }
};

// -----------------------------------------------------------------------------
class LessOrEqual : public DataCriterium
{
  /// Upper data value threshold
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  LessOrEqual(double threshold) : _Threshold(threshold) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return value <= _Threshold;
  }
};

// -----------------------------------------------------------------------------
class GreaterThan : public DataCriterium
{
  /// Lower data value threshold
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  GreaterThan(double threshold) : _Threshold(threshold) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return value > _Threshold;
  }
};

// -----------------------------------------------------------------------------
class GreaterOrEqual : public DataCriterium
{
  /// Lower data value threshold
  mirtkPublicAttributeMacro(double, Threshold);

public:

  /// Constructor
  GreaterOrEqual(double threshold) : _Threshold(threshold) {}

  /// Evaluate criterium for given data value
  virtual bool Select(double value) const
  {
    return value >= _Threshold;
  }
};

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
  bool output_surface = false;
  bool recalc_normals = false;

  double                     value;
  Array<double>              values;
  SharedPtr<LogicalOperator> selector(new LogicalAnd());

  for (ALL_OPTIONS) {
    if (OPTION("-select")) {
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
    else HANDLE_BOOLEAN_OPTION("surface", output_surface);
    else HANDLE_BOOLEAN_OPTION("normals", recalc_normals);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
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
