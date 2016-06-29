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
#include "mirtk/PointSetUtils.h"

#include "vtkSmartPointer.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkDataSetAttributes.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <source> <target> [<output>] [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Copies point and/or cell data from a source point set to a target point set" << endl;
  cout << "  and writes the resulting amended point set to the specified output file." << endl;
  cout << "  When no separate output file is specified, the target point set is overwritten." << endl;
  cout << endl;
  cout << "  If the point sets have differing number of points or cells, respectively," << endl;
  cout << "  zero entries are either added to the target arrays or only the first n tuples" << endl;
  cout << "  of the source arrays are copied. In case of the -points option, the remaining" << endl;
  cout << "  target points are kept unchanged." << endl;
  cout << endl;
  cout << "  If the <attr> argument is given to any of the option below, the copied data" << endl;
  cout << "  array is set as the new \"scalars\", \"vectors\", \"tcoords\", or \"normals\"" << endl;
  cout << "  attribute, respectively, of the output point set." << endl;
  cout << endl;
  cout << "  When no options are given, the :option:`-scalars` are copied over." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  source   Point set from which to copy the point/cell data." << endl;
  cout << "  target   Point set to which to add the point/cell data." << endl;
  cout << "  output   Output point set with copied point/cell data." << endl;
  cout << "           If not specified, the target file is overwritten." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -scalars" << endl;
  cout << "      Copy \"scalars\" data set attribute. (default)" << endl;
  cout << endl;
  cout << "  -vectors" << endl;
  cout << "      Copy \"vectors\" data set attribute." << endl;
  cout << endl;
  cout << "  -vectors-as-points" << endl;
  cout << "      Set points of target to \"vectors\" data set attribute of source." << endl;
  cout << endl;
  cout << "  -normals" << endl;
  cout << "      Copy \"normals\" data set attribute." << endl;
  cout << endl;
  cout << "  -tcoords" << endl;
  cout << "      Copy \"tcoords\" data set attribute." << endl;
  cout << endl;
  cout << "  -tcoords-as-points" << endl;
  cout << "      Set points of target to \"tcoords\" data set attribute of source." << endl;
  cout << endl;
  cout << "  -points [<target_name> [<attr>]]" << endl;
  cout << "      Add points of source as point data of target." << endl;
  cout << "      When no <target_name> (and <attr>) specified, the points of the" << endl;
  cout << "      target point set are replaced by those of the source point set." << endl;
  cout << "      If the target point set has more points than the source point set," << endl;
  cout << "      the remaining points of the target point set are unchanged." << endl;
  cout << endl;
  cout << "  -pointdata <name|index> [<target_name> [<attr>]]" << endl;
  cout << "      Name or index of source point data to copy." << endl;
  cout << endl;
  cout << "  -pointdata-as-points <name|index>" << endl;
  cout << "      Replaces the points of the target point set by the first three" << endl;
  cout << "      components of the specified source point data array." << endl;
  cout << endl;
  cout << "  -pointmask <name>" << endl;
  cout << "      Add point data array which indicates copied point data." << endl;
  cout << endl;
  cout << "  -celldata  <name|index> [<target_name> [<attr>]]" << endl;
  cout << "      Name or index of source cell data to copy." << endl;
  cout << endl;
  cout << "  -cellmask <name>" << endl;
  cout << "      Add cell data array which indicates copied cell data." << endl;
  cout << endl;
  cout << "  -case-[in]sensitive" << endl;
  cout << "      Lookup source arrays by case [in]sensitive name. (default: insensitive)" << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
struct ArrayInfo
{
  int                                  _SourceIndex;
  const char                          *_SourceName;
  vtkDataSetAttributes::AttributeTypes _SourceAttribute;
  int                                  _TargetIndex;
  const char                          *_TargetName;
  vtkDataSetAttributes::AttributeTypes _TargetAttribute;

  ArrayInfo()
  :
    _SourceIndex(-1),
    _SourceName(NULL),
    _SourceAttribute(vtkDataSetAttributes::NUM_ATTRIBUTES),
    _TargetIndex(-1),
    _TargetName(NULL),
    _TargetAttribute(vtkDataSetAttributes::NUM_ATTRIBUTES)
  {}
};

// -----------------------------------------------------------------------------
/// Copy tuples from the input data array
vtkSmartPointer<vtkDataArray> Copy(vtkDataArray *array, vtkIdType n)
{
  vtkSmartPointer<vtkDataArray> copy;
  copy.TakeReference(array->NewInstance());
  copy->SetName(array->GetName());
  copy->SetNumberOfComponents(array->GetNumberOfComponents());
  copy->SetNumberOfTuples(n);
  vtkIdType ptId = 0;
  while (ptId < n && ptId < array->GetNumberOfTuples()) {
    copy->SetTuple(ptId, array->GetTuple(ptId));
    ++ptId;
  }
  while (ptId < n) {
    for (int j = 0; j < copy->GetNumberOfComponents(); ++j) {
      copy->SetComponent(ptId, j, .0);
    }
    ++ptId;
  }
  return copy;
}

// -----------------------------------------------------------------------------
/// Add data array to data set attributes
void AddArray(vtkDataSetAttributes *data, vtkDataArray *array, vtkDataSetAttributes::AttributeTypes type)
{
  switch (type) {
    case vtkDataSetAttributes::SCALARS:     data->SetScalars(array); break;
    case vtkDataSetAttributes::VECTORS:     data->SetVectors(array); break;
    case vtkDataSetAttributes::NORMALS:     data->SetNormals(array); break;
    case vtkDataSetAttributes::TCOORDS:     data->SetTCoords(array); break;
    case vtkDataSetAttributes::TENSORS:     data->SetTensors(array); break;
    case vtkDataSetAttributes::GLOBALIDS:   data->SetGlobalIds(array); break;
    case vtkDataSetAttributes::PEDIGREEIDS: data->SetPedigreeIds(array); break;
    default: data->AddArray(array); break;
  }
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  vtkSmartPointer<vtkDataArray> array;
  vtkSmartPointer<vtkDataArray> copy;

  // Positional arguments
  REQUIRES_POSARGS(2);

  if (NUM_POSARGS > 3) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  const char *source_name = POSARG(1);
  const char *target_name = POSARG(2);
  const char *output_name = (NUM_POSARGS == 3 ? POSARG(3) : target_name);

  // Optional arguments
  Array<ArrayInfo> pd, cd;
  const char *point_mask_name = NULL;
  const char *cell_mask_name  = NULL;
  bool case_sensitive = false;

  for (ALL_OPTIONS) {
    if (OPTION("-points")) {
      ArrayInfo info;
      info._SourceIndex = -2;
      if (HAS_ARGUMENT) info._TargetName = ARGUMENT;
      else              info._TargetName = NULL;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(info._TargetAttribute);
      pd.push_back(info);
    }
    else if (OPTION("-scalars") || OPTION("-vectors") || OPTION("-normals") || OPTION("-tcoords")) {
      ArrayInfo info;
      FromString(OPTNAME + 1, info._SourceAttribute);
      pd.push_back(info);
      cd.push_back(info);
    }
    else if (OPTION("-vectors-as-points") || OPTION("-tcoords-as-points")) {
      ArrayInfo info;
      FromString(OPTNAME + 1, info._SourceAttribute);
      info._TargetIndex = -2;
      pd.push_back(info);
    }
    else if (OPTION("-pointdata") || OPTION("-celldata")) {
      ArrayInfo info;
      info._SourceName = ARGUMENT;
      if (FromString(info._SourceName, info._SourceIndex)) {
        info._SourceName = NULL;
      } else {
        info._SourceIndex = -1;
      }
      if (HAS_ARGUMENT) {
        info._TargetName = ARGUMENT;
      }
      if (HAS_ARGUMENT) PARSE_ARGUMENT(info._TargetAttribute);
      if (strcmp(OPTNAME, "-pointdata") == 0) pd.push_back(info);
      else                                    cd.push_back(info);
    }
    else if (OPTION("-pointdata-as-points")) {
      ArrayInfo info;
      info._SourceName = ARGUMENT;
      if (FromString(info._SourceName, info._SourceIndex)) {
        info._SourceName = NULL;
      } else {
        info._SourceIndex = -1;
      }
      info._TargetIndex = -2;
      pd.push_back(info);
    }
    else if (OPTION("-pointmask")) point_mask_name = ARGUMENT;
    else if (OPTION("-cellmask"))  cell_mask_name  = ARGUMENT;
    else if (OPTION("-case-sensitive"))   case_sensitive = true;
    else if (OPTION("-case-insensitive")) case_sensitive = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input data sets
  vtkSmartPointer<vtkPointSet> source = ReadPointSet(source_name);
  vtkSmartPointer<vtkPointSet> target = ReadPointSet(target_name);

  // Add point data
  const vtkIdType npoints = target->GetNumberOfPoints();
  vtkPointData *sourcePD = source->GetPointData();
  vtkPointData *targetPD = target->GetPointData();
  for (size_t i = 0; i < pd.size(); ++i) {
    if (pd[i]._SourceIndex == -2) {
      vtkIdType end = min(npoints, source->GetNumberOfPoints());
      if (pd[i]._TargetName) {
        copy = vtkSmartPointer<vtkFloatArray>::New();
        copy->SetNumberOfComponents(3);
        copy->SetNumberOfTuples(npoints);
        for (vtkIdType ptId = 0; ptId < end; ++ptId) {
          copy->SetTuple(ptId, source->GetPoint(ptId));
        }
        const double zero[3] = {.0, .0, .0};
        for (vtkIdType ptId = end; ptId < npoints; ++ptId) {
          copy->SetTuple(ptId, zero);
        }
      } else {
        vtkPoints *points = target->GetPoints();
        for (vtkIdType ptId = 0; ptId < end; ++ptId) {
          points->SetPoint(ptId, source->GetPoint(ptId));
        }
        copy = NULL;
      }
    } else {
      if (pd[i]._SourceAttribute < vtkDataSetAttributes::NUM_ATTRIBUTES) {
        array = sourcePD->GetAttribute(pd[i]._SourceAttribute);
        if (!array) continue;
      } else if (pd[i]._SourceName) {
        if (case_sensitive) {
          array = sourcePD->GetArray(pd[i]._SourceName);
        } else {
          array = GetArrayByCaseInsensitiveName(sourcePD, pd[i]._SourceName);
        }
        if (array == NULL) {
          FatalError("Source has no point data array named " << pd[i]._SourceName);
        }
      } else {
        if (pd[i]._SourceIndex < 0 || pd[i]._SourceIndex > sourcePD->GetNumberOfArrays()) {
          FatalError("Source has no point data array with index " << pd[i]._SourceIndex);
        }
        array = sourcePD->GetArray(pd[i]._SourceIndex);
      }
      if (pd[i]._TargetIndex == -2) {
        double p[3] = {.0};
        vtkPoints *points = target->GetPoints();
        vtkIdType end = min(npoints, source->GetNumberOfPoints());
        int       dim = min(3,       array->GetNumberOfComponents());
        for (vtkIdType ptId = 0; ptId < end; ++ptId) {
          for (int i = 0; i < dim; ++i) p[i] = array->GetComponent(ptId, i);
          points->SetPoint(ptId, p);
        }
        copy = NULL;
      } else {
        copy = Copy(array, npoints);
      }
    }
    if (copy) {
      if (pd[i]._TargetName) copy->SetName(pd[i]._TargetName);
      AddArray(targetPD, copy, pd[i]._TargetAttribute);
    }
  }
  if (point_mask_name) {
    vtkSmartPointer<vtkDataArray> mask = NewVTKDataArray(VTK_UNSIGNED_CHAR);
    mask->SetName(point_mask_name);
    mask->SetNumberOfComponents(1);
    mask->SetNumberOfTuples(npoints);
    mask->FillComponent(0, .0);
    for (vtkIdType ptId = 0; ptId < npoints && ptId < source->GetNumberOfPoints(); ++ptId) {
      mask->SetComponent(ptId, 0, 1.0);
    }
    targetPD->AddArray(mask);
  }

  // Add cell data
  const vtkIdType ncells = target->GetNumberOfCells();
  vtkCellData *sourceCD = source->GetCellData();
  vtkCellData *targetCD = target->GetCellData();
  for (size_t i = 0; i < cd.size(); ++i) {
    if (cd[i]._SourceAttribute < vtkDataSetAttributes::NUM_ATTRIBUTES) {
      array = sourceCD->GetAttribute(cd[i]._SourceAttribute);
      if (!array) continue;
    } else if (cd[i]._SourceName) {
      if (case_sensitive) {
        array = sourceCD->GetArray(cd[i]._SourceName);
      } else {
        array = GetArrayByCaseInsensitiveName(sourceCD, cd[i]._SourceName);
      }
      if (array == NULL) {
        FatalError("Source has no cell data array named " << cd[i]._SourceName);
      }
    } else {
      if (cd[i]._SourceIndex < 0 || cd[i]._SourceIndex > sourceCD->GetNumberOfArrays()) {
        FatalError("Source has no cell data array with index " << cd[i]._SourceIndex);
      }
      array = sourceCD->GetArray(cd[i]._SourceIndex);
    }
    copy = Copy(array, ncells);
    if (cd[i]._TargetName) copy->SetName(cd[i]._TargetName);
    AddArray(targetCD, copy, cd[i]._TargetAttribute);
  }
  if (cell_mask_name) {
    vtkSmartPointer<vtkDataArray> mask = NewVTKDataArray(VTK_UNSIGNED_CHAR);
    mask->SetName(cell_mask_name);
    mask->SetNumberOfComponents(1);
    mask->SetNumberOfTuples(ncells);
    mask->FillComponent(0, .0);
    for (vtkIdType cellId = 0; cellId < ncells && cellId < source->GetNumberOfCells(); ++cellId) {
      mask->SetComponent(cellId, 0, 1.0);
    }
    targetCD->AddArray(mask);
  }

  // Write resulting data set
  WritePointSet(output_name, target);
}
