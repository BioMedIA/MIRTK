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
#include "vtkPointDataToCellData.h"
#include "vtkCellDataToPointData.h"

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
  cout << "  This command can also convert from point data to cell data and vice versa." << endl;
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
  #if MIRTK_IO_WITH_GIFTI
  cout << "           Can also be a GIFTI file with only point data arrays." << endl;
  cout << "           In this case, the coordinates and topology of the target" << endl;
  cout << "           surface are used to define the source surface." << endl;
  #endif
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
  cout << "  -pointdata-as-celldata <name|index>" << endl;
  cout << "      Converts the specified point data array of the source point set to cell" << endl;
  cout << "      data and adds a corresponding cell data array to the target point set." << endl;
  cout << endl;
  cout << "  -pointmask <name>" << endl;
  cout << "      Add point data array which indicates copied point data." << endl;
  cout << endl;
  cout << "  -celldata  <name|index> [<target_name> [<attr>]]" << endl;
  cout << "      Name or index of source cell data to copy." << endl;
  cout << endl;
  cout << "  -celldata-as-pointdata <name|index>" << endl;
  cout << "      Converts the specified cell data array of the source point set to point" << endl;
  cout << "      data and adds a corresponding point data array to the target point set." << endl;
  cout << endl;
  cout << "  -cellmask <name>" << endl;
  cout << "      Add cell data array which indicates copied cell data." << endl;
  cout << endl;
  cout << "  -case-[in]sensitive" << endl;
  cout << "      Lookup source arrays by case [in]sensitive name. (default: insensitive)" << endl;
  cout << endl;
  cout << "  -majority" << endl;
  cout << "      Use majority vote for resampling categorical data. (default)" << endl;
  cout << endl;
  cout << "  -unanimous" << endl;
  cout << "      Use unanimous vote for resampling categorical data." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
using AttributeType                = vtkDataSetAttributes::AttributeTypes;
const AttributeType NUM_ATTRIBUTES = vtkDataSetAttributes::NUM_ATTRIBUTES;

// -----------------------------------------------------------------------------
struct ArrayInfo
{
  int           _SourceIndex;
  const char   *_SourceName;
  AttributeType _SourceAttribute;
  int           _TargetIndex;
  const char   *_TargetName;
  AttributeType _TargetAttribute;
  int           _TargetDataType;
  bool          _PointCellConversion;

  ArrayInfo()
  :
    _SourceIndex(-1),
    _SourceName(nullptr),
    _SourceAttribute(NUM_ATTRIBUTES),
    _TargetIndex(-1),
    _TargetName(nullptr),
    _TargetAttribute(NUM_ATTRIBUTES),
    _TargetDataType(VTK_VOID),
    _PointCellConversion(false)
  {}
};

// -----------------------------------------------------------------------------
/// Get source point data array
vtkDataArray *GetSourceArray(const char *type, vtkDataSetAttributes *attr, const ArrayInfo &info, bool case_sensitive)
{
  vtkDataArray *array;
  if (info._SourceAttribute < NUM_ATTRIBUTES) {
    array = attr->GetAttribute(info._SourceAttribute);
  } else if (info._SourceName) {
    if (case_sensitive) {
      array = attr->GetArray(info._SourceName);
    } else {
      array = GetArrayByCaseInsensitiveName(attr, info._SourceName);
    }
    if (array == nullptr) {
      FatalError("Source has no " << type << " data array named " << info._SourceName);
    }
  } else {
    if (info._SourceIndex < 0 || info._SourceIndex >= attr->GetNumberOfArrays()) {
      FatalError("Source has no " << type << " data array with index " << info._SourceIndex);
    }
    array = attr->GetArray(info._SourceIndex);
  }
  return array;
}

// -----------------------------------------------------------------------------
/// Parse VTK data type string
int ParseVtkDataType(const char *arg)
{
  const string lstr = ToLower(Trim(arg));
  if (lstr == "char")   return VTK_CHAR;
  if (lstr == "uchar")  return VTK_UNSIGNED_CHAR;
  if (lstr == "binary") return VTK_UNSIGNED_CHAR;
  if (lstr == "short")  return VTK_SHORT;
  if (lstr == "grey")   return VTK_SHORT;
  if (lstr == "ushort") return VTK_UNSIGNED_SHORT;
  if (lstr == "int")    return VTK_INT;
  if (lstr == "uint")   return VTK_UNSIGNED_INT;
  if (lstr == "long")   return VTK_LONG_LONG;
  if (lstr == "ulong")  return VTK_UNSIGNED_LONG_LONG;
  if (lstr == "float")  return VTK_FLOAT;
  if (lstr == "double") return VTK_DOUBLE;
  if (lstr == "real") {
    #if MIRTK_USE_FLOAT_BY_DEFAULT
      return VTK_FLOAT;
    #else
      return VTK_DOUBLE;
    #endif
  }
  FatalError("Unknown data type: " << arg);
}

/// Enumeration of categorical data interpolation methods
enum LabelInterpolationMode
{
  MajorityVote,
  UnanimousVote
};

// -----------------------------------------------------------------------------
/// Convert point data labels to cell data
vtkSmartPointer<vtkDataArray>
ConvertPointToCellLabels(vtkPointSet *pset, vtkDataArray *labels,
                         LabelInterpolationMode mode = MajorityVote)
{
  vtkSmartPointer<vtkDataArray> output;
  output.TakeReference(labels->NewInstance());
  output->SetName(labels->GetName());
  for (int j = 0; j < labels->GetNumberOfComponents(); ++j) {
    output->SetComponentName(j, labels->GetComponentName(j));
  }
  output->SetNumberOfComponents(labels->GetNumberOfComponents());
  output->SetNumberOfTuples(pset->GetNumberOfCells());

  switch (mode) {
    case MajorityVote: {
      int    count;
      double label;
      OrderedMap<double, int> hist;
      OrderedMap<double, int>::iterator bin;
      vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
      for (vtkIdType cellId = 0; cellId < pset->GetNumberOfCells(); ++cellId) {
        pset->GetCellPoints(cellId, ptIds);
        for (int j = 0; j < labels->GetNumberOfComponents(); ++j) {
          hist.clear();
          for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
            label = labels->GetComponent(ptIds->GetId(i), j);
            bin = hist.find(label);
            if (bin == hist.end()) hist[label] = 1;
            else bin->second += 1;
          }
          label = 0., count = 0;
          for (bin = hist.begin(); bin != hist.end(); ++bin) {
            if (bin->second > count) {
              label = bin->first;
              count = bin->second;
            }
          }
          output->SetComponent(cellId, j, label);
        }
      }
    } break;
    case UnanimousVote: {
      const double invalid = NaN;
      double       label;
      vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
      for (vtkIdType cellId = 0; cellId < pset->GetNumberOfCells(); ++cellId) {
        pset->GetCellPoints(cellId, ptIds);
        for (int j = 0; j < labels->GetNumberOfComponents(); ++j) {
          if (ptIds->GetNumberOfIds() == 0) {
            label = invalid;
          } else {
            label = labels->GetComponent(ptIds->GetId(0), j);
            for (vtkIdType i = 1; i < ptIds->GetNumberOfIds(); ++i) {
              if (labels->GetComponent(ptIds->GetId(i), j) != label) {
                label = invalid;
                break;
              }
            }
          }
          output->SetComponent(cellId, j, label);
        }
      }
    } break;
  }

  return output;
}

// -----------------------------------------------------------------------------
/// Convert cell data labels to point data
vtkSmartPointer<vtkDataArray>
ConvertCellToPointLabels(vtkPointSet *pset, vtkDataArray *labels,
                         LabelInterpolationMode mode = MajorityVote)
{
  vtkSmartPointer<vtkDataArray> output;
  output.TakeReference(labels->NewInstance());
  output->SetName(labels->GetName());
  for (int j = 0; j < labels->GetNumberOfComponents(); ++j) {
    output->SetComponentName(j, labels->GetComponentName(j));
  }
  output->SetNumberOfComponents(labels->GetNumberOfComponents());
  output->SetNumberOfTuples(pset->GetNumberOfPoints());

  switch (mode) {
    case MajorityVote: {
      int    count;
      double label;
      OrderedMap<double, int> hist;
      OrderedMap<double, int>::iterator bin;
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      for (vtkIdType ptId = 0; ptId < pset->GetNumberOfPoints(); ++ptId) {
        pset->GetPointCells(ptId, cellIds);
        for (int j = 0; j < labels->GetNumberOfComponents(); ++j) {
          hist.clear();
          for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
            label = labels->GetComponent(cellIds->GetId(i), j);
            bin = hist.find(label);
            if (bin == hist.end()) hist[label] = 1;
            else bin->second += 1;
          }
          label = 0., count = 0;
          for (bin = hist.begin(); bin != hist.end(); ++bin) {
            if (bin->second > count) {
              label = bin->first;
              count = bin->second;
            }
          }
          output->SetComponent(ptId, j, label);
        }
      }
    } break;
    case UnanimousVote: {
      double       label;
      const double invalid = NaN;
      vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
      for (vtkIdType ptId = 0; ptId < pset->GetNumberOfPoints(); ++ptId) {
        pset->GetPointCells(ptId, cellIds);
        for (int j = 0; j < labels->GetNumberOfComponents(); ++j) {
          if (cellIds->GetNumberOfIds() == 0) {
            label = invalid;
          } else {
            label = labels->GetComponent(cellIds->GetId(0), j);
            for (vtkIdType i = 1; i < cellIds->GetNumberOfIds(); ++i) {
              if (labels->GetComponent(cellIds->GetId(i), j) != label) {
                label = invalid;
                break;
              }
            }
          }
          output->SetComponent(ptId, j, label);
        }
      }
    } break;
  }

  return output;
}

// -----------------------------------------------------------------------------
/// Copy tuples from the input data array
vtkSmartPointer<vtkDataArray> Copy(vtkDataArray *array, vtkIdType n, int dtype = VTK_VOID)
{
  vtkSmartPointer<vtkDataArray> copy;
  if (dtype == VTK_VOID) {
    copy.TakeReference(array->NewInstance());
  } else {
    copy = NewVtkDataArray(dtype);
  }
  copy->SetName(array->GetName());
  copy->SetNumberOfComponents(array->GetNumberOfComponents());
  copy->SetNumberOfTuples(n);
  for (int j = 0; j < array->GetNumberOfComponents(); ++j) {
    copy->SetComponentName(j, array->GetComponentName(j));
  }
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
void AddArray(vtkDataSetAttributes *data, vtkDataArray *array, AttributeType type)
{
  switch (type) {
    case AttributeType::SCALARS:     data->SetScalars(array); break;
    case AttributeType::VECTORS:     data->SetVectors(array); break;
    case AttributeType::NORMALS:     data->SetNormals(array); break;
    case AttributeType::TCOORDS:     data->SetTCoords(array); break;
    case AttributeType::TENSORS:     data->SetTensors(array); break;
    case AttributeType::GLOBALIDS:   data->SetGlobalIds(array); break;
    case AttributeType::PEDIGREEIDS: data->SetPedigreeIds(array); break;
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
  vtkDataSetAttributes         *sourceAttr, *targetAttr;
  vtkIdType                     ntuples;

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
  const char *point_mask_name = nullptr;
  const char *cell_mask_name  = nullptr;
  bool        case_sensitive  = false;
  LabelInterpolationMode label_mode = MajorityVote;

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
    else if (OPTION("-pd") || OPTION("-pointdata") || OPTION("-pd-as-cd") || OPTION("-pointdata-as-celldata") ||
             OPTION("-cd") || OPTION("-celldata")  || OPTION("-cd-as-pd") || OPTION("-celldata-as-pointdata")) {
      ArrayInfo info;
      if (strcmp(OPTNAME, "-pointdata-as-celldata") == 0 || strcmp(OPTNAME, "-pd-as-cd") == 0 ||
          strcmp(OPTNAME, "-celldata-as-pointdata") == 0 || strcmp(OPTNAME, "-cd-as-pd") == 0) {
        info._PointCellConversion = true;
      } else {
        info._PointCellConversion = false;
      }
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
      if (HAS_ARGUMENT) info._TargetDataType = ParseVtkDataType(ARGUMENT);
      if (strncmp(OPTNAME, "-pointdata", 10) == 0 || strncmp(OPTNAME, "-pd", 3) == 0) {
        pd.push_back(info);
      } else {
        cd.push_back(info);
      }
    }
    else if (OPTION("-pointdata-as-points") || OPTION("-pd-as-points")) {
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
    else if (OPTION("-majority"))  label_mode = MajorityVote;
    else if (OPTION("-unanimous")) label_mode = UnanimousVote;
    else if (OPTION("-pointmask")) point_mask_name = ARGUMENT;
    else if (OPTION("-cellmask"))  cell_mask_name  = ARGUMENT;
    else if (OPTION("-case-sensitive"))   case_sensitive = true;
    else if (OPTION("-case-insensitive")) case_sensitive = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Read input data sets
  vtkSmartPointer<vtkPointSet> source, target;
  target = ReadPointSet(target_name);
  #if MIRTK_IO_WITH_GIFTI
    const string source_ext = Extension(source_name, EXT_Last);
    if (source_ext == ".gii") {
      source = ReadGIFTI(source_name, vtkPolyData::SafeDownCast(target));
    } else if (source_ext == ".csv" || source_ext == ".tsv" || source_ext == ".txt") {
      char sep = ',';
      if (source_ext == ".tsv") sep = '\t';
      source = ReadPointSetTable(source_name, sep, target);
    }
  #endif
  if (source == nullptr) {
    source = ReadPointSet(source_name);
  }

  const vtkIdType npoints = target->GetNumberOfPoints();
  const vtkIdType ncells  = target->GetNumberOfCells();

  vtkPointData * const sourcePD = source->GetPointData();
  vtkCellData  * const sourceCD = source->GetCellData();

  vtkPointData * const targetPD = target->GetPointData();
  vtkCellData  * const targetCD = target->GetCellData();

  // Copy all point set attributes by default
  if (pd.empty() && cd.empty()) {
    int attr;
    for (int i = 0; i < sourcePD->GetNumberOfArrays(); ++i) {
      attr = sourcePD->IsArrayAnAttribute(i);
      if (attr < 0) attr = NUM_ATTRIBUTES;
      ArrayInfo info;
      info._SourceIndex     = i;
      info._TargetAttribute = static_cast<AttributeType>(attr);
      pd.push_back(info);
    }
    for (int i = 0; i < sourceCD->GetNumberOfArrays(); ++i) {
      attr = sourceCD->IsArrayAnAttribute(i);
      if (attr < 0) attr = NUM_ATTRIBUTES;
      ArrayInfo info;
      info._SourceIndex     = i;
      info._TargetAttribute = static_cast<AttributeType>(attr);
      cd.push_back(info);
    }
  }

  // Convert point data to cell data if necessary
  vtkSmartPointer<vtkCellData> pd_as_cd;
  for (size_t i = 0; i < pd.size(); ++i) {
    if (pd[i]._PointCellConversion) {
      vtkNew<vtkPointDataToCellData> pd_to_cd;
      SetVTKInput(pd_to_cd, source);
      pd_to_cd->PassPointDataOff();
      pd_to_cd->Update();
      pd_as_cd = pd_to_cd->GetOutput()->GetCellData();
      break;
    }
  }

  // Convert cell data to point data if necessary
  vtkSmartPointer<vtkPointData> cd_as_pd;
  for (size_t i = 0; i < cd.size(); ++i) {
    if (cd[i]._PointCellConversion) {
      vtkNew<vtkCellDataToPointData> cd_to_pd;
      SetVTKInput(cd_to_pd, source);
      cd_to_pd->PassCellDataOff();
      cd_to_pd->Update();
      cd_as_pd = cd_to_pd->GetOutput()->GetPointData();
      break;
    }
  }

  // Add point data
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
        copy = nullptr;
      }
    } else {
      if (pd[i]._PointCellConversion && pd[i]._TargetIndex != -2) {
        sourceAttr = pd_as_cd;
        targetAttr = targetCD;
        ntuples    = ncells;
      } else {
        sourceAttr = sourcePD;
        targetAttr = targetPD;
        ntuples    = npoints;
      }
      array = GetSourceArray("point", sourceAttr, pd[i], case_sensitive);
      if (pd[i]._TargetIndex == -2) {
        double p[3] = {.0};
        vtkPoints *points = target->GetPoints();
        vtkIdType end = min(npoints, source->GetNumberOfPoints());
        int       dim = min(3,       array->GetNumberOfComponents());
        for (vtkIdType ptId = 0; ptId < end; ++ptId) {
          for (int i = 0; i < dim; ++i) p[i] = array->GetComponent(ptId, i);
          points->SetPoint(ptId, p);
        }
        copy = nullptr;
      } else {
        if (sourceAttr == pd_as_cd.GetPointer()) {
          const auto sourceName = (array->GetName()  ? array->GetName()  : string());
          const auto targetName = (pd[i]._TargetName ? pd[i]._TargetName : string());
          if (IsCategoricalArrayName(sourceName) || IsCategoricalArrayName(targetName)) {
            vtkDataArray *labels = GetSourceArray("point", sourcePD, pd[i], case_sensitive);
            array = ConvertPointToCellLabels(source, labels, label_mode);
          }
        }
        copy = Copy(array, ntuples, pd[i]._TargetDataType);
      }
    }
    if (copy) {
      if (pd[i]._TargetName) copy->SetName(pd[i]._TargetName);
      AddArray(targetAttr, copy, pd[i]._TargetAttribute);
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
  for (size_t i = 0; i < cd.size(); ++i) {
    if (cd[i]._PointCellConversion) {
      sourceAttr = cd_as_pd;
      targetAttr = targetPD;
      ntuples    = npoints;
    } else {
      sourceAttr = sourceCD;
      targetAttr = targetCD;
      ntuples    = ncells;
    }
    array = GetSourceArray("cell", sourceAttr, cd[i], case_sensitive);
    if (sourceAttr == cd_as_pd.GetPointer()) {
      const auto sourceName = (array->GetName()  ? array->GetName()  : string());
      const auto targetName = (cd[i]._TargetName ? cd[i]._TargetName : string());
      if (IsCategoricalArrayName(sourceName) || IsCategoricalArrayName(targetName)) {
        vtkDataArray *labels = GetSourceArray("cell", sourceCD, cd[i], case_sensitive);
        array = ConvertCellToPointLabels(source, labels, label_mode);
      }
    }
    copy = Copy(array, ntuples, cd[i]._TargetDataType);
    if (cd[i]._TargetName) copy->SetName(cd[i]._TargetName);
    AddArray(targetAttr, copy, cd[i]._TargetAttribute);
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
