/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2017 Imperial College London
 * Copyright 2013-2017 Andreas Schuh
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
#include "mirtk/DilatePointData.h"
#include "mirtk/ErodePointData.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/EdgeConnectivity.h"

#include "mirtk/Vtk.h"
#include "mirtk/VtkMath.h"

#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkFloatArray.h"
#include "vtkCellLocator.h"
#include "vtkPolyDataNormals.h"
#include "vtkModifiedBSPTree.h"
#include "vtkKdTreePointLocator.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <target> <source> [<output>] [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Evaluate distance between two given point sets. With increased verbosity (see :option:`-v`),\n";
  cout << "  the mean and standard deviation of the measured distances  (verbosity level >=1) and the\n";
  cout << "  individual distance for each target point is reported (verbosity level >=2).\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  target   Target point set/surface.\n";
  cout << "  source   Source point set/surface.\n";
  cout << "  output   Output point set/surface with point data array of measured distances. (default: none)\n";
  cout << "\n";
  cout << "Optional arguments:\n";
  cout << "  -mask-name <name>\n";
  cout << "      Name of input mask. (default: none)\n";
  cout << "  -mask-erosion <n>\n";
  cout << "      Number of iterations of mask erosions. (default: 0)\n";
  cout << "  -mask-dilation <n>\n";
  cout << "      Number of iterations of mask dilations. (default: 0)\n";
  cout << "  -pad-value, -padding-value <value>\n";
  cout << "      Floating point value for points with zero mask value. (default: NaN)\n";
  cout << "  -undef-value, -undefined-value <value>\n";
  cout << "      Value for undefined distance measures, e.g., normal does not intersect source. (default: NaN)\n";
  cout << "  -nan-value <value>\n";
  cout << "      Set both :option:`-pad-value` and :option:`-undef-value` to the specified value.\n";
  cout << "  -dist-name, -name, -array <name>\n";
  cout << "      Name of output point data array or column header. (default: Distance)\n";
  cout << "  -digits, -precision <n>\n";
  cout << "      Number of digits after the decimal point to print/write. (default: 5)\n";
  cout << "  -table [<name>|stdout|cout|print]\n";
  cout << "      Write statistics of distance measure in table format to named output file.\n";
  cout << "      When <name> is 'stdout', 'cout',  or 'print', print table to standard output stream. (default: off)\n";
  cout << "  -[no]append\n";
  cout << "      Whether to append row to existing :option:`-table` file.\n";
  cout << "      If table file exists, skip header and write measurements only. (default: off)\n";
  cout << "  -[no]header\n";
  cout << "      When :option:`-table` given, write table header before row of measurements. (default: on)\n";
  cout << "  -delim, -delimiter, -sep, -separator <str>\n";
  cout << "      Column separator used for table of measured point distances. (default: \\t)\n";
  cout << "  -point\n";
  cout << "      Evaluate minimum distance to closest point. (default for point clouds)\n";
  cout << "  -cell\n";
  cout << "      Evaluate minimum distance to closest surface point. (default for surface meshes)\n";
  cout << "  -normal\n";
  cout << "      Evaluate minimum distance along surface normal.\n";
  cout << "  -index\n";
  cout << "      Evaluate distance between points with identical index, e.g., pairs of fiducial markers.\n";
  cout << "  -hausdorff\n";
  cout << "      Report Hausdorff distance. (default: off)\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Evaluate distance of target to closest source points
void EvaluateClosestPointDistance(vtkPointSet *target, vtkPointSet *source, vtkDataArray *distance,
                                  vtkDataArray *mask = nullptr, double pad_value = NaN, double nan_value = NaN)
{
  vtkIdType otherPtId;
  double    a[3], b[3], dist;

  vtkNew<vtkKdTreePointLocator> locator;
  locator->SetDataSet(source);
  locator->BuildLocator();

  for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
    if (!mask || mask->GetComponent(ptId, 0) != 0) {
      target->GetPoint(ptId, a);
      otherPtId = locator->FindClosestPoint(a);
      source->GetPoint(otherPtId, b);
      dist = sqrt(vtkMath::Distance2BetweenPoints(a, b));
      distance->SetComponent(ptId, 0, IsNaN(dist) ? nan_value : dist);
    } else {
      distance->SetComponent(ptId, 0, pad_value);
    }
  }
}

// -----------------------------------------------------------------------------
/// Evaluate distance of target to closest source cell point
void EvaluateClosestCellDistance(vtkPointSet *target, vtkPointSet *source, vtkDataArray *distance,
                                 vtkDataArray *mask = nullptr, double pad_value = NaN, double nan_value = NaN)
{
  int       subId;
  vtkIdType cellId;
  double    a[3], b[3], dist;

  vtkNew<vtkCellLocator> locator;
  locator->SetDataSet(source);
  locator->BuildLocator();

  for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
    if (!mask || mask->GetComponent(ptId, 0) != 0) {
      target->GetPoint(ptId, a);
      locator->FindClosestPoint(a, b, cellId, subId, dist);
      dist = sqrt(dist);
      distance->SetComponent(ptId, 0, IsNaN(dist) ? nan_value : dist);
    } else {
      distance->SetComponent(ptId, 0, pad_value);
    }
  }
}

// -----------------------------------------------------------------------------
/// Evaluate distance of corresponding points
void EvaluateCorrespondingPointDistance(vtkPointSet *target, vtkPointSet *source, vtkDataArray *distance,
                                        vtkDataArray *mask = nullptr, double pad_value = NaN, double nan_value = NaN)
{
  double a[3], b[3], dist;

  if (target->GetNumberOfPoints() != source->GetNumberOfPoints()) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Point sets must have equal number of points!");
  }
  if (mask && mask->GetNumberOfTuples() != target->GetNumberOfPoints()) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Mask array must have one tuple for each target point!");
  }

  for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
    if (!mask || mask->GetComponent(ptId, 0) != 0) {
      target->GetPoint(ptId, a);
      source->GetPoint(ptId, b);
      dist = sqrt(vtkMath::Distance2BetweenPoints(a, b));
      distance->SetComponent(ptId, 0, IsNaN(dist) ? nan_value : dist);
    } else {
      distance->SetComponent(ptId, 0, pad_value);
    }
  }
}

// -----------------------------------------------------------------------------
/// Evaluate distance of target to closest source cell point in normal direction
void EvaluateNormalDistance(vtkPointSet *target, vtkPointSet *source, vtkDataArray *distance,
                            vtkDataArray *mask = nullptr, double pad_value = NaN, double nan_value = NaN)
{
  vtkPolyData *surface = vtkPolyData::SafeDownCast(target);
  if (surface == nullptr) {
    Throw(ERR_InvalidArgument, __FUNCTION__, "Cannot compute distance along normal for non-polygonal point set!");
  }

  vtkNew<vtkPolyDataNormals> calc_normals;
  calc_normals->ComputeCellNormalsOff();
  calc_normals->ComputePointNormalsOn();
  calc_normals->AutoOrientNormalsOn();
  calc_normals->SplittingOff();
  SetVTKInput(calc_normals, surface);
  calc_normals->Update();

  vtkSmartPointer<vtkDataArray> normals;
  normals = calc_normals->GetOutput()->GetPointData()->GetNormals();

  int    subId;
  double a[3], b[3], n[3], x[3], t, pcoords[3], dist;

  const double max_dist = 1000.0;

  vtkNew<vtkModifiedBSPTree> locator;
  locator->SetDataSet(source);
  locator->BuildLocator();

  for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
    if (!mask || mask->GetComponent(ptId, 0) != 0) {
      dist = inf;
      target->GetPoint(ptId, a);
      normals->GetTuple(ptId, n);
      b[0] = a[0] + max_dist * n[0];
      b[1] = a[1] + max_dist * n[1];
      b[2] = a[2] + max_dist * n[2];
      if (locator->IntersectWithLine(a, b, .0, t, x, pcoords, subId)) {
        dist = sqrt(vtkMath::Distance2BetweenPoints(a, x));
      }
      b[0] = a[0] - max_dist * n[0];
      b[1] = a[1] - max_dist * n[1];
      b[2] = a[2] - max_dist * n[2];
      if (locator->IntersectWithLine(a, b, .0, t, x, pcoords, subId)) {
        dist = min(dist, sqrt(vtkMath::Distance2BetweenPoints(a, x)));
      }
      if (IsInf(dist)) dist = nan_value;
    } else {
      dist = pad_value;
    }
    distance->SetComponent(ptId, 0, dist);
  }
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Parse arguments
  REQUIRES_POSARGS(2);

  const char *target_name = POSARG(1);
  const char *source_name = POSARG(2);
  const char *output_name = nullptr;

  if (NUM_POSARGS == 3) {
    output_name = POSARG(3);
  } else if (NUM_POSARGS > 3) {
    PrintHelp(EXECNAME);
    exit(1);
  }

  enum DistanceType
  {
    Default,
    ClosestPoint,
    ClosestCell,
    AlongNormal,
    CorrespondingPoint
  };

  enum MaskOperation
  {
    DilateMask,
    ErodeMask
  };

  const char  *array_name     = "Distance";
  const char  *separator      = nullptr;
  const char  *table_name     = nullptr;
  bool         table_header   = true;
  bool         table_append   = false;
  int          table_digits   = 5;
  DistanceType dist_type      = Default;
  bool         eval_hausdorff = false;
  FileOption   fopt           = FO_Default;
  double       nan_value      = NaN;
  double       pad_value      = NaN;
  const char  *mask_name      = nullptr;
  Array<Pair<MaskOperation, int> > mask_ops;

  for (ALL_OPTIONS) {
    if (OPTION("-dist-name") || OPTION("-array") || OPTION("-name")) {
      array_name = ARGUMENT;
    }
    else if (OPTION("-mask-name")) {
      mask_name = ARGUMENT;
    }
    else if (OPTION("-mask-erosion")) {
      int n;
      PARSE_ARGUMENT(n);
      if (n > 0) {
        mask_ops.push_back(MakePair(ErodeMask, n));
      }
    }
    else if (OPTION("-mask-dilation")) {
      int n;
      PARSE_ARGUMENT(n);
      if (n > 0) {
        mask_ops.push_back(MakePair(DilateMask, n));
      }
    }
    else if (OPTION("-undefined-value") || OPTION("-undef-value")) {
      PARSE_ARGUMENT(nan_value);
    }
    else if (OPTION("-padding-value") || OPTION("-pad-value")) {
      PARSE_ARGUMENT(pad_value);
    }
    else if (OPTION("-nan-value")) {
      PARSE_ARGUMENT(nan_value);
      pad_value = nan_value;
    }
    else if (OPTION("-point"))  dist_type = ClosestPoint;
    else if (OPTION("-cell"))   dist_type = ClosestCell;
    else if (OPTION("-normal")) dist_type = AlongNormal;
    else if (OPTION("-index"))  dist_type = CorrespondingPoint;
    else if (OPTION("-hausdorff")) eval_hausdorff = true;
    else if (OPTION("-sep") || OPTION("-separator") || OPTION("-delimiter") || OPTION("-delim")) {
      separator = ARGUMENT;
      if (strcmp(separator, "\\t") == 0) separator = "\t";
    }
    else if (OPTION("-table")) {
      if (HAS_ARGUMENT) {
        table_name = ARGUMENT;
      } else {
        table_name = "stdout";
      }
    }
    else if (OPTION("-digits") || OPTION("-precision")) {
      PARSE_ARGUMENT(table_digits);
    }
    else HANDLE_BOOLEAN_OPTION("header", table_header);
    else HANDLE_BOOLEAN_OPTION("append", table_append);
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }
  if (!mask_name && !mask_ops.empty()) {
    Warning("Ignoring -mask-* options because no -mask-name is specified.");
  }
  if (table_name) {
    const auto table_lname = ToLower(table_name);
    if (table_lname == "" || table_lname == "cout" || table_lname == "stdout" || table_lname == "print") {
      table_name = "stdout";
    } else if (!separator) {
      const auto ext = Extension(table_lname, EXT_Last);
      if      (ext == ".csv") separator = ",";
      else if (ext == ".tsv") separator = "\t";
      else                    separator = " ";
    }
  }
  if (!separator) separator = "\t";

  // read input point sets
  FileOption source_fopt;
  vtkSmartPointer<vtkPointSet> target, source;
  target = ReadPointSet(target_name);
  source = ReadPointSet(source_name, source_fopt);
  if (fopt == FO_Default) fopt = source_fopt;

  if (dist_type == Default) {
    if (target->GetNumberOfCells() > 0) dist_type = ClosestCell;
    else                                dist_type = ClosestPoint;
  }
  if (dist_type == ClosestCell) {
    if (source->GetNumberOfCells() == 0) {
      FatalError("Cannot evaluate -cell distance when source has no cells");
    }
    if (eval_hausdorff && target->GetNumberOfCells() == 0) {
      FatalError("Cannot evaluate -cell -hausdorff distance when target has no cells");
    }
  }
  if (eval_hausdorff && dist_type != ClosestPoint && dist_type != ClosestCell) {
    FatalError("Hausdorff distance currenly only implemented for -point or -cell distance");
  }

  // check mask array
  vtkSmartPointer<vtkDataArray> mask;
  if (mask_name) {
    mask = GetArrayByCaseInsensitiveName(target->GetPointData(), mask_name);
    if (!mask) {
      FatalError("Input target point set does not have a point data array named '" << mask_name << "'");
    }
    if (!mask_ops.empty()) {
      vtkPolyData *surface = vtkPolyData::SafeDownCast(target);
      if (!surface) {
        FatalError("Morphological mask operations only supported for surfaces");
      }
      SharedPtr<EdgeTable> edgeTable(new EdgeTable(surface));
      SharedPtr<EdgeConnectivity> neighbors(new EdgeConnectivity(surface, 1, edgeTable.get()));
      for (const auto op : mask_ops) {
        if (op.second > 0) {
          if (op.first == DilateMask) {
            DilatePointData dilate;
            dilate.Input(surface);
            dilate.InputData(mask);
            dilate.EdgeTable(edgeTable);
            dilate.Neighbors(neighbors);
            dilate.Iterations(op.second);
            dilate.Run();
            mask = dilate.OutputData();
          } else if (op.first == ErodeMask) {
            ErodePointData erode;
            erode.Input(surface);
            erode.InputData(mask);
            erode.EdgeTable(edgeTable);
            erode.Neighbors(neighbors);
            erode.Iterations(op.second);
            erode.Run();
            mask = erode.OutputData();
          } else {
            Throw(ERR_LogicError, EXECNAME, "Unknown mask operation: ", op.first);
          }
        }
      }
    }
  }

  // find closest source points
  int       subId;
  vtkIdType otherPtId, cellId;
  double    a[3], b[3], dist;

  vtkSmartPointer<vtkFloatArray> distance;
  distance = vtkSmartPointer<vtkFloatArray>::New();
  distance->SetName(array_name);
  distance->SetNumberOfComponents(1);
  distance->SetNumberOfTuples(target->GetNumberOfPoints());
  target->GetPointData()->AddArray(distance);

  switch (dist_type) {
    case Default: // never the case
    case ClosestPoint:       EvaluateClosestPointDistance      (target, source, distance, mask, pad_value, nan_value); break;
    case ClosestCell:        EvaluateClosestCellDistance       (target, source, distance, mask, pad_value, nan_value); break;
    case AlongNormal:        EvaluateNormalDistance            (target, source, distance, mask, pad_value, nan_value); break;
    case CorrespondingPoint: EvaluateCorrespondingPointDistance(target, source, distance, mask, pad_value, nan_value); break;
  }

  // calculate min/max and average distance
  Array<double> dists;
  dists.reserve(target->GetNumberOfPoints());
  for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
    if (!mask || mask->GetComponent(ptId, 0) != 0.) {
      dist = distance->GetComponent(ptId, 0);
      if (!IsNaN(dist) && !IsInf(dist)) {
        dists.push_back(dist);
      }
    }
  }
  if (dists.empty()) {
    Warning("No. of valid target point distance measures = 0");
  }
  const auto range = MinMaxElement(dists);

  // calculate Hausdorff distance
  double hausdorff_dist = range.second;
  if (eval_hausdorff) {
    if (dist_type == ClosestCell) {

      vtkNew<vtkCellLocator> locator;
      locator->SetDataSet(target);
      locator->BuildLocator();
      vtkSmartPointer<vtkIdList> ptIds;
      ptIds = vtkSmartPointer<vtkIdList>::New();

      for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
        source->GetPoint(ptId, a);
        locator->FindClosestPoint(a, b, cellId, subId, dist);
        dist = sqrt(dist);
        if (mask) {
          target->GetCellPoints(cellId, ptIds);
          for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
            if (mask->GetComponent(ptIds->GetId(i), 0) != 0.) {
              if (dist > hausdorff_dist) hausdorff_dist = dist;
              break;
            }
          }
        } else {
          if (dist > hausdorff_dist) hausdorff_dist = dist;
        }
      }

    } else {

      vtkNew<vtkKdTreePointLocator> locator;
      locator->SetDataSet(target);
      locator->BuildLocator();

      for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
        source->GetPoint(ptId, a);
        otherPtId = locator->FindClosestPoint(a);
        if (!mask || mask->GetComponent(otherPtId, 0) != 0.) {
          source->GetPoint(otherPtId, b);
          dist = sqrt(vtkMath::Distance2BetweenPoints(a, b));
          if (dist > hausdorff_dist) hausdorff_dist = dist;
        }
      }

    }
  }

  // write output
  if (output_name) {
    if (!WritePointSet(output_name, target, fopt)) {
      FatalError("Failed to write output point set to " << output_name);
    }
  }

  // print distances and summary
  if (verbose > 1) {
    cout << "ID" << separator << array_name << "\n";
    for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
      cout << (ptId + 1) << separator << distance->GetComponent(ptId, 0) << "\n";
    }
    cout << "\n";
  }
  if (verbose || !output_name || table_name) {
    const auto mean = Accumulate(dists) / dists.size();
    const auto median = MedianValue(dists);
    if (table_name) {
      ofstream ofs;
      ostream *os = nullptr;
      if (strcmp(table_name, "stdout") == 0) {
        os = &cout;
      } else {
        if (table_append) {
          ifstream ifs(table_name);
          if (ifs.is_open()) {
            table_header = false;
            ifs.close();
          } else {
            table_append = false;
          }
        }
        ofs.open(table_name, table_append ? ostream::app : ostream::out);
        if (!ofs.is_open()) {
          FatalError("Failed to create/open output table file: " << table_name);
        }
        os = &ofs;
      }
      if (table_header) {
        (*os) << "min" << separator << "max" << separator << "mean" << separator << "median";
        if (eval_hausdorff) {
          (*os) << separator << "hausdorff";
        }
        (*os) << "\n";
      }
      fixed(*os);
      os->precision(table_digits);
      (*os) << range.first << separator << range.second << separator << mean << separator << median;
      if (eval_hausdorff) {
        (*os) << separator << hausdorff_dist;
      }
      (*os) << "\n";
      ofs.close();
    }
    if (!table_name && (verbose || !output_name)) {
      fixed(cout);
      cout.precision(table_digits);
      cout << "Minimum distance   = " << range.first << "\n";
      cout << "Maximum distance   = " << range.second << "\n";
      cout << "Average distance   = " << mean << "\n";
      cout << "Median distance    = " << median << "\n";
    }
  }
  if (!table_name && eval_hausdorff) {
    cout << "Hausdorff distance = " << hausdorff_dist << "\n";
  }

  return 0;
}
