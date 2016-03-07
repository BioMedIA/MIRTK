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

#include <mirtkCommon.h>
#include <mirtkOptions.h>

#include <mirtkPointSetUtils.h>
#include <mirtkVtkMath.h>

#include <vtkPointSet.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkFloatArray.h>
#include <vtkCellLocator.h>
#include <vtkPolyDataNormals.h>
#include <vtkModifiedBSPTree.h>
#include <vtkKdTreePointLocator.h>
#include <vtkDataReader.h> // VTK_ASCII, VTK_BINARY defines

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <target> <source> [<output>] [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Evaluate distance between two given point sets. With increased verbosity (see :option:`-v`)," << endl;
  cout << "  the mean and standard deviation of the measured distances  (verbosity level >=1) and the" << endl;
  cout << "  individual distance for each target point is reported (verbosity level >=2)." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  target   Target point set/surface." << endl;
  cout << "  source   Source point set/surface." << endl;
  cout << "  output   Output point set/surface with point data array of measured distances. (default: none)" << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -name <name>       Name of output point data array or column header. (default: Distance)" << endl;
  cout << "  -separator <sep>   Column separator used for table of measured point distances. (default: \\t)" << endl;
  cout << "  -point             Evaluate minimum distance to closest point. (default for point clouds)" << endl;
  cout << "  -cell              Evaluate minimum distance to closest surface point. (default for surface meshes)" << endl;
  cout << "  -normal            Evaluate minimum distance along surface normal." << endl;
  cout << "  -index             Evaluate distance between points with identical index, e.g., pairs of fiducial markers." << endl;
  cout << "  -hausdorff         Report Hausdorff distance. (default: off)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Evaluate distance of target to closest source points
void EvaluateClosestPointDistance(vtkPointSet *target, vtkPointSet *source, vtkDataArray *distance)
{
  vtkIdType otherPtId;
  double    a[3], b[3];

  vtkNew<vtkKdTreePointLocator> locator;
  locator->SetDataSet(source);
  locator->BuildLocator();

  for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
    target->GetPoint(ptId, a);
    otherPtId = locator->FindClosestPoint(a);
    source->GetPoint(otherPtId, b);
    distance->SetComponent(ptId, 0, sqrt(vtkMath::Distance2BetweenPoints(a, b)));
  }
}

// -----------------------------------------------------------------------------
/// Evaluate distance of target to closest source cell point
void EvaluateClosestCellDistance(vtkPointSet *target, vtkPointSet *source, vtkDataArray *distance)
{
  int       subId;
  vtkIdType cellId;
  double    a[3], b[3], dist2;

  vtkNew<vtkCellLocator> locator;
  locator->SetDataSet(source);
  locator->BuildLocator();

  for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
    target->GetPoint(ptId, a);
    locator->FindClosestPoint(a, b, cellId, subId, dist2);
    distance->SetComponent(ptId, 0, sqrt(dist2));
  }
}

// -----------------------------------------------------------------------------
/// Evaluate distance of corresponding points
void EvaluateCorrespondingPointDistance(vtkPointSet *target, vtkPointSet *source, vtkDataArray *distance)
{
  if (target->GetNumberOfPoints() != source->GetNumberOfPoints()) {
    FatalError("Point sets must have equal number of points to evaluate -corr point distance");
  }

  double a[3], b[3];

  for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
    target->GetPoint(ptId, a);
    source->GetPoint(ptId, b);
    distance->SetComponent(ptId, 0, sqrt(vtkMath::Distance2BetweenPoints(a, b)));
  }
}

// -----------------------------------------------------------------------------
/// Evaluate distance of target to closest source cell point in normal direction
void EvaluateNormalDistance(vtkPointSet *target, vtkPointSet *source, vtkDataArray *distance,
                            double default_value = numeric_limits<double>::quiet_NaN())
{
  vtkPolyData *surface = vtkPolyData::SafeDownCast(target);
  if (surface == NULL) {
    FatalError("Cannot compute distance along -normal for non-polygonal point set");
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
    dist = numeric_limits<double>::infinity();
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
    if (IsInf(dist)) dist = default_value;
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

  const char  *array_name     = "Distance";
  const char  *separator      = "\t";
  DistanceType dist_type      = Default;
  bool         eval_hausdorff = false;

  for (ALL_OPTIONS) {
    if      (OPTION("-array") || OPTION("-name")) array_name = ARGUMENT;
    else if (OPTION("-point"))     dist_type      = ClosestPoint;
    else if (OPTION("-cell"))      dist_type      = ClosestCell;
    else if (OPTION("-normal"))    dist_type      = AlongNormal;
    else if (OPTION("-index"))     dist_type      = CorrespondingPoint;
    else if (OPTION("-hausdorff")) eval_hausdorff = true;
    else if (OPTION("-sep") || OPTION("-separator") || OPTION("-delimiter")) {
      separator = ARGUMENT;
      if (strcmp(separator, "\\t") == 0) separator = "\t";
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // read input point sets
  int output_type = VTK_BINARY;
  vtkSmartPointer<vtkPointSet> target, source;
  target = ReadPointSet(target_name);
  source = ReadPointSet(source_name, &output_type);

  if (dist_type == Default) {
    if (target->GetNumberOfCells() > 0) dist_type = ClosestCell;
    else                                dist_type = ClosestPoint;
  }
  if (dist_type == ClosestCell) {
    if (source->GetNumberOfCells() == 0) {
      FatalError("Cannot evaluation -cell distance when source has no cells");
    }
    if (eval_hausdorff && target->GetNumberOfCells() == 0) {
      FatalError("Cannot evaluation -cell -hausdorff distance when target has no cells");
    }
  }

  if (eval_hausdorff && dist_type != ClosestPoint && dist_type != ClosestCell) {
    FatalError("Hausdorff distance currenly only implemented for -point or -cell distance");
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
    case ClosestPoint:       EvaluateClosestPointDistance      (target, source, distance); break;
    case ClosestCell:        EvaluateClosestCellDistance       (target, source, distance); break;
    case AlongNormal:        EvaluateNormalDistance            (target, source, distance); break;
    case CorrespondingPoint: EvaluateCorrespondingPointDistance(target, source, distance); break;
  }

  // calculate min/max and average distance
  double sum_dist = .0;
  double min_dist = +numeric_limits<double>::infinity();
  double max_dist = -numeric_limits<double>::infinity();

  for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
    dist = distance->GetComponent(ptId, 0);
    sum_dist += dist;
    if (dist < min_dist) min_dist = dist;
    if (dist > max_dist) max_dist = dist;
  }

  // calculate Hausdorff distance
  double hausdorff_dist = max_dist;
  if (eval_hausdorff) {
    if (dist_type == ClosestCell) {

      vtkNew<vtkCellLocator> locator;
      locator->SetDataSet(target);
      locator->BuildLocator();

      for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
        source->GetPoint(ptId, a);
        locator->FindClosestPoint(a, b, cellId, subId, dist);
        dist = sqrt(dist);
        if (dist > hausdorff_dist) hausdorff_dist = dist;
      }

    } else {

      vtkNew<vtkKdTreePointLocator> locator;
      locator->SetDataSet(target);
      locator->BuildLocator();

      for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
        source->GetPoint(ptId, a);
        otherPtId = locator->FindClosestPoint(a);
        source->GetPoint(otherPtId, b);
        dist = sqrt(vtkMath::Distance2BetweenPoints(a, b));
        if (dist > hausdorff_dist) hausdorff_dist = dist;
      }

    }
  }

  // write output
  if (output_name) {
    const bool compress = true;
    if (!WritePointSet(output_name, target, compress, output_type == VTK_ASCII)) {
      FatalError("Failed to write output point set to " << output_name);
    }
  }

  // print distances and summary
  if (verbose > 1) {
    cout << "ID" << separator << array_name << "\n";
    for (vtkIdType ptId = 0; ptId < target->GetNumberOfPoints(); ++ptId) {
      cout << (ptId+1) << separator << distance->GetComponent(ptId, 0) << "\n";
    }
    cout << "\n";
  }
  if (verbose || !output_name) {
    cout << "Minimum distance   = " << min_dist << "\n";
    cout << "Maximum distance   = " << max_dist << "\n";
    cout << "Average distance   = " << sum_dist / target->GetNumberOfPoints() << "\n";
  }
  if (eval_hausdorff) cout << "Hausdorff distance = " << hausdorff_dist << "\n";

  return 0;
}
