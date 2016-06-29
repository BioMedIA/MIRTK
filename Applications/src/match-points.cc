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

#include "mirtk/Transformation.h"
#include "mirtk/PointCorrespondence.h"
#include "mirtk/RegisteredPointSet.h"
#include "mirtk/PointCorrespondence.h"
#include "mirtk/PointSetIO.h"

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkCellArray.h"
#include "vtkIdTypeArray.h"
#include "vtkPolyData.h"
#include "vtkCellData.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <target> <source> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Writes for each point in a target point set its corresponding point in" << endl;
  cout << "  the source point set. The found correspondences can be either written" << endl;
  cout << "  to an output point set with indices and difference vectors stored as" << endl;
  cout << "  point data, or a text file listing for each target point index the" << endl;
  cout << "  index of the corresponding source point." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  target   Target point set." << endl;
  cout << "  source   Source point set." << endl;
  cout << endl;
  cout << "Optional arguments:" << endl;
  cout << "  -dofin <file>                Target point set transformation. (default: none)" << endl;
  cout << "  -closest-point               Find the closest point in the source data set. (default)" << endl;
  cout << "  -closest-cell                Find the closest cell point in the source data set." << endl;
  cout << "  -cor <name>                  Name of correspondence type, e.g., 'closest point', 'closest cell'." << endl;
  cout << "  -corpar <name> <value>       Set parameter of chosen correspondence type." << endl;
  cout << "  -feature <name> [<weight>]   Use named point data as feature. (default: spatial coordinates)" << endl;
  cout << "  -corin <file>                Read correspondences from text file which lists on each line" << endl;
  cout << "                               the indices of corresponding target and source points." << endl;
  cout << "  -corout <file>               Write indices of corresponding points to text file." << endl;
  cout << "  -output <file>               Write point set of corresponding points." << endl;
  cout << "  -vertices                    Add a vertex cell for each point." << endl;
  cout << "  -vectors                     Add displacement vectors from corresponding output points to" << endl;
  cout << "                               input target points. (default: off)" << endl;
  cout << "  -indices                     Store indices of corresponding points as point data." << endl;
  cout << "                               Indices are cell indices in case of :option:`-closest-cell`." << endl;
  cout << "  -copy-pointdata              Copy input point data to output point set. (default: off)" << endl;
  cout << "  -copy-celldata               Copy input cell data to output point set. (default: off)" << endl;
  cout << "  -[no]compress                Whether to compress XML VTK file data. (default: on)" << endl;
  cout << "  -binary                      Write legacy VTK file in binary format. (default)" << endl;
  cout << "  -ascii                       Write legacy VTK file in ASCII  format." << endl;
  PrintCommonOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char **argv)
{
  // Positional arguments
  EXPECTS_POSARGS(2);

  const char *target_name = POSARG(1);
  const char *source_name = POSARG(2);

  // Parse options
  const char *dofin_name   = nullptr;
  const char *corrout_name = nullptr;
  const char *output_name  = nullptr;
  FileOption  output_fopt  = FO_Default;
  bool        add_vertices = false;
  bool        add_indices  = false;
  bool        add_vectors  = false;
  bool        copy_pdata   = false;
  bool        copy_cdata   = false;

  PointCorrespondence::TypeId ctype = PointCorrespondence::ClosestPoint;
  ParameterList param;
  Array<string> feature_name;
  Array<double> feature_weight;

  for (ALL_OPTIONS) {
    if      (OPTION("-dofin"))          dofin_name = ARGUMENT;
    else if (OPTION("-closest-point"))  ctype      = PointCorrespondence::ClosestPoint;
    else if (OPTION("-closest-cell"))   ctype      = PointCorrespondence::ClosestCell;
    else if (OPTION("-corin")) {
      ctype = PointCorrespondence::FiducialMatch;
      Insert(param, "Correspondence map", ARGUMENT);
    }
    else if (OPTION("-c") || OPTION("-cor") || OPTION("-corr") || OPTION("-correspondence")) {
      const char *arg = ARGUMENT;
      if (!FromString(arg, ctype)) FatalError("Invalid -correspondence argument: " << arg);
    }
    else if (OPTION("-p") || OPTION("-corpar") || OPTION("-corrpar")) {
      Insert(param, ARGUMENT, ARGUMENT);
    }
    else if (OPTION("-feature") || OPTION("-f")) {
      feature_name.push_back(ARGUMENT);
      if (HAS_ARGUMENT) feature_weight.push_back(atof(ARGUMENT));
      else              feature_weight.push_back(1.0);
    }
    else if (OPTION("-corout") || OPTION("-corrout")) corrout_name = ARGUMENT;
    else if (OPTION("-o") || OPTION("-output")) output_name = ARGUMENT;
    else if (OPTION("-vertices"))   add_vertices = true;
    else if (OPTION("-novertices")) add_vertices = false;
    else if (OPTION("-indices"))    add_indices  = true;
    else if (OPTION("-noindices"))  add_indices  = false;
    else if (OPTION("-vectors")   || OPTION("-disp"))   add_vectors = true;
    else if (OPTION("-novectors") || OPTION("-nodisp")) add_vectors = false;
    else if (OPTION("-copy-pd") || OPTION("-copy-pointdata")) copy_pdata = true;
    else if (OPTION("-copy-cd") || OPTION("-copy-celldata"))  copy_cdata = true;
    else HANDLE_POINTSETIO_OPTION(output_fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  if (feature_name.empty()) {
    feature_name  .push_back("spatial coordinates");
    feature_weight.push_back(1.0);
  }

  // Read target and source data sets
  if (verbose) cout << "Read input point sets...", cout.flush();
  FileOption target_fopt;
  RegisteredPointSet target, source;
  target.InputPointSet(ReadPointSet(target_name, target_fopt));
  source.InputPointSet(ReadPointSet(source_name));
  if (output_fopt == FO_Default) output_fopt = target_fopt;
  if (verbose) cout << " done" << endl;

  // Read target transformation
  UniquePtr<Transformation> dof;
  if (dofin_name) {
    if (verbose) cout << "Reading transformation...", cout.flush();
    dof.reset(Transformation::New(dofin_name));
    if (verbose) cout << " done" << endl;
  }
  target.Transformation(dof.get());

  // Initialize (transformed) data sets
  source.Initialize(), source.Update();
  target.Initialize(), target.Update();
  const vtkIdType n = target.NumberOfPoints();

  if (target.NumberOfPoints() == 0) {
    FatalError("Failed to open target point set or point set contains no points");
  }
  if (source.NumberOfPoints() == 0) {
    FatalError("Failed to open source point set or point set contains no points");
  }
  if (verbose > 1) {
    cout << "Number of target points: " << target.NumberOfPoints() << endl;
    cout << "Number of source points: " << source.NumberOfPoints() << endl;
  }

  // Initialize point correspondence map
  if (verbose) cout << "Initialize correspondence map...", cout.flush();
  UniquePtr<PointCorrespondence> cmap(PointCorrespondence::New(ctype));
  cmap->FromTargetToSource(true);
  cmap->FromSourceToTarget(false);
  cmap->Parameter(param);
  cmap->Target(&target);
  cmap->Source(&source);
  for (size_t i = 0; i < feature_name.size(); ++i) {
    cmap->AddFeature(feature_name[i].c_str(), feature_weight[i]);
  }
  cmap->Initialize();
  if (verbose) cout << " done" << endl;

  // Find corresponding points
  if (verbose) cout << "Find corresponding points...", cout.flush();
  cmap->Update();
  if (verbose) cout << " done" << endl;

  // Write index map to text file
  if (corrout_name) {
    if (verbose) cout << "Write index map to text file...", cout.flush();
    Array<int> corr(n, -1);
    Point p;
    for (vtkIdType i = 0; i < n; ++i) {
      corr[i] = cmap->GetIndex(i);
      if (corr[i] < 0) {
        FatalError("Chosen correspondence is not a point to point correspondence map! Option -corout cannot be used.");
      }
    }
    ofstream ofs(corrout_name);
    for (vtkIdType i = 0; i < n; ++i) {
      ofs << i << "\t" << corr[i] << endl;
    }
    ofs.close();
    if (verbose) cout << " done" << endl;
  }

  // Write point set of corresponding points
  if (output_name) {
    if (verbose) cout << "Write output point set...", cout.flush();
    vtkSmartPointer<vtkPoints>      points;
    vtkSmartPointer<vtkIdTypeArray> indices;
    vtkSmartPointer<vtkFloatArray>  vectors;
    vtkSmartPointer<vtkCellArray>   vertices;
    vtkSmartPointer<vtkPointSet>    output;

    points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(n);

    if (add_vertices) {
      vertices = vtkSmartPointer<vtkCellArray>::New();
      vertices->Allocate(n);
    }
    if (add_indices) {
      indices = vtkSmartPointer<vtkIdTypeArray>::New();
      indices->SetName("indices");
      indices->SetNumberOfComponents(1);
      indices->SetNumberOfTuples(n);
    }
    if (add_vectors) {
      vectors = vtkSmartPointer<vtkFloatArray>::New();
      vectors->SetName("displacements");
      vectors->SetNumberOfComponents(3);
      vectors->SetNumberOfTuples(n);
    }

    Point p1, p2;
    int j;

    for (vtkIdType i = 0; i < n; ++i) {
      cmap->GetPoint(i, p1);
      points->SetPoint(i, p1._x, p1._y, p1._z);
      if (add_vertices) vertices->InsertNextCell(1, &i);
      if (add_indices) {
        j = cmap->GetIndex(i);
        if (j < 0) {
          FatalError("Chosen correspondence is not a point to point correspondence map! Option -indices cannot be used.");
        }
        indices->SetTuple1(i, j);
      }
      if (add_vectors) {
        target.GetPoint(i, p2);
        vectors->SetTuple3(i, p2._x - p1._x, p2._y - p1._y, p2._z - p1._z);
      }
    }

    output.TakeReference(target.InputPointSet()->NewInstance());
    output->DeepCopy(target.InputPointSet());
    if (!copy_pdata) output->GetPointData()->Initialize();
    if (!copy_cdata) output->GetCellData()->Initialize();
    output->SetPoints(points);
    if (add_vertices) {
      vtkPolyData *polydata = vtkPolyData::SafeDownCast(output);
      if (polydata) polydata->SetVerts(vertices);
    }
    if (add_indices)  output->GetPointData()->SetScalars(indices);
    if (add_vectors)  output->GetPointData()->SetVectors(vectors);
    WritePointSet(output_name, output, output_fopt);
    if (verbose) cout << " done" << endl;
  }

  return 0;
}
