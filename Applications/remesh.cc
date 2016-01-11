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

#include <mirtkEdgeTable.h>
#include <mirtkPointSetUtils.h>
#include <mirtkPolyDataRemeshing.h>

#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkCellLocator.h>
#include <vtkMath.h>

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Remeshes a surface mesh such that the resulting mesh has an average" << endl;
  cout << "  edge length within the specified limits." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -write-all                         Write also intermediate meshes when more than one" << endl;
  cout << "                                     desired average edge length was specified. Output file" << endl;
  cout << "                                     names contain the level number as suffix. (default: off)" << endl;
  cout << "  -target <file>                     Find closest cell points on target surface" << endl;
  cout << "                                     and use these as new point locations." << endl;
  cout << "  -edgelength    <float>...          Average edge length." << endl;
  cout << "  -minedgelength <float>...          Minimum edge length. (default: 0)" << endl;
  cout << "  -maxedgelength <float>...          Maximum edge length. (default: inf)" << endl;
  cout << "  -adaptiveedgelength <name>         Name of point data array to use for adapting the edge length range. (default: none)" << endl;
  cout << "  -meltorder <index|area|shortest>   Order in which to process cells in melting pass. (default: area)" << endl;
  cout << "  -[no]meltnodes                     Whether to allow removal of nodes with connectivity 3. (default: off)" << endl;
  cout << "  -[no]melttriangles                 Whether to melt triangles when all three edges are too short. (default: off)" << endl;
  cout << "  -ascii/-binary                     Write legacy VTK in ASCII or binary format. (default: binary)" << endl;
  cout << "  -[no]compress                      Write XML VTK file with or without compression. (default: compress)" << endl;
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

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  string dir  = Directory(output_name);
  string name = FileName(output_name);
  string ext  = Extension(output_name);

  Array<double> minlength;
  Array<double> maxlength;
  Array<double> minangle;
  Array<double> maxangle;

  bool        melt_nodes     = false;
  bool        melt_triangles = false;
  const char *adaptive_name  = NULL;
  const char *target_name    = NULL;
  bool        ascii          = false;
  bool        compress       = true;
  bool        write_all      = false;

  PolyDataRemeshing::Order melt_order = PolyDataRemeshing::AREA;

  for (ALL_OPTIONS) {
    if (OPTION("-minedgelength")) {
      do {
        minlength.push_back(atof(ARGUMENT));
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-maxedgelength")) {
      do {
        maxlength.push_back(atof(ARGUMENT));
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-edgelength")) {
      if (minlength.size() < maxlength.size()) minlength.resize(maxlength.size(), minlength.back());
      if (minlength.size() > maxlength.size()) maxlength.resize(minlength.size(), maxlength.back());
      do {
        double l = atof(ARGUMENT);
        minlength.push_back(l);
        maxlength.push_back(l);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-adaptive") || OPTION("-adaptiveedgelength")) {
      adaptive_name = ARGUMENT;
    }
    else if (OPTION("-minangle")) {
      do {
        minangle.push_back(atof(ARGUMENT));
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-maxangle")) {
      do {
        maxangle.push_back(atof(ARGUMENT));
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-meltorder")) {
      string arg = ARGUMENT;
      transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
      if (arg == "area" || arg == "smallest area" || arg == "smallestarea")
        melt_order = PolyDataRemeshing::AREA;
      else if (arg == "edge" ||
               arg == "shortest" ||
               arg == "shortest edge" ||
               arg == "shortestedge")
        melt_order = PolyDataRemeshing::SHORTEST_EDGE;
      else if (arg == "index" || arg == "none")
        melt_order = PolyDataRemeshing::INDEX;
      else {
        cerr << "Invalid -meltorder: " << arg << endl;
        exit(1);
      }
    }
    else if (OPTION("-write-all")) write_all = true;
    else if (OPTION("-meltnodes")) melt_nodes = true;
    else if (OPTION("-nomeltnodes")) melt_nodes = false;
    else if (OPTION("-melttriangles")) melt_triangles = true;
    else if (OPTION("-nomelttriangles")) melt_triangles = false;
    else if (OPTION("-target")) target_name = ARGUMENT;
    else if (OPTION("-ascii" ) || OPTION("-nobinary")) ascii = true;
    else if (OPTION("-binary") || OPTION("-noascii" )) ascii = false;
    else if (OPTION("-compress"))   compress = true;
    else if (OPTION("-nocompress")) compress = false;
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  const size_t nlevels = max(max(minlength.size(), maxlength.size()), maxangle.size());
  const double &inf = numeric_limits<double>::infinity();
  minlength.resize(nlevels, minlength.empty() ?   .0  : minlength.back());
  maxlength.resize(nlevels, maxlength.empty() ?  inf  : maxlength.back());
  minangle .resize(nlevels, minangle .empty() ? 180.0 : minangle .back());
  maxangle .resize(nlevels, maxangle .empty() ? 180.0 : maxangle .back());
  for (size_t i = 0; i < nlevels; ++i) {
    if (minangle[i] > maxangle[i]) minangle[i] = maxangle[i];
  }

  // Read input surface
  vtkSmartPointer<vtkPolyData> mesh = ReadPolyData(input_name);

  const vtkIdType npoints = mesh->GetNumberOfPoints();
  const vtkIdType ncells  = mesh->GetNumberOfCells();

  if (verbose) {
    cout << "No. of input vertices    = " << npoints << endl;
    cout << "No. of input triangles   = " << ncells  << endl;
  }

  // Use point locations of closest target surface as initial vertex positions
  if (target_name) {
    vtkSmartPointer<vtkPolyData>    target  = ReadPolyData(target_name);
    vtkSmartPointer<vtkCellLocator> locator = vtkSmartPointer<vtkCellLocator>::New();
    locator->SetDataSet(target);
    locator->BuildLocator();
    double p[3], dist2;
    vtkIdType cellId;
    int       subId;
    for (vtkIdType i = 0; i < mesh->GetNumberOfPoints(); ++i) {
      locator->FindClosestPoint(mesh->GetPoint(i), p, cellId, subId, dist2);
      mesh->GetPoints()->SetPoint(i, p);
    }
  }

  // Remesh to desired average edge length
  PolyDataRemeshing remesher;
  remesher.MeltingOrder(melt_order);
  remesher.MeltNodes(melt_nodes);
  remesher.MeltTriangles(melt_triangles);
  if (adaptive_name) {
    vtkDataArray *array = mesh->GetPointData()->GetArray(adaptive_name);
    if (!array) {
      cerr << "Error: Input mesh has no point data array named: " << adaptive_name << endl;
      exit(1);
    }
    remesher.AdaptiveEdgeLengthArray(array);
  }

  for (size_t i = 0; i < nlevels; ++i) {
    if (verbose) {
      cout << "\nRemeshing surface with edge length range [" << minlength[i] << ", " << maxlength[i] << "]";
      if (minangle[i] < 180.0) cout << " and feature angle in [" << minangle[i] << ", " << maxangle[i] << "] degrees";
      cout << endl;
    }
    remesher.Input(mesh);
    remesher.MinEdgeLength  (minlength[i]);
    remesher.MaxEdgeLength  (maxlength[i]);
    remesher.MinFeatureAngle(minangle [i]);
    remesher.MaxFeatureAngle(maxangle [i]);
    remesher.Run();
    if (verbose > 1) {
      cout << "No. of node-meltings     = " << remesher.NumberOfMeltedNodes() << endl;
      cout << "No. of edge-meltings     = " << remesher.NumberOfMeltedEdges() << endl;
      cout << "No. of triangle-meltings = " << remesher.NumberOfMeltedCells() << endl;
      cout << "No. of inversions        = " << remesher.NumberOfInversions() << endl;
      cout << "No. of bisections        = " << remesher.NumberOfBisections() << endl;
      cout << "No. of trisections       = " << remesher.NumberOfTrisections() << endl;
      cout << "No. of quadsections      = " << remesher.NumberOfQuadsections() << endl;
    }
    mesh = remesher.Output();
    if (verbose) {
      cout << "No. of output vertices   = " << mesh->GetNumberOfPoints()
           << " (" << (100.0 * mesh->GetNumberOfPoints() / npoints) << "%)" << endl;
      cout << "No. of output triangles  = " << mesh->GetNumberOfCells()
           << " (" << (100.0 * mesh->GetNumberOfCells() / ncells) << "%)" << endl;
    }
    if (write_all && i < minlength.size() - 1) {
      ostringstream os;
      if (!dir.empty()) os << dir << PATHSEP;
      os << name << "_" << (i + 1) << ext;
      if (verbose > 1) cout << "Write " << os.str() << endl;
      WritePolyData(os.str().c_str(), mesh);
    }
  }

  // Report properties of resulting surface mesh
  if (verbose) {
    vtkIdType ptId1, ptId2;
    double p1[3], p2[3], l, sum = .0, max = .0;
    double min = numeric_limits<double>::infinity();

    EdgeTable edgeTable(mesh);
    EdgeIterator it(edgeTable);
    for (it.InitTraversal(); it.GetNextEdge(ptId1, ptId2) != -1;) {
      mesh->GetPoint(ptId1, p1);
      mesh->GetPoint(ptId2, p2);
      l = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      if (l < min) min = l;
      if (l > max) max = l;
      sum += l;
    }

    cout << endl;
    cout << "Minimum edge length      = " << min << endl;
    cout << "Maximum edge length      = " << max << endl;
    cout << "Average edge length      = " << sum / edgeTable.NumberOfEdges() << endl;
  }

  // Write surface mesh
  return WritePolyData(output_name, mesh, compress, ascii) ? 0 : 1;
}
