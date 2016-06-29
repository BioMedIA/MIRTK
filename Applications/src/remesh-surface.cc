/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2013-2016 Imperial College London
 * Copyright 2013-2016 Andreas Schuh
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

#include "mirtk/EdgeTable.h"
#include "mirtk/SurfaceRemeshing.h"
#include "mirtk/PointSetIO.h"

#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellArray.h"
#include "vtkCellLocator.h"

using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  SurfaceRemeshing remesher; // with default settings
  cout << endl;
  cout << "Usage: " << name << " <input> <output> [options]" << endl;
  cout << endl;
  cout << "Description:" << endl;
  cout << "  Remeshes a surface mesh such that the resulting mesh has an average" << endl;
  cout << "  edge length within the specified limits. The input surface mesh is" << endl;
  cout << "  triangulated when necessary before the local remeshing passes." << endl;
  cout << endl;
  cout << "Arguments:" << endl;
  cout << "  input    Input surface mesh." << endl;
  cout << "  output   Output surface mesh." << endl;
  cout << endl;
  cout << "Options:" << endl;
  cout << "  -target <file>" << endl;
  cout << "      Find closest cell points on target surface and use these as new point locations." << endl;
  cout << "      (default: use points of input mesh)" << endl;
  cout << "  -edgelength <float>..." << endl;
  cout << "      Average edge length. (default: [0, inf)" << endl;
  cout << "  -min-edgelength <float>..." << endl;
  cout << "      Minimum edge length. (default: 0)" << endl;
  cout << "  -max-edgelength <float>..." << endl;
  cout << "      Maximum edge length. (default: inf)" << endl;
  cout << "  -adaptive-edgelength <name>" << endl;
  cout << "      Name of point data array to use for adapting the edge length range. (default: none)" << endl;
  cout << "  -melting-order <none|area|shortest>" << endl;
  cout << "      Order in which to process cells in melting pass. (default: ";
  switch (remesher.MeltingOrder()) {
    case SurfaceRemeshing::INDEX:         cout << "none"; break;
    case SurfaceRemeshing::AREA:          cout << "area"; break;
    case SurfaceRemeshing::SHORTEST_EDGE: cout << "shortest"; break;
  }
  cout << ")" << endl;
  cout << "  -[no]melt-nodes" << endl;
  cout << "      Whether to allow removal of adjacent nodes with connectivity three" << endl;
  cout << "      during melting pass. (default: " << ToLower(ToString(remesher.MeltNodes())) << ")" << endl;
  cout << "  -[no]melt-triangles" << endl;
  cout << "      Whether to melt triangles when all three edges are too short. (default: "
       << ToLower(ToString(remesher.MeltTriangles())) << ")" << endl;
  cout << "  -[no]invert-long-edges" << endl;
  cout << "      Enable/disable inversion of triangles sharing one too long edge. (default: "
       << ToLower(ToString(remesher.InvertTrianglesSharingOneLongEdge())) << ")" << endl;
  cout << "  -[no]invert-min-height" << endl;
  cout << "      Enable/disable inversion of triangles when it increases the minimum height. (default: "
       << ToLower(ToString(remesher.InvertTrianglesToIncreaseMinHeight())) << ")" << endl;
  cout << "  -noinversion" << endl;
  cout << "      Disable :option:`-invert-long-edges` and :option:`-invert-min-height`." << endl;
  cout << endl;
  cout << "Output options:" << endl;
  cout << "  -write-all" << endl;
  cout << "      Write also intermediate meshes when more than one edge length range" << endl;
  cout << "      was specified. Output file names contain the level number as suffix. (default: off)" << endl;
  cout << "  -[no]ascii" << endl;
  cout << "      Write legacy VTK in ASCII or binary format. (default: input type)" << endl;
  cout << "  -[no]binary" << endl;
  cout << "      Write legacy VTK in ASCII or binary format. (default: input type)" << endl;
  cout << "  -[no]compress" << endl;
  cout << "      Write XML VTK file with or without compression. (default: on)" << endl;
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  SurfaceRemeshing remesher;
  FileOption fopt = FO_Default;

  // Parse positional arguments
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  const string dir  = Directory(output_name);
  const string name = FileName (output_name);
  const string ext  = Extension(output_name);

  // Read input mesh
  vtkSmartPointer<vtkPolyData> mesh = ReadPolyData(input_name, fopt);
  const vtkIdType npoints = mesh->GetNumberOfPoints();
  const vtkIdType ncells  = mesh->GetNumberOfCells();

  // Parse optional arguments
  Array<double> minlength;
  Array<double> maxlength;
  Array<double> minangle;
  Array<double> maxangle;

  const char *target_name = nullptr;
  bool        write_all   = false;

  for (ALL_OPTIONS) {
    if (OPTION("-min-edge-length") || OPTION("-min-edgelength") || OPTION("-minedgelength")) {
      do {
        minlength.push_back(atof(ARGUMENT));
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-max-edge-length") || OPTION("-max-edgelength") || OPTION("-maxedgelength")) {
      do {
        maxlength.push_back(atof(ARGUMENT));
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-edge-length") || OPTION("-edgelength")) {
      if (minlength.size() < maxlength.size()) minlength.resize(maxlength.size(), minlength.back());
      if (minlength.size() > maxlength.size()) maxlength.resize(minlength.size(), maxlength.back());
      do {
        double l = atof(ARGUMENT);
        minlength.push_back(l);
        maxlength.push_back(l);
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-adaptive") || OPTION("-adaptive-edge-length") ||
             OPTION("-adaptive-edgelength") || OPTION("-adaptiveedgelength")) {
      const char *adaptive_name = ARGUMENT;
      vtkDataArray *array = mesh->GetPointData()->GetArray(adaptive_name);
      if (!array) {
        cerr << "Error: Input mesh has no point data array named: " << adaptive_name << endl;
        exit(1);
      }
      remesher.AdaptiveEdgeLengthArray(array);
    }
    else if (OPTION("-min-angle") || OPTION("-minangle")) {
      do {
        minangle.push_back(atof(ARGUMENT));
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-max-angle") || OPTION("-maxangle")) {
      do {
        maxangle.push_back(atof(ARGUMENT));
      } while (HAS_ARGUMENT);
    }
    else if (OPTION("-melting-order") || OPTION("-melt-order") || OPTION("-meltorder")) {
      string arg = ARGUMENT;
      transform(arg.begin(), arg.end(), arg.begin(), ::tolower);
      if (arg == "area" || arg == "smallest area" || arg == "smallestarea") {
        remesher.MeltingOrder(SurfaceRemeshing::AREA);
      } else if (arg == "edge" ||
               arg == "shortest" ||
               arg == "shortest edge" ||
               arg == "shortestedge") {
        remesher.MeltingOrder(SurfaceRemeshing::SHORTEST_EDGE);
      } else if (arg == "index" || arg == "none") {
        remesher.MeltingOrder(SurfaceRemeshing::INDEX);
      } else {
        cerr << "Invalid -meltorder: " << arg << endl;
        exit(1);
      }
    }
    else if (OPTION("-write-all")) write_all = true;
    else if (OPTION("-melt-nodes") || OPTION("-meltnodes")) {
      bool melt_nodes = true;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(melt_nodes);
      remesher.MeltNodes(melt_nodes);
    }
    else if (OPTION("-nomelt-nodes") || OPTION("-nomeltnodes")) {
      remesher.MeltNodesOff();
    }
    else if (OPTION("-melt-triangles") || OPTION("-melttriangles")) {
      bool melt_triangles = true;;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(melt_triangles);
      remesher.MeltTriangles(melt_triangles);
    }
    else if (OPTION("-nomelt-triangles") || OPTION("-nomelttriangles")) {
      remesher.MeltTrianglesOff();
    }
    else if (OPTION("-invert-long-edges")) {
      bool invert = true;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(invert);
      remesher.InvertTrianglesSharingOneLongEdge(invert);
    }
    else if (OPTION("-noinvert-long-edges")) {
      remesher.InvertTrianglesSharingOneLongEdgeOff();
    }
    else if (OPTION("-invert-min-height")) {
      bool invert = true;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(invert);
      remesher.InvertTrianglesToIncreaseMinHeight(invert);
    }
    else if (OPTION("-noinvert-min-height")) {
      remesher.InvertTrianglesToIncreaseMinHeightOff();
    }
    else if (OPTION("-noinversion")) {
      remesher.InvertTrianglesSharingOneLongEdgeOff();
      remesher.InvertTrianglesToIncreaseMinHeightOff();
    }
    else if (OPTION("-target")) target_name = ARGUMENT;
    else HANDLE_POINTSETIO_OPTION(fopt);
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
  SharedPtr<EdgeTable> edgeTable = NewShared<EdgeTable>(mesh);
  vtkIdType nedges = edgeTable->NumberOfEdges();
  if (verbose) {
    vtkIdType euler = mesh->GetNumberOfPoints() - nedges + mesh->GetNumberOfCells();
    cout << "No. of input vertices  = " << npoints << "\n";
    cout << "No. of input triangles = " << ncells  << "\n";
    cout << "No. of input edges     = " << nedges  << "\n";
    cout << "Euler characteristic   = " << euler << "\n";
    cout.flush();
  }
  for (size_t i = 0; i < nlevels; ++i) {
    if (verbose) {
      cout << "\nRemeshing surface with edge length range [" << minlength[i] << ", " << maxlength[i] << "]";
      if (minangle[i] < 180.0) cout << " and feature angle in [" << minangle[i] << ", " << maxangle[i] << "] degrees";
      cout << endl;
    }
    remesher.Input(mesh);
    remesher.EdgeTable(edgeTable);
    remesher.MinEdgeLength  (minlength[i]);
    remesher.MaxEdgeLength  (maxlength[i]);
    remesher.MinFeatureAngle(minangle [i]);
    remesher.MaxFeatureAngle(maxangle [i]);
    remesher.Run();
    mesh = remesher.Output();
    edgeTable->Initialize(mesh);
    nedges = edgeTable->NumberOfEdges();
    if (verbose > 1) {
      if (debug_time) cout << endl;
      cout << "\n";
      cout << "  No. of node-meltings     = " << remesher.NumberOfMeltedNodes()  << "\n";
      cout << "  No. of edge-meltings     = " << remesher.NumberOfMeltedEdges()  << "\n";
      cout << "  No. of triangle-meltings = " << remesher.NumberOfMeltedCells()  << "\n";
      cout << "  No. of inversions        = " << remesher.NumberOfInversions()   << "\n";
      cout << "  No. of bisections        = " << remesher.NumberOfBisections()   << "\n";
      cout << "  No. of trisections       = " << remesher.NumberOfTrisections()  << "\n";
      cout << "  No. of quadsections      = " << remesher.NumberOfQuadsections() << "\n";
    }
    if (verbose) {
      vtkIdType euler = mesh->GetNumberOfPoints() - nedges + mesh->GetNumberOfCells();
      cout << "\n";
      cout << "  No. of output vertices   = " << mesh->GetNumberOfPoints()
           << " (" << (100.0 * mesh->GetNumberOfPoints() / npoints) << "%)" << "\n";
      cout << "  No. of output triangles  = " << mesh->GetNumberOfCells()
           << " (" << (100.0 * mesh->GetNumberOfCells() / ncells) << "%)" << "\n";
      cout << "  No. of output edges      = " << nedges
           << " (" << (100.0 * nedges / nedges) << "%)" << "\n";
      cout << "  Euler characteristic     = " << euler << "\n";
      cout.flush();
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
    EdgeIterator it(*edgeTable);
    for (it.InitTraversal(); it.GetNextEdge(ptId1, ptId2) != -1;) {
      mesh->GetPoint(ptId1, p1);
      mesh->GetPoint(ptId2, p2);
      l = sqrt(vtkMath::Distance2BetweenPoints(p1, p2));
      if (l < min) min = l;
      if (l > max) max = l;
      sum += l;
    }
    cout << "\n";
    cout << "Minimum edge length = " << min << "\n";
    cout << "Maximum edge length = " << max << "\n";
    cout << "Average edge length = " << sum / nedges << "\n";
    cout.flush();
  }

  // Write surface mesh
  if (!WritePolyData(output_name, mesh, fopt)) {
    FatalError("Failed to write output surface to file " << output_name);
  }

  return 0;
}
