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

#include "mirtk/Algorithm.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/SurfaceCollisions.h"

#include "mirtk/GenericImage.h"
#include "mirtk/InterpolateImageFunction.h"

#include "mirtk/Vtk.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkIdTypeArray.h"
#include "vtkPolyDataConnectivityFilter.h"


using namespace mirtk;


// =============================================================================
// Help
// =============================================================================

// -----------------------------------------------------------------------------
void PrintHelp(const char *name)
{
  cout << "\n";
  cout << "Usage: " << name << " <surface> [options]\n";
  cout << "\n";
  cout << "Description:\n";
  cout << "  Prints surface mesh quality measures and topology information.\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  surface   Surface mesh file.\n";
  cout << "\n";
  cout << "Options:\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliaries
// =============================================================================

// -----------------------------------------------------------------------------
/// Enumeration of mesh quality measures
enum MeshProperty
{
  // General
  MESH_Attributes,
  // Topology
  MESH_EulerCharacteristic,
  MESH_Genus,
  // Self-intersections
  MESH_Collisions
};

// -----------------------------------------------------------------------------
/// Number of cells of a given type
int CountCellsOfType(vtkDataSet *dataset, int type)
{
  int n = 0;
  for (vtkIdType cellId = 0; cellId < dataset->GetNumberOfCells(); ++cellId) {
    if (dataset->GetCellType(cellId) == type) {
      ++n;
    }
  }
  return n;
}

// =============================================================================
// WM/cGM interface, i.e., cortical white surface mesh
// =============================================================================

// =============================================================================
// cGM/CSF interface, i.e., cortical pial surface mesh
// =============================================================================

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  REQUIRES_POSARGS(1);
  if (NUM_POSARGS > 2) {
    FatalError("Too many positional arguments!");
  }

  const char *input_name  = POSARG(1);
  const char *output_name = (NUM_POSARGS == 2 ? POSARG(2) : nullptr);

  vtkSmartPointer<vtkPolyData> surface = ReadPolyData(input_name);
  surface->BuildLinks();
  EdgeTable edgeTable(surface);

  Array<MeshProperty> measures;
  double              min_frontface_dist = 1e-2;
  double              min_backface_dist  = 1e-2;

  for (ALL_OPTIONS) {
    if (OPTION("-attributes") || OPTION("-attr")) {
      measures.push_back(MESH_Attributes);
    }
    else if (OPTION("-genus")) {
      measures.push_back(MESH_Genus);
    }
    else if (OPTION("-chi") || OPTION("-euler") ||
             OPTION("-euler-characteristic") || OPTION("-euler-number")) {
      measures.push_back(MESH_EulerCharacteristic);
    }
    else if (OPTION("-topology") || OPTION("-topo")) {
      measures.push_back(MESH_EulerCharacteristic);
      measures.push_back(MESH_Genus);
    }
    else if (OPTION("-collisions") || OPTION("-self-intersections") || OPTION("-coll")) {
      measures.push_back(MESH_Collisions);
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(min_frontface_dist);
        if (HAS_ARGUMENT) PARSE_ARGUMENT(min_backface_dist);
        else min_backface_dist = min_frontface_dist;
      } else {
        min_frontface_dist = min_backface_dist = 1e-2;
      }
    }
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }
  if (measures.empty()) {
    measures.push_back(MESH_Attributes);
    measures.push_back(MESH_EulerCharacteristic);
    measures.push_back(MESH_Genus);
  }

  // Group measures into sections
  const char *section = "";
  Sort(measures);

  // Print file name in verbose mode or when more than just a scalar measure is printed
  if (verbose > 0 || measures.size() > 1) {
    cout << "Surface mesh file = " << input_name << endl;
  }

  for (auto measure : measures) {
    switch (measure) {

      // -----------------------------------------------------------------------
      // Print surface mesh attributes
      case MESH_Attributes: {
        section = "";

        // Components
        vtkSmartPointer<vtkIdTypeArray> components;
        {
          vtkNew<vtkPolyDataConnectivityFilter> connectivity;
          connectivity->SetExtractionModeToAllRegions();
          connectivity->ScalarConnectivityOff();
          SetVTKInput(connectivity, surface);
          connectivity->Update();
          components = connectivity->GetRegionSizes();
        }
        cout << "\nComponents:\n";
        cout << "  No. of components    = " << components->GetNumberOfTuples() << "\n";
        if (verbose > 0) {
          for (vtkIdType i = 0; i < components->GetNumberOfTuples(); ++i) {
            cout << "  Component " << i+1 << " has " << components->GetComponent(i, 0) << " faces\n";
          }
        }

        // Nodes
        cout << "\nNodes:\n";
        int npoints = NumberOfPoints(surface);
        cout << "  No. of points        = " << surface->GetNumberOfPoints() << "\n";
        cout << "  No. of used points   = " << npoints << "\n";
        cout << "  No. of unused points = " << surface->GetNumberOfPoints() - npoints << "\n";

        // Faces
        cout << "\nFaces:\n";
        int nfaces = static_cast<int>(surface->GetNumberOfCells());
        int nempty = CountCellsOfType(surface, VTK_EMPTY_CELL);
        int ntri   = CountCellsOfType(surface, VTK_TRIANGLE);
        int nquad  = CountCellsOfType(surface, VTK_QUAD);
        int nmisc  = surface->GetNumberOfCells() - ntri - nquad - nempty;
        cout << "  No. of faces         = " << nfaces << "\n";
        cout << "  No. of triangles     = " << ntri   << "\n";
        cout << "  No. of quadrangles   = " << nquad  << "\n";
        cout << "  No. of other faces   = " << nmisc  << "\n";
        cout << "  No. of empty faces   = " << nempty << "\n";
        cout << "  Is triangular mesh   = " << (nfaces == ntri  + nempty ? "yes" : "no") << "\n";
        cout << "  Is quadrangular mesh = " << (nfaces == nquad + nempty ? "yes" : "no") << "\n";

        // Edges
        double l_min, l_max;
        GetMinMaxEdgeLength(surface->GetPoints(), edgeTable, l_min, l_max);
        cout << "\nEdges:\n";
        cout << "  No. of edges         = " << edgeTable.NumberOfEdges() << "\n";
        cout << "  Average edge length  = " << AverageEdgeLength(surface->GetPoints(), edgeTable) << "\n";
        cout << "  Minimum edge length  = " << l_min << "\n";
        cout << "  Maximum edge length  = " << l_max << "\n";

        // Boundaries
        cout << "\nBoundaries:\n";
        cout << "  No. of boundary segments = " << NumberOfBoundarySegments(surface) << "\n";
      } break;

      // -----------------------------------------------------------------------
      // Topology
      case MESH_EulerCharacteristic: {
        int chi = EulerCharacteristic(surface, edgeTable);
        if (verbose <= 0 && measures.size() == 1) {
          cout << chi << "\n";
        } else {
          const char * const section_name = "Topology";
          if (strcmp(section, section_name) != 0) {
            section = section_name;
            cout << "\n" << section << ":\n";
          }
          cout << "  Euler characteristic = " << chi << "\n";
        }
      } break;
      case MESH_Genus: {
        double genus = Genus(surface, edgeTable);
        if (verbose <= 0 && measures.size() == 1) {
          cout << genus << "\n";
        } else {
          const char * const section_name = "Topology";
          if (strcmp(section, section_name) != 0) {
            section = section_name;
            cout << "\n" << section << ":\n";
          }
          cout << "  Genus                = " << genus << "\n";
        }
      } break;

      // -----------------------------------------------------------------------
      // Collisions and self-intersections
      case MESH_Collisions: {
        const char * const section_name = "Collisions";
        if (strcmp(section, section_name) != 0) {
          section = section_name;
          cout << "\n" << section << ":\n";
        }

        SurfaceCollisions collisions;
        collisions.Input(surface);
        collisions.MinFrontfaceDistance(min_frontface_dist);
        collisions.MinBackfaceDistance(min_backface_dist);
        collisions.MaxAngle(20.0);
        collisions.AdjacentIntersectionTestOn();
        collisions.NonAdjacentIntersectionTestOn();
        collisions.FrontfaceCollisionTestOn();
        collisions.BackfaceCollisionTestOn();
        collisions.StoreIntersectionDetailsOn();
        collisions.StoreCollisionDetailsOn();
        collisions.Run();
        surface = collisions.Output();

        int    nintersections     = 0;
        int    ncollisions        = 0;
        int    ncollisions_front  = 0;
        int    ncollisions_back   = 0;
        double min_distance       = numeric_limits<double>::max();
        double min_distance_front = numeric_limits<double>::max();
        double min_distance_back  = numeric_limits<double>::max();
        double intersection_area  = 0.;
        double collision_area     = 0.;
        double total_area         = 0.;

        for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
          auto type = collisions.GetCollisionType(cellId);
          switch (type) {
            case SurfaceCollisions::Ambiguous:
            case SurfaceCollisions::Intersection:
            case SurfaceCollisions::SelfIntersection:
            case SurfaceCollisions::AdjacentIntersection: {
              ++nintersections;
            } break;
            case SurfaceCollisions::Collision: {
              ++ncollisions;
              ++ncollisions_front;
              ++ncollisions_back;
            } break;
            case SurfaceCollisions::FrontfaceCollision: {
              ++ncollisions;
              ++ncollisions_front;
            } break;
            case SurfaceCollisions::BackfaceCollision: {
              ++ncollisions;
              ++ncollisions_back;
            } break;
            case SurfaceCollisions::NoCollision: {
            } break;
          }
          double area = ComputeArea(surface->GetCell(cellId));
          if (type != SurfaceCollisions::NoCollision) {
            if (collisions.IsCollision(type)) {
              collision_area += area;
              for (const auto &info : collisions.Collisions(cellId)) {
                if (verbose > 1) {
                  cout << "  Triangle " << cellId << " collides with triangle " << info._CellId << "\n";
                }
                if (info._Type == SurfaceCollisions::Collision) {
                  if (info._Distance < min_distance_front) min_distance_front = info._Distance;
                  if (info._Distance < min_distance_back)  min_distance_back  = info._Distance;
                } else if (info._Type == SurfaceCollisions::FrontfaceCollision) {
                  if (info._Distance < min_distance_front) min_distance_front = info._Distance;
                } else if (info._Type == SurfaceCollisions::BackfaceCollision) {
                  if (info._Distance < min_distance_back) min_distance_back = info._Distance;
                }
                if (info._Distance < min_distance) min_distance = info._Distance;
              }
            }
            if (collisions.IsIntersection(type)) {
              intersection_area += area;
              if (verbose > 0) {
                for (const auto &info : collisions.Intersections(cellId)) {
                  cout << "  Triangle " << cellId << " intersects with triangle " << info._CellId << "\n";
                }
              }
            }
          }
          total_area += area;
        }
        if (ncollisions       == 0) min_distance       = mirtk::nan;
        if (ncollisions_front == 0) min_distance_front = mirtk::nan;
        if (ncollisions_back  == 0) min_distance_back  = mirtk::nan;
        if ((verbose > 0 && nintersections > 0) || (verbose > 1 && ncollisions > 0)) {
          cout << "\n";
        }
        cout << "  Total no. of checked faces   = " << surface->GetNumberOfCells() << "\n";
        cout << "  No. of self-intersections    = " << nintersections << "\n";
        cout << "  No. of near-miss collisions  = " << ncollisions    << "\n";
        cout << "  No. of front-face collisions = " << ncollisions_front << "\n";
        cout << "  No. of back-face collisions  = " << ncollisions_back << "\n";
        cout << "\n";
        cout << "  Minimum inter-face distance  = " << min_distance << "\n";
        cout << "  Minimum front-face distance  = " << min_distance_front << "\n";
        cout << "  Minimum back-face distance   = " << min_distance_back << "\n";
        cout << "\n";
        cout << "  Area of intersecting faces   = " << 100. * intersection_area / total_area << "%\n";
        cout << "  Area of collided faces       = " << 100. * collision_area / total_area << "%\n";
      } break;
    }
    cout.flush();
  }

  if (output_name) {
    if (!WritePolyData(output_name, surface)) {
      FatalError("Failed to write surface mesh to file " << output_name);
    }
  }

  return 0;
}
