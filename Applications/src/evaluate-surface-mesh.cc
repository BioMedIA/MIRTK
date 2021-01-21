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
#include "mirtk/DataSelection.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/SurfaceCollisions.h"

#include "mirtk/GenericImage.h"
#include "mirtk/InterpolateImageFunction.h"

#include "mirtk/Vtk.h"
#include "vtkSmartPointer.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkCellData.h"
#include "vtkIdTypeArray.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkExtractCells.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkPointLocator.h"


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

// -----------------------------------------------------------------------------
bool AddCellData(vtkPolyData *input, vtkPolyData *surface, vtkDataArray *arr)
{
  if (arr == nullptr) return false;
  vtkDataArray * const ids = surface->GetCellData()->GetArray("OriginalIds");
  if (ids == nullptr) {
    if (input == surface) {
      vtkCellData * const cd = input->GetCellData();
      for (int i = 0; i < cd->GetNumberOfArrays(); ++i) {
        if (cd->GetArray(i) == arr) return true;
      }
      if (arr->GetName()) cd->RemoveArray(arr->GetName());
      cd->AddArray(arr);
      return true;
    }
    return false;
  }
  vtkSmartPointer<vtkDataArray> out;
  out.TakeReference(arr->NewInstance());
  out->SetNumberOfComponents(arr->GetNumberOfComponents());
  out->SetNumberOfTuples(input->GetNumberOfCells());
  if (arr->GetName()) out->SetName(arr->GetName());
  for (int j = 0; j < arr->GetNumberOfComponents(); ++j) {
    out->SetComponentName(j, arr->GetComponentName(j));
  }
  out->FillComponent(0, 0.);
  for (vtkIdType cellId = 0, origId; cellId < surface->GetNumberOfCells(); ++cellId) {
    origId = static_cast<vtkIdType>(ids->GetComponent(cellId, 0));
    out->SetTuple(origId, arr->GetTuple(cellId));
  }
  input->GetCellData()->AddArray(out);
  return true;
}

// -----------------------------------------------------------------------------
bool AddCellData(vtkPolyData *input, vtkPolyData *surface, const char *name)
{
  return AddCellData(input, surface, surface->GetCellData()->GetArray(name));
}

// -----------------------------------------------------------------------------
/// Number of redundant cells, i.e., cells with same shared points
int NumberOfRedundantCells(vtkPolyData *dataset, const char *mask_name = nullptr)
{
  vtkSmartPointer<vtkDataArray> mask;
  if (mask_name) {
    mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, dataset->GetNumberOfCells(), 1, mask_name);
    mask->SetName(mask_name);
    mask->FillComponent(0, 0.);
    dataset->GetCellData()->RemoveArray(mask->GetName());
    dataset->GetCellData()->AddArray(mask);
  }
  int n = 0;
  vtkSmartPointer<vtkPolyData> surface;
  surface.TakeReference(dataset->NewInstance());
  surface->ShallowCopy(dataset);
  surface->DeleteCells();
  surface->BuildLinks();
  vtkNew<vtkIdList> ptIds1, ptIds2, cellIds;
  for (vtkIdType cellId = 0; cellId < surface->GetNumberOfCells(); ++cellId) {
    if (surface->GetCellType(cellId) != VTK_EMPTY_CELL) {
      GetCellPoints(surface, cellId, ptIds1.GetPointer());
      for (vtkIdType i = 0; i < ptIds1->GetNumberOfIds(); ++i) {
        surface->GetPointCells(ptIds1->GetId(i), cellIds.GetPointer());
        for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
          if (cellIds->GetId(j) > cellId && surface->GetCellType(cellIds->GetId(j)) != VTK_EMPTY_CELL) {
            GetCellPoints(surface, cellIds->GetId(j), ptIds2.GetPointer());
            if (ptIds1->GetNumberOfIds() == ptIds2->GetNumberOfIds()) {
              ptIds2->IntersectWith(ptIds1.GetPointer());
              if (ptIds1->GetNumberOfIds() == ptIds2->GetNumberOfIds()) {
                surface->DeleteCell(cellIds->GetId(j));
                if (mask) {
                  mask->SetComponent(cellId,            0, 1.);
                  mask->SetComponent(cellIds->GetId(j), 0, 1.);
                }
                ++n;
              }
            }
          }
        }
      }
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
  verbose = (output_name ? 0 : 1);

  const char *redundant_cells_mask = nullptr;
  const char *boundary_point_mask  = nullptr;
  const char *boundary_cell_mask   = nullptr;

  if (output_name) {
    redundant_cells_mask = "DuplicatedMask";
    boundary_point_mask  = "BoundaryMask";
    boundary_cell_mask   = "BoundaryMask";
  }

  SharedPtr<LogicalOp> selector(new LogicalAnd());
  const char *select_name = nullptr;
  int         select_comp = 0;

  Array<MeshProperty> measures;
  double min_frontface_dist  = 1e-2;
  double min_backface_dist   = 1e-2;
  double max_collision_angle = 20.0;
  bool adjacent_collision_test = true;
  bool fast_collision_test = false;

  double value;
  for (ALL_OPTIONS) {
    if (OPTION("-select") || OPTION("-where")) {
      select_name = ARGUMENT;
      if (HAS_ARGUMENT) PARSE_ARGUMENT(select_comp);
      else select_comp = 0;
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
    else if (OPTION("-attributes") || OPTION("-attr")) {
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
      max_collision_angle = 20.0;
      if (HAS_ARGUMENT) {
        PARSE_ARGUMENT(min_frontface_dist);
        if (HAS_ARGUMENT) {
          PARSE_ARGUMENT(min_backface_dist);
          if (HAS_ARGUMENT) PARSE_ARGUMENT(max_collision_angle);
        } else {
          min_backface_dist = min_frontface_dist;
        }
      } else {
        min_frontface_dist = min_backface_dist = 1e-2;
      }
    }
    else HANDLE_BOOLEAN_OPTION("adjacent-collision-test", adjacent_collision_test);
    else HANDLE_BOOLEAN_OPTION("fast-collision-test", fast_collision_test);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }
  if (measures.empty()) {
    measures.push_back(MESH_Attributes);
    if (!output_name || verbose > 0) {
      measures.push_back(MESH_EulerCharacteristic);
      measures.push_back(MESH_Genus);
    }
  } else if (output_name) {
    for (auto m : measures) {
      if (m == MESH_EulerCharacteristic || m == MESH_Genus) {
        verbose = 1;
        break;
      }
    }
  }

  // Group measures into sections
  const char *section = "";
  Sort(measures);

  // Print file name in verbose mode or when more than just a scalar measure is printed
  if (verbose > 1 || (verbose > 0 && measures.size() > 1)) {
    cout << "Surface mesh file = " << input_name << endl;
  }

  // Read input surface mesh
  vtkSmartPointer<vtkPolyData> input   = ReadPolyData(input_name);
  vtkSmartPointer<vtkPolyData> surface = input;
  input->BuildLinks();

  // Extract sub-mesh to evaluate
  if (selector->NumberOfCriteria() > 0) {
    vtkDataArray *scalars = nullptr;
    if (select_name) {
      scalars = GetArrayByCaseInsensitiveName(input->GetCellData(), select_name);
      if (scalars == nullptr) {
        FatalError("Input surface has no cell data array named: " << select_name);
      }
    } else {
      FatalError("Cell selection options require the -where option.");
    }
    if (select_comp < 0 || select_comp >= scalars->GetNumberOfComponents()) {
      FatalError("Cell data array " << select_name << " has only " << scalars->GetNumberOfComponents() << " component(s)");
    }
    Array<double> values(scalars->GetNumberOfTuples());
    for (vtkIdType cellId = 0; cellId < scalars->GetNumberOfTuples(); ++cellId) {
      values[cellId] = scalars->GetComponent(cellId, select_comp);
    }
    const auto selection = selector->Evaluate(values);
    surface.TakeReference(input->NewInstance());
    surface->SetPoints(input->GetPoints());
    surface->Allocate(input->GetNumberOfCells());
    vtkIdType cellId;
    vtkNew<vtkIdList> ptIds;
    vtkSmartPointer<vtkIdTypeArray> origCellIds;
    origCellIds = vtkSmartPointer<vtkIdTypeArray>::New();
    origCellIds->SetName("OriginalIds");
    origCellIds->SetNumberOfComponents(1);
    origCellIds->SetNumberOfTuples(selection.size());
    surface->GetCellData()->AddArray(origCellIds);
    for (auto origCellId : selection) {
      GetCellPoints(input, origCellId, ptIds.GetPointer());
      cellId = surface->InsertNextCell(input->GetCellType(origCellId), ptIds.GetPointer());
      origCellIds->SetValue(cellId, origCellId);
    }
    surface->BuildLinks();
    surface->Squeeze();
  }

  const EdgeTable edgeTable(surface);

  for (auto measure : measures) {
    switch (measure) {

      // -----------------------------------------------------------------------
      // Print surface mesh attributes
      case MESH_Attributes: {
        section = "";

        // Components
        if (verbose > 0) {
          cout << "\nComponents:\n";
          cout.flush();
        }
        vtkSmartPointer<vtkIdTypeArray> components;
        {
          vtkNew<vtkPolyDataConnectivityFilter> connectivity;
          connectivity->SetExtractionModeToAllRegions();
          connectivity->ScalarConnectivityOff();
          connectivity->SetColorRegions(output_name != nullptr);
          connectivity->SetColorRegions(1);
          SetVTKInput(connectivity, surface);
          connectivity->Update();
          components = connectivity->GetRegionSizes();
          if (output_name) {
            int compId;
            vtkPolyData * const output = connectivity->GetOutput();
            const char * const COMPONENT_ID = "ComponentId";
            // Re-assign components IDs from largest to smallest
            vtkSmartPointer<vtkDataArray> regionIds, compIds, pointIds, cellIds;
            const int ncomp = static_cast<int>(connectivity->GetNumberOfExtractedRegions());
            regionIds = output->GetPointData()->GetArray("RegionId");
            compIds   = NewVtkDataArray(ncomp < 256 ? VTK_UNSIGNED_CHAR :
                                          (ncomp < 65535 ? VTK_UNSIGNED_SHORT : VTK_INT),
                                        regionIds->GetNumberOfTuples(), 1, COMPONENT_ID);
            Array<int> compPts(ncomp);
            for (int i = 0; i < ncomp; ++i) {
              compPts[i] = static_cast<int>(components->GetComponent(i, 0));
            }
            compId = 0;
            for (auto i : DecreasingOrder(compPts)) {
              ++compId;
              for (vtkIdType ptId = 0; ptId < regionIds->GetNumberOfTuples(); ++ptId) {
                if (static_cast<int>(regionIds->GetComponent(ptId, 0)) == i) {
                  compIds->SetComponent(ptId, 0, static_cast<double>(compId));
                }
              }
            }
            // Add surface component IDs point data
            Point p, q;
            pointIds.TakeReference(compIds->NewInstance());
            pointIds->SetName(compIds->GetName());
            pointIds->SetNumberOfComponents(1);
            pointIds->SetNumberOfTuples(input->GetNumberOfPoints());
            pointIds->FillComponent(0, 0.);
            vtkNew<vtkPointLocator> locator;
            locator->SetDataSet(connectivity->GetOutput());
            locator->BuildLocator();
            for (vtkIdType ptId = 0, tupleId; ptId < input->GetNumberOfPoints(); ++ptId) {
              input->GetPoint(ptId, p);
              tupleId = locator->FindClosestPoint(p);
              connectivity->GetOutput()->GetPoint(tupleId, q);
              if (p.SquaredDistance(q) < 1e-12) {
                pointIds->SetComponent(ptId, 0, compIds->GetComponent(tupleId, 0));
              }
            }
            input->GetPointData()->RemoveArray(pointIds->GetName());
            input->GetPointData()->AddArray(pointIds);
            // Add surface component IDs cell data
            vtkNew<vtkIdList> ptIds;
            output->BuildLinks();
            cellIds = NewVtkDataArray(ncomp < 256 ? VTK_UNSIGNED_CHAR :
                                        (ncomp < 65535 ? VTK_UNSIGNED_SHORT : VTK_INT),
                                      surface->GetNumberOfCells(), 1, COMPONENT_ID);
            for (vtkIdType cellId = 0; cellId < output->GetNumberOfCells(); ++cellId) {
              GetCellPoints(output, cellId, ptIds.GetPointer());
              cellIds->SetComponent(cellId, 0, (ptIds->GetNumberOfIds() == 0 ? 0. : compIds->GetComponent(ptIds->GetId(0), 0)));
            }
            AddCellData(input, surface, cellIds);
          }
        }
        if (verbose > 0) {
          cout << "  No. of components    = " << components->GetNumberOfTuples() << "\n";
          if (verbose > 1) {
            for (vtkIdType i = 0; i < components->GetNumberOfTuples(); ++i) {
              cout << "  Component " << i+1 << " has " << components->GetComponent(i, 0) << " faces\n";
            }
          }
        }

        // Nodes
        if (verbose > 0) {
          cout << "\nNodes:\n";
          int npoints = NumberOfPoints(surface);
          cout << "  No. of points        = " << surface->GetNumberOfPoints() << "\n";
          cout << "  No. of used points   = " << npoints << "\n";
          cout << "  No. of unused points = " << surface->GetNumberOfPoints() - npoints << "\n";
        }

        // Edges
        if (verbose > 0) {
          double l_min, l_max;
          GetMinMaxEdgeLength(surface->GetPoints(), edgeTable, l_min, l_max);
          cout << "\nEdges:\n";
          cout << "  No. of edges         = " << edgeTable.NumberOfEdges() << "\n";
          cout << "  Average edge length  = " << AverageEdgeLength(surface->GetPoints(), edgeTable) << "\n";
          cout << "  Minimum edge length  = " << l_min << "\n";
          cout << "  Maximum edge length  = " << l_max << "\n";
        }

        // Faces
        int ndup = NumberOfRedundantCells(surface, redundant_cells_mask);
        if (verbose > 0) {
          cout << "\nFaces:\n";
          int nfaces = static_cast<int>(surface->GetNumberOfCells());
          int nempty = CountCellsOfType(surface, VTK_EMPTY_CELL);
          int ntri   = CountCellsOfType(surface, VTK_TRIANGLE);
          int nquad  = CountCellsOfType(surface, VTK_QUAD);
          int nmisc  = surface->GetNumberOfCells() - ntri - nquad - nempty;
          cout << "  No. of faces           = " << nfaces << "\n";
          cout << "  No. of triangles       = " << ntri   << "\n";
          cout << "  No. of quadrangles     = " << nquad  << "\n";
          cout << "  No. of other faces     = " << nmisc  << "\n";
          cout << "  No. of empty faces     = " << nempty << "\n";
          cout << "  No. of redundant faces = " << ndup   << "\n";
          cout << "  Is triangular mesh     = " << (nfaces == ntri  + nempty ? "yes" : "no") << "\n";
          cout << "  Is quadrangular mesh   = " << (nfaces == nquad + nempty ? "yes" : "no") << "\n";
        }
        if (output_name) {
          AddCellData(input, surface, redundant_cells_mask);
        }

        // Boundaries
        const UnorderedSet<int> boundaryPtIds = BoundaryPoints(surface, &edgeTable);
        const EdgeList          boundaryEdges = BoundaryEdges(surface, edgeTable);
        if (verbose > 0) {
          cout << "\nBoundaries:\n";
          cout << "  No. of boundary segments = " << NumberOfBoundarySegments(surface) << "\n";
          cout << "  No. of boundary points   = " << boundaryPtIds.size() << "\n";
          cout << "  No. of boundary edges    = " << boundaryEdges.size() << "\n";
        }
        if (output_name) {
          if (boundary_point_mask) {
            vtkSmartPointer<vtkDataArray> mask;
            mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, surface->GetNumberOfPoints(), 1, boundary_point_mask);
            mask->FillComponent(0, 0.);
            for (auto ptId : boundaryPtIds) mask->SetComponent(ptId, 0, 1.);
            input->GetPointData()->RemoveArray(mask->GetName());
            input->GetPointData()->AddArray(mask);
          }
          if (boundary_cell_mask) {
            vtkSmartPointer<vtkDataArray> mask;
            mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, surface->GetNumberOfCells(), 1, boundary_cell_mask);
            mask->FillComponent(0, 0.);
            vtkNew<vtkIdList> cellIds1, cellIds2;
            for (auto edge : boundaryEdges) {
              surface->GetPointCells(edge.first,  cellIds1.GetPointer());
              surface->GetPointCells(edge.second, cellIds2.GetPointer());
              cellIds1->IntersectWith(cellIds2.GetPointer());
              for (vtkIdType i = 0; i < cellIds1->GetNumberOfIds(); ++i) {
                mask->SetComponent(cellIds1->GetId(i), 0, 1.);
              }
            }
            AddCellData(input, surface, mask);
          }
        }
      } break;

      // -----------------------------------------------------------------------
      // Topology
      case MESH_EulerCharacteristic: {
        const auto chi = EulerCharacteristic(surface, edgeTable);
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
        const auto genus = Genus(surface, edgeTable);
        if (verbose <= 0 && measures.size() == 1) {
          cout << genus << "\n";
        } else {
          const char * const section_name = "Topology";
          if (strcmp(section, section_name) != 0) {
            section = section_name;
            cout << "\n" << section << ":\n";
          }
          cout << "  Genus                = " << genus << "\n";
          if (BoundaryPoints(surface, &edgeTable).empty()) {
            const auto chi = EulerCharacteristic(surface, edgeTable);
            cout << "  Demigenus            = " << 2 - chi << "\n";
          } else {
            cout << "  Demigenus            = NA\n";
          }
        }
      } break;

      // -----------------------------------------------------------------------
      // Collisions and self-intersections
      case MESH_Collisions: {
        if (verbose > 0) {
          const char * const section_name = "Collisions";
          if (strcmp(section, section_name) != 0) {
            section = section_name;
            cout << "\n" << section << ":\n";
          }
        }

        SurfaceCollisions collisions;
        collisions.Input(surface);
        collisions.AdjacentCollisionTest(adjacent_collision_test);
        collisions.FastCollisionTest(fast_collision_test);
        collisions.MinFrontfaceDistance(min_frontface_dist);
        collisions.MinBackfaceDistance(min_backface_dist);
        collisions.MaxAngle(max_collision_angle);
        collisions.AdjacentIntersectionTestOn();
        collisions.NonAdjacentIntersectionTestOn();
        collisions.FrontfaceCollisionTestOn();
        collisions.BackfaceCollisionTestOn();
        collisions.StoreIntersectionDetailsOn();
        collisions.StoreCollisionDetailsOn();
        collisions.Run();

        if (output_name) {
          vtkNew<vtkIdList> cellIds;
          vtkSmartPointer<vtkDataArray> mask;
          mask = NewVtkDataArray(VTK_UNSIGNED_CHAR, surface->GetNumberOfPoints(), 1, "CollisionMask");
          for (vtkIdType ptId = 0; ptId < surface->GetNumberOfPoints(); ++ptId) {
            mask->SetComponent(ptId, 0, 0.);
            surface->GetPointCells(ptId, cellIds.GetPointer());
            for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
              if (collisions.GetCollisionType(cellIds->GetId(i)) != SurfaceCollisions::NoCollision) {
                mask->SetComponent(ptId, 0, 1.);
                break;
              }
            }
          }
          input->GetPointData()->AddArray(mask);
          AddCellData(input, surface, collisions.GetCollisionTypeArray());
        }
        if (verbose > 0) {
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
                  if (verbose > 2) {
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
                if (verbose > 1) {
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
          if ((verbose > 1 && nintersections > 0) || (verbose > 2 && ncollisions > 0)) {
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
        }
      } break;
    }
    cout.flush();
  }

  if (output_name && !WritePolyData(output_name, input)) {
    FatalError("Failed to write surface mesh to file " << output_name);
  }

  return 0;
}
