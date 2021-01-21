/*
 * Medical Image Registration ToolKit (MIRTK)
 *
 * Copyright 2020 Andreas Schuh
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
#include "mirtk/PointSetIO.h"
#include "mirtk/PointSetUtils.h"
#include "mirtk/Queue.h"

#include "mirtk/Vtk.h"
#include "mirtk/VtkMath.h"

#include "vtkPolyData.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkSmartPointer.h"

#include "vtkCleanPolyData.h"
#include "vtkFillHolesFilter.h"
#include "vtkPolyDataNormals.h"

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
  cout << "  Clean surface mesh geometry and/or topology.\n";
  cout << "\n";
  cout << "  Order of mesh clean up steps, where disabled steps are simply skipped:\n";
  cout << "  1. Remove pointy spikes\n";
  cout << "  2. Merge coincident points\n";
  cout << "  3. Remove redundant faces\n";
  cout << "  4. Remove boundary faces\n";
  cout << "\n";
  cout << "  An example usage is to fix up the geometry and/or topology of a closed cortical\n";
  cout << "  surface mesh obtained using a deformable surface mesh (cf. deform-mesh).\n";
  cout << "\n";
  cout << "Arguments:\n";
  cout << "  input    Input surface mesh.\n";
  cout << "  output   Output surface mesh.\n";
  cout << "\n";
  cout << "Options:\n";
  cout << "  -all\n";
  cout << "      Enable all steps. Individual steps can be disabled again after this option.\n";
  cout << "  -[no]merge-points [on|off]\n";
  cout << "      Merge coincident points. (default: off)\n";
  cout << "  -[no]remove-duplicates [on|off]\n";
  cout << "      Delete duplicate cells. (default: off)\n";
  cout << "  -[no]remove-boundaries [on|off]\n";
  cout << "      Delete cells with at least one edge that is not shared with another cell. (default: off)\n";
  cout << "  -[no]remove-spikes [on|off]\n";
  cout << "      Whether to identify and remove spikes in surface mesh.\n";
  cout << "      This option assumes a closed triangular surface mesh. It uses vtkFillHolesFilter\n";
  cout << "      with maximum hole size 5 is used to fill holes created by removal of spikes\n";
  cout << "      The tip of a spike is identified as a mesh vertex with incoming edges which\n";
  cout << "      all point more or less in the same direction. Starting from this tip, triangles\n";
  cout << "      with smallest angle less than 40 degrees are removed.\n";
  cout << "\n";
  cout << "Output options:\n";
  cout << "  -[no]ascii\n";
  cout << "      Write legacy VTK in ASCII or binary format. (default: input type)\n";
  cout << "  -[no]binary\n";
  cout << "      Write legacy VTK in ASCII or binary format. (default: input type)\n";
  cout << "  -[no]compress\n";
  cout << "      Write XML VTK file with or without compression. (default: on)\n";
  PrintStandardOptions(cout);
  cout << endl;
}

// =============================================================================
// Auxiliary functions
// =============================================================================

// -----------------------------------------------------------------------------
void RemoveUnusedPoinst(vtkPolyData *mesh)
{
  vtkNew<vtkCleanPolyData> cleaner;
  SetVTKInput(cleaner, mesh);
  cleaner->PointMergingOff();
  cleaner->ConvertStripsToPolysOff();
  cleaner->ConvertPolysToLinesOff();
  cleaner->ConvertLinesToPointsOff();
  cleaner->Update();
  mesh->ShallowCopy(cleaner->GetOutput());
}

// -----------------------------------------------------------------------------
bool IsSpikeCell(vtkPolyData *mesh, vtkIdType cellId, const Vector3 &direction, double tol)
{
  Vector3 p, q, u;
  vtkIdType ptId1, ptId2;
  vtkCell *edge;
  double dp;
  auto cell = mesh->GetCell(cellId);
  int num_parallel = 0;
  for (int edgeId = 0; edgeId < cell->GetNumberOfEdges(); ++edgeId) {
    edge = cell->GetEdge(edgeId);
    mirtkAssert(edge->GetNumberOfPoints() == 2, "edges have two points");
    ptId1 = edge->GetPointId(0);
    ptId2 = edge->GetPointId(1);
    mesh->GetPoint(ptId1, p);
    mesh->GetPoint(ptId2, q);
    u = p - q;
    u.Normalize();
    dp = abs(abs(direction.Dot(u)) - 1);
    if (dp <= tol) num_parallel += 1;
  }
  return num_parallel > 1;
}

// -----------------------------------------------------------------------------
int RemoveSpike(vtkPolyData *mesh, vtkIdType ptId, const Vector3 &direction, double tol)
{
  int num_deleted = 0;

  vtkIdType cellId;
  Queue<vtkIdType> active;
  vtkNew<vtkIdList> cellIds, ptIds;

  mesh->GetPointCells(ptId, cellIds.GetPointer());
  for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
    active.push(cellIds->GetId(i));
  }
  while (!active.empty()) {
    cellId = active.front();
    active.pop();
    if (IsSpikeCell(mesh, cellId, direction, tol)) {
      GetCellPoints(mesh, cellId, ptIds.GetPointer());
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        mesh->GetPointCells(ptIds->GetId(i), cellIds.GetPointer());
        for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
          if (cellIds->GetId(i) != cellId) {
            active.push(cellIds->GetId(i));
          }
        }
      }
      mesh->DeleteCell(cellId);
      num_deleted += 1;
    }
  }

  return num_deleted;
}

// -----------------------------------------------------------------------------
// Remove spikes from surface
int RemoveSpikes(vtkPolyData *mesh)
{
  mesh->BuildLinks();
  const EdgeTable edges(mesh);

  int numAdjPts;
  const int *adjPtIds;
  Vector3 p, q, u, v;

  int num_spikes = 0;
  for (vtkIdType ptId = 0; ptId < mesh->GetNumberOfPoints(); ++ptId) {
    edges.GetAdjacentPoints(static_cast<int>(ptId), numAdjPts, adjPtIds);
    if (numAdjPts > 2) {
      bool is_tip_of_spike = true;
      mesh->GetPoint(ptId, p);
      mesh->GetPoint(adjPtIds[0], q);
      u = q - p;
      u.Normalize();
      for (int i = 1; i < numAdjPts; ++i) {
        mesh->GetPoint(adjPtIds[i], q);
        v = q - p;
        v.Normalize();
        if (1 - u.Dot(v) > 0.001) {
          is_tip_of_spike = false;
          break;
        }
      }
      if (is_tip_of_spike) {
        RemoveSpike(mesh, ptId, u, 0.234);
        num_spikes += 1;
      }
    }
  }

  mesh->RemoveDeletedCells();
  RemoveUnusedPoinst(mesh);

  vtkNew<vtkFillHolesFilter> filler;
  SetVTKInput(filler, mesh);
  filler->SetHoleSize(5);
  filler->Update();

  vtkNew<vtkPolyDataNormals> normals;
  SetVTKConnection(normals, filler);
  normals->SplittingOff();
  normals->ConsistencyOn();
  normals->NonManifoldTraversalOn();
  normals->AutoOrientNormalsOn();
  normals->SetComputeCellNormals(mesh->GetCellData()->GetNormals() != nullptr);
  normals->SetComputePointNormals(mesh->GetPointData()->GetNormals() != nullptr);
  normals->Update();

  mesh->ShallowCopy(normals->GetOutput());
  return num_spikes;
}

// -----------------------------------------------------------------------------
void MergePoints(vtkPolyData *mesh)
{
  vtkNew<vtkCleanPolyData> cleaner;
  SetVTKInput(cleaner, mesh);
  cleaner->PointMergingOn();
  cleaner->SetAbsoluteTolerance(1e-12);
  cleaner->ToleranceIsAbsoluteOn();
  cleaner->ConvertStripsToPolysOff();
  cleaner->ConvertPolysToLinesOff();
  cleaner->ConvertLinesToPointsOff();
  cleaner->Update();
  mesh->ShallowCopy(cleaner->GetOutput());
}

// -----------------------------------------------------------------------------
void DeleteBoundaryCells(vtkPolyData *mesh)
{
  mesh->BuildLinks();
  vtkNew<vtkIdList> cellIds, otherCellIds;
  while (mesh->GetNumberOfCells() > 0) {
    const auto edges = BoundaryEdges(mesh);
    if (edges.empty()) {
      break;
    }
    for (const auto &edge : edges) {
      mesh->GetPointCells(edge.first, cellIds.GetPointer());
      mesh->GetPointCells(edge.second, otherCellIds.GetPointer());
      cellIds->IntersectWith(otherCellIds.GetPointer());
      for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
        mesh->DeleteCell(cellIds->GetId(i));
      }
    }
  }
  mesh->RemoveDeletedCells();
  RemoveUnusedPoinst(mesh);
}

// -----------------------------------------------------------------------------
void DeleteRedundantCells(vtkPolyData *mesh)
{
  mesh->DeleteCells();
  mesh->BuildLinks();
  vtkNew<vtkIdList> ptIds1, ptIds2, cellIds;
  for (vtkIdType cellId = 0; cellId < mesh->GetNumberOfCells(); ++cellId) {
    if (mesh->GetCellType(cellId) != VTK_EMPTY_CELL) {
      GetCellPoints(mesh, cellId, ptIds1.GetPointer());
      for (vtkIdType i = 0; i < ptIds1->GetNumberOfIds(); ++i) {
        mesh->GetPointCells(ptIds1->GetId(i), cellIds.GetPointer());
        for (vtkIdType j = 0; j < cellIds->GetNumberOfIds(); ++j) {
          if (cellIds->GetId(j) > cellId && mesh->GetCellType(cellIds->GetId(j)) != VTK_EMPTY_CELL) {
            GetCellPoints(mesh, cellIds->GetId(j), ptIds2.GetPointer());
            if (ptIds1->GetNumberOfIds() == ptIds2->GetNumberOfIds()) {
              ptIds2->IntersectWith(ptIds1.GetPointer());
              if (ptIds1->GetNumberOfIds() == ptIds2->GetNumberOfIds()) {
                mesh->DeleteCell(cellIds->GetId(j));
              }
            }
          }
        }
      }
    }
  }
  mesh->RemoveDeletedCells();
}

// =============================================================================
// Main
// =============================================================================

// -----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  // Parse positional arguments
  EXPECTS_POSARGS(2);

  const char *input_name  = POSARG(1);
  const char *output_name = POSARG(2);

  // Parse optional arguments
  FileOption fopt = FO_Default;
  bool merge_points = false;
  bool remove_boundaries = false;
  bool remove_duplicates = false;
  bool remove_spikes = false;

  for (ALL_OPTIONS) {
    if (OPTION("-all")) {
      merge_points = true;
      remove_boundaries = true;
      remove_duplicates = true;
      remove_spikes = true;
    }
    else HANDLE_BOOLEAN_OPTION("merge-points", merge_points);
    else HANDLE_BOOLEAN_OPTION("remove-duplicates", remove_duplicates);
    else HANDLE_BOOLEAN_OPTION("remove-boundaries", remove_boundaries);
    else HANDLE_BOOLEAN_OPTION("remove-spikes", remove_spikes);
    else HANDLE_POINTSETIO_OPTION(fopt);
    else HANDLE_COMMON_OR_UNKNOWN_OPTION();
  }

  // Load input mesh
  vtkSmartPointer<vtkPolyData> mesh = ReadPolyData(input_name, fopt);

  // Remove pointy spikes from surface mesh
  // MUST be done before point merging!
  if (remove_spikes) {
    RemoveSpikes(mesh);
  }
  // Merge points after other operations
  if (merge_points) {
    MergePoints(mesh);
  }
  // Remove redundant faces (cf. evaluate-surface-mesh)
  if (remove_duplicates) {
    DeleteRedundantCells(mesh);
  }
  // Remove boundary triangles (e.g., may be introduced by point merging)
  // AFTER removal of redundant faces which may introduce boundary cells.
  if (remove_boundaries) {
    DeleteBoundaryCells(mesh);
  }

  // Write surface mesh
  if (!WritePolyData(output_name, mesh, fopt)) {
    FatalError("Failed to write output surface to file " << output_name);
  }

  return 0;
}
