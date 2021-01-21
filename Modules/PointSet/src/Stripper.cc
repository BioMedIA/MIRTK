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

#include "mirtk/Stripper.h"

#include "mirtk/List.h"
#include "mirtk/Vtk.h"

#include "vtkStripper.h"
#include "vtkCellArray.h"


namespace mirtk {


// =============================================================================
// Auxiliares
// =============================================================================

namespace {


// -----------------------------------------------------------------------------
void GrowLine(vtkPolyData *output, List<vtkIdType> &line)
{
  vtkIdType cellId;
  vtkNew<vtkIdList> cellIds, ptIds;
  while (true) {
    cellId = -1;
    output->GetPointCells(line.back(), cellIds.GetPointer());
    for (vtkIdType i = 0; i < cellIds->GetNumberOfIds(); ++i) {
      if (output->GetCellType(cellIds->GetId(i)) == VTK_LINE) {
        if (cellId == -1) {
          cellId = cellIds->GetId(i);
        } else {
          cellId = -1;
          break;
        }
      }
    }
    if (cellId == -1) break;
    GetCellPoints(output, cellId, ptIds.GetPointer());
    output->RemoveCellReference(cellId);
    output->DeleteCell(cellId);
    if (ptIds->GetNumberOfIds() == 0) break;
    if (ptIds->GetId(0) == line.back()) {
      for (vtkIdType i = 1; i < ptIds->GetNumberOfIds(); ++i) {
        line.push_back(ptIds->GetId(i));
      }
    } else if (ptIds->GetId(ptIds->GetNumberOfIds() - 1) == line.back()) {
      for (vtkIdType i = ptIds->GetNumberOfIds() - 2; i >= 0; --i) {
        line.push_back(ptIds->GetId(i));
      }
    } else {
      break;
    }
  }
}


} // namespace

// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void Stripper::CopyAttributes(const Stripper &other)
{
  _StripLines     = other._StripLines;
  _StripTriangles = other._StripTriangles;
}

// -----------------------------------------------------------------------------
Stripper::Stripper()
:
  _StripLines(true),
  _StripTriangles(true)
{
}

// -----------------------------------------------------------------------------
Stripper::Stripper(const Stripper &other)
:
  MeshFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
Stripper &Stripper::operator =(const Stripper &other)
{
  if (this != &other) {
    MeshFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
Stripper::~Stripper()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void Stripper::Execute()
{
  // Strip lines
  if (_StripLines && _Output->GetNumberOfLines() > 0) {
    List<List<vtkIdType>> lines;
    vtkNew<vtkIdList> ptIds;

    // Pick each line segment as possible seed for contiguous line
    for (vtkIdType seedId = 0; seedId < _Output->GetNumberOfCells(); ++seedId) {
      // Skip non-line cells and those that are already marked as deleted
      if (_Output->GetCellType(seedId) != VTK_LINE) continue;
      // Start new line in reverse order
      List<vtkIdType> line;
      GetCellPoints(_Output, seedId, ptIds.GetPointer());
      for (vtkIdType i = 0; i < ptIds->GetNumberOfIds(); ++i) {
        line.push_front(ptIds->GetId(i));
      }
      _Output->RemoveCellReference(seedId);
      _Output->DeleteCell(seedId);
      // Append line at the front
      GrowLine(_Output, line);
      // Reverse line
      line.reverse();
      // Append line at the back
      GrowLine(_Output, line);
      // Add to list of joined lines
      lines.push_back(move(line));
    }

    // Remove lines which are being joined
    _Output->RemoveDeletedCells();

    // Add new joined line (**after** RemoveDeletedCells)
    vtkCellArray * const arr = _Output->GetLines();
    for (const auto &line : lines) {
      arr->InsertNextCell(static_cast<int>(line.size()));
      for (const auto &ptId : line) {
        arr->InsertCellPoint(ptId);
      }
    }
  }

  // Strip triangles
  if (_StripTriangles && _Output->GetNumberOfPolys() > 0) {
    vtkNew<vtkStripper> stripper;
    stripper->SetInputData(_Output);
    stripper->Update();
    _Output = stripper->GetOutput();
  }
}


} // namespace mirtk
