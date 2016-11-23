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

#include "vtkNew.h"
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
  unsigned short ncells;
  vtkIdType      npts, *pts, *cells, cellId;

  while (true) {
    cellId = -1;
    output->GetPointCells(line.back(), ncells, cells);
    for (unsigned short i = 0; i < ncells; ++i) {
      if (output->GetCellType(cells[i]) == VTK_LINE) {
        if (cellId == -1) {
          cellId = cells[i];
        } else {
          cellId = -1;
          break;
        }
      }
    }
    if (cellId == -1) break;
    output->GetCellPoints(cellId, npts, pts);
    output->RemoveCellReference(cellId);
    output->DeleteCell(cellId);
    if (npts == 0) break;
    if (pts[0] == line.back()) {
      for (vtkIdType i = 1; i < npts; ++i) {
        line.push_back(pts[i]);
      }
    } else if (pts[npts-1] == line.back()) {
      for (vtkIdType i = npts-2; i >= 0; --i) {
        line.push_back(pts[i]);
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
    vtkIdType             npts, *pts;

    // Pick each line segment as possible seed for contiguous line
    for (vtkIdType seedId = 0; seedId < _Output->GetNumberOfCells(); ++seedId) {
      // Skip non-line cells and those that are already marked as deleted
      if (_Output->GetCellType(seedId) != VTK_LINE) continue;
      // Start new line in reverse order
      List<vtkIdType> line;
      _Output->GetCellPoints(seedId, npts, pts);
      for (vtkIdType i = 0; i < npts; ++i) {
        line.push_front(pts[i]);
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
