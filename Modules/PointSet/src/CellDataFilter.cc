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

#include "mirtk/CellDataFilter.h"

#include "mirtk/PointSetUtils.h"

#include "vtkNew.h"
#include "vtkGenericCell.h"


namespace mirtk {


// =============================================================================
// Construction/destruction
// =============================================================================

// -----------------------------------------------------------------------------
void CellDataFilter::CopyAttributes(const CellDataFilter &other)
{
  _DataName   = other._DataName;
  _InputData  = other._InputData;
  _OutputData = _Output->GetCellData()->GetArray(_DataName.c_str());
}

// -----------------------------------------------------------------------------
CellDataFilter::CellDataFilter()
{
}

// -----------------------------------------------------------------------------
CellDataFilter::CellDataFilter(const CellDataFilter &other)
:
  MeshFilter(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
CellDataFilter &CellDataFilter::operator =(const CellDataFilter &other)
{
  if (this != &other) {
    MeshFilter::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
CellDataFilter::~CellDataFilter()
{
}

// =============================================================================
// Execution
// =============================================================================

// -----------------------------------------------------------------------------
void CellDataFilter::Initialize()
{
  // Initialize base class
  MeshFilter::Initialize();

  // Get input cell data array
  if (_InputData) {
    if (_DataName.empty()) {
      if (_InputData->GetName() == nullptr) {
        Throw(ERR_RuntimeError, __FUNCTION__, "Input data array must have a name or specify output data array name");
      }
      _DataName = _InputData->GetName();
    }
  } else if (_DataName.empty()) {
    _InputData = _Input->GetCellData()->GetScalars();
    Throw(ERR_RuntimeError, __FUNCTION__, "Input mesh has no active scalars cell data array, data array name required");
  } else {
    _InputData = _Input->GetCellData()->GetArray(_DataName.c_str());
    if (_InputData == nullptr) {
      _InputData = GetArrayByCaseInsensitiveName(_Input->GetCellData(), _DataName.c_str());
      if (_InputData == nullptr) {
        Throw(ERR_RuntimeError, __FUNCTION__, "Input mesh has no cell data array named ", _DataName);
      }
    }
  }

  // Add new output cell data array
  int attr = -1;
  vtkCellData * const outputCD = _Output->GetCellData();
  for (int i = 0; i < outputCD->GetNumberOfArrays(); ++i) {
    if (outputCD->GetArray(i) == _InputData) {
      attr = outputCD->IsArrayAnAttribute(i);
      break;
    }
  }
  _OutputData = NewCellArray(_DataName.c_str(), _InputData->GetNumberOfComponents(), _InputData->GetDataType());
  _OutputData->CopyComponentNames(_InputData);
  outputCD->RemoveArray(_DataName.c_str());
  const int idx = outputCD->AddArray(_OutputData);
  if (attr >= 0) outputCD->SetActiveAttribute(idx, attr);
}

// -----------------------------------------------------------------------------
void CellDataFilter::GetNodeNeighbors(int cellId, UnorderedSet<int> &cellIds) const
{
  cellIds.clear();
  unsigned short ncells;
  vtkIdType npts, *pts, *cells;
  _Output->GetCellPoints(cellId, npts, pts);
  for (vtkIdType i = 0; i < npts; ++i) {
    _Output->GetPointCells(pts[i], ncells, cells);
    for (unsigned short j = 0; j < ncells; ++j) {
      if (cells[j] != cellId) {
        cellIds.insert(cells[j]);
      }
    }
  }
}

// -----------------------------------------------------------------------------
void CellDataFilter::GetEdgeNeighbors(int cellId, UnorderedSet<int> &cellIds) const
{
  cellIds.clear();
  vtkIdType ptId1, ptId2;
  vtkNew<vtkGenericCell> cell;
  vtkNew<vtkIdList> nbrIds;
  nbrIds->Allocate(10);
  _Output->GetCell(cellId, cell.GetPointer());
  const int nedges = cell->GetNumberOfEdges();
  for (int edgeId = 0; edgeId < nedges; ++edgeId) {
    vtkCell * const edge = cell->GetEdge(edgeId);
    if (edge->IsLinear()) {
      const int npts = edge->GetNumberOfPoints();
      if (npts > 1) {
        ptId1 = edge->PointIds->GetId(0);
        for (int i = 1; i < npts; ++i, ptId1 = ptId2) {
          ptId2 = edge->PointIds->GetId(i);;
          _Output->GetCellEdgeNeighbors(cellId, ptId1, ptId2, nbrIds.GetPointer());
          for (vtkIdType j = 0; j < nbrIds->GetNumberOfIds(); ++j) {
            cellIds.insert(nbrIds->GetId(j));
          }
        }
      }
    }
  }
}


} // namespace mirtk
