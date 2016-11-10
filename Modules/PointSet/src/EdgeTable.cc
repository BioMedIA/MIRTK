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

#include "mirtk/EdgeTable.h"

#include "mirtk/Assert.h"
#include "mirtk/Memory.h"
#include "mirtk/Pair.h"
#include "mirtk/Algorithm.h" // Sort

#include "vtkGenericCell.h"


namespace mirtk {


// -----------------------------------------------------------------------------
void EdgeTable::CopyAttributes(const EdgeTable &other)
{
  _Mesh          = other._Mesh;
  _NumberOfEdges = other._NumberOfEdges;
}

// -----------------------------------------------------------------------------
EdgeTable::EdgeTable(vtkDataSet *mesh)
:
  _NumberOfEdges(0)
{
  Initialize(mesh);
}

// -----------------------------------------------------------------------------
EdgeTable::EdgeTable(const EdgeTable &other)
:
  GenericSparseMatrix(other)
{
  CopyAttributes(other);
}

// -----------------------------------------------------------------------------
EdgeTable &EdgeTable::operator =(const EdgeTable &other)
{
  if (this != &other) {
    GenericSparseMatrix::operator =(other);
    CopyAttributes(other);
  }
  return *this;
}

// -----------------------------------------------------------------------------
EdgeTable::~EdgeTable()
{
}

// -----------------------------------------------------------------------------
void EdgeTable::Initialize(vtkDataSet *mesh)
{
  if (mesh == nullptr) {
    this->Clear();
    return;
  }
  MIRTK_START_TIMING();

  const int       numPts   = static_cast<int>(mesh->GetNumberOfPoints());
  const vtkIdType numCells = mesh->GetNumberOfCells();
  Array<Entries>  entries(numPts);

  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

  Entries::const_iterator entry;
  vtkIdType cellId, edgeId, ptId1, ptId2;
  int numCellEdges, numEdgePts;
  bool new_edge;
  vtkCell *edge;

  _Mesh = mesh;
  _NumberOfEdges = 0;

  for (cellId = 0; cellId < numCells; ++cellId) {
    mesh->GetCell(cellId, cell);
    numCellEdges = cell->GetNumberOfEdges();
    for (edgeId = 0; edgeId < numCellEdges; ++edgeId) {
      edge = cell->GetEdge(edgeId);
      if (edge->IsLinear()) {
        numEdgePts = edge->GetNumberOfPoints();
        if (numEdgePts > 1) {
          ptId1 = edge->PointIds->GetId(0);
          for (int i = 1; i < numEdgePts; ++i, ptId1 = ptId2) {
            ptId2 = edge->PointIds->GetId(i);
            new_edge = true;
            for (entry = entries[ptId1].begin(); entry != entries[ptId1].end(); ++entry) {
              if (entry->first == ptId2) {
                new_edge = false;
                break;
              }
            }
            if (new_edge) {
              // Symmetric entries such that AdjacentPoints is efficient
              ++_NumberOfEdges; // edgeId + 1 such that entries are non-zero
              if (ptId1 == ptId2) {
                entries[ptId1].push_back(MakePair(static_cast<int>(ptId2), _NumberOfEdges));
              } else {
                entries[ptId1].push_back(MakePair(static_cast<int>(ptId2), _NumberOfEdges));
                entries[ptId2].push_back(MakePair(static_cast<int>(ptId1), _NumberOfEdges));
              }
            }
          }
        }
      } else {
        cerr << this->NameOfType() << "::Initialize: Only linear edges supported" << endl;
        exit(1);
      }
    }
  }

  for (auto &&col : entries) Sort(col);
  GenericSparseMatrix::Initialize(numPts, numPts, entries, true);

  // Reassign edge IDs -- same order as edges are visited by EdgeIterator!
  int i1, i2, j1, j2;
  const int *row = _Row;
  const int *col = _Col;
  if (_Layout == CRS) swap(row, col);

  int edgeIdPlusOne = 0;
  for (ptId2 = 0; ptId2 < numPts; ++ptId2) {
    for (i1 = col[ptId2], j1 = col[ptId2 + 1]; i1 < j1; ++i1) {
      ptId1 = row[i1];
      if (ptId1 > ptId2) break;
      i2 = col[ptId1];
      j2 = col[ptId1 + 1];
      while (i2 < j2 && row[i2] < ptId2) ++i2;
      mirtkAssert(i2 < j2 && row[i2] == ptId2, "matrix is symmetric");
      // one-based IDs to have non-zero entries in sparse matrix
      _Data[i1] = _Data[i2] = ++edgeIdPlusOne;
    }
  }
  mirtkAssert(edgeIdPlusOne == _NumberOfEdges, "edge ID reassigned is consistent");

  MIRTK_DEBUG_TIMING(5, "initialization of edge table");
}

// -----------------------------------------------------------------------------
void EdgeTable::Clear()
{
  GenericSparseMatrix::Clear();
  _NumberOfEdges = 0;
  _Mesh          = nullptr;
}

// -----------------------------------------------------------------------------
int EdgeTable::MaxNumberOfAdjacentPoints() const
{
  int n = 0;
  for (int i = 0; i < NumberOfPoints(); ++i) {
    n = max(n, NumberOfAdjacentPoints(i));
  }
  return n;
}


} // namespace mirtk
