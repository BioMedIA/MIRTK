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

#include <mirtkEdgeTable.h>

#include <mirtkAssert.h>
#include <mirtkMemory.h>
#include <mirtkPair.h>
#include <mirtkAlgorithm.h> // sort

#include <vtkSmartPointer.h>
#include <vtkGenericCell.h>
#include <vtkDataSet.h>


namespace mirtk {


// -----------------------------------------------------------------------------
EdgeTable::EdgeTable(vtkDataSet *mesh)
{
  if (mesh) Initialize(mesh);
}

// -----------------------------------------------------------------------------
EdgeTable::EdgeTable(const EdgeTable &other)
:
  GenericSparseMatrix(other)
{
}

// -----------------------------------------------------------------------------
EdgeTable &EdgeTable::operator =(const EdgeTable &other)
{
  GenericSparseMatrix::operator =(other);
  return *this;
}

// -----------------------------------------------------------------------------
EdgeTable::~EdgeTable()
{
}

// -----------------------------------------------------------------------------
void EdgeTable::Initialize(vtkDataSet *mesh)
{
  MIRTK_START_TIMING();

  _NumberOfEdges = 0;

  const int       numPts   = static_cast<int>(mesh->GetNumberOfPoints());
  const vtkIdType numCells = mesh->GetNumberOfCells();

  Entries *entries = new Entries[numPts];

  vtkSmartPointer<vtkGenericCell> cell = vtkSmartPointer<vtkGenericCell>::New();

  Entries::const_iterator entry;
  vtkIdType cellId, edgeId, ptId1, ptId2;
  int numCellEdges, numEdgePts;
  bool new_edge;
  vtkCell *edge;

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
                entries[ptId1].push_back(MakePair(ptId2, _NumberOfEdges));
              } else {
                entries[ptId1].push_back(MakePair(ptId2, _NumberOfEdges));
                entries[ptId2].push_back(MakePair(ptId1, _NumberOfEdges));
              }
            }
          }
        }
      } else {
        cerr << "WARNING: EdgeTable::Initialize: Only linear edges supported" << endl;
      }
    }
  }

  for (int i = 0; i < numPts; ++i) {
    sort(entries[i].begin(), entries[i].end());
  }

  GenericSparseMatrix::Initialize(numPts, numPts, entries, true);
  delete[] entries;

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


} // namespace mirtk
