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

#include "mirtk/EdgeConnectivity.h"

#include "mirtk/Pair.h"
#include "mirtk/UnorderedSet.h"
#include "mirtk/Memory.h"
#include "mirtk/EdgeTable.h"
#include "mirtk/Parallel.h"
#include "mirtk/Profiling.h"
#include "mirtk/Algorithm.h" // sort
#include "mirtk/VtkMath.h"

#include "vtkSmartPointer.h"
#include "vtkGenericCell.h"


namespace mirtk {


// =============================================================================
// Auxiliary functors
// =============================================================================

namespace EdgeConnectivityUtils {


// -----------------------------------------------------------------------------
/// Determine edge-connectivity of nodes
struct ComputeEdgeConnectivity
{
  typedef GenericSparseMatrix<int>::Entries Entries;

  const EdgeTable *_EdgeTable;
  Entries         *_Connectivity;
  int              _Maximum;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    const int        *adjPts;
    int               numAdjPts, adjPtId;
    UnorderedSet<int> visited, set1, set2, *cur, *nxt;
    UnorderedSet<int>::const_iterator nbrPtId;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      cur = &set1;
      nxt = &set2;
      cur->clear();
      cur->insert(ptId);
      visited.clear();
      visited.insert(ptId);
      for (int c = 1; c <= _Maximum || _Maximum == 0; ++c) {
        nxt->clear();
        for (nbrPtId = cur->begin(); nbrPtId != cur->end(); ++nbrPtId) {
          _EdgeTable->GetAdjacentPoints(*nbrPtId, numAdjPts, adjPts);
          for (int i = 0; i < numAdjPts; ++i) {
            adjPtId = adjPts[i];
            if (visited.find(adjPtId) == visited.end()) {
              visited.insert(adjPtId);
              _Connectivity[ptId].push_back(MakePair(adjPtId, c));
              if (c < _Maximum || _Maximum == 0) {
                nxt->insert(adjPtId);
              }
            }
          }
        }
        if (nxt->empty()) break;
        swap(cur, nxt);
      }
    }
  }
};

// -----------------------------------------------------------------------------
/// Determine edge-connectivity of nodes
///
/// FIXME: Use geodesic distance instead of Euclidean distance.
struct ComputeEdgeConnectivityWithinRadius
{
  typedef GenericSparseMatrix<int>::Entries Entries;

  vtkDataSet      *_DataSet;
  const EdgeTable *_EdgeTable;
  Entries         *_Connectivity;
  double           _Radius2;

  void operator ()(const blocked_range<int> &ptIds) const
  {
    const int        *adjPts;
    int               numAdjPts, adjPtId;
    double            p0[3], p[3], dist2;
    UnorderedSet<int> visited, set1, set2, *cur, *nxt;
    UnorderedSet<int>::const_iterator nbrPtId;

    for (int ptId = ptIds.begin(); ptId != ptIds.end(); ++ptId) {
      cur = &set1;
      nxt = &set2;
      cur->clear();
      cur->insert(ptId);
      visited.clear();
      visited.insert(ptId);
      _DataSet->GetPoint(ptId, p0);
      for (int c = 1; !cur->empty(); ++c) {
        nxt->clear();
        for (nbrPtId = cur->begin(); nbrPtId != cur->end(); ++nbrPtId) {
          _EdgeTable->GetAdjacentPoints(*nbrPtId, numAdjPts, adjPts);
          for (int i = 0; i < numAdjPts; ++i) {
            adjPtId = adjPts[i];
            if (visited.find(adjPtId) == visited.end()) {
              visited.insert(adjPtId);
              _DataSet->GetPoint(adjPtId, p);
              dist2 = vtkMath::Distance2BetweenPoints(p0, p);
              if (dist2 <= _Radius2) {
                _Connectivity[ptId].push_back(MakePair(adjPtId, c));
                nxt->insert(adjPtId);
              }
            }
          }
        }
        swap(cur, nxt);
      }
    }
  }
};


} // namespace EdgeConnectivityUtils
using namespace EdgeConnectivityUtils;

// =============================================================================
// Construction/Destruction
// =============================================================================

// -----------------------------------------------------------------------------
EdgeConnectivity::EdgeConnectivity(vtkDataSet *mesh, int n, const EdgeTable *edgeTable)
:
  _Maximum(-1)
{
  if (mesh) Initialize(mesh, n, edgeTable);
}

// -----------------------------------------------------------------------------
EdgeConnectivity::EdgeConnectivity(vtkDataSet *mesh, double r, const EdgeTable *edgeTable)
  :
  _Maximum(-1)
{
  if (mesh) Initialize(mesh, r, edgeTable);
}

// -----------------------------------------------------------------------------
EdgeConnectivity::EdgeConnectivity(const EdgeConnectivity &other)
:
  GenericSparseMatrix(other),
  _Maximum(other._Maximum)
{
}

// -----------------------------------------------------------------------------
EdgeConnectivity &EdgeConnectivity::operator =(const EdgeConnectivity &other)
{
  if (this != &other) {
    GenericSparseMatrix::operator =(other);
    _Maximum = other._Maximum;
  }
  return *this;
}

// -----------------------------------------------------------------------------
EdgeConnectivity::~EdgeConnectivity()
{
}

// -----------------------------------------------------------------------------
void EdgeConnectivity::Initialize(vtkDataSet *mesh, int n, const EdgeTable *edgeTable)
{
  MIRTK_START_TIMING();
  if (n < 0) n = 0;

  const int numPts = static_cast<int>(mesh->GetNumberOfPoints());

  _Maximum = n;
  if (_Maximum == 0) {
    GenericSparseMatrix::Initialize(numPts, numPts);
    return;
  }

  EdgeTable _edgeTable;
  if (edgeTable == NULL || edgeTable->Rows() == 0) {
    _edgeTable.Initialize(mesh);
    edgeTable = &_edgeTable;
  }
  Array<Entries> entries(numPts);

  ComputeEdgeConnectivity eval;
  eval._EdgeTable    = edgeTable;
  eval._Maximum      = _Maximum;
  eval._Connectivity = entries.data();
  parallel_for(blocked_range<int>(0, numPts), eval);

  for (auto &&col : entries) Sort(col);
  GenericSparseMatrix::Initialize(numPts, numPts, entries, true);

  MIRTK_DEBUG_TIMING(5, "initialization of edge-connectivity table");
}

// -----------------------------------------------------------------------------
void EdgeConnectivity::Initialize(vtkDataSet *mesh, double r, const EdgeTable *edgeTable)
{
  MIRTK_START_TIMING();

  const int numPts = static_cast<int>(mesh->GetNumberOfPoints());

  if (r <= .0) {
    _Maximum = 0;
    GenericSparseMatrix::Initialize(numPts, numPts);
    return;
  }

  EdgeTable _edgeTable;
  if (edgeTable == NULL || edgeTable->Rows() == 0) {
    _edgeTable.Initialize(mesh);
    edgeTable = &_edgeTable;
  }
  Array<Entries> entries(numPts);

  ComputeEdgeConnectivityWithinRadius eval;
  eval._DataSet      = mesh;
  eval._EdgeTable    = edgeTable;
  eval._Radius2      = r * r;
  eval._Connectivity = entries.data();
  parallel_for(blocked_range<int>(0, numPts), eval);

  for (auto &&col : entries) Sort(col);
  GenericSparseMatrix::Initialize(numPts, numPts, entries, true);

  _Maximum = 0;
  for (int i = 0; i < _NNZ; ++i) {
    if (_Data[i] > _Maximum) _Maximum = _Data[i];
  }

  MIRTK_DEBUG_TIMING(5, "initialization of edge-connectivity table (r = " << r << ")");
}

// -----------------------------------------------------------------------------
void EdgeConnectivity::Clear()
{
  GenericSparseMatrix::Clear();
  _Maximum = 0;
}


} // namespace mirtk
