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

#ifndef MIRTK_EdgeConnectivity_H
#define MIRTK_EdgeConnectivity_H

#include "mirtk/SparseMatrix.h"

#include "vtkDataSet.h"


namespace mirtk {


// Table of adjacent nodes
class EdgeTable;


/**
 * Table of edge-connectivities / n-connected neighbors
 *
 * The entries of this sparse matrix represent the minimum number of edges
 * that connect two given nodes. Unlike an EdgeTable, which stores a
 * unique ID for each edge connecting adjacent nodes (1-connected), this
 * table can be used to store a larger neighborhood of not only adjacent
 * nodes. It can, however, not be used to iterate the edges/paths in a
 * certain order. An irtkNeighborsTable is instead used to simply iterate
 * over the n-connected neighbors of a given node, e.g., to compute the
 * irtkMetricDistortion.
 *
 * Alternatively, this table can be used to query the edge-connectivity
 * for every pair of mesh nodes. If two nodes are not connected,
 # the edge-connectivity is zero.
 */
class EdgeConnectivity : public GenericSparseMatrix<int>
{
  mirtkObjectMacro(EdgeConnectivity);

  /// Maximum (considered) edge-connectivity
  mirtkReadOnlyAttributeMacro(int, Maximum);

public:

  /// Construct edge-connectivity table for given dataset
  EdgeConnectivity(vtkDataSet * = nullptr, int n = 3, const EdgeTable * = nullptr);

  /// Construct edge-connectivity table for given dataset
  EdgeConnectivity(vtkDataSet *, double r, const EdgeTable * = nullptr);

  /// Copy constructor
  EdgeConnectivity(const EdgeConnectivity &);

  /// Assignment operator
  EdgeConnectivity &operator =(const EdgeConnectivity &);

  /// Destructor
  virtual ~EdgeConnectivity();

  /// Initialize edge-connectivity table from given dataset
  void Initialize(vtkDataSet *, int n = 3, const EdgeTable * = nullptr);

  /// Initialize edge-connectivity table from given dataset
  void Initialize(vtkDataSet *, double r, const EdgeTable * = nullptr);

  /// Clear edge-connectivity table
  virtual void Clear();

  /// Number of nodes
  int NumberOfPoints() const;

  /// Get number of nodes with edge-connectivity less or equal to n
  int NumberOfConnectedPoints(int, int n = -1) const;

  /// Get number of adjacent nodes, i.e., nodes with edge-connectivity equal one
  int NumberOfAdjacentPoints(int) const;

  /// Access list of nodes with edge-connectivity less or equal to n (thread-safe)
  void GetConnectedPoints(int, int &, const int *&, int n = -1) const;

  /// Get start and end pointer into list of nodes with edge-connectivity
  /// less or equal to n (thread-safe)
  void GetConnectedPoints(int, const int *&, const int *&, int n = -1) const;

  /// Access list of adjacent (i.e. 1-connected) nodes (thread-safe)
  void GetAdjacentPoints(int, int &, const int *&) const;

  /// Get start and end pointer into list of adjacent (i.e. 1-connected) nodes (thread-safe)
  void GetAdjacentPoints(int, const int *&, const int *&) const;
};

////////////////////////////////////////////////////////////////////////////////
// Inline definitions
////////////////////////////////////////////////////////////////////////////////

// -----------------------------------------------------------------------------
inline int EdgeConnectivity::NumberOfPoints() const
{
  return Rows();
}

// -----------------------------------------------------------------------------
inline int EdgeConnectivity::NumberOfConnectedPoints(int ptId, int n) const
{
  if (n == 0) return 0;
  if (n <  0 || n == _Maximum) {
    if (_Layout == CCS) {
      return _Col[ptId+1] - _Col[ptId];
    } else {
      return _Row[ptId+1] - _Row[ptId];
    }
  } else {
    int count = 0;
    if (_Layout == CCS) {
      for (int i = _Col[ptId], j = _Col[ptId+1]; i != j; ++i) {
        if (_Data[i] > n) break;
        ++count;
      }
    } else {
      for (int i = _Row[ptId], j = _Row[ptId+1]; i != j; ++i) {
        if (_Data[i] > n) break;
        ++count;
      }
    }
    return count;
  }
}

// -----------------------------------------------------------------------------
inline int EdgeConnectivity::NumberOfAdjacentPoints(int ptId) const
{
  return NumberOfConnectedPoints(ptId, 1);
}

// -----------------------------------------------------------------------------
inline void EdgeConnectivity
::GetConnectedPoints(int ptId, const int *&begin, const int *&end, int n) const
{
  if (n == 0) {
    if (_Layout == CCS) {
      begin = end = _Row;
    } else {
      begin = end = _Col;
    }
  } else if (n < 0 || n == _Maximum) {
    if (_Layout == CCS) {
      begin = _Row + _Col[ptId];
      end   = _Row + _Col[ptId + 1];
    } else {
      begin = _Col + _Row[ptId];
      end   = _Col + _Row[ptId + 1];
    }
  } else {
    if (_Layout == CCS) {
      int i = _Col[ptId  ];
      int j = _Col[ptId+1];
      begin = end = _Row + i;
      while (i != j && *end <= n) ++i, ++end;
    } else {
      int i = _Row[ptId  ];
      int j = _Row[ptId+1];
      begin = end = _Col + i;
      while (i != j && *end <= n) ++i, ++end;
    }
  }
}

// -----------------------------------------------------------------------------
inline void EdgeConnectivity
::GetConnectedPoints(int ptId, int &numNbrPts, const int *&nbrPtIds, int n) const
{
  const int *end;
  GetConnectedPoints(ptId, nbrPtIds, end, n);
  numNbrPts = end - nbrPtIds;
}

// -----------------------------------------------------------------------------
inline void EdgeConnectivity
::GetAdjacentPoints(int ptId, int &numAdjPts, const int *&adjPtIds) const
{
  return GetConnectedPoints(ptId, numAdjPts, adjPtIds, 1);
}

// -----------------------------------------------------------------------------
inline void EdgeConnectivity
::GetAdjacentPoints(int ptId, const int *&begin, const int *&end) const
{
  return GetConnectedPoints(ptId, begin, end, 1);
}


} // namespace mirtk

#endif // MIRTK_EdgeConnectivity_H
